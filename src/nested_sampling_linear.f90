module nested_sampling_linear_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSamplingL(loglikelihood,priors,M,settings)
        use priors_module,     only: prior
        use model_module,      only: model
        use utils_module,      only: logzero,loginf,DBL_FMT,read_resume_unit,stdout_unit,write_dead_unit
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp
        use read_write_module, only: write_resume_file,write_posterior_file
        use feedback_module
        use random_module,     only: random_integer,random_direction

        implicit none

        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings
        type(prior), dimension(:), intent(in) :: priors



        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal,settings%nlive) :: live_data

        double precision, allocatable, dimension(:,:) :: posterior_array
        double precision, dimension(settings%nDims+settings%nDerived+2) :: posterior_point
        integer :: nposterior
        integer :: insertion_index(1)
        integer :: late_index(1)

        logical :: more_samples_needed

        ! The new-born baby point
        double precision,    dimension(settings%nTotal)   :: baby_point
        double precision :: baby_likelihood

        ! The recently dead point
        double precision,    dimension(settings%nTotal)   :: late_point
        double precision :: late_likelihood
        double precision :: late_logweight

        ! Point to seed a new one from
        double precision,    dimension(settings%nTotal)   :: seed_point


        ! Evidence info
        double precision, dimension(6) :: evidence_vec


        logical :: resume=.false.
        ! Means to be calculated
        double precision :: mean_likelihood_calls
        integer :: total_likelihood_calls

        integer :: ndead

        double precision :: lognlive 
        double precision :: lognlivep1 
        double precision :: logminimumweight


        double precision :: nhats(M%nDims,settings%num_chords)



        call write_opening_statement(M,settings) 

        ! Check to see whether there's a resume file present, and record in the
        ! variable 'resume'
        inquire(file=trim(settings%file_root)//'.resume',exist=resume)

        ! Check if we actually want to resume
        resume = settings%read_resume .and. resume

        if(resume .and. settings%feedback>=0) write(stdout_unit,'("Resuming from previous run")')


        !======= 1) Initialisation =====================================
        ! (i)   generate initial live points by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables

        !~~~ (i) Generate Live Points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(resume) then
            ! If there is a resume file present, then load the live points from that
            open(read_resume_unit,file=trim(settings%file_root)//'.resume',action='read')
            ! Read the live data
            read(read_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
        else !(not resume)
            call write_started_generating(settings%feedback)

            ! Otherwise generate them anew:
            live_data = GenerateLivePoints(loglikelihood,priors,M,settings%nlive)

            call write_finished_generating(settings%feedback) !Flag to note that we're done generating
        end if !(resume)






        !~~~ (ii) Initialise all variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! There are several variables used throughout the rest of the
        ! algorithm that need to be initialised here
        !  (a) evidence_vec           | Vector containing the evidence, its error, and any other 
        !                             |  things that need to be accumulated over the run.
        !                             |  we need to initialise its sixth argument.
        !  (b) ndead                  | Number of iterations/number of dead points
        !  (c) mean_likelihood_calls  | Mean number of likelihood calls over the past nlive iterations
        !  (d) posterior_array        | Array of weighted posterior points

        ! (a)
        if(resume) then
            ! If resuming, get the accumulated stats to calculate the
            ! evidence from the resume file
            read(read_resume_unit,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
        else !(not resume) 
            ! Otherwise compute the average loglikelihood and initialise the evidence vector accordingly
            evidence_vec = logzero
            evidence_vec(4) = logsumexp(live_data(M%l0,:)) - log(settings%nlive+0d0)
        end if !(resume) 

        ! (b) get number of dead points
        if(resume) then
            ! If resuming, then get the number of dead points from the resume file
            read(read_resume_unit,'(I)') ndead
        else !(not resume) 
            ! Otherwise no dead points originally
            ndead = 0
        end if !(resume) 

        ! (c) initialise the mean and total number of likelihood calls
        if(resume) then
            ! If resuming, then get the mean likelihood calls from the resume file
            read(read_resume_unit,'(E<DBL_FMT(1)>.<DBL_FMT(2)>)') mean_likelihood_calls
            ! Also get the total likelihood calls
            read(read_resume_unit,'(I)') total_likelihood_calls
        else
            mean_likelihood_calls = 1d0
            total_likelihood_calls = settings%nlive
        end if


        ! (d) Posterior array

        allocate(posterior_array(settings%nDims+settings%nDerived+2,settings%nmax_posterior))
        nposterior = 0
        ! set all of the loglikelihoods and logweights to be zero initially
        posterior_array(1:2,:) = logzero

        ! set the posterior coordinates to be zero initially
        posterior_array(3:,:) = 0d0

        if(resume) then
            ! Read the actual number we've used so far
            read(read_resume_unit,'(I)') nposterior
            !...followed by the posterior array itself
            read(read_resume_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)
        end if !(resume) 

        ! Close the resume file if we've openend it
        if(resume) close(read_resume_unit)

        ! Calculate these global variables so we don't need to again
        lognlive         = log(settings%nlive+0d0)
        lognlivep1       = log(settings%nlive+1d0)
        logminimumweight = log(settings%minimum_weight)



        ! Write a resume file before we start
        if(settings%write_resume) call write_resume_file(settings,M,live_data,evidence_vec,ndead,mean_likelihood_calls,total_likelihood_calls,nposterior,posterior_array) 

        ! Open a dead points file if desired
        if(settings%save_all) open(write_dead_unit,file=trim(settings%file_root)//'_dead.dat',action='write') 


        !======= 2) Main loop body =====================================

        call write_started_sampling(settings%feedback)

        ! definitely more samples needed than this
        more_samples_needed = .true.

        do while ( more_samples_needed )

            ! (1) Get the late point

            ! Find the point with the lowest likelihood...
            late_index = minloc(live_data(settings%l0,:))
            ! ...and save it.
            late_point = live_data(:,late_index(1))
            ! Get the likelihood contour
            late_likelihood = late_point(settings%l0)
            ! Calculate the late logweight
            late_logweight = (ndead-1)*lognlive - ndead*lognlivep1 

            ! (2) Generate a new baby point
            ! Select a seed point for the generator
            !  -excluding the points which have likelihoods equal to the
            !   loglikelihood bound
            seed_point(settings%l0)=late_likelihood
            do while (seed_point(settings%l0)<=late_likelihood)
                ! get a random number in [1,nlive]
                ! get this point from live_data 
                seed_point = live_data(:,random_integer(settings%nlive))
            end do

            ! Record the likelihood bound which this seed will generate from
            seed_point(settings%l1) = late_likelihood

            ! Generate a set of directions
            ! do i_dim=1,M%nDims
            !     nhat2 = 0
            !     do while(nhat2==0)
            !         transformation_matrix(:,i_dim) = live_data(M%h0:M%h1,random_integer(settings%nlive)) - live_data(M%h0:M%h1,random_integer(settings%nlive))
            !         nhat2 = dot_product(transformation_matrix(:,i_dim),transformation_matrix(:,i_dim))
            !     end do
            !     transformation_matrix(:,i_dim) = transformation_matrix(:,i_dim) / sqrt(nhat2)
            ! end do

            call settings%generate_directions(M,live_data,nhats)

            ! Generate a new point within the likelihood bound of the late point
            baby_point = settings%sampler(loglikelihood,priors,seed_point, M)
            baby_likelihood  = baby_point(settings%l0)


            ! (3) Insert the baby point into the set of live points (over the
            !     old position of the dead points

            ! Insert the baby point over the late point
            live_data(:,late_index(1)) = baby_point

            ! record that we have a new dead point
            ndead = ndead + 1

            ! If we've put a limit on the maximum number of iterations, then
            ! check to see if we've reached this
            if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

            ! (4) Calculate the new evidence (and check to see if we're accurate enough)
            call settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed,evidence_vec)




            ! (5) Update the set of weighted posteriors
            if( settings%calculate_posterior .and. late_point(settings%l0) + late_logweight - evidence_vec(1) > logminimumweight ) then
                ! If the late point has a sufficiently large weighting, then we
                ! should add it to the set of saved posterior points

                ! calculate a new point for insertion
                posterior_point(1)  = late_point(M%l0) + late_logweight
                posterior_point(2)  = late_point(M%l0)
                posterior_point(2+1:2+M%nDims) = late_point(M%p0:M%p1)
                posterior_point(2+M%nDims+1:2+M%nDims+M%nDerived) = late_point(M%d0:M%d1)

                if(nposterior<settings%nmax_posterior) then
                    ! If we're still able to use a restricted array,

                    ! Find the closest point in the array which is beneath the minimum weight
                    insertion_index = minloc(posterior_array(1,:nposterior),mask=posterior_array(1,:nposterior)<logminimumweight+evidence_vec(1))

                    if(insertion_index(1)==0) then
                        ! If there are no points to overwrite, then we should
                        ! expand the available storage array
                        nposterior=nposterior+1
                        posterior_array(:,nposterior) = posterior_point
                    else
                        ! Otherwise overwrite the 
                        posterior_array(:,insertion_index(1)) = posterior_point
                    end if

                else
                    ! Otherwise we have to overwrite the smallest element
                    insertion_index = minloc(posterior_array(1,:nposterior))
                    posterior_array(:,insertion_index(1)) = posterior_point
                end if

            end if

            live_data(M%last_chord,:) = live_data(M%last_chord,:)/  (1d0+1d0/(M%nDims*settings%nlive) )


            ! (6) Command line feedback

            ! update the mean number of likelihood calls
            mean_likelihood_calls = mean_likelihood_calls + (baby_point(settings%nlike) - late_point(settings%nlike) ) / (settings%nlive + 0d0)

            ! update the total number of likelihood calls
            total_likelihood_calls = total_likelihood_calls + baby_point(settings%nlike)

            ! update ndead file if we're that way inclined
            if(settings%save_all) write(write_dead_unit,'(<2*M%nDims+M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>,I8,3E<DBL_FMT(1)>.<DBL_FMT(2)>)')  late_point(M%h0:M%h1), late_point(M%p0:M%p1), late_point(M%d0:M%d1), late_point(M%nlike),late_point(M%last_chord), late_point(M%l0:M%l1)


            ! Feedback to command line every nlive iterations
            if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                write(stdout_unit,'("ndead     = ", I20                  )') ndead
                write(stdout_unit,'("efficiency= ", F20.2                )') mean_likelihood_calls
                write(stdout_unit,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                write(stdout_unit,'("")')
            end if

            ! (7) Update the resume and posterior files every update_resume iterations, or at program termination
            if (mod(ndead,settings%update_resume) .eq. 0 .or.  more_samples_needed==.false.)  then
                if(settings%write_resume) call write_resume_file(settings,M,live_data,evidence_vec,ndead,mean_likelihood_calls,total_likelihood_calls,nposterior,posterior_array) 
                if(settings%calculate_posterior) call write_posterior_file(settings,M,posterior_array,evidence_vec(1),nposterior)  
            end if

        end do ! End main loop


        if(settings%save_all) close(write_dead_unit) 

        call write_final_results(M,evidence_vec,ndead,total_likelihood_calls,settings%feedback)

    end subroutine NestedSamplingL




    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePoints(loglikelihood,priors,M,nlive) result(live_data)
        use priors_module,   only: prior
        use model_module,    only: model, calculate_point
        use random_module,   only: random_reals
        use utils_module,    only: logzero

        implicit none

        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface


        !> The model details (loglikelihood, priors, ndims etc...)
        type(model), intent(in) :: M

        !> The prior information
        type(prior), dimension(:), intent(in) :: priors

        !> The number of points to be generated
        integer, intent(in) :: nlive

        !live_data(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal,nlive) :: live_data

        ! Loop variable
        integer i_live

        ! initialise live points at zero
        live_data = 0d0

        do i_live=1,nlive

            ! Generate a random coordinate
            live_data(:,i_live) = random_reals(settings%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( loglikelihood, priors,M, live_data(:,i_live) )

        end do

        ! Set the number of likelihood calls for each point to 1
        live_data(settings%nlike,:) = 1

        ! Set the initial trial values of the chords as the diagonal of the hypercube
        live_data(settings%last_chord,:) = sqrt(settings%nDims+0d0)


    end function GenerateLivePoints




end module nested_sampling_linear_module
