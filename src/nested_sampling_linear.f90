module nested_sampling_linear_module
    use utils_module,      only: flag_blank,flag_gestating,flag_waiting  
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSamplingL(loglikelihood,priors,settings) result(output_info)
        use priors_module,     only: prior,prior_log_volume
        use utils_module,      only: logzero,loginf,DBL_FMT,read_resume_unit,stdout_unit,write_dead_unit,calc_cholesky,calc_covmat
        use settings_module
        use utils_module,      only: logsumexp
        use read_write_module, only: write_resume_file,write_posterior_file,write_phys_live_points
        use feedback_module
        use evidence_module,   only: infer_evidence,KeetonEvidence
        use chordal_module,    only: SliceSampling,GradedSliceSampling,AdaptiveParallelSliceSampling
        use random_module,     only: random_integer

        use grades_module,     only: calc_graded_choleskys

        implicit none

        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        type(prior), dimension(:), intent(in) :: priors
        type(program_settings), intent(in) :: settings

        ! Output of the program
        ! 1) log(evidence)
        ! 2) error(log(evidence))
        ! 3) ndead
        ! 4) number of likelihood calls
        ! 5) log(evidence) + log(prior volume)
        double precision, dimension(5) :: output_info



        !> This is a very important array. live_points(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal,settings%nstack) :: live_points

        double precision, dimension(settings%nDims,settings%nDims) :: covmat
        double precision, dimension(settings%nDims,settings%nDims) :: cholesky
        double precision, dimension(settings%nDims,settings%nDims,settings%grades%min_grade:settings%grades%max_grade) :: choleskys

        double precision, dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior) :: posterior_array
        integer :: nposterior

        logical :: more_samples_needed

        ! The new-born baby points
        double precision,    dimension(settings%nTotal,settings%chain_length)   :: baby_points
        double precision :: baby_likelihood

        ! The recently dead point
        double precision :: late_likelihood

        ! Point to seed a new one from
        double precision,    dimension(settings%nTotal)   :: seed_point


        ! Evidence info
        double precision, allocatable, dimension(:) :: evidence_vec


        logical :: resume=.false.
        ! Means to be calculated
        double precision :: mean_likelihood_calls
        integer :: total_likelihood_calls

        integer :: ndead

        double precision, dimension(settings%max_ndead) :: dead_likes

        integer :: stack_size
        logical :: first_loop


        call write_opening_statement(settings) 

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

        ! Allocate the evidence vector 
        more_samples_needed = KeetonEvidence(settings,baby_likelihood,late_likelihood,ndead,evidence_vec)

        !~~~ (i) Generate Live Points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(resume) then
            ! If there is a resume file present, then load the live points from that
            open(read_resume_unit,file=trim(settings%file_root)//'.resume',action='read')

            read(read_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_points
            read(read_resume_unit,'(<size(evidence_vec)>E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
            read(read_resume_unit,'(I)') ndead
            read(read_resume_unit,'(E<DBL_FMT(1)>.<DBL_FMT(2)>)') mean_likelihood_calls
            read(read_resume_unit,'(I)') total_likelihood_calls
            read(read_resume_unit,'(I)') nposterior
            read(read_resume_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)

            close(read_resume_unit)

        else !(not resume)
            call write_started_generating(settings%feedback)

            ! Zero the array
            live_points = 0

            ! Set them all to blank initially
            live_points(settings%point_type,:) = blank_type

            ! Otherwise generate them anew:
            live_points = GenerateLivePoints(loglikelihood,priors,settings,settings%nlive)

            stack_size=settings%nlive

            call write_finished_generating(settings%feedback) !Flag to note that we're done generating

            ! Otherwise compute the average loglikelihood and initialise the evidence vector accordingly
            evidence_vec = logzero
            evidence_vec(4) = logsumexp(live_points(settings%l0,:)) - log(settings%nlive+0d0)

            mean_likelihood_calls = 1d0
            total_likelihood_calls = settings%nlive

            ! Otherwise no dead points originally
            ndead = 0

            nposterior = 0
            ! set all of the loglikelihoods and logweights to be zero initially
            posterior_array(1:2,:) = logzero

            ! set the posterior coordinates to be zero initially
            posterior_array(3:,:) = 0d0
        end if !(resume)



        ! Initialise the late likelihood
        late_likelihood = minval(live_points(settings%l0,:stack_size), mask=nint(live_points(settings%point_type,:stack_size))==live_type) 


        ! Write a resume file before we start
        if(settings%write_resume) call write_resume_file(settings,live_points,evidence_vec,ndead,mean_likelihood_calls,total_likelihood_calls,nposterior,posterior_array) 


        !======= 2) Main loop body =====================================

        call write_started_sampling(settings%feedback)

        ! definitely more samples needed than this
        more_samples_needed = .true.

        do while ( more_samples_needed )

            ! (1) Update the covariance matrix of the distribution of live points
            if(mod(ndead,settings%nlive) .eq.0) then
                select case(settings%sampler)

                case(sampler_covariance)
                    ! Calculate the covariance matrix
                    covmat = calc_covmat( live_points(settings%h0:settings%h1,:stack_size), settings%nDims,stack_size )
                    ! Calculate the cholesky decomposition
                    cholesky = calc_cholesky(covmat,settings%nDims)

                case(sampler_graded_covariance)
                    ! Calculate the covariance matrix
                    covmat = calc_covmat( live_points(settings%h0:settings%h1,:stack_size), settings%nDims,stack_size )
                    ! Calculate the graded cholesky matrices
                    choleskys = calc_graded_choleskys(covmat,settings%nDims,settings%grades)

                end select
            end if

            ! (2) Generate a new set of baby points
            ! Select a seed point for the generator
            first_loop = .true.
            do while (seed_point(settings%l0)<late_likelihood .or. first_loop)
                seed_point = live_points(:,random_integer(stack_size))
                ! Record the likelihood bound which this seed will generate from
                seed_point(settings%l1) = late_likelihood
                first_loop=.false.
            end do

            ! Generate a new set of points within the likelihood bound of the late point
            select case(settings%sampler)

            case(sampler_covariance)
                baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)

            case(sampler_graded_covariance)
                baby_points = GradedSliceSampling(loglikelihood,priors,settings,choleskys,seed_point)

            case(sampler_adaptive_parallel)
                baby_points = AdaptiveParallelSliceSampling(loglikelihood,priors,settings,live_points(:,:stack_size),seed_point)

            end select

            ! The new likelihood is the last point
            baby_likelihood  = baby_points(settings%l0,settings%chain_length)

            if(baby_likelihood>late_likelihood) then

                ! (3) Calculate the new evidence (and check to see if we're accurate enough)
                more_samples_needed = KeetonEvidence(settings,baby_likelihood,late_likelihood,ndead,evidence_vec)

                ! Record the loglikelihoods if we're inferring the evidence
                if (settings%infer_evidence) dead_likes(ndead) = late_likelihood

                ! (4) Update the stack of live points and the posterior array
                !     This function does multiple things:
                !     1) Insert baby_points into live_points
                !     2) Remove points from live_point that have died this round
                !     3) Add any of these which are at a high enough likelihood to the posterior_array
                !     4) re-calculate stack_size and nposterior
                !     5) update the late_likelihood
                !     6) Update ndead
                more_samples_needed = more_samples_needed .or. update_stacks(settings,baby_points,live_points,stack_size,posterior_array,nposterior,late_likelihood,evidence_vec(1),ndead)

                ! (5) Feedback to command line every nlive iterations
                if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                    write(stdout_unit,'("ndead     = ", I20                  )') ndead
                    write(stdout_unit,'("stack size= ", I20, "/", I20        )') stack_size, settings%nstack
                    write(stdout_unit,'("nposterior= ", I20                  )') nposterior
                    !write(stdout_unit,'("efficiency= ", F20.2                )') mean_likelihood_calls
                    write(stdout_unit,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                    write(stdout_unit,'("")')
                end if

                ! (6) Update the resume and posterior files every update_resume iterations, or at program termination
                if (mod(ndead,settings%update_resume) .eq. 0 .or.  more_samples_needed==.false.)  then
                    if(settings%write_resume) call write_resume_file(settings,live_points,evidence_vec,ndead,mean_likelihood_calls,total_likelihood_calls,nposterior,posterior_array) 
                    if(settings%calculate_posterior) call write_posterior_file(settings,posterior_array,evidence_vec(1),nposterior)  
                    if(settings%write_live) call write_phys_live_points(settings,live_points(:,:stack_size),stack_size)
                end if

                ! If we've put a limit on the maximum number of iterations, then
                ! check to see if we've reached this
                if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

            end if




        end do ! End main loop


        if(settings%infer_evidence) call infer_evidence(settings,dead_likes(:ndead))

        ! Create the output array
        ! (1) log evidence
        ! (2) Error in the log evidence
        ! (3) Number of dead points
        ! (4) Number of likelihood calls
        ! (5) log(evidence * prior volume)
        output_info(1) = evidence_vec(1) - 0.5d0*log(1+exp(evidence_vec(2)-2*evidence_vec(1)))
        output_info(2) = sqrt(log(1+exp(evidence_vec(2)-2*evidence_vec(1))))
        output_info(3) = ndead
        output_info(4) = total_likelihood_calls
        output_info(5) = output_info(1)+prior_log_volume(priors)


        call write_final_results(output_info,settings%feedback,priors)


    end function NestedSamplingL




    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePoints(loglikelihood,priors,settings,nlive) result(live_points)
        use priors_module,    only: prior
        use settings_module,  only: program_settings,live_type
        use random_module,    only: random_reals
        use utils_module,     only: logzero
        use calculate_module, only: calculate_point

        implicit none

        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        !> The prior information
        type(prior), dimension(:), intent(in) :: priors

        !> Program settings
        type(program_settings), intent(in) :: settings

        !> The number of points to be generated
        integer, intent(in) :: nlive

        !live_points(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal,nlive) :: live_points

        ! Loop variable
        integer i_live

        ! initialise live points at zero
        live_points = 0d0

        do i_live=1,nlive

            ! Generate a random coordinate
            live_points(:,i_live) = random_reals(settings%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( loglikelihood, priors, live_points(:,i_live), settings )

        end do

        ! Set the number of likelihood calls for each point to 1
        live_points(settings%nlike,:) = 1

        ! Set the initial trial values of the chords as the diagonal of the hypercube
        live_points(settings%last_chord,:) = sqrt(settings%nDims+0d0)

        ! Set the likelihood contours to logzero for now
        live_points(settings%l1,:) = logzero

        ! These are all true live points
        live_points(settings%point_type,:) = live_type


    end function GenerateLivePoints



    function GenerateSeed(settings,live_points,seed_pos) result(seed_point)
        use settings_module,   only: program_settings
        use random_module,     only: random_integer
        implicit none
        type(program_settings), intent(in) :: settings
        double precision, intent(inout), dimension(settings%nTotal,settings%nstack) :: live_points

        integer, intent(out) :: seed_pos

        ! Point to seed a new one from
        double precision,    dimension(settings%nTotal)   :: seed_point

        seed_pos =random_integer(settings%nlive)
        seed_point = live_points(:,seed_pos)

    end function GenerateSeed



    function update_stacks(settings,baby_points,live_points,stack_size,posterior_array,nposterior,late_likelihood,evidence,ndead) result(more_samples_needed)
        use settings_module,   only: program_settings,blank_type,live_type
        implicit none
        type(program_settings), intent(in)                                                                           :: settings
        double precision,       intent(in),    dimension(settings%nTotal,settings%chain_length)                      :: baby_points
        double precision,       intent(inout), dimension(settings%nTotal,settings%nstack)                            :: live_points
        integer,                intent(inout)                                                                        :: stack_size
        double precision,       intent(inout), dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior) :: posterior_array
        integer,                intent(inout)                                                                        :: nposterior
        double precision,       intent(inout)                                                                        :: late_likelihood
        double precision,       intent(in)                                                                           :: evidence
        integer,                intent(inout)                                                                        :: ndead

        logical :: more_samples_needed

        integer :: late_index(1)
        integer :: insertion_index(1)

        integer :: i_live

        double precision :: late_logweight

        double precision, dimension(settings%nDims+settings%nDerived+2) :: posterior_point

        double precision :: lognmax_posterior
        double precision :: max_logweight


        late_logweight = (ndead-1)*log(settings%nlive+0d0) - ndead*log(settings%nlive+1d0)                

        ! Start by finding the original lowest likelihood live point (about to be deleted)
        late_index = minloc(live_points(settings%l0,:stack_size), mask=nint(live_points(settings%point_type,:stack_size))==live_type)

        ! Update the late likelihood
        late_likelihood = live_points(settings%l0,late_index(1))

        ! Add the discarded point to the posterior array
        posterior_point(1)  = live_points(settings%l0,late_index(1)) + late_logweight
        posterior_point(2)  = live_points(settings%l0,late_index(1))
        posterior_point(2+1:2+settings%nDims) = live_points(settings%p0:settings%p1,late_index(1))
        posterior_point(2+settings%nDims+1:2+settings%nDims+settings%nDerived) = live_points(settings%d0:settings%d1,late_index(1))

        ! Replace the late point with the new baby point
        live_points(:,late_index(1)) = baby_points(:,settings%chain_length)


        nposterior=nposterior+1
        posterior_array(:,nposterior) = posterior_point
 
        ! Add the remaining baby points to the end of the array, and update the stack size
        stack_size=stack_size+settings%chain_length-1
        live_points(:,stack_size-settings%chain_length+2:stack_size) = baby_points(:,:settings%chain_length-1)

        ! Now run through the stack and strip out any points that are less
        ! than the new late_likelihood, replacing them with points drawn from
        ! the end 

        i_live=1
        do while(i_live<=stack_size)
            if( live_points(settings%l0,i_live) < late_likelihood ) then

                ! Add the discarded point to the posterior array
                posterior_point(1)  = live_points(settings%l0,i_live) + late_logweight
                posterior_point(2)  = live_points(settings%l0,i_live)
                posterior_point(2+1:2+settings%nDims) = live_points(settings%p0:settings%p1,i_live)
                posterior_point(2+settings%nDims+1:2+settings%nDims+settings%nDerived) = live_points(settings%d0:settings%d1,i_live)

                nposterior=nposterior+1
                posterior_array(:,nposterior) = posterior_point

                ! Overwrite the discarded point with a point from the end...
                live_points(:,i_live) = live_points(:,stack_size)
                ! ...and reduce the stack size
                stack_size=stack_size-1
            else
                i_live=i_live+1
            end if
        end do

        if(nposterior>settings%nmax_posterior) write(*,*) 'over the top'

        ! Clean out the posterior array

        ! Find the maximum weighted posterior point
        max_logweight = maxval(posterior_array(1,:nposterior))

        lognmax_posterior = log(settings%nmax_posterior+0d0)

        i_live=1
        do while(i_live<=nposterior)
            if( posterior_array(1,i_live) - max_logweight + lognmax_posterior < 0 ) then
                ! Overwrite the discarded point with a point from the end...
                posterior_array(:,i_live) = posterior_array(:,nposterior)
                ! ...and reduce the stack size
                nposterior=nposterior-1
            else
                i_live=i_live+1
            end if
        end do

        live_points(settings%last_chord,:) = live_points(settings%last_chord,:)/  (1d0+1d0/(settings%nDims*settings%nlive) )


        ! Find the new late likelihood
        late_likelihood = minval(live_points(settings%l0,:stack_size), mask=nint(live_points(settings%point_type,:stack_size))==live_type) 

        ! Increment the number of dead points
        ndead=ndead+1

        ! Test to see if we need more samples
        !more_samples_needed = any(posterior_points(1,:) - max_logweight+lognmax_posterior > 0 )
        more_samples_needed = .false.

    end function update_stacks





end module nested_sampling_linear_module
