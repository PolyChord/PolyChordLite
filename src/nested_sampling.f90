module nested_sampling_module
    use mpi

    implicit none

    integer, parameter :: RUNTAG=0
    integer, parameter :: ENDTAG=1


    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood,priors,settings,mpi_communicator) result(output_info)
        use priors_module,     only: prior,prior_log_volume
        use utils_module,      only: logzero,loginf,DBL_FMT,read_resume_unit,stdout_unit,write_dead_unit,calc_cholesky,calc_covmat,TwoPi
        use settings_module
        use utils_module,      only: logsumexp
        use read_write_module, only: write_resume_file,write_posterior_file,write_phys_live_points
        use feedback_module
        use evidence_module,   only: KeetonEvidence
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: SNN_clustering,NN_clustering

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

        integer, intent(in) :: mpi_communicator

        ! Output of the program
        ! 1) log(evidence)
        ! 2) error(log(evidence))
        ! 3) ndead
        ! 4) number of likelihood calls
        ! 5) log(evidence) + log(prior volume)
        double precision, dimension(5) :: output_info



        !> This is a very important array.
        double precision, dimension(settings%nTotal,settings%nlive) :: live_points
        double precision, dimension(settings%nTotal,settings%nstack) :: phantom_points

        double precision, dimension(settings%nDims,settings%nDims) :: covmat
        double precision, dimension(settings%nDims,settings%nDims) :: cholesky

        double precision, dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior) :: posterior_array
        integer :: nposterior

        logical :: more_samples_needed

        ! The new-born baby points
        double precision, dimension(settings%nTotal,settings%num_babies)   :: baby_points
        double precision :: baby_likelihood

        ! The recently dead point
        double precision :: late_likelihood

        ! Point to seed a new one from
        double precision,    dimension(settings%nTotal)   :: seed_point


        ! Evidence info
        double precision, allocatable, dimension(:) :: evidence_vec


        logical :: baby_made
        logical :: resume=.false.
        ! Means to be calculated
        double precision :: mean_likelihood_calls
        integer :: total_likelihood_calls

        integer :: ndead

        integer :: stack_size
        logical :: first_loop

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        integer :: send_start
        integer :: nprocs
        integer :: myrank
        integer :: root
        logical :: linear_mode
        integer :: mpierror

        integer :: i
        double precision :: mu1(2), mu2(2), radius,sigma


        integer :: clusters(settings%nlive)
        double precision, dimension(settings%nlive,settings%nlive) :: similarity_matrix












        ! Get the number of MPI procedures
        call MPI_COMM_SIZE(mpi_communicator, nprocs, mpierror)
        send_start=nprocs-1
        linear_mode = nprocs==1

        ! Get the MPI label of the current processor
        call MPI_COMM_RANK(mpi_communicator, myrank, mpierror)

        ! Assign the root
        call MPI_ALLREDUCE(myrank,root,1,MPI_INTEGER,MPI_MIN,mpi_communicator,mpierror)

        if(myrank==root) then
            call write_opening_statement(settings) 

            ! Allocate the evidence vector 
            more_samples_needed = KeetonEvidence(settings,baby_likelihood,late_likelihood,ndead,evidence_vec)
        end if




        !======= 1) Initialisation =====================================
        ! (i)   generate initial live points by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables

        !~~~ (i) Generate Live Points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Check if we actually want to resume
        inquire(file=trim(settings%file_root)//'.resume',exist=resume)
        resume = settings%read_resume .and. resume

        if(resume) then
            if(myrank==root) then
                ! Check to see whether there's a resume file present, and record in the
                ! variable 'resume'
                if(settings%feedback>=0) write(stdout_unit,'("Resuming from previous run")')

                ! If there is a resume file present, then load the live points from that
                open(read_resume_unit,file=trim(settings%file_root)//'.resume',action='read')

                read(read_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_points
                read(read_resume_unit,'(I)') stack_size
                read(read_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') phantom_points(:,:stack_size)
                read(read_resume_unit,'(<size(evidence_vec)>E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
                read(read_resume_unit,'(I)') ndead
                read(read_resume_unit,'(I)') total_likelihood_calls
                read(read_resume_unit,'(I)') nposterior
                read(read_resume_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)

                close(read_resume_unit)
            endif ! only root

        else !(not resume)
            ! Otherwise generate them anew:
            if(linear_mode) then
                live_points = GenerateLivePointsL(loglikelihood,priors,settings)
            else
                live_points = GenerateLivePointsP(loglikelihood,priors,settings,mpi_communicator,root)
            end if

            if(myrank==root) then

                stack_size=settings%nlive

                call write_finished_generating(settings%feedback) !Flag to note that we're done generating

                ! Otherwise compute the average loglikelihood and initialise the evidence vector accordingly
                evidence_vec = logzero
                evidence_vec(4) = logsumexp(live_points(settings%l0,:settings%nlive)) - log(settings%nlive+0d0)

                total_likelihood_calls = sum(live_points(settings%nlike,:stack_size))

                ! Otherwise no dead points originally
                ndead = 0

                nposterior = 0
                ! set all of the loglikelihoods and logweights to be zero initially
                posterior_array(1:2,:) = logzero

                ! set the posterior coordinates to be zero initially
                posterior_array(3:,:) = 0d0
            endif ! only root

        end if !(resume/not resume)


        if(myrank==root) then


            ! Initialise the late likelihood
            late_likelihood = minval(live_points(settings%l0,:))


            ! Write a resume file before we start
            if(settings%write_resume) call write_resume_file(settings,live_points,stack_size,phantom_points,evidence_vec,ndead,total_likelihood_calls,nposterior,posterior_array) 

            ! Calculate the covariance matrix
            covmat = calc_covmat( (/ live_points(settings%h0:settings%h1,:),phantom_points(settings%h0:settings%h1,:stack_size)/), settings%nDims,stack_size+settings%nlive )
            ! Calculate the cholesky decomposition
            cholesky = calc_cholesky(covmat,settings%nDims)


            !======= 2) Main loop body =====================================

            call write_started_sampling(settings%feedback)

            ! definitely more samples needed than this
            more_samples_needed = .true.

            do while ( more_samples_needed )

                ! (1) Update the covariance matrix of the distribution of live points
                if(mod(ndead,settings%nlive) .eq.0) then

                    ! Calculate the covariance matrix
                    covmat = calc_covmat( (/ live_points(settings%h0:settings%h1,:),phantom_points(settings%h0:settings%h1,:stack_size)/), settings%nDims,stack_size+settings%nlive )
                    ! Calculate the cholesky decomposition
                    cholesky = calc_cholesky(covmat,settings%nDims)


                end if


                ! (2) Generate a new set of baby points
                ! Select a seed point for the generator
                first_loop = .true.
                do while (seed_point(settings%l0)<late_likelihood .or. nint(seed_point(settings%point_type)) /= live_type .or. first_loop)
                    seed_point = live_points(:,random_integer(settings%nlive))
                    ! Record the likelihood bound which this seed will generate from
                    seed_point(settings%l1) = late_likelihood
                    first_loop=.false.
                end do


                if(linear_mode) then

                    ! Generate a new set of points within the likelihood bound of the late point
                    baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)

                    baby_made=.true.
                else
                    if(send_start==0) then
                        ! (2) Recieve newly generated baby point from any slave
                        call MPI_RECV(baby_points,settings%nTotal*settings%num_babies,&
                            MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)
                        baby_made=.true.
                    else
                        mpi_status(MPI_SOURCE)=send_start
                        send_start=send_start-1
                        baby_made=.false.
                    end if

                    ! Send a seed point back to that slave
                    call MPI_SEND(seed_point,settings%nTotal,&
                        MPI_DOUBLE_PRECISION,mpi_status(MPI_SOURCE),RUNTAG,mpi_communicator,mpierror)

                    ! Send the information needed
                    call MPI_SEND(cholesky,settings%nDims*settings%nDims,&
                        MPI_DOUBLE_PRECISION,mpi_status(MPI_SOURCE),RUNTAG,mpi_communicator,mpierror)

                end if

                ! The new likelihood is the last point
                baby_likelihood  = baby_points(settings%l0,settings%num_babies)

                if(baby_likelihood>late_likelihood .and. baby_made) then

                    ! (3) Calculate the new evidence (and check to see if we're accurate enough)
                    more_samples_needed = KeetonEvidence(settings,baby_likelihood,late_likelihood,ndead,evidence_vec)

                    ! (4) Update the stack of live points and the posterior array
                    !     This function does multiple things:
                    !     1) Insert baby_points into live_points
                    !     2) Remove points from live_point that have died this round
                    !     3) Add any of these which are at a high enough likelihood to the posterior_array
                    !     4) re-calculate stack_size and nposterior
                    !     5) update the late_likelihood
                    !     6) Update ndead
                    more_samples_needed = more_samples_needed .or. update_stacks(settings,baby_points,live_points,stack_size,phantom_points,posterior_array,nposterior,late_likelihood,ndead,total_likelihood_calls,mpi_communicator)

                    ! (5) Feedback to command line every nlive iterations
                    if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                        mean_likelihood_calls = sum( (/live_points(settings%nlike,:),phantom_points(settings%nlike,:stack_size) /) )/(settings%nlive+0d0)
                        write(stdout_unit,'("ndead     = ", I20                  )') ndead
                        write(stdout_unit,'("stack size= ", I20, "/", I20        )') stack_size, settings%nstack
                        if(settings%calculate_posterior) &
                        write(stdout_unit,'("nposterior= ", I20                  )') nposterior
                        write(stdout_unit,'("efficiency= ", F20.2                )') mean_likelihood_calls
                        write(stdout_unit,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                        write(stdout_unit,'("")')
                    end if

                    ! (6) Update the resume and posterior files every update_resume iterations, or at program termination
                    if (mod(ndead,settings%update_resume) .eq. 0 .or.  more_samples_needed==.false.)  then
                        if(settings%do_clustering) then 

                            ! Compute the similarity matrix
                            similarity_matrix = spread( [( dot_product(live_points(settings%h0:settings%h1,i),live_points(settings%h0:settings%h1,i)),i=1,settings%nlive )], dim=2,ncopies=settings%nlive )
                            similarity_matrix = similarity_matrix + transpose(similarity_matrix) - 2d0 * matmul( transpose(live_points(settings%h0:settings%h1,:)),live_points(settings%h0:settings%h1,:) )

                            clusters = NN_clustering(similarity_matrix,settings%SNN_k)
                        end if

                        if(settings%write_resume) call write_resume_file(settings,live_points,stack_size,phantom_points,evidence_vec,ndead,total_likelihood_calls,nposterior,posterior_array) 
                        if(settings%calculate_posterior) call write_posterior_file(settings,posterior_array,evidence_vec(1),nposterior)  
                        if(settings%write_live) call write_phys_live_points(settings,live_points)

                        open(100, file='chains/circle.dat') 
                        sigma=0.01
                        mu1 = (/5d-1 - 10*sigma,5d-1/)
                        mu2 = (/5d-1 + 10*sigma,5d-1/)
                        radius = sqrt(-2 *(log(2d0) + late_likelihood + settings%nDims*log( sqrt(TwoPi)* sigma)))*sigma
                        do i=1,500
                            write(100,'(2E17.8)') mu1 + random_direction(2)*radius
                            write(100,'(2E17.8)') mu2 + random_direction(2)*radius
                        end do
                        close(100) 
                    end if

                    ! If we've put a limit on the maximum number of iterations, then
                    ! check to see if we've reached this
                    if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

                end if




            end do ! End main loop

            if(.not.linear_mode) then
                ! Kill off the final slaves
                ! If we're done, then clean up by receiving the last piece of
                ! data from each node (and throw it away) and then send a kill signal back to it
                do send_start=1,nprocs-1
                    call MPI_RECV(baby_points,settings%nTotal*settings%num_babies, &
                        MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)
                    call MPI_SEND(seed_point,settings%nTotal, &
                        MPI_DOUBLE_PRECISION,mpi_status(MPI_SOURCE),ENDTAG,mpi_communicator,mpierror)
                end do
            end if


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





        else !(myrank/=root)
            
            do while(.true.)
                ! Listen for a signal from the master
                call MPI_RECV(seed_point,settings%nTotal, &
                    MPI_DOUBLE_PRECISION,root,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)

                ! If we receive a kill signal, then exit the loop
                if(mpi_status(MPI_TAG)==ENDTAG) exit

                call MPI_RECV(cholesky,settings%nDims*settings%nDims, &
                    MPI_DOUBLE_PRECISION,root,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)

                baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)


                ! Send the baby points back
                call MPI_SEND(baby_points,settings%nTotal*settings%num_babies, &
                    MPI_DOUBLE_PRECISION,root,RUNTAG,mpi_communicator,mpierror)

            end do

        end if

    end function NestedSampling

    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePointsP(loglikelihood,priors,settings,mpi_communicator,root) result(live_points)
        use priors_module,    only: prior
        use settings_module,  only: program_settings,live_type,blank_type
        use random_module,   only: random_reals
        use utils_module,    only: logzero
        use calculate_module, only: calculate_point
        use read_write_module, only: write_phys_live_points
        use feedback_module,  only: write_started_generating

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


        integer, intent(in) :: mpi_communicator
        integer, intent(in) :: root

        !> The rank of the processor
        integer :: myrank
        integer :: nprocs
        integer :: active_procs
        integer :: mpierror

        double precision, dimension(settings%nTotal,settings%nlive) :: live_points

        !live_points(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal) :: live_point

        ! Loop variable
        integer i_live

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        integer :: empty_buffer(0)

        integer :: tag

        integer :: nlike

        ! Get the number of MPI procedures
        call MPI_COMM_SIZE(mpi_communicator, nprocs, mpierror)
        ! Get the MPI label of the current processor
        call MPI_COMM_RANK(mpi_communicator, myrank, mpierror)

        ! initialise live points at zero
        live_points = 0d0

        if(myrank==root) then

            call write_started_generating(settings%feedback)

            ! The root node just recieves data from all other processors
            active_procs=nprocs-1
            i_live=0
            nlike=0
            do while(active_procs>0) 

                ! Recieve a point from any slave
                call MPI_RECV(live_point,settings%nTotal, &
                    MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)

                ! If its valid, and we need more points, add it to the array
                if(live_point(settings%l0)>logzero .and. i_live<settings%nlive) then
                    i_live=i_live+1
                    live_points(:,i_live) = live_point
                    live_points(settings%nlike,i_live) = live_points(settings%nlike,i_live) + nlike
                    if(settings%write_live) call write_phys_live_points(settings,live_points)
                    nlike=0
                else
                    nlike = nlike+live_point(settings%nlike)
                end if

                ! If we still need more points, send a signal to have another go
                if(i_live<settings%nlive) then
                    tag=RUNTAG
                else
                    tag=ENDTAG
                    active_procs=active_procs-1
                end if

                call MPI_SEND(empty_buffer,0,MPI_INT,mpi_status(MPI_SOURCE),tag,mpi_communicator,mpierror)

            end do


            ! Set the initial trial values of the chords as the diagonal of the hypercube
            live_points(settings%last_chord,:) = sqrt(settings%nDims+0d0)


            ! Set the likelihood contours to logzero for now
            live_points(settings%l1,:) = logzero

            ! These are all true live points
            live_points(settings%point_type,:settings%nlive) = live_type
            live_points(settings%point_type,settings%nlive+1:) = blank_type


        else

            generating_loop: do while(.true.)

                ! Zero the likelihood calls 
                live_point(settings%nlike) = 0

                ! Generate a random hypercube coordinate
                live_point(settings%h0:settings%h1) = random_reals(settings%nDims)

                ! Compute physical coordinates, likelihoods and derived parameters
                call calculate_point( loglikelihood, priors, live_point, settings )

                ! Send it to the root node
                call MPI_SEND(live_point,settings%nTotal, &
                    MPI_DOUBLE_PRECISION,root,0,mpi_communicator,mpierror)

                ! Recieve signal as to whether we should keep generating
                call MPI_RECV(empty_buffer,0,MPI_INT,root,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)

                if(mpi_status(MPI_TAG) == ENDTAG ) exit generating_loop

            end do generating_loop

        end if





    end function GenerateLivePointsP


    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePointsL(loglikelihood,priors,settings) result(live_points)
        use priors_module,    only: prior
        use settings_module,  only: program_settings,live_type,blank_type
        use random_module,    only: random_reals
        use utils_module,     only: logzero
        use calculate_module, only: calculate_point
        use feedback_module,  only: write_started_generating

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

        !live_points(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal,settings%nlive) :: live_points

        ! Loop variable
        integer i_live

        call write_started_generating(settings%feedback)

        ! initialise live points at zero
        live_points = 0d0

        do i_live=1,settings%nlive

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
        live_points(settings%point_type,:settings%nlive) = live_type
        live_points(settings%point_type,settings%nlive+1:) = blank_type


    end function GenerateLivePointsL



    function update_stacks(settings,baby_points,live_points,stack_size,phantom_points,posterior_array,nposterior,late_likelihood,ndead,total_likelihood_calls,mpi_communicator) result(more_samples_needed)
        use settings_module,   only: program_settings,live_type
        use random_module, only: random_real
        use utils_module, only: stdout_unit
        implicit none
        type(program_settings), intent(in)                                                                           :: settings
        double precision,       intent(in),    dimension(settings%nTotal,settings%num_babies)                        :: baby_points
        double precision,       intent(inout), dimension(settings%nTotal,settings%nlive)                             :: live_points
        integer,                intent(inout)                                                                        :: stack_size
        double precision,       intent(inout), dimension(settings%nTotal,settings%nstack)                            :: phantom_points
        double precision,       intent(inout), dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior) :: posterior_array
        integer,                intent(inout)                                                                        :: nposterior
        double precision,       intent(inout)                                                                        :: late_likelihood
        integer,                intent(inout)                                                                        :: ndead
        integer,                intent(inout)                                                                        :: total_likelihood_calls
        integer,                intent(in)                                                                           :: mpi_communicator

        logical :: more_samples_needed

        integer :: late_index(1)

        integer :: i_live

        double precision :: late_logweight

        double precision, dimension(settings%nDims+settings%nDerived+2) :: posterior_point

        double precision :: lognmax_posterior
        double precision :: max_logweight
        integer :: errorcode
        integer :: mpierror


        late_logweight = (ndead-1)*log(settings%nlive+0d0) - ndead*log(settings%nlive+1d0)                

        ! Start by finding the original lowest likelihood live point (about to be deleted)
        late_index = minloc(live_points(settings%l0,:))

        ! Update the late likelihood
        late_likelihood = live_points(settings%l0,late_index(1))

        if(settings%calculate_posterior) then
            ! Add the discarded point to the posterior array
            posterior_point(1)  = live_points(settings%l0,late_index(1)) + late_logweight
            posterior_point(2)  = live_points(settings%l0,late_index(1))
            posterior_point(2+1:2+settings%nDims) = live_points(settings%p0:settings%p1,late_index(1))
            posterior_point(2+settings%nDims+1:2+settings%nDims+settings%nDerived) = live_points(settings%d0:settings%d1,late_index(1))

            nposterior=nposterior+1
            if(nposterior>settings%nmax_posterior) then
                write(stdout_unit,'(" Too many posterior points. Consider increasing nmax_posterior ")')
                call MPI_ABORT(mpi_communicator,errorcode,mpierror)
            end if
            posterior_array(:,nposterior) = posterior_point
        end if

        ! Replace the late point with the new baby point
        live_points(:,late_index(1)) = baby_points(:,settings%num_babies)


        ! Add the remaining baby points to the end of the array, and update the stack size
        stack_size=stack_size+settings%num_babies-1
        if(stack_size>settings%nstack) then
            write(stdout_unit,'(" Stack size too small, increase nstack ")')
            call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
        end if
        phantom_points(:,stack_size-settings%num_babies+2:stack_size) = baby_points(:,:settings%num_babies-1)

        ! Now run through the stack and strip out any points that are less
        ! than the new late_likelihood, replacing them with points drawn from
        ! the end 

        i_live=1
        do while(i_live<=stack_size)
            if( phantom_points(settings%l0,i_live) < late_likelihood ) then

                if(settings%calculate_posterior .and. random_real() < settings%thin_posterior) then
                    ! Add the discarded point to the posterior array
                    posterior_point(1)  = phantom_points(settings%l0,i_live) + late_logweight
                    posterior_point(2)  = phantom_points(settings%l0,i_live)
                    posterior_point(2+1:2+settings%nDims) = phantom_points(settings%p0:settings%p1,i_live)
                    posterior_point(2+settings%nDims+1:2+settings%nDims+settings%nDerived) = phantom_points(settings%d0:settings%d1,i_live)

                    nposterior=nposterior+1
                    if(nposterior>settings%nmax_posterior) then
                        write(stdout_unit,'(" Too many posterior points. Consider increasing nmax_posterior ")')
                        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
                    end if
                    posterior_array(:,nposterior) = posterior_point
                end if

                ! Update the total likelihood calls
                total_likelihood_calls = total_likelihood_calls + phantom_points(settings%nlike,i_live)

                ! Overwrite the discarded point with a point from the end...
                phantom_points(:,i_live) = phantom_points(:,stack_size)
                ! ...and reduce the stack size
                stack_size=stack_size-1
            else
                i_live=i_live+1
            end if
        end do

        if(settings%calculate_posterior) then

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
        end if

        live_points(settings%last_chord,:) = live_points(settings%last_chord,:)/  (1d0+1d0/(settings%nDims*settings%nlive) )


        ! Find the new late likelihood
        late_likelihood = minval(live_points(settings%l0,:))

        ! Increment the number of dead points
        ndead=ndead+1

        ! Test to see if we need more samples
        !more_samples_needed = any(posterior_points(1,:) - max_logweight+lognmax_posterior > 0 )
        more_samples_needed = .false.

    end function update_stacks





end module nested_sampling_module
