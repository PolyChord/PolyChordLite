module nested_sampling_module
    use utils_module, only: dp

#ifdef MPI
    use mpi_module, only: get_mpi_information,mpi_bundle,is_root,linear_mode,catch_babies,throw_babies,throw_seed,catch_seed,broadcast_integers,mpi_synchronise
#else
    use mpi_module, only: get_mpi_information,mpi_bundle,is_root,linear_mode
#endif

    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood, prior, dumper, cluster, settings, mpi_communicator) result(output_info)
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp,calc_similarity_matrix,swap_integers,cyc,time
        use read_write_module
        use feedback_module
        use run_time_module,   only: run_time_info,replace_point,calculate_logZ_estimate,calculate_covmats,delete_cluster,update_posteriors,delete_outermost_point
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: do_clustering
        use generate_module,   only: GenerateSeed,GenerateLivePoints,GenerateLivePointsFromSeed
        use maximise_module,   only: maximise
#ifdef MPI
        use utils_module, only: normal_fb,stdout_unit
#else
        use utils_module, only: stdout_unit
#endif

        implicit none

        ! Program Inputs
        ! ==============      

        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in), dimension(:)  :: theta
                real(dp), intent(out), dimension(:) :: phi
                real(dp) :: loglikelihood
            end function
        end interface
        interface
            function prior(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function
        end interface
        interface
            subroutine dumper(live, dead, logweights, logZ, logZerr)
                import :: dp
                real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
                real(dp), intent(in) :: logZ, logZerr
            end subroutine dumper
        end interface
        interface
            function cluster(points) result(cluster_list)
                import :: dp
                real(dp), intent(in), dimension(:,:) :: points
                integer, dimension(size(points,2)) :: cluster_list
            end function
        end interface

        !> Program settings
        type(program_settings), intent(in) :: settings

        !> MPI handle
        integer, intent(in) :: mpi_communicator
        
        ! Program Outputs
        ! ===============      
        ! 1) log(evidence)
        ! 2) error(log(evidence))
        ! 3) ndead
        ! 4) number of likelihood calls
        ! 5) log(evidence) + log(prior volume)
        real(dp), dimension(4) :: output_info



        ! Local variables
        ! ===============

        ! The run time info
        ! ----------------
        ! very important, see src/run_time_info.f90
        type(run_time_info) :: RTI

        ! The number of repeats within each parameter speed to do
        integer, dimension(size(settings%grade_dims)) :: num_repeats
        integer, dimension(size(settings%grade_dims)) :: nlike      ! Temporary storage for number of likelihood calls
        integer, dimension(size(settings%grade_dims)) :: nlikesum   ! rolling average of number of loglikelihood calls


        ! Temporary variables
        ! -------------------
        ! Seed point to generate babies from
        real(dp), dimension(settings%nTotal) :: seed_point
        ! Cholesky matrix to generate babies from
        real(dp), dimension(settings%nDims,settings%nDims)   :: cholesky
        ! Loglikelihood bound to generate babies from
        real(dp) :: logL
        ! New-born baby points, created by slice sampling routine
        real(dp), allocatable, dimension(:,:) :: baby_points

        integer :: cluster_id ! Cluster identifier

        integer :: failures, nfail


        ! Logical Switches
        ! ----------------
        logical :: temp_logical
        logical :: update

        ! MPI process variable
        ! --------------------
        type(mpi_bundle) :: mpi_information


#ifdef MPI
        ! MPI specific variables
        ! ----------------------
        integer                            :: i_worker       ! Worker iterator
        integer                            :: worker_id      ! Worker identifier
        integer, dimension(:), allocatable :: worker_cluster ! The cluster the worker is currently working on
        real(dp) :: time0,time1,slice_time,wait_time

        ! Worker switch
        ! -------------
        ! This prevents workers delivering points to incorrect clusters after clustering
        ! has reorganised the cluster indices
        integer ::  worker_epoch
        integer ::  administrator_epoch

        ! Nursary for storing babies in synchronous parallel mode
        real(dp), allocatable, dimension(:,:,:) :: nursary
        integer :: i_nursary
        integer, allocatable, dimension(:) :: worker_epochs
        integer, allocatable, dimension(:,:) :: nlikes
#endif


        ! A) Initialisation
        !    ==============
        ! MPI initialisation
        mpi_information = get_mpi_information(mpi_communicator)

#ifdef MPI
        allocate(worker_cluster(mpi_information%nprocs-1)) ! Allocate the worker arrays
        worker_cluster = 1                          ! initialise with 1

        ! worker switch
        worker_epoch=0
        administrator_epoch=0
        i_nursary=0
#endif

        ! Rolling loglikelihood calculation
        nlikesum=0

        ! Number of failed spawns
        if (settings%nfail <= 0) then
            nfail = settings%nlive
        else
            nfail = settings%nfail
        end if
        failures = 0


        !-------------------------------------------------------!
        if(is_root(mpi_information)) call check_directories(settings)
        if(is_root(mpi_information)) call write_opening_statement(settings) !
        !-------------------------------------------------------!



        ! Check if we actually want to resume
        if ( settings%read_resume .and. resume_file_exists(settings) ) then 

            ! Read the resume file on root
            if(is_root(mpi_information)) then
                call read_resume_file(settings,RTI) 
                ! -------------------------------------------- !
                call write_resuming(settings%feedback)
                ! -------------------------------------------- !
            end if


        else 
            
            ! Delete any existing files if we're going to be producing our own new ones
            if(is_root(mpi_information)) then
                if(settings%write_resume) call delete_files(settings) 
            end if

            ! Intialise the run by setting all of the relevant run time info, and generating live points
            if (settings%generate_from_seed) then
                call GenerateLivePointsFromSeed(loglikelihood,prior,settings,RTI,mpi_information)
            else
                call GenerateLivePoints(loglikelihood,prior,settings,RTI,mpi_information)
            end if

            if(is_root(mpi_information)) then
                if(settings%write_prior) call write_prior_file(settings,RTI) 
            end if

            if (is_root(mpi_information)) then
                do while(RTI%nlive(1) > settings%nlive )
                    call delete_outermost_point(settings,RTI)
                end do
            end if

            ! Write a resume file (as the generation of live points can be intensive)
            if(is_root(mpi_information)) then
                if(settings%write_resume) then
                    call write_resume_file(settings,RTI) 
                    call rename_files(settings,RTI)
                end if
            end if


        end if 

        if(is_root(mpi_information)) then
            num_repeats = RTI%num_repeats
            call write_num_repeats(num_repeats,settings%feedback)
        end if
#ifdef MPI
        call broadcast_integers(num_repeats,mpi_information)
        allocate(nursary(settings%nTotal,sum(num_repeats), mpi_information%nprocs-1))
        allocate(worker_epochs(mpi_information%nprocs-1), nlikes(size(settings%grade_dims),mpi_information%nprocs-1))
#endif
        allocate(baby_points(settings%nTotal,sum(num_repeats)))



        ! B) Main loop body
        !    ==============

        if(is_root(mpi_information)) then

            ! -------------------------------------------- !
            call write_started_sampling(settings%feedback) !
            ! -------------------------------------------- !

            do while ( more_samples_needed(settings,RTI) .and. failures <= nfail )


                ! 1) Generate a new live point
                !    -------------------------
                ! Generate a seed point --- update this
                seed_point = GenerateSeed(settings,RTI,cluster_id)

                ! Choose the cholesky decomposition for the cluster
                cholesky = RTI%cholesky(:,:,cluster_id)

                ! Get the loglikelihood contour we're generating from
                logL = RTI%logLp(cluster_id)


                if(linear_mode(mpi_information)) then
                    ! Linear Mode
                    ! -----------

                    ! Generate a new set of points within the likelihood bound of the late point
                    baby_points = SliceSampling(loglikelihood,prior,settings,logL,seed_point,cholesky,nlike,num_repeats)
                    baby_points(settings%b0,:) = logL ! Note the moment it is born at
#ifdef MPI
                else if(settings%synchronous) then
                    ! Parallel synchronous mode
                    ! -------------------------

                    if(i_nursary == 0) then
                        do worker_id=1,mpi_information%nprocs-1
                            seed_point = GenerateSeed(settings,RTI,cluster_id)
                            cholesky = RTI%cholesky(:,:,cluster_id)
                            logL = RTI%logLp(cluster_id)
                            call throw_seed(seed_point,cholesky,logL,mpi_information,worker_id,administrator_epoch,.true.)
                            worker_cluster(worker_id) = cluster_id
                        end do
                        do i_worker=1,mpi_information%nprocs-1
                            i_nursary = catch_babies(baby_points,nlike,worker_epoch,mpi_information)
                            nursary(:,:,i_nursary) = baby_points
                            worker_epochs(i_nursary) = worker_epoch
                            nlikes(:,i_nursary) = nlike
                        end do
                        i_nursary = mpi_information%nprocs-1
                    end if
                    baby_points = nursary(:,:,i_nursary)
                    cluster_id = worker_cluster(i_nursary)
                    worker_epoch = worker_epochs(i_nursary)
                    nlike = nlikes(:,i_nursary)
                    i_nursary = i_nursary-1

                else
                    ! Parallel mode
                    ! -------------

                    ! Recieve any new baby points from any worker currently sending
                    worker_id = catch_babies(baby_points,nlike,worker_epoch,mpi_information)

                    ! and throw seeding information back to worker (true => keep going)
                    call throw_seed(seed_point,cholesky,logL,mpi_information,worker_id,administrator_epoch,.true.)

                    ! set cluster_id to be the cluster identity of the babies just recieved 
                    ! (saved in worker_cluster from the last send) and set worker_cluster to 
                    ! be the bound just sent off.
                    call swap_integers(cluster_id,worker_cluster(worker_id))

#endif
                end if

                ! Add the likelihood calls to our counter
                RTI%nlike = RTI%nlike + nlike
                nlikesum  = nlikesum  + nlike


                ! See if this point is suitable to be added to the arrays
#ifdef MPI
                if( linear_mode(mpi_information) .or. administrator_epoch==worker_epoch ) then
#endif
                    if(replace_point(settings,RTI,baby_points,cluster_id)) then
                        failures = 0
                    else
                        failures = failures + 1
                    end if

                    update = logsumexp(RTI%logXp) <= RTI%logX_last_update + log(settings%compression_factor) 
                    if (update) then
                        RTI%logX_last_update = logsumexp(RTI%logXp)

                        ! Update the posterior array
                        call update_posteriors(settings,RTI)  

                        ! Update the resume files
                        if(settings%write_resume)                  call write_resume_file(settings,RTI)
                        if(settings%write_live)                    call write_phys_live_points(settings,RTI)
                        if(settings%write_dead)                    call write_dead_points(settings,RTI)   
                        if(settings%write_stats)                   call write_stats_file(settings,RTI,nlikesum)
                        if(settings%equals.or.settings%posteriors) call write_posterior_file(settings,RTI)   
                        call rename_files(settings,RTI)
                        call dump(dumper,settings,RTI)

                    end if

                    if(delete_cluster(settings,RTI)) then
#ifdef MPI
                        administrator_epoch = administrator_epoch+1
#endif
                    end if! Delete any clusters as necessary
                    if (RTI%ncluster == 0) exit

                    if(update) then
                        !--------------------------------------------!
                        call write_intermediate_results(settings,RTI,nlikesum)
                        nlikesum=0
                        !--------------------------------------------!
                        if(settings%do_clustering) then

                            ! If we want to cluster on sub dimensions, then do this first
                            if(allocated(settings%sub_clustering_dimensions)) then
                                if( do_clustering(settings,RTI,cluster,settings%sub_clustering_dimensions) )  then
#ifdef MPI
                                    administrator_epoch = administrator_epoch+1
#endif
                                end if
                            end if

                            if( do_clustering(settings,RTI,cluster) )  then
#ifdef MPI
                                administrator_epoch = administrator_epoch+1
#endif
                            end if
                        end if
                        call calculate_covmats(settings,RTI)
                    end if
#ifdef MPI
                end if
#endif

            end do ! End of main loop body

            if(settings%write_resume)                  call write_resume_file(settings,RTI)

            ! Do maximisation if required
            if(is_root(mpi_information) .and. settings%maximise) call maximise(loglikelihood,prior,settings,RTI)

            do while(RTI%ncluster > 0)
                call delete_outermost_point(settings,RTI)
                temp_logical = delete_cluster(settings,RTI)
            end do

            call update_posteriors(settings,RTI) 
            if(settings%write_live)                    call write_phys_live_points(settings,RTI)
            if(settings%equals.or.settings%posteriors) call write_posterior_file(settings,RTI)   
            if(settings%write_dead)                    call write_dead_points(settings,RTI)   
            if(settings%write_stats)                   call write_stats_file(settings,RTI,nlikesum)
            call rename_files(settings,RTI)
            call dump(dumper,settings,RTI)

            ! Create the output array
            ! (1) log evidence
            ! (2) variance in the log evidence
            ! (3) Number of dead points
            ! (4) Number of slow likelihood calls
            ! (5) log(evidence * prior volume)
            call calculate_logZ_estimate(RTI,output_info(1),output_info(2))
            output_info(3) = RTI%ndead
            output_info(4) = RTI%nlike(1)

            ! ------------------------------------------------------------ !
            call write_final_results(output_info,settings%feedback)        !
            ! ------------------------------------------------------------ !
            if (failures > nfail) then
                write(stdout_unit,'("Warning, unable to proceed after ",I6,": failed spawn events")') failures
            end if







            ! C) Clean up
            !    ========

#ifdef MPI
            ! MPI cleanup
            ! -----------
            ! Kill off the final workers.
            ! If we're done, then clean up by receiving the last piece of
            ! data from each node (and throw it away) and then send a kill signal back to it
            if (settings%synchronous) then
                do worker_id=mpi_information%nprocs-1,1,-1
                    call throw_seed(seed_point,cholesky,logL,mpi_information,worker_id,administrator_epoch,.false.) 
                end do
            else
                do i_worker=mpi_information%nprocs-1,1,-1

                    ! Recieve baby point from worker worker_id
                    worker_id = catch_babies(baby_points,nlike,worker_epoch,mpi_information)

                    ! Add the likelihood calls to our counter
                    RTI%nlike = RTI%nlike + nlike

                    ! Send kill signal to worker worker_id (note that we no longer care about seed_point, so we'll just use the last one
                    call throw_seed(seed_point,cholesky,logL,mpi_information,worker_id,administrator_epoch,.false.) 
                end do
            end if


        else !(myrank/=root)

            ! These are the worker tasks
            ! --------------------------
            !
            ! This is considerably simpler than that of the administrator.
            ! All workers do is:
            ! 1) recieve a seed point,cholesky decomposition and loglikelihood
            !    contour from the administrator
            ! 2) using the above, generate a new set of baby points
            ! 3) send the baby points and nlike back to the administrator.


            ! On the first loop, send a nonsense set of baby_points
            ! to indicate that we're ready to start receiving

            if (.not. settings%synchronous) then
                baby_points = 0d0                              ! Avoid sending nonsense
                baby_points(settings%l0,:) = settings%logzero  ! zero contour to ensure these are all thrown away
                baby_points(settings%b0,:) = settings%logzero  ! zero birth contour
                nlike = 0                                      ! no likelihood calls in this round
                call throw_babies(baby_points,nlike,worker_epoch,mpi_information)
            end if
            wait_time = 0
            slice_time = 0
            time1 = time()




            ! 1) Listen for a seed point being sent by the administrator
            !    Note that this also tests for a kill signal sent by the administrator
            do while(catch_seed(seed_point,cholesky,logL,worker_epoch,mpi_information))
                time0 = time()
                ! 2) Generate a new set of baby points
                baby_points = SliceSampling(loglikelihood,prior,settings,logL,seed_point,cholesky,nlike,num_repeats)
                baby_points(settings%b0,:) = logL ! Note the moment it is born at


                wait_time = wait_time + time0-time1
                time1 = time()
                slice_time = slice_time + time1-time0


                ! 3) Send the baby points back
                call throw_babies(baby_points,nlike,worker_epoch,mpi_information)

            end do

            if(slice_time<wait_time) then
                if(settings%feedback>=normal_fb) write(stdout_unit,'("Worker",I3,": Inefficient MPI parallelisation, I spend more time waiting than slicing ", E17.8, ">", E17.8 )') mpi_information%rank, wait_time,slice_time
            else
                if(settings%feedback>=normal_fb) write(stdout_unit,'("Worker",I3,": efficient MPI parallelisation; wait_time/slice_time= ", E17.8 )') mpi_information%rank, wait_time/slice_time 
            end if

#endif
        end if !(myrank==root / myrank/=root) 

#ifdef MPI
        call mpi_synchronise(mpi_information)
#endif




    end function NestedSampling


    !> This function checks whether we've done enough to stop
    function more_samples_needed(settings,RTI)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info,live_logZ
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: RTI

        logical :: more_samples_needed

        ! Set it to default to true
        more_samples_needed = .true.

        ! If we've put a maximum number of iterations on the algorithm, then
        ! we'll stop if we've reached that number. 
        ! If we don't want a maximum number of iterations, then max_ndead should
        ! be set negative
        if(settings%max_ndead==0) then
            more_samples_needed = .false. 
        else if(settings%max_ndead>0 .and. RTI%ndead >= settings%max_ndead) then
            more_samples_needed = .false. 

            ! If the evidence in the live points is less than precision_criterion %
            ! of the total accumulated evidence, then stop.
        else if( settings%precision_criterion > 0 .and. live_logZ(settings,RTI) < log(settings%precision_criterion) + RTI%logZ )  then
            more_samples_needed = .false.
        end if


    end function more_samples_needed


    subroutine dump(dumper, settings, RTI)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info, calculate_logZ_estimate
        use utils_module,      only: logsumexp
        implicit none
        type(run_time_info), intent(in) :: RTI
        type(program_settings), intent(in) :: settings

        interface
            subroutine dumper(live, dead, logweights, logZ, logZerr)
                import :: dp
                real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
                real(dp), intent(in) :: logZ, logZerr
            end subroutine dumper
        end interface

        real(dp), dimension(settings%nDims+settings%nDerived+2,sum(RTI%nlive)) :: live
        real(dp), dimension(settings%nDims+settings%nDerived+2,RTI%ndead) :: dead
        real(dp), dimension(RTI%ndead) :: logweights
        real(dp) :: logZ, varlogZ

        integer i_cluster, nlive,n0,n1

        dead(1:settings%nDims,:) = RTI%dead(settings%p0:settings%p1,:RTI%ndead)
        dead(settings%nDims+1:settings%nDims+settings%nDerived,:) = RTI%dead(settings%d0:settings%d1,:RTI%ndead)
        dead(settings%nDims+settings%nDerived+1,:) = RTI%dead(settings%b0,:RTI%ndead)
        dead(settings%nDims+settings%nDerived+2,:) = RTI%dead(settings%l0,:RTI%ndead)

        logweights = RTI%logweights(:RTI%ndead) + RTI%dead(settings%l0,:RTI%ndead)
        logweights = logweights - logsumexp(logweights)

        do i_cluster=1,RTI%ncluster
            nlive = RTI%nlive(i_cluster)
            n0 = 1+sum(RTI%nlive(:i_cluster-1))
            n1 = sum(RTI%nlive(:i_cluster))
            live(1:settings%nDims,n0:n1) = RTI%live(settings%p0:settings%p1,:nlive,i_cluster)
            live(settings%nDims+1:settings%nDims+settings%nDerived,n0:n1) = RTI%live(settings%d0:settings%d1,:nlive,i_cluster)
            live(settings%nDims+settings%nDerived+1,n0:n1) = RTI%live(settings%b0,:nlive,i_cluster)
            live(settings%nDims+settings%nDerived+2,n0:n1) = RTI%live(settings%l0,:nlive,i_cluster)
        end do
        call calculate_logZ_estimate(RTI,logZ,varlogZ)
        call dumper(live, dead, logweights, logZ, sqrt(varlogZ))


    end subroutine dump





end module nested_sampling_module
