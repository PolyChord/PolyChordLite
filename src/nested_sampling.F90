module nested_sampling_module

#ifdef MPI
    use mpi_module, only: get_rank,get_nprocs,get_root,catch_babies,throw_babies,throw_seed,catch_seed,broadcast_integers
#endif

    implicit none


    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood,priors,settings,mpi_communicator) result(output_info)
        use priors_module,     only: prior,prior_log_volume
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp,calc_similarity_matrix,swap_integers,logzero
        use read_write_module
        use feedback_module
        use run_time_module,   only: run_time_info,replace_point,calculate_logZ_estimate,calculate_covmats,delete_cluster,update_posteriors
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: do_clustering
        use generate_module,   only: GenerateSeed,GenerateLivePoints

        implicit none

        ! Program Inputs
        ! ==============      

        !> The loglikelihood function
        interface
            function loglikelihood(theta,phi)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                double precision :: loglikelihood
            end function
        end interface

        !> Prior information
        type(prior), dimension(:), intent(in) :: priors

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
        double precision, dimension(5) :: output_info



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
        double precision, dimension(settings%nTotal) :: seed_point
        ! Cholesky matrix to generate babies from
        double precision, dimension(settings%nDims,settings%nDims)   :: cholesky
        ! Loglikelihood bound to generate babies from
        double precision :: logL
        ! New-born baby points, created by slice sampling routine
        double precision, allocatable, dimension(:,:) :: baby_points

        integer :: cluster_id ! Cluster identifier


        ! Logical Switches
        ! ----------------
        logical :: need_more_samples

        ! MPI process variables
        ! ---------------------
        integer :: nprocs ! number of MPI processes
        integer :: myrank ! rank of this MPI process
        integer :: root   ! root MPI process (defines 'master')


#ifdef MPI
        ! MPI specific variables
        ! ----------------------
        integer                            :: i_slave       ! Slave iterator
        integer                            :: slave_id      ! Slave identifier
        integer, dimension(:), allocatable :: slave_cluster ! The cluster the slave is currently working on
#endif










        ! A) Initialisation
        !    ==============
#ifdef MPI
        ! MPI initialisation
        nprocs = get_nprocs(mpi_communicator)         ! Get the number of MPI processes
        myrank = get_rank(mpi_communicator)           ! Get the identity of this process
        root   = get_root(myrank,mpi_communicator)    ! Assign the root process as the minimum integer

        allocate(slave_cluster(nprocs-1)) ! Allocate the slave arrays
        slave_cluster = 1                 ! initialise with 1

#else 
        ! non-MPI initialisation
        root=0      ! Define root node
        myrank=root ! set the only node's rank to be the root
        nprocs=1    ! there's only one process

#endif

        ! Rolling loglikelihood calculation
        nlikesum=0


        !-------------------------------------------------------!
        if(myrank==root) call write_opening_statement(settings) !
        !-------------------------------------------------------!



        ! Check if we actually want to resume
        if ( settings%read_resume .and. resume_file_exists(settings) ) then 

            ! Read the resume file on root
            if(myrank==root) call read_resume_file(settings,RTI) 

        else 
            
            ! Delete any existing files if we're going to be producing our own new ones
            if(myrank==root.and.settings%write_resume) call delete_files(settings) 

            ! Intialise the run by setting all of the relevant run time info, and generating live points
            call GenerateLivePoints(loglikelihood,priors,settings,RTI,mpi_communicator,nprocs,myrank,root)

            ! Write a resume file (as the generation of live points can be intensive)
            if(myrank==root.and.settings%write_resume) call write_resume_file(settings,RTI) 

        end if 

        if(myrank==root) then
            num_repeats = RTI%num_repeats
            call write_num_repeats(num_repeats,settings%feedback)
        end if
#ifdef MPI
        call broadcast_integers(num_repeats,mpi_communicator,root)
#endif
        allocate(baby_points(settings%nTotal,sum(num_repeats)))



        ! B) Main loop body
        !    ==============

        if(myrank==root) then

            ! -------------------------------------------- !
            call write_started_sampling(settings%feedback) !
            ! -------------------------------------------- !

            ! Definitely need more samples than this
            need_more_samples = .true.

            do while ( need_more_samples )


                ! 1) Generate a new live point
                !    -------------------------
                ! Generate a seed point --- update this
                seed_point = GenerateSeed(settings,RTI,cluster_id)

                ! Choose the cholesky decomposition for the cluster
                cholesky = RTI%cholesky(:,:,cluster_id)

                ! Get the loglikelihood contour we're generating from
                logL = RTI%logLp(cluster_id)


                if(nprocs==1) then
                    ! Linear Mode
                    ! -----------

                    ! Generate a new set of points within the likelihood bound of the late point
                    baby_points = SliceSampling(loglikelihood,priors,settings,logL,seed_point,cholesky,nlike,num_repeats)
#ifdef MPI
                else
                    ! Parallel mode
                    ! -------------

                    ! Recieve any new baby points from any slave currently sending
                    slave_id = catch_babies(baby_points,nlike,mpi_communicator)

                    ! and throw seeding information back to slave (true => keep going)
                    call throw_seed(seed_point,cholesky,logL,mpi_communicator,slave_id,.true.)

                    ! set cluster_id to be the cluster identity of the babies just recieved 
                    ! (saved in slave_cluster from the last send) and set slave_cluster to 
                    ! be the bound just sent off.
                    call swap_integers(cluster_id,slave_cluster(slave_id))

#endif
                end if

                ! Add the likelihood calls to our counter
                RTI%nlike = RTI%nlike + nlike
                nlikesum  = nlikesum  + nlike


                if( replace_point(settings,RTI,baby_points,cluster_id) ) then

                    need_more_samples = more_samples_needed(settings,RTI) 

                    ! Update the resume files every settings%update_resume iterations,
                    ! or at the end of the run
                    if( mod(RTI%ndead,settings%update_resume)==0 .or. .not. need_more_samples ) then

                        call update_posteriors(settings,RTI) 

                        if(settings%write_resume)        call write_resume_file(settings,RTI)
                        call write_posterior_file(settings,RTI)  
                        if(settings%write_live)          call write_phys_live_points(settings,RTI)
                        if(settings%write_stats)         call write_stats_file(settings,RTI)
                    end if

                    call delete_cluster(settings,RTI) ! Delete any clusters as necessary

                    if( mod(RTI%ndead,settings%nlive)==0 ) then
                        !--------------------------------------------!
                        call write_intermediate_results(settings,RTI,nlikesum,num_repeats)
                        !--------------------------------------------!
                        nlikesum=0
                        if(settings%do_clustering) call do_clustering(settings,RTI)
                        call calculate_covmats(settings,RTI)
                    end if


                end if

            end do ! End of main loop body


            ! Create the output array
            ! (1) log evidence
            ! (2) Error in the log evidence
            ! (3) Number of dead points
            ! (4) Number of slow likelihood calls
            ! (5) log(evidence * prior volume)
            call calculate_logZ_estimate(RTI,output_info(1),output_info(2))
            output_info(3) = RTI%ndead
            output_info(4) = RTI%nlike(1)
            output_info(5) = output_info(1)+prior_log_volume(priors)

            ! ------------------------------------------------------------ !
            call write_final_results(output_info,settings%feedback,priors) !
            ! ------------------------------------------------------------ !







            ! C) Clean up
            !    ========

#ifdef MPI
            ! MPI cleanup
            ! -----------
            ! Kill off the final slaves.
            ! If we're done, then clean up by receiving the last piece of
            ! data from each node (and throw it away) and then send a kill signal back to it
            do i_slave=nprocs-1,1,-1

                ! Recieve baby point from slave slave_id
                slave_id = catch_babies(baby_points,nlike,mpi_communicator)

                ! Add the likelihood calls to our counter
                RTI%nlike = RTI%nlike + nlike

                ! Send kill signal to slave slave_id (note that we no longer care about seed_point, so we'll just use the last one
                call throw_seed(seed_point,cholesky,logL,mpi_communicator,slave_id,.false.) 

            end do


        else !(myrank/=root)

            ! These are the slave tasks
            ! -------------------------
            !
            ! This is considerably simpler than that of the master.
            ! All slaves do is:
            ! 1) recieve a seed point,cholesky decomposition and loglikelihood
            !    contour from the master
            ! 2) using the above, generate a new set of baby points
            ! 3) send the baby points and nlike back to the master.


            ! On the first loop, send a nonsense set of baby_points
            ! to indicate that we're ready to start receiving

            baby_points = 0d0                     ! Avoid sending nonsense
            baby_points(settings%l0,:) = logzero  ! zero contour to ensure these are all thrown away
            nlike = 0                             ! no likelihood calls in this round
            call throw_babies(baby_points,nlike,mpi_communicator,root)



            ! 1) Listen for a seed point being sent by the master
            !    Note that this also tests for a kill signal sent by the master
            do while(catch_seed(seed_point,cholesky,logL,mpi_communicator,root))

                ! 2) Generate a new set of baby points
                baby_points = SliceSampling(loglikelihood,priors,settings,logL,seed_point,cholesky,nlike,num_repeats)

                ! 3) Send the baby points back
                call throw_babies(baby_points,nlike,mpi_communicator,root)

            end do

#endif
        end if !(myrank==root / myrank/=root) 




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
        ! be set negative or 0
        if(settings%max_ndead>0 .and. RTI%ndead >= settings%max_ndead) then
            more_samples_needed = .false. 
            return
        end if

        ! If the evidence in the live points is less than precision_criterion %
        ! of the total accumulated evidence, then stop.
        if( live_logZ(settings,RTI) < log(settings%precision_criterion) + RTI%logZ )  then
            more_samples_needed = .false.
            return
        end if


    end function more_samples_needed







end module nested_sampling_module
