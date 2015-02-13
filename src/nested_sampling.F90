module nested_sampling_module

#ifdef MPI
    use mpi_module, only: get_rank,get_nprocs,get_root,catch_babies,throw_babies,throw_seed,catch_seed
#endif

    implicit none


    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood,priors,settings,mpi_communicator) result(output_info)
        use priors_module,     only: prior,prior_log_volume
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp,calc_similarity_matrix,stdout_unit,swap_integers,logzero
        use read_write_module
        use feedback_module
        use run_time_module,   only: run_time_info,replace_point,calculate_logZ_estimate
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: NN_clustering
        use generate_module,   only: GenerateSeed,GenerateLivePoints

        implicit none

        ! Program Inputs
        ! ==============      

        !> The loglikelihood function
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
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

        ! Temporary variables
        ! -------------------
        ! Seed point to generate babies from
        double precision, dimension(settings%nTotal) :: seed_point
        ! Cholesky matrix to generate babies from
        double precision, dimension(settings%nDims,settings%nDims)   :: cholesky
        ! Loglikelihood bound to generate babies from
        double precision :: logL
        ! New-born baby points, created by slice sampling routine
        double precision, dimension(settings%nTotal,settings%num_babies) :: baby_points
        integer :: cluster_id ! Cluster identifier
        integer :: nlike      ! Temporary storage for number of likelihood calls

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

#else 
        ! non-MPI initialisation
        root=0      ! Define root node
        myrank=root ! set the only node's rank to be the root
        nprocs=1    ! there's only one process

#endif



        !-------------------------------------------------------!
        if(myrank==root) call write_opening_statement(settings) !
        !-------------------------------------------------------!



        ! Check if we actually want to resume
        if ( settings%read_resume .and. resume_file_exists(settings) ) then 

            ! Read the resume file on root
            if(myrank==root) call read_resume_file(settings,RTI) 

        else 

            ! Intialise the run by setting all of the relevant run time info, and generating live points
            call GenerateLivePoints(loglikelihood,priors,settings,RTI,mpi_communicator,nprocs,myrank,root)

            ! Write a resume file (as the generation of live points can be intensive)
            if(myrank==root.and.settings%write_resume) then
                call delete_files(settings)     ! Delete any existing files
                if(settings%write_resume) call write_resume_file(settings,RTI) ! Write a new resume file
            end if

        end if 




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
                    baby_points = SliceSampling(loglikelihood,priors,settings,logL,seed_point,cholesky,nlike)
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


                if( replace_point(settings,RTI,baby_points,cluster_id) ) then

                    need_more_samples = more_samples_needed(settings,RTI) 

                    ! Update the resume files every settings%update_resume iterations,
                    ! or at the end of the run
                    if( mod(RTI%ndead,settings%update_resume)==0 .or. need_more_samples ) then
                        ! Write the resume file if desired
                        if(settings%write_resume)        call write_resume_file(settings,RTI)
                        if(settings%calculate_posterior) call write_unnormalised_posterior_file(settings,RTI)  
                        if(settings%write_live)          call write_phys_live_points(settings,RTI)
                        if(settings%feedback>=1)         call write_intermediate_results(settings,RTI)
                        if(settings%write_stats)         call write_stats_file(settings,RTI)
                    end if

                    !call do_clustering



                end if

            end do ! End of main loop body




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
#endif


            ! Create the output array
            ! (1) log evidence
            ! (2) Error in the log evidence
            ! (3) Number of dead points
            ! (4) Number of likelihood calls
            ! (5) log(evidence * prior volume)
            call calculate_logZ_estimate(RTI,output_info(1),output_info(2))
            output_info(3) = RTI%ndead
            output_info(4) = RTI%nlike
            output_info(5) = output_info(1)+prior_log_volume(priors)

            ! ------------------------------------------------------------ !
            call write_final_results(output_info,settings%feedback,priors) !
            ! ------------------------------------------------------------ !





#ifdef MPI
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
                baby_points = SliceSampling(loglikelihood,priors,settings,logL,seed_point,cholesky,nlike)

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




    function identify_cluster(settings,RTI,point) result(cluster)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info
        use utils_module,      only: loginf,distance2
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: RTI

        double precision, dimension(settings%nTotal),intent(in)   :: point

        integer :: cluster

        integer :: i_cluster
        integer :: i_live

        double precision :: temp_distance2
        double precision :: closest_distance2

        closest_distance2=loginf

        ! Find the cluster this point is nearest to
        do i_cluster=1,RTI%ncluster
            do i_live=1,RTI%nlive(i_cluster)
                temp_distance2 = distance2(point(settings%h0:settings%h1),RTI%live(settings%h0:settings%h1,i_live,i_cluster) )
                if(temp_distance2 < closest_distance2) then
                    cluster = i_cluster
                    closest_distance2 = temp_distance2
                end if
            end do
        end do

    end function identify_cluster






    !> Calculate a posterior point from a live/phantom point, suitable for
    !! adding to a .txt file
    function calc_posterior_point(settings,point,logweight,evidence) result(posterior_point)
        use settings_module,   only: program_settings
        use utils_module, only: logincexp
        implicit none

        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nTotal),intent(in) :: point
        double precision,intent(in) :: logweight
        double precision,intent(in) :: evidence
        double precision, dimension(settings%nposterior) :: posterior_point


        ! Un-normalised weighting (needs to be unnormalised since the evidence is only correct at the end)
        posterior_point(settings%pos_w)  = point(settings%l0) + logweight
        ! un-normalise cumulative weighting
        posterior_point(settings%pos_Z)  = evidence
        ! Likelihood
        posterior_point(settings%pos_l)  = point(settings%l0)
        ! Physical parameters
        posterior_point(settings%pos_p0:settings%pos_p1) = point(settings%p0:settings%p1)
        ! Derived parameters
        posterior_point(settings%pos_d0:settings%pos_d1) = point(settings%d0:settings%d1)

    end function calc_posterior_point

















end module nested_sampling_module
