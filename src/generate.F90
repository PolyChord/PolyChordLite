!> This module contains 'generating tools', namely:
!!
!! * GenerateSeed
!! ** Generates a seed point
!! * GenerateLivePointsP
!! * GenerateLivePointsL
module generate_module


    implicit none

    contains


    !> This function generates a seed point for slice sampling.
    !!
    !! It uses the existing set of live points and the details of the clustering
    !! in order to select a seed point to generate from in proportion to the
    !! estimates of the prior volume of each cluster.
    function GenerateSeed(settings,RTI,seed_cluster) result(seed_point)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info
        use random_module,     only: random_integer,random_integer_P,bernoulli_trial
        use utils_module,      only: logsumexp
        implicit none

        !> Program settings
        type(program_settings), intent(in) :: settings
        !> The evidence storage
        type(run_time_info), intent(in) :: RTI
        !> The cluster number chosen
        integer,intent(out) :: seed_cluster


        ! The seed point to be produced
        double precision, dimension(settings%nTotal) :: seed_point

        integer :: seed_choice

        double precision, dimension(RTI%ncluster) :: probs


        ! 0) Calculate an array proportional to the volumes
        probs = RTI%logXp(:RTI%ncluster)  ! prob_p = log( X_p )
        probs = probs - logsumexp(probs)  ! prob_p = log( X_p/(sum_q X_q) )
        probs = exp(probs)                ! prob_p = X_p/(sum_q X_q)

        ! 1) Pick cluster in proportion to the set of volume estimates of the active clusters
        seed_cluster = random_integer_P(probs)

        ! 2) Pick whether to draw from phantom or live points
        if(bernoulli_trial(RTI%nlive(seed_cluster)+0d0,RTI%nphantom(seed_cluster)+0d0)) then

            ! 3a) Pick a random integer in between 1 and the number of live points in the cluster 'p'
            seed_choice = random_integer(RTI%nlive(seed_cluster))

            ! 4a) Select the live point at index 'seed_choice' in cluster 'p' for the seed point
            seed_point = RTI%live(:,seed_choice,seed_cluster)
        else

            ! 3b) Pick a random integer in between 1 and the number of phantom points in cluster 'p'
            seed_choice = random_integer(RTI%nphantom(seed_cluster))

            ! 4b) Select the phantom point at index 'seed_choice' in cluster 'p' for the seed point
            seed_point = RTI%phantom(:,seed_choice,seed_cluster)
        end if

    end function GenerateSeed




    !> Generate an initial set of live points distributed uniformly in the unit hypercube in parallel
    subroutine GenerateLivePoints(loglikelihood,priors,settings,RTI,mpi_communicator,nprocs,myrank,root)
        use priors_module,    only: prior
        use settings_module,  only: program_settings
        use random_module,   only: random_reals
        use utils_module,    only: logzero,write_phys_unit,DB_FMT,fmt_len,minpos
        use calculate_module, only: calculate_point
        use read_write_module, only: phys_live_file
        use feedback_module,  only: write_started_generating,write_finished_generating
        use run_time_module,   only: run_time_info,initialise_run_time_info
        use abort_module
#ifdef MPI
        use mpi_module, only: throw_point,catch_point,more_points_needed,sum_nlike,request_point,no_more_points
#endif

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

        ! The run time info (very important, see src/run_time_info.f90)
        type(run_time_info) :: RTI

        integer, intent(in) :: mpi_communicator !> MPI handle
        integer, intent(in) :: nprocs           !> The number of processes
        integer, intent(in) :: myrank           !> The rank of the processor
        integer, intent(in) :: root             !> The root process
#ifdef MPI
        integer             :: active_slaves    !  Number of currently working slaves
        integer             :: slave_id         !  Slave identifier to signal who to throw back to
#endif

        double precision, dimension(settings%nTotal) :: live_point ! Temporary live point array

        integer i_live ! Loop variable

        character(len=fmt_len) :: fmt_dbl ! writing format variable

        integer :: nlike ! number of likelihood calls


        ! Initialise number of likelihood calls to zero here
        nlike = 0


        if(myrank==root) then
            ! ---------------------------------------------- !
            call write_started_generating(settings%feedback)
            ! ---------------------------------------------- !

            ! Initialise the format
            write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT

            ! Open the live points file to sequentially add live points
            open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')

            ! Allocate the run time arrays, and set the default values for the variables
            call initialise_run_time_info(settings,RTI)

        end if



        if(nprocs==1) then
            !===================== LINEAR MODE =========================

            i_live=0 ! Initially no live points have been generated
            do while(i_live<settings%nlive)

                ! Generate a random coordinate
                live_point(settings%h0:settings%h1) = random_reals(settings%nDims)

                ! Compute physical coordinates, likelihoods and derived parameters
                call calculate_point( loglikelihood, priors, live_point, settings, nlike)

                ! If its valid, and we need more points, add it to the array
                if(live_point(settings%l0)>logzero) then

                    i_live=i_live+1                      ! Increase the live point counter
                    RTI%live(:,i_live,1) = live_point    ! Add the new live point to the array
                    RTI%nlive(1) = i_live                ! note down the number of live point in RTI

                    if(settings%write_live) then
                        ! Write the live points to the live_points file
                        write(write_phys_unit,fmt_dbl) live_point(settings%p0:settings%d1), live_point(settings%l0)
                        ! flush the unit to force write
                        call flush(write_phys_unit)
                    end if

                end if

            end do


        else if(nprocs>1) then
#ifdef MPI
            !===================== PARALLEL MODE =======================

            if(myrank==root) then
                ! The root node just recieves data from all other processors


                active_slaves=nprocs-1 ! Set the number of active processors to the number of slaves
                i_live=0               ! No live points initially

                do while(active_slaves>0) 

                    ! Recieve a point from any slave
                    slave_id = catch_point(live_point,mpi_communicator)

                    ! If its valid, and we need more points, add it to the array
                    if(live_point(settings%l0)>logzero .and. i_live<settings%nlive) then

                        i_live=i_live+1                      ! Increase the live point counter
                        RTI%live(:,i_live,1) = live_point    ! Add the new live point to the array
                        RTI%nlive(1) = i_live                ! note down the number of live point in RTI

                        if(settings%write_live) then
                            ! Write the live points to the live_points file
                            write(write_phys_unit,fmt_dbl) live_point(settings%p0:settings%d1), live_point(settings%l0)
                            ! flush the unit to force write
                            call flush(write_phys_unit)
                        end if

                    end if


                    if(i_live<settings%nlive) then
                        ! If we still need more points, send a signal to have another go
                        call request_point(mpi_communicator,slave_id)
                    else
                        ! Otherwise, send a signal to stop
                        call no_more_points(mpi_communicator,slave_id)

                        active_slaves=active_slaves-1 ! decrease the active slave counter
                    end if

                end do




            else

                ! The slaves simply generate and send points until they're told to stop by
                ! the master

                do while(.true.)

                    ! Generate a random hypercube coordinate
                    live_point(settings%h0:settings%h1) = random_reals(settings%nDims)

                    ! Compute physical coordinates, likelihoods and derived parameters
                    call calculate_point( loglikelihood, priors, live_point, settings,nlike)

                    ! Send it to the root node
                    call throw_point(live_point,mpi_communicator,root)

                    ! If we've recieved a kill signal, then exit this loop
                    if(.not. more_points_needed(mpi_communicator,root)) exit

                end do

            end if




#else
            ! If we don't have MPI configured, we can't generate in parallel
            call halt_program('generate error: cannot have nprocs>1 without MPI')
#endif
        else !(nprocs<0)
            call halt_program('generate error: nprocs<0')
        end if !(nprocs case)

#ifdef MPI
        nlike = sum_nlike(nlike,mpi_communicator) ! Gather the likelihood calls onto one node
#endif


        if(myrank==root) then

            ! Pass over the number of likelihood calls
            RTI%nlike = nlike

            ! Set the local and global loglikelihood bounds
            RTI%p     = 1                                 ! Only one cluster at this point
            RTI%i     = minpos(RTI%live(settings%l0,:,1)) ! Find the position of the minimum loglikelihood
            RTI%logL  = RTI%live(settings%l0,RTI%i,1)     ! Store the value of the minimum loglikelihood 
            RTI%logLp = RTI%logL                          ! global = local for one cluster

            ! Close the file
            close(write_phys_unit)

            ! ----------------------------------------------- !
            call write_finished_generating(settings%feedback)  
            ! ----------------------------------------------- !
        end if






    end subroutine GenerateLivePoints


end module generate_module
