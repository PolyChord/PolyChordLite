!> This module contains 'generating tools', namely:
!!
!! * GenerateSeed
!! * GenerateLivePoints
module generate_module
    use utils_module, only: dp


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
        real(dp), dimension(settings%nTotal) :: seed_point

        integer :: seed_choice

        real(dp), dimension(RTI%ncluster) :: probs

        ! 0) Calculate an array proportional to the volumes
        probs = RTI%logXp                 ! prob_p = log( X_p )
        probs = probs - logsumexp(probs)  ! prob_p = log( X_p/(sum_q X_q) )
        probs = exp(probs)                ! prob_p = X_p/(sum_q X_q)

        ! 1) Pick cluster in proportion to the set of volume estimates of the active clusters
        seed_cluster = random_integer_P(probs)

        ! 3) Pick a random integer in between 1 and the number of live points in the cluster 'p'
        seed_choice = random_integer(RTI%nlive(seed_cluster))

        ! 4) Select the live point at index 'seed_choice' in cluster 'p' for the seed point
        seed_point = RTI%live(:,seed_choice,seed_cluster)

    end function GenerateSeed




    !> Generate an initial set of live points distributed uniformly in the unit hypercube in parallel
    subroutine GenerateLivePoints(loglikelihood,prior,settings,RTI,mpi_information)
        use settings_module,  only: program_settings
        use random_module,   only: random_reals
        use utils_module,    only: write_phys_unit,DB_FMT,fmt_len,minpos,time
        use calculate_module, only: calculate_point
        use read_write_module, only: phys_live_file
        use feedback_module,  only: write_started_generating,write_finished_generating,write_generating_live_points
        use run_time_module,   only: run_time_info,initialise_run_time_info
        use array_module,     only: add_point
        use abort_module
#ifdef MPI
        use mpi_module, only: mpi_bundle,is_root,linear_mode,throw_point,catch_point,more_points_needed,sum_integers,sum_doubles,request_point,no_more_points
#else
        use mpi_module, only: mpi_bundle,is_root,linear_mode
#endif

        implicit none

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


        !> Program settings
        type(program_settings), intent(in) :: settings

        ! The run time info (very important, see src/run_time_info.f90)
        type(run_time_info) :: RTI

        type(mpi_bundle),intent(in) :: mpi_information
#ifdef MPI
        integer             :: active_slaves    !  Number of currently working slaves
        integer             :: slave_id         !  Slave identifier to signal who to throw back to
#endif

        real(dp), dimension(settings%nTotal) :: live_point ! Temporary live point array


        character(len=fmt_len) :: fmt_dbl ! writing format variable

        integer :: nlike ! number of likelihood calls
        integer :: nprior

        real(dp) :: time0,time1,total_time
        real(dp),dimension(size(settings%grade_dims)) :: speed


        ! Initialise number of likelihood calls to zero here
        nlike = 0


        if(is_root(mpi_information)) then
            ! ---------------------------------------------- !
            call write_started_generating(settings%feedback)
            ! ---------------------------------------------- !

            ! Initialise the format
            write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT

            ! Open the live points file to sequentially add live points
            if(settings%write_live) open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')

            ! Allocate the run time arrays, and set the default values for the variables
            call initialise_run_time_info(settings,RTI)

        end if

        if (settings%nprior <= 0) then
            nprior = settings%nlive
        else
            nprior = settings%nprior
        end if

        total_time=0
        if(linear_mode(mpi_information)) then
            !===================== LINEAR MODE =========================

            do while(RTI%nlive(1)<nprior)

                ! Generate a random coordinate
                live_point(settings%h0:settings%h1) = random_reals(settings%nDims)

                ! Compute physical coordinates, likelihoods and derived parameters
                time0 = time()
                call calculate_point( loglikelihood, prior, live_point, settings, nlike)
                time1 = time()
                live_point(settings%b0) = settings%logzero

                ! If its valid, and we need more points, add it to the array
                if(live_point(settings%l0)>settings%logzero) then
                    total_time =total_time+ time1-time0

                    call add_point(live_point,RTI%live,RTI%nlive,1) ! Add this point to the array

                    !-------------------------------------------------------------------------------!
                    call write_generating_live_points(settings%feedback,RTI%nlive(1),nprior)
                    !-------------------------------------------------------------------------------!

                    if(settings%write_live) then
                        ! Write the live points to the live_points file
                        write(write_phys_unit,fmt_dbl) live_point(settings%p0:settings%d1), live_point(settings%l0)
                        flush(write_phys_unit) ! flush the unit to force write
                    end if

                end if

            end do


#ifdef MPI
        else 
            !===================== PARALLEL MODE =======================

            if(is_root(mpi_information)) then
                ! The root node just recieves data from all other processors


                active_slaves=mpi_information%nprocs-1 ! Set the number of active processors to the number of slaves

                do while(active_slaves>0) 

                    ! Recieve a point from any slave
                    slave_id = catch_point(live_point,mpi_information)

                    ! If its valid, and we need more points, add it to the array
                    if(live_point(settings%l0)>settings%logzero .and. RTI%nlive(1)<nprior) then

                        call add_point(live_point,RTI%live,RTI%nlive,1) ! Add this point to the array

                        !-------------------------------------------------------------------------------!
                        call write_generating_live_points(settings%feedback,RTI%nlive(1),nprior)
                        !-------------------------------------------------------------------------------!

                        if(settings%write_live) then
                            ! Write the live points to the live_points file
                            write(write_phys_unit,fmt_dbl) live_point(settings%p0:settings%d1), live_point(settings%l0)
                            flush(write_phys_unit) ! flush the unit to force write
                        end if

                    end if


                    if(RTI%nlive(1)<nprior) then
                        call request_point(mpi_information,slave_id)  ! If we still need more points, send a signal to have another go
                    else
                        call no_more_points(mpi_information,slave_id) ! Otherwise, send a signal to stop
                        active_slaves=active_slaves-1                  ! decrease the active slave counter
                    end if

                end do




            else

                ! The slaves simply generate and send points until they're told to stop by the master
                do while(.true.)
        
                    live_point(settings%h0:settings%h1) = random_reals(settings%nDims)       ! Generate a random hypercube coordinate
                    time0 = time()
                    call calculate_point( loglikelihood, prior, live_point, settings,nlike) ! Compute physical coordinates, likelihoods and derived parameters
                    time1 = time()
                    live_point(settings%b0) = settings%logzero
                    if(live_point(settings%l0)>settings%logzero) total_time = total_time + time1-time0
                    call throw_point(live_point,mpi_information)                                    ! Send it to the root node
                    if(.not. more_points_needed(mpi_information)) exit                              ! If we've recieved a kill signal, then exit this loop

                end do
            end if
#endif
        end if !(nprocs case)

#ifdef MPI
        nlike = sum_integers(nlike,mpi_information) ! Gather the likelihood calls onto one node
        total_time = sum_doubles(total_time,mpi_information) ! Sum up the total time taken
#endif


        ! ----------------------------------------------- !
        if(is_root(mpi_information)) call write_finished_generating(settings%feedback)  
        ! ----------------------------------------------- !
        ! Find the average time taken
        speed(1) = total_time/nprior
        call time_speeds(loglikelihood,prior,settings,RTI,speed,mpi_information) 



        if(is_root(mpi_information)) then

            ! Pass over the number of likelihood calls
            RTI%nlike(1) = nlike

            ! Set the local and global loglikelihood bounds
            RTI%i(1)  = minpos(RTI%live(settings%l0,:,1)) ! Find the position of the minimum loglikelihood
            RTI%logLp = RTI%live(settings%l0,RTI%i(1),1)  ! Store the value of the minimum loglikelihood 

            ! Close the file
            if(settings%write_live) close(write_phys_unit)

            if(.not. allocated(RTI%num_repeats) ) allocate(RTI%num_repeats(size(settings%grade_dims)))
            RTI%num_repeats(1) = settings%num_repeats
            RTI%num_repeats(2:) = nint(settings%grade_frac(2:)/(settings%grade_frac(1)+0d0)*RTI%num_repeats(1)*speed(1)/speed(2:))

            ! Set the posterior thinning factor
            if(settings%boost_posterior<0d0) then
                RTI%thin_posterior = 1d0
            else
                RTI%thin_posterior = (settings%boost_posterior+0d0)/(sum(RTI%num_repeats)+0d0)
            end if

        end if





    end subroutine GenerateLivePoints



    subroutine time_speeds(loglikelihood,prior,settings,RTI,speed,mpi_information)
        use settings_module,  only: program_settings
        use run_time_module,   only: run_time_info
        use random_module,   only: random_reals
        use utils_module,    only: normal_fb,stdout_unit,fancy_fb,time
        use calculate_module, only: calculate_point
        use abort_module
#ifdef MPI
        use mpi_module, only: mpi_bundle,is_root,sum_doubles,sum_integers
#else
        use mpi_module, only: mpi_bundle,is_root
#endif

        implicit none

        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
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


        !> Program settings
        type(program_settings), intent(in) :: settings
       
        ! The run time info (very important, see src/run_time_info.f90)
        type(run_time_info) :: RTI

        real(dp),dimension(size(settings%grade_dims)) :: speed

        type(mpi_bundle), intent(in) :: mpi_information

        integer :: i_speed
        integer :: i_live
        integer :: nlike
        integer :: h0,h1
        integer :: failed

        real(dp) :: time0,time1,total_time

        real(dp), dimension(settings%nTotal) :: live_point ! Temporary live point array

        nlike=0

        ! Calculate a slow likelihood
        do 
            live_point(settings%h0:settings%h1) = random_reals(settings%nDims)
            call calculate_point( loglikelihood, prior, live_point, settings, nlike)
            live_point(settings%b0) = settings%logzero
            if (live_point(settings%l0)> settings%logzero) exit
        end do

        if (is_root(mpi_information)) then
            if(settings%feedback>=normal_fb) write(stdout_unit,'(A1,"Speed ",I2," = ",E10.3, " seconds")') char(13), 1, speed(1)
        end if
        do i_speed=2,size(speed)

            h0 = settings%h0+sum(settings%grade_dims(:i_speed-1))
            h1 = settings%h1

            i_live=0
            total_time=0
            
            if (is_root(mpi_information)) then
                if(settings%feedback>=fancy_fb) then
                    write(stdout_unit,'(A1,"Speed ",I2, " = ? (calculating)")',advance='no') char(13), i_speed
                    flush(stdout_unit)
                end if
            end if

            failed = 0
            do while( total_time/settings%grade_frac(i_speed)< speed(1)/settings%grade_frac(1) *settings%nlive / dble(mpi_information%nprocs ))
                live_point(h0:h1) = random_reals(h1-h0+1)

                time0 = time()
                call calculate_point( loglikelihood, prior, live_point, settings, nlike)
                time1 = time()
                live_point(settings%b0) = settings%logzero

                if(live_point(settings%l0)>settings%logzero) then
                    total_time=total_time+time1-time0
                    i_live=i_live+1
                else
                    failed = failed + 1
                    if (failed > 10*(i_live+1)) then
                        ! Do another slow likelihood calculation
                        do 
                            live_point(settings%h0:settings%h1) = random_reals(settings%nDims)
                            call calculate_point( loglikelihood, prior, live_point, settings, nlike)
                            live_point(settings%b0) = settings%logzero
                            if (live_point(settings%l0)> settings%logzero) exit
                        end do
                        failed = 0
                    end if
                end if

            end do
#ifdef MPI
            total_time=sum_doubles(total_time,mpi_information)
            i_live = sum_integers(i_live,mpi_information)
            nlike = sum_integers(nlike,mpi_information)
#endif
            speed(i_speed) = total_time/i_live
            if (is_root(mpi_information)) then
                RTI%nlike(i_speed) =  RTI%nlike(i_speed) + nlike 
                if(settings%feedback>=fancy_fb) then
                    write(stdout_unit,'(A1,"Speed ",I2," = ",E10.3, " seconds     ")') char(13), i_speed, speed(i_speed)
                else if(settings%feedback>=normal_fb) then
                    write(stdout_unit,'("Speed ",I2," = ",E10.3, " seconds     ")') i_speed, speed(i_speed)
                end if
            end if


        end do

    end subroutine time_speeds


    subroutine GenerateLivePointsFromSeed(loglikelihood,prior,settings,RTI,mpi_information)
        use settings_module,  only: program_settings
        use utils_module,    only: write_phys_unit,DB_FMT,fmt_len,minpos,time,stdout_unit
        use calculate_module, only: calculate_point
        use read_write_module, only: phys_live_file
        use feedback_module,  only: write_started_generating,write_finished_generating,write_generating_live_points
        use run_time_module,   only: run_time_info,initialise_run_time_info
        use array_module,     only: add_point
        use abort_module
        use chordal_module, only: slice_sample
#ifdef MPI
        use mpi_module, only: mpi_bundle,is_root,linear_mode,throw_point,catch_point,more_points_needed,sum_integers,sum_doubles,request_point,no_more_points
#else
        use mpi_module, only: mpi_bundle,is_root,linear_mode
#endif

        implicit none

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


        !> Program settings
        type(program_settings), intent(in) :: settings

        ! The run time info (very important, see src/run_time_info.f90)
        type(run_time_info) :: RTI

        type(mpi_bundle),intent(in) :: mpi_information
#ifdef MPI
        integer             :: active_slaves    !  Number of currently working slaves
        integer             :: slave_id         !  Slave identifier to signal who to throw back to
#endif

        real(dp), dimension(settings%nTotal) :: live_point ! Temporary live point array


        character(len=fmt_len) :: fmt_dbl ! writing format variable

        integer :: nprior

        real(dp) :: time0,time1
        real(dp),dimension(size(settings%grade_dims)) :: speed
        real(dp),dimension(size(settings%grade_dims)) :: times
        integer, dimension(size(settings%grade_dims)) :: nlikes

        integer :: i_dim,i_repeat,i_grade
        integer, dimension(settings%nDims) :: grade
        real(dp), dimension(settings%nDims) :: nhat

        do i_grade = size(settings%grade_dims),1,-1
            grade(:sum(settings%grade_dims(:i_grade))) = i_grade
        end do
        nlikes = 0
        times = 0d0

        if(is_root(mpi_information)) then
            ! ---------------------------------------------- !
            call write_started_generating(settings%feedback)
            ! ---------------------------------------------- !

            ! Initialise the format
            write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT

            ! Open the live points file to sequentially add live points
            if(settings%write_live) open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')

            ! Allocate the run time arrays, and set the default values for the variables
            call initialise_run_time_info(settings,RTI)

        end if

        if (settings%nprior <= 0) then
            nprior = settings%nlive
        else
            nprior = settings%nprior
        end if

        if(linear_mode(mpi_information)) then
            !===================== LINEAR MODE =========================

            do while(RTI%nlive(1)<nprior)
                live_point = settings%seed_point
                do i_repeat = 1,settings%nprior_repeat
                    i_dim = 1+mod(i_repeat,settings%nDims)
                    i_grade = grade(i_dim)
                    nhat = 0d0
                    nhat(i_dim) = 1d0

                    time0 = time()
                    live_point = slice_sample(loglikelihood,prior,settings%logzero,nhat,live_point,1d0,settings,nlikes(i_grade))
                    time1 = time()
                    times(i_grade) = times(i_grade) + time1 - time0
                end do

                call add_point(live_point,RTI%live,RTI%nlive,1) ! Add this point to the array
                call write_generating_live_points(settings%feedback,RTI%nlive(1),nprior)
                if(settings%write_live) then
                    ! Write the live points to the live_points file
                    write(write_phys_unit,fmt_dbl) live_point(settings%p0:settings%d1), live_point(settings%l0)
                    flush(write_phys_unit) ! flush the unit to force write
                end if

            end do

#ifdef MPI
        else 
            !===================== PARALLEL MODE =======================

            if(is_root(mpi_information)) then
                ! The root node just recieves data from all other processors

                active_slaves=mpi_information%nprocs-1 ! Set the number of active processors to the number of slaves

                do while(active_slaves>0) 

                    ! Recieve a point from any slave
                    slave_id = catch_point(live_point,mpi_information)

                    ! If its valid, and we need more points, add it to the array
                    if(RTI%nlive(1)<nprior) then

                        call add_point(live_point,RTI%live,RTI%nlive,1) ! Add this point to the array

                        !-------------------------------------------------------------------------------!
                        call write_generating_live_points(settings%feedback,RTI%nlive(1),nprior)
                        !-------------------------------------------------------------------------------!

                        if(settings%write_live) then
                            ! Write the live points to the live_points file
                            write(write_phys_unit,fmt_dbl) live_point(settings%p0:settings%d1), live_point(settings%l0)
                            flush(write_phys_unit) ! flush the unit to force write
                        end if

                    end if


                    if(RTI%nlive(1)<nprior) then
                        call request_point(mpi_information,slave_id)  ! If we still need more points, send a signal to have another go
                    else
                        call no_more_points(mpi_information,slave_id) ! Otherwise, send a signal to stop
                        active_slaves=active_slaves-1                  ! decrease the active slave counter
                    end if

                end do




            else

                ! The slaves simply generate and send points until they're told to stop by the master
                live_point = settings%seed_point
                do while(.true.)
                    do i_repeat = 1,5
                        do i_dim=1,settings%nDims
                            i_grade = grade(i_dim)
                            nhat = 0d0
                            nhat(i_dim) = 1d0

                            time0 = time()
                            live_point = slice_sample(loglikelihood,prior,settings%logzero,nhat,live_point,1d0,settings,nlikes(i_grade))
                            time1 = time()
                            times(i_grade) = times(i_grade) + time1 - time0
                        end do
                    end do
        
                    call throw_point(live_point,mpi_information)                                    ! Send it to the root node
                    if(.not. more_points_needed(mpi_information)) exit                              ! If we've recieved a kill signal, then exit this loop

                end do
            end if
#endif
        end if !(nprocs case)

#ifdef MPI
        do i_grade=1,size(settings%grade_dims)
            nlikes(i_grade) = sum_integers(nlikes(i_grade),mpi_information) ! Gather the likelihood calls onto one node
            times(i_grade) = sum_doubles(times(i_grade),mpi_information) ! Sum up the total time taken
        end do
#endif


        ! ----------------------------------------------- !
        if(is_root(mpi_information)) call write_finished_generating(settings%feedback)  
        ! ----------------------------------------------- !

        if(is_root(mpi_information)) then
            ! Find the average time taken
            speed = times/nlikes

            do i_grade = 1,size(settings%grade_dims)
                write(stdout_unit,'("Speed ",I2," = ",E10.3, " seconds     ")') i_grade, speed(i_grade)
            end do


            ! Pass over the number of likelihood calls
            RTI%nlike = nlikes

            ! Set the local and global loglikelihood bounds
            RTI%i(1)  = minpos(RTI%live(settings%l0,:,1)) ! Find the position of the minimum loglikelihood
            RTI%logLp = RTI%live(settings%l0,RTI%i(1),1)  ! Store the value of the minimum loglikelihood 

            ! Close the file
            if(settings%write_live) close(write_phys_unit)

            if(.not. allocated(RTI%num_repeats) ) allocate(RTI%num_repeats(size(settings%grade_dims)))
            RTI%num_repeats(1) = settings%num_repeats
            RTI%num_repeats(2:) = nint(settings%grade_frac(2:)/(settings%grade_frac(1)+0d0)*RTI%num_repeats(1)*speed(1)/speed(2:))

            ! Set the posterior thinning factor
            if(settings%boost_posterior<0d0) then
                RTI%thin_posterior = 1d0
            else
                RTI%thin_posterior = (settings%boost_posterior+0d0)/(sum(RTI%num_repeats)+0d0)
            end if

        end if





    end subroutine GenerateLivePointsFromSeed







end module generate_module
