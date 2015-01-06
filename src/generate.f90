!> This module contains 'generating tools', namely:
!!
!! * GenerateSeed
!! ** Generates a seed point
!! * GenerateLivePointsP
!! * GenerateLivePointsL
module generate_module

#ifdef MPI
    use mpi_module
#endif

    implicit none

    contains


    !> This function generates a seed point for slice sampling.
    !!
    !! It uses the existing set of live points and the details of the clustering
    !! in order to select a seed point to generate from in proportion to the
    !! estimates of the prior volume of each cluster.
    function GenerateSeed(settings,info,live_points,cluster_choice) result(seed_point)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use random_module,     only: random_integer,random_integer_P
        use utils_module,      only: logsumexp
        implicit none

        !> Program settings
        type(program_settings), intent(in) :: settings
        !> The evidence storage
        type(run_time_info), intent(in) :: info
        !> The live points
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(in) :: live_points
        !> The cluster number chosen
        integer,intent(out) :: cluster_choice


        ! The seed point to be produced
        double precision, dimension(settings%nTotal) :: seed_point

        integer :: seed_choice

        double precision, dimension(info%ncluster_A) :: probs


        ! 0) Calculate an array proportional to the volumes
        probs = info%logX(:info%ncluster_A) ! prob_i = log( X_i )
        probs = probs - logsumexp(probs)    ! prob_i = log( X_i/(sum_j X_j) )
        probs = exp(probs)                  ! prob_i = X_i/(sum_j X_j)

        ! 1) Pick cluster in proportion to the set of volume estimates of the active clusters
        cluster_choice = random_integer_P(probs)

        ! 2) Pick a random integer in between 1 and the number of live points in the chosen cluster
        seed_choice = random_integer(info%n(cluster_choice))

        ! 3) Select the live point at index 'seed_choice' in cluster 'cluster_choice' for the seed point
        seed_point = live_points(:,seed_choice,cluster_choice)
        
        ! 4) Give the seed point the likelihood contour of cluster ! 'cluster_choice'
        seed_point(settings%l1) = info%logL(cluster_choice)

    end function GenerateSeed







#ifdef MPI

    function GenerateLivePointsFromSeedP(loglikelihood,priors,settings,mpi_communicator,root) result(live_points) 
        use priors_module,    only: prior
        use settings_module,  only: program_settings
        use random_module,   only: random_reals,random_distinct_integers
        use utils_module,    only: write_phys_unit,DB_FMT,fmt_len
        use calculate_module, only: calculate_point
        use read_write_module, only: phys_live_file
        use feedback_module,  only: write_started_generating,write_finished_generating
        use utils_module, only: calc_cholesky,calc_covmat

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
        integer :: mpierror

        double precision, dimension(settings%nTotal,settings%nlive) :: live_points
        double precision, dimension(settings%nTotal,settings%ngenerate) :: live_stack
        double precision, dimension(settings%nTotal,settings%num_babies)   :: baby_points
        double precision,    dimension(settings%nTotal)   :: seed_point


        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        integer :: ngenerate

        integer :: i_dims
        integer :: i_baby
        integer :: i_slave
        logical :: slave_sending

        integer, dimension(:), allocatable :: burn_in

        double precision, dimension(settings%nTotal,0) :: blank_stack
        double precision, dimension(settings%nDims,settings%nDims) :: cholesky
        double precision, dimension(settings%nDims,settings%nDims) :: covmat

        character(len=fmt_len) :: fmt_dbl_nphys

        ! Initialise the formats
        write(fmt_dbl_nphys,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT

        ! Get the number of MPI procedures
        call MPI_COMM_SIZE(mpi_communicator, nprocs, mpierror)
        ! Get the MPI label of the current processor
        call MPI_COMM_RANK(mpi_communicator, myrank, mpierror)

        ! initialise live points at zero
        live_points = 0d0

        ngenerate = 0

        if(myrank==root) then


            allocate(burn_in(nprocs-1))
            burn_in = 0

            ! Initialise the covmats at the identity
            covmat = 0d0
            cholesky=0d0
            do i_dims=1,settings%nDims
                covmat(i_dims,i_dims) = 1d0
                cholesky(i_dims,i_dims) = 1d0
            end do

            !======= 2) Main loop body =====================================

            ! -------------------------------------------- !
            call write_started_generating(settings%feedback)
            ! -------------------------------------------- !
            open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')

            do while ( ngenerate<=settings%ngenerate )

                    ! Recieve any new baby points from any slave currently sending

                    ! Loop through all the slaves
                    do i_slave=1,nprocs-1
                        ! Use MPI_IPROBE to see if the slave at i_slave is sending
                        ! If it is, the logical variable 'slave_sending' is true
                        call MPI_IPROBE(i_slave,MPI_ANY_TAG,mpi_communicator,slave_sending,mpi_status,mpierror)

                        if(slave_sending) then
                            ! If this slave is sending, then recieve the newly generated baby points
                            call MPI_RECV(baby_points,settings%nTotal*settings%num_babies,&
                                MPI_DOUBLE_PRECISION,i_slave,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)

                            ! If these baby points aren't nonsense (i.e. the first send) ...
                            if(mpi_status(MPI_TAG)/=tag_run_no_points) then

                                burn_in(i_slave) = burn_in(i_slave)+1

                                if(burn_in(i_slave)>settings%generate_burn_in) then

                                    live_stack(:,ngenerate+1:min(settings%ngenerate,ngenerate+settings%num_babies))= baby_points
                                    ngenerate = ngenerate+settings%num_babies

                                    do i_baby=1,settings%num_babies
                                        write(write_phys_unit,fmt_dbl_nphys) baby_points(settings%p0:settings%d1,i_baby), baby_points(settings%l0,i_baby)
                                    end do
                                    call flush(write_phys_unit)
                                end if

                                ! choose the seed point as the last baby point
                                seed_point = baby_points(:,settings%num_babies)
                            else
                                ! Otherwise initialise at the start seed
                                seed_point = settings%seed_point
                            end if

                            ! Send the seed point back to this slave
                            call MPI_SEND(seed_point,settings%nTotal,MPI_DOUBLE_PRECISION,i_slave,tag_run_new_seed,mpi_communicator,mpierror)

                            ! Send the cholesky decomposition
                            call MPI_SEND(cholesky,settings%nDims*settings%nDims,MPI_DOUBLE_PRECISION,i_slave,tag_run_new_cholesky,mpi_communicator,mpierror)
                        end if !(slave_sending)

                    end do !(i_slave=1,nprocs-1)

                    if(ngenerate>settings%nlive) then
                        covmat   = calc_covmat(live_stack(settings%h0:settings%h1,:ngenerate),blank_stack)
                        cholesky = calc_cholesky(covmat)
                    end if

            end do ! End main loop
            close(write_phys_unit)

            ! Now thin the stack
            live_points = live_stack(:,random_distinct_integers(settings%nstack,settings%nlive))
        end if


    end function GenerateLivePointsFromSeedP










    !> Generate an initial set of live points distributed uniformly in the unit hypercube in parallel
    function GenerateLivePointsP(loglikelihood,priors,settings,mpi_communicator,root) result(live_points)
        use priors_module,    only: prior
        use settings_module,  only: program_settings
        use random_module,   only: random_reals
        use utils_module,    only: logzero,write_phys_unit,DB_FMT,fmt_len
        use calculate_module, only: calculate_point
        use read_write_module, only: phys_live_file
        use feedback_module,  only: write_started_generating,write_finished_generating

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
        integer :: active_slaves
        integer :: mpierror

        double precision, dimension(settings%nTotal,settings%nlive) :: live_points

        !live_points(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(settings%nTotal) :: live_point

        ! Loop variable
        integer i_live

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        integer :: empty_buffer(0)

        integer :: nlike
        character(len=fmt_len) :: fmt_dbl_nphys

        ! Initialise the formats
        write(fmt_dbl_nphys,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT

        ! Get the number of MPI procedures
        call MPI_COMM_SIZE(mpi_communicator, nprocs, mpierror)
        ! Get the MPI label of the current processor
        call MPI_COMM_RANK(mpi_communicator, myrank, mpierror)

        ! initialise live points at zero
        live_points = 0d0

        if(myrank==root) then
            ! The root node just recieves data from all other processors

            ! ---------------------------------------------- !
            call write_started_generating(settings%feedback)
            ! ---------------------------------------------- !

            
            active_slaves=nprocs-1 ! Set the number of active processors to the number of slaves
            i_live=0               ! No live points initially
            nlike=0                ! No wasted likelihood calls initially

            ! Open the live points file to sequentially add live points
            open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')

            do while(active_slaves>0) 

                ! Recieve a point from any slave
                call MPI_RECV(live_point,settings%nTotal, &
                    MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag_gen_new_point,mpi_communicator,mpi_status,mpierror)

                ! If its valid, and we need more points, add it to the array
                if(live_point(settings%l0)>logzero .and. i_live<settings%nlive) then

                    ! Increase the live point counter
                    i_live=i_live+1
                    ! Add the new live point to the live point array
                    live_points(:,i_live) = live_point
                    ! Write the live points to the live_points file
                    write(write_phys_unit,fmt_dbl_nphys) live_point(settings%p0:settings%d1), live_point(settings%l0)
                    call flush(write_phys_unit)

                else
                    ! If it failed for whatever reason, then record that we've
                    ! had a 'lost' likelihood evaluation
                    nlike = nlike+1
                end if


                if(i_live<settings%nlive) then
                    ! If we still need more points, send a signal to have another go
                    call MPI_SEND(empty_buffer,0,MPI_INT,mpi_status(MPI_SOURCE),tag_gen_continue,mpi_communicator,mpierror)
                else
                    ! Otherwise, send a signal to stop
                    call MPI_SEND(empty_buffer,0,MPI_INT,mpi_status(MPI_SOURCE),tag_gen_stop,mpi_communicator,mpierror)
                    ! and decrease the counter for the number of active slaves
                    active_slaves=active_slaves-1
                end if

            end do


            ! Set the initial trial values of the chords as 1
            live_points(settings%last_chord,:) = 1

            ! Set the likelihood contours to logzero for now
            live_points(settings%l1,:) = logzero

            ! Add the the number of 'wasted' likelihood calls to the first live
            ! point
            live_points(settings%nlike,1) = nlike

            ! Close the file
            close(write_phys_unit)

            ! ----------------------------------------------- !
            call write_finished_generating(settings%feedback)  
            ! ----------------------------------------------- !


        else

            ! The slaves simply generate and send points until they're told to stop by
            ! the master

            do while(.true.)

                ! Zero the likelihood calls 
                live_point(settings%nlike) = 0

                ! Generate a random hypercube coordinate
                live_point(settings%h0:settings%h1) = random_reals(settings%nDims)

                ! Compute physical coordinates, likelihoods and derived parameters
                call calculate_point( loglikelihood, priors, live_point, settings )

                ! Send it to the root node
                call MPI_SEND(live_point,settings%nTotal, &
                    MPI_DOUBLE_PRECISION,root,tag_gen_new_point,mpi_communicator,mpierror)

                ! Recieve signal as to whether we should keep generating
                call MPI_RECV(empty_buffer,0,MPI_INT,root,MPI_ANY_TAG,mpi_communicator,mpi_status,mpierror)

                ! If we've recieved a kill signal, then exit this loop
                if(mpi_status(MPI_TAG) == tag_gen_stop ) exit

            end do

        end if

    end function GenerateLivePointsP
#endif


    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePointsL(loglikelihood,priors,settings) result(live_points)
        use priors_module,    only: prior
        use settings_module,  only: program_settings
        use random_module,    only: random_reals
        use utils_module,     only: logzero
        use calculate_module, only: calculate_point
        use feedback_module,  only: write_started_generating,write_finished_generating

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

        ! ---------------------------------------------- !
        call write_started_generating(settings%feedback)
        ! ---------------------------------------------- !

        ! initialise live points at zero
        live_points = 0d0

        do i_live=1,settings%nlive

            ! Generate a random coordinate
            live_points(settings%h0:settings%h1,i_live) = random_reals(settings%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( loglikelihood, priors, live_points(:,i_live), settings )

        end do

        ! Set the number of likelihood calls for each point to 1
        live_points(settings%nlike,:) = 1

        ! Set the initial trial values of the chords as the diagonal of the hypercube
        live_points(settings%last_chord,:) = sqrt(settings%nDims+0d0)

        ! Set the likelihood contours to logzero for now
        live_points(settings%l1,:) = logzero

        ! ----------------------------------------------- !
        call write_finished_generating(settings%feedback)  
        ! ----------------------------------------------- !
    end function GenerateLivePointsL


end module generate_module
