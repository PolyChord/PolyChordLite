module mpi_module

    use mpi
    implicit none

    integer :: mpierror

    integer, parameter :: tag_gen_point=1
    integer, parameter :: tag_gen_request=2
    integer, parameter :: tag_gen_stop=3

    integer, parameter :: tag_run_baby=4
    integer, parameter :: tag_run_seed=5
    integer, parameter :: tag_run_cholesky=6
    integer, parameter :: tag_run_logL=7
    integer, parameter :: tag_run_nlike=8
    integer, parameter :: tag_run_stop=9



    contains


    !> Procedure to get the number of mpi processors
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Comm_size.html
    function get_nprocs(mpi_communicator) result(nprocs)
        implicit none
        integer, intent(in) :: mpi_communicator
        integer :: nprocs

        call MPI_COMM_SIZE(   & 
            mpi_communicator, &!handle
            nprocs,           &!return number of processors
            mpierror          &!error flag
            )

    end function get_nprocs

    !> Procedure to get the number of mpi processors
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Comm_rank.html
    function get_rank(mpi_communicator) result(myrank)
        implicit none
        integer, intent(in) :: mpi_communicator
        integer :: myrank

        call MPI_COMM_RANK(   & 
            mpi_communicator, &!handle
            myrank,           &!return rank of calling processor 
            mpierror          &!error flag
            )

    end function get_rank


    !> This defines the root node for this mpi communicator
    !!
    !! All processes call this simultaneously, and the root is defined
    !! as the process with the lowest rank. This value is returned
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Allreduce.html
    !!
    function get_root(myrank,mpi_communicator) result(root)
        implicit none
        integer, intent(in) :: myrank
        integer, intent(in) :: mpi_communicator
        integer :: root


        call MPI_ALLREDUCE(    &
            myrank,            &!send buffer 
            root,              &!recieve buffer
            1,                 &!number of elements sent
            MPI_INTEGER,       &!type of element sent
            MPI_MIN,           &!reduce by finding the minimum
            mpi_communicator,  &!handle
            mpierror           &!error flag
            )

    end function get_root


    !> This sums a whole set of integers across all processes
    !!
    !! All processes call this simultaneously
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Allreduce.html
    !!
    function sum_nlike(nlike_local,mpi_communicator) result(nlike)
        implicit none
        integer, intent(in) :: nlike_local
        integer, intent(in) :: mpi_communicator
        integer :: nlike

        call MPI_ALLREDUCE(    &
            nlike_local,       &!send buffer 
            nlike,             &!recieve buffer
            1,                 &!number of elements sent
            MPI_INTEGER,       &!type of element sent
            MPI_SUM,           &!reduce by finding the minimum
            mpi_communicator,  &!handle
            mpierror           &!error flag
            )

    end function sum_nlike

    !============== Throwing and catching routines ====================
    ! There are four points in the algorithm where arrays need to be transferred 
    ! between master and slaves
    !
    ! 1) live point during initial generation
    !     any slave    ---->  root
    !     throw_point         catch_point
    !
    ! 2) baby points for use in nested sampling
    !     any slave   ----> root
    !     throw_babies      catch_babies
    !
    ! 3) seed information for generation of babies
    !     root      ----> specific slave
    !     throw_seed      catch_seed


    !> Root catch
    !!
    !! This a process by which the root node 'catches' a thrown point from 
    !! any slave, and returns the slave identifier for use to throw back
    function catch_point(live_point,mpi_communicator) result(slave_id)
        implicit none

        double precision,intent(out),dimension(:) :: live_point !> The caught live point
        integer, intent(in) :: mpi_communicator                 !> The mpi communicator

        integer :: slave_id ! slave identifier

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status ! status identifier

        call MPI_RECV(             &!
            live_point,            &!
            size(live_point),      &!
            MPI_DOUBLE_PRECISION,  &!
            MPI_ANY_SOURCE,        &!
            tag_gen_point,         &!
            mpi_communicator,      &!
            mpi_status,            &!
            mpierror               &!
            )

        slave_id = mpi_status(MPI_SOURCE) ! Pass on the slave id

    end function catch_point


    !> Slave throws a point to the master
    !!
    !! This a process by which a slave node 'throws' a point to the root
    subroutine throw_point(live_point,mpi_communicator,root)
        implicit none

        double precision,intent(in),dimension(:) :: live_point !> live point to throw
        integer, intent(in) :: mpi_communicator                !> mpi communicator
        integer, intent(in) :: root                            !> 

        call MPI_SEND(             &!
            live_point,            &!
            size(live_point),      &!
            MPI_DOUBLE_PRECISION,  &!
            root,                  &!
            tag_gen_point,         &!
            mpi_communicator,      &!
            mpierror               &!
            )

    end subroutine throw_point







    !> Master catches babies thrown by any slave, and returns the slave identity that did the throwing
    function catch_babies(baby_points,nlike,mpi_communicator) result(slave_id)
        implicit none

        double precision,intent(out),dimension(:,:) :: baby_points !> The babies to be caught
        integer, intent(out) :: nlike                              !> The number of likelihood evaluations to be caught
        integer, intent(in) :: mpi_communicator                    !> The mpi communicator

        integer :: slave_id ! slave identifier

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status ! status identifier

        call MPI_RECV(                               &!
            baby_points,                             &!
            size(baby_points,1)*size(baby_points,2), &!
            MPI_DOUBLE_PRECISION,                    &!
            MPI_ANY_SOURCE,                          &!
            tag_run_baby,                            &!
            mpi_communicator,                        &!
            mpi_status,                              &!
            mpierror                                 &!
            )

        ! Pass on the slave id
        slave_id = mpi_status(MPI_SOURCE)

        call MPI_RECV(         &! 
            nlike,             &! 
            1,                 &! 
            MPI_INT,           &! 
            slave_id,          &! 
            tag_run_nlike,     &! 
            mpi_communicator,  &! 
            mpi_status,        &! 
            mpierror           &! 
            )

    end function catch_babies

    !> Slave throws babies to the master
    subroutine throw_babies(baby_points,nlike,mpi_communicator,root)
        implicit none

        double precision,intent(in),dimension(:,:) :: baby_points !> The babies to be thrown
        integer, intent(in) :: nlike                              !> The number of likelihood evaluations to be thrown
        integer, intent(in) :: mpi_communicator                   !> The mpi communicator
        integer, intent(in) :: root                               !> root node to throw to

        call MPI_SEND(                               &! 
            baby_points,                             &! 
            size(baby_points,1)*size(baby_points,2), &! 
            MPI_DOUBLE_PRECISION,                    &! 
            root,                                    &! 
            tag_run_baby,                            &! 
            mpi_communicator,                        &! 
            mpierror                                 &! 
            )                                        
        call MPI_SEND(         &!  
            nlike,             &!  
            1,                 &!  
            MPI_INT,           &!  
            root,              &!  
            tag_run_nlike,     &!  
            mpi_communicator,  &!  
            mpierror           &!  
            )

    end subroutine throw_babies






    !> slave catches seed thrown by master
    function catch_seed(seed_point,cholesky,logL,mpi_communicator,root) result(more_points_needed)
        use abort_module, only: halt_program
        implicit none


        double precision,intent(out),dimension(:) :: seed_point  !> The seed point to be caught
        double precision,intent(out),dimension(:,:) :: cholesky  !> Cholesky matrix to be caught
        double precision,intent(in)                :: logL       !> loglikelihood contour to be caught
        integer, intent(in) :: mpi_communicator                  !> The mpi communicator
        integer, intent(in) :: root                              !> The root node

        logical :: more_points_needed ! whether or not we need more points

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status ! status identifier


        call MPI_RECV(                      &!
            seed_point,                     &!
            size(seed_point),               &!
            MPI_DOUBLE_PRECISION,           &!
            root,                           &!
            MPI_ANY_TAG,                    &!
            mpi_communicator,               &!
            mpi_status,                     &!
            mpierror                        &!
            )
        if(mpi_status(MPI_TAG) == tag_run_stop ) then
            more_points_needed = .false.
            return
        else if(mpi_status(MPI_TAG) == tag_run_seed) then
            more_points_needed = .true.
        else
            call halt_program('slave error: unrecognised tag')
        end if

        call MPI_RECV(                        &!
            cholesky,                         &!
            size(cholesky,1)*size(cholesky,1),&!
            MPI_DOUBLE_PRECISION,             &!
            root,                             &!
            tag_run_cholesky,                 &!
            mpi_communicator,                 &!
            mpi_status,                       &!
            mpierror                          &!
            )
        call MPI_RECV(                 &!
            logL,                      &!
            1,                         &!
            MPI_DOUBLE_PRECISION,      &!
            root,                      &!
            tag_run_logL,              &!
            mpi_communicator,          &!
            mpi_status,                &!
            mpierror                   &!
            )

    end function catch_seed



    !> root throws seed to slave
    subroutine throw_seed(seed_point,cholesky,logL,mpi_communicator,slave_id,keep_going)
        implicit none

        double precision,intent(in),dimension(:) :: seed_point   !> seed to be thrown
        double precision,intent(in),dimension(:,:) :: cholesky   !> cholesky to be thrown
        double precision,intent(in)                :: logL       !> loglikelihood contour to be thrown
        integer, intent(in) :: mpi_communicator                  !> mpi handle
        integer, intent(in) :: slave_id                          !> identity of target slave
        logical, intent(in) :: keep_going                        !> Further signal whether to keep going 

        integer :: tag ! tag variable to

        tag = tag_run_stop                 ! Default tag is stop tag
        if(keep_going) tag = tag_run_seed  ! If we want to keep going then change this to the seed tag


        call MPI_SEND(             &!  
            seed_point,            &!  
            size(seed_point),      &!  
            MPI_DOUBLE_PRECISION,  &!  
            slave_id,              &!  
            tag,                   &!  
            mpi_communicator,      &!  
            mpierror               &!  
            )

        if(.not. keep_going) return ! Stop here if we're wrapping up

        call MPI_SEND(                          &!  
            cholesky,                           &!  
            size(cholesky,1)*size(cholesky,2),  &!  
            MPI_DOUBLE_PRECISION,               &!  
            slave_id,                           &!  
            tag_run_cholesky,                   &!  
            mpi_communicator,                   &!  
            mpierror                            &!  
            )
        call MPI_SEND(              &!  
            logL,                   &!  
            1,                      &!  
            MPI_DOUBLE_PRECISION,   &!  
            slave_id,               &!  
            tag_run_logL,           &!  
            mpi_communicator,       &!  
            mpierror                &!  
            )


    end subroutine throw_seed




    !============== Pure messaging routines ===========================
    ! During initial live point generation, the master needs to signal to the slaves
    ! whether or not to keep generating live points, or whether to stop
    !
    ! master         ----> slave
    ! request_point        more_points_needed -> true
    ! no_more_points       more_points_needed -> false
    ! 

    !> Request point
    !!
    !! This subroutine is used by the root node to request a new live point
    subroutine request_point(mpi_communicator,slave_id)
        implicit none
        integer, intent(in) :: mpi_communicator !> mpi handle
        integer, intent(in) :: slave_id         !> Slave to request a new point from


        integer :: empty_buffer(0) ! empty buffer to send

        call MPI_SEND(         &
            empty_buffer,      &! not sending anything
            0,                 &! size of nothing
            MPI_INT,           &! sending no integers
            slave_id,          &! process id to send to
            tag_gen_request,   &! continuation tag
            mpi_communicator,  &! mpi handle
            mpierror           &! error flag
            )

    end subroutine request_point


    !> No more points please
    !!
    !! This subroutine is used by the root node to signal that no more points are required
    subroutine no_more_points(mpi_communicator,slave_id)
        implicit none
        integer, intent(in) :: mpi_communicator !> mpi handle
        integer, intent(in) :: slave_id         !> Slave to request a new point from


        integer :: empty_buffer(0) ! empty buffer to send

        call MPI_SEND(         &
            empty_buffer,      &! not sending anything
            0,                 &! size of nothing
            MPI_INT,           &! sending no integers
            slave_id,          &! process id to send to
            tag_gen_stop,      &! continuation tag
            mpi_communicator,  &! mpi handle
            mpierror           &! error flag
            )

    end subroutine no_more_points

    !> See if more points are needed
    !!
    !! This subroutine is used by the root node to request a new live point
    function more_points_needed(mpi_communicator,root)
        use abort_module
        implicit none
        integer, intent(in) :: mpi_communicator !> mpi handle
        integer, intent(in) :: root             !> root to recieve request

        integer :: empty_buffer(0) ! empty buffer to send
        integer, dimension(MPI_STATUS_SIZE) :: mpi_status  ! status identifier

        logical :: more_points_needed !> Whether we need more points or not

        call MPI_RECV(       &!
            empty_buffer,    &!
            0,               &!
            MPI_INT,         &!
            root,            &!
            MPI_ANY_TAG,     &!
            mpi_communicator,&!
            mpi_status,      &!
            mpierror         &!
            )

        ! If we've recieved a kill signal, then exit this loop
        if(mpi_status(MPI_TAG) == tag_gen_stop ) then
            more_points_needed = .false.
        else if(mpi_status(MPI_TAG) == tag_gen_request) then
            more_points_needed = .true.
        else
            call halt_program('generate error: unrecognised tag')
        end if

    end function more_points_needed








end module mpi_module
