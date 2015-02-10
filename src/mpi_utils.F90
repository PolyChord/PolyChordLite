module mpi_module

    use mpi
    implicit none

    integer :: mpierror

    integer, parameter :: STARTTAG=0
    integer, parameter :: RUNTAG=1
    integer, parameter :: ENDTAG=2

    integer, parameter :: tag_run_new_points=0
    integer, parameter :: tag_run_no_points=2
    integer, parameter :: tag_run_new_seed=3
    integer, parameter :: tag_run_new_cholesky=4
    integer, parameter :: tag_run_end=5



    integer, parameter :: tag_gen_continue=6
    integer, parameter :: tag_gen_stop=7
    integer, parameter :: tag_gen_new_point=8



    contains

    subroutine abort_all(message)
        use utils_module, only: stdout_unit
        implicit none

        character(LEN=*), intent(in), optional :: message

        integer :: errorcode
        integer :: mpierror

        integer :: strlen


        if (present(message)) then
            strlen = len(trim(adjustl(message)))
            write(stdout_unit,'( 20("=") )')
            write(stdout_unit,'(A)') trim(adjustl(message))
            write(stdout_unit,'( 20("=") )')
        end if 


        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)

    end subroutine abort_all


    !> Procedure to get the number of mpi processors
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Comm_size.html
    function get_nprocs(mpi_communicator) result(nprocs)
        implicit none
        integer, intent(in) :: mpi_communicator
        integer :: nprocs

        call MPI_COMM_SIZE(       & 
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

        call MPI_COMM_RANK(       & 
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


        call MPI_ALLREDUCE(        &
                myrank,            &!send buffer 
                root,              &!recieve buffer
                1,                 &!number of elements sent
                MPI_INTEGER,       &!type of element sent
                MPI_MIN,           &!reduce by finding the minimum
                mpi_communicator,  &!handle
                mpierror           &!error flag
                )

    end function assign_root


    !> Root catch
    !!
    !! This a process by which the root node 'catches' a thrown point from 
    !! any slave, and returns the slave identifier for use to throw back
    function catch_point(live_point,mpi_communicator) result(slave_id)
        implicit none

        !> The caught live point
        double precision,intent(out),dimension(:) :: live_point

        !> The mpi communicator
        integer, intent(in) :: mpi_communicator

        ! slave identifier
        integer :: slave_id

        ! status identifier
        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        call MPI_RECV(&
                live_point,&
                size(live_point), &
                MPI_DOUBLE_PRECISION,&
                MPI_ANY_SOURCE,&
                tag_gen_new_point,&
                mpi_communicator,&
                mpi_status,&
                mpierror&
                )

        ! Pass on the slave id
        slave_id = mpi_status(MPI_SOURCE)

    end function catch_point


    !> Slave throw
    !!
    !! This a process by which a slave node 'throws' a point to the root
    subroutine throw_point(live_point,mpi_communicator,root) result(slave_id)
        implicit none

        !> The caught live point
        double precision,intent(out),dimension(:) :: live_point

        !> The mpi communicator
        integer, intent(in) :: mpi_communicator
        integer, intent(in) :: root

        ! slave identifier
        integer :: slave_id

        call MPI_SEND(&
                live_point,&
                size(live_point),&
                MPI_DOUBLE_PRECISION,&
                root,&
                tag_gen_new_point,&
                mpi_communicator,&
                mpierror&
                )

    end subroutine throw_point

    !> Request point
    !!
    !! This subroutine is used by the root node to request a new live point
    subroutine request_point(slave_id,mpi_communicator)
        implicit none
        integer, intent(in) :: mpi_communicator !> mpi handle
        integer, intent(in) :: slave_id         !> Slave to request a new point from

        
        integer :: empty_buffer(0) ! empty buffer to send

        call MPI_SEND(&
                empty_buffer,&       ! not sending anything
                0,&                  ! size of nothing
                MPI_INT,&            ! sending no integers
                slave_id,&           ! process id to send to
                tag_gen_continue,&   ! continuation tag
                mpi_communicator,&   ! mpi handle
                mpierror&            ! error flag
                )

    end subroutine request_point


    !> No more points please
    !!
    !! This subroutine is used by the root node to signal that no more points are required
    subroutine no_more_points(slave_id,mpi_communicator)
        implicit none
        integer, intent(in) :: mpi_communicator !> mpi handle
        integer, intent(in) :: slave_id         !> Slave to request a new point from

        
        integer :: empty_buffer(0) ! empty buffer to send

        call MPI_SEND(&
                empty_buffer,&       ! not sending anything
                0,&                  ! size of nothing
                MPI_INT,&            ! sending no integers
                slave_id,&           ! process id to send to
                tag_gen_continue,&   ! continuation tag
                mpi_communicator,&   ! mpi handle
                mpierror&            ! error flag
                )

    end subroutine no_more_points

    !> See if more points are needed
    !!
    !! This subroutine is used by the root node to request a new live point
    function more_points_needed(mpi_communicator,root)
        implicit none
        integer, intent(in) :: mpi_communicator !> mpi handle
        integer, intent(in) :: root             !> root to recieve request

        integer :: empty_buffer(0) ! empty buffer to send
        integer, dimension(MPI_STATUS_SIZE) :: mpi_status  ! status identifier

        logical :: more_points_needed !> Whether we need more points or not

        call MPI_RECV(&
            empty_buffer,&
            0,&
            MPI_INT,&
            root,&
            MPI_ANY_TAG,&
            mpi_communicator,&
            mpi_status,&
            mpierror&
            )

        ! If we've recieved a kill signal, then exit this loop
        if(mpi_status(MPI_TAG) == tag_gen_stop ) then
            more_points_needed = .false.
        else if(mpi_status(MPI_TAG) == tag_gen_continue) then
            more_points_needed = .true.
        else
            call abort('generate error: unrecognised tag')
        end if




    end function more_points_needed


end module mpi_module
