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


end module mpi_module
