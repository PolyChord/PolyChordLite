module mpi_module

    use mpi
    implicit none

    integer :: mpierror

    contains

    subroutine abort_all(message)
        use utils_module, only: stdout_unit
        implicit none

        character(LEN=*), intent(in), optional :: message

        integer :: errorcode
        integer :: mpierror

        if (present(message)) write(stdout_unit,trim(message))


        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)

    end subroutine abort_all




end module mpi_module
