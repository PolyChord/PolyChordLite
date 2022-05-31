module abort_module
    use utils_module, only: dp
    contains

    subroutine halt_program(message)
        use utils_module, only: stderr_unit
        implicit none
#ifdef USE_MPI
        include 'mpif.h'
#endif

        character(LEN=*), intent(in), optional :: message

#ifdef USE_MPI
        integer :: errorcode=1
        integer :: mpierror
#endif

        if (present(message)) then
            write(stderr_unit,'( 20("=") )')
            write(stderr_unit,'(A)') trim(adjustl(message))
            write(stderr_unit,'( 20("=") )')
        end if 

#ifdef USE_MPI
        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
#else
        stop 1
#endif

    end subroutine halt_program


end module abort_module
