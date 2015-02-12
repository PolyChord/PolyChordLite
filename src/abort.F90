module abort_module
#ifdef MPI
    use mpi
#endif
    contains

    subroutine halt_program(message)
        use utils_module, only: stdout_unit
        implicit none

        character(LEN=*), intent(in), optional :: message

#ifdef MPI
        integer :: errorcode
        integer :: mpierror
#endif

        if (present(message)) then
            write(stdout_unit,'( 20("=") )')
            write(stdout_unit,'(A)') trim(adjustl(message))
            write(stdout_unit,'( 20("=") )')
        end if 


#ifdef MPI
        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
#else
        stop 1
#endif

    end subroutine halt_program


end module abort_module
