module mpi_module

    use mpi
    
    implicit none
    integer :: mpierror

    contains

    subroutine abort()
        implicit none
        integer :: i
        call MPI_ABORT(MPI_COMM_WORLD,i,mpierror)
    end subroutine abort

    subroutine mpi_initialise
        implicit none
        call MPI_INIT(mpierror)
    end subroutine mpi_initialise

    subroutine mpi_finalise
        implicit none
        call MPI_FINALIZE(mpierror)
    end subroutine mpi_finalise





end module mpi_module
