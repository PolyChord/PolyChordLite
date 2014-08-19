module mpi_module

    use mpi
    
    implicit none
    integer :: mpierror

    contains



    function mpi_rank
        implicit none
        integer :: mpi_rank
        call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpierror)
    end function mpi_rank

    function mpi_size 
        integer :: mpi_size
        call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpierror)
    end function mpi_size

    subroutine mpi_synchronise
        implicit none
        call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
    end subroutine mpi_synchronise

    subroutine mpi_initialise
        implicit none
        call MPI_INIT(mpierror)
    end subroutine mpi_initialise

    subroutine mpi_finalise
        implicit none
        call MPI_FINALIZE(mpierror)
    end subroutine mpi_finalise





end module mpi_module
