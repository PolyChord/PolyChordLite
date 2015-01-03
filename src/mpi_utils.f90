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




end module mpi_module
