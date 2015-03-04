module loglikelihood_module
    !> This module is where your likelihood code should be placed.
    !!
    !! * The loglikelihood is called by the subroutine loglikelihood.
    !! * The likelihood is set up by setup_loglikelihood.
    !! * You can store any global/saved variables in the module.


    !============================================================
    ! insert likelihood variables here
    !
    !
    !============================================================

    contains

    !> Main loglikelihood function
    !!
    !! Either write your likelihood code directly into this function, or call an
    !! external library from it. This should return the logarithm of the
    !! likelihood, i.e:
    !!
    !! loglikelihood = log_e ( P (data | parameters, model ) )
    !!
    !! theta are the values of the input parameters (in the physical space 
    !! NB: not the hypercube space).
    !!
    !! phi are any derived parameters that you would like to save with your
    !! likelihood.
    !!
    !! the variable loglikelihood should contain the return value.
    !!
    !! If you would like to see some example likelihood functions, look in:
    !! likelihoods/examples/
    !!
    function loglikelihood(theta,phi)
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision                            :: loglikelihood ! loglikelihood value to output

        !============================================================
        ! insert likelihood code here
        !
        !
        !============================================================



    end function loglikelihood



    !> Setup of the loglikelihood
    !!
    !! This is called before nested sampling, but after the priors and settings
    !! have been set up.
    !!
    !! This is the time at which you should load any files that the likelihoods
    !! need, and do any initial calculations.
    !!
    !! This module can be used to save variables in between calls
    !! (at the top of the file).
    !!
    !! All MPI threads will call this function simultaneously, but you may need
    !! to use mpi utilities to synchronise them. This should be done through the
    !! integer mpi_communicator (which is normally MPI_COMM_WORLD).
    !!
    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

        !============================================================
        ! insert likelihood setup here
        !
        !
        !============================================================

    end subroutine setup_loglikelihood

end module loglikelihood_module
