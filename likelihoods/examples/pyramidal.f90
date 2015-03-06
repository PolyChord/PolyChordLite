module loglikelihood_module

    contains

    !> Pyramidal likelihood centered on 0.5.
    !!
    !! This is effectively a gaussian loglikelihood with an L_\infty norm
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    function loglikelihood(theta,phi)
        use utils_module, only: logTwoPi,Vn
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi

        double precision :: loglikelihood

        double precision, dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(size(theta)) :: mu    ! Mean

        mu= 5d-1   ! mean in the center
        sigma = 1d-1 

        ! normalisation
        loglikelihood =   - log(gamma(1d0+size(theta)/2d0)) -sum(log(2d0*sigma))

        ! theta dependence
        loglikelihood = loglikelihood - maxval(abs(theta-mu)/sigma)**2


    end function loglikelihood

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

    end subroutine setup_loglikelihood

end module loglikelihood_module
