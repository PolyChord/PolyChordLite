module loglikelihood_module
    use utils_module, only: dp

    contains

    !> Basic Gaussian likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, and all sigmas at 0.1
    function loglikelihood(theta,phi)
        use utils_module, only: logTwoPi,Vn
        implicit none
        !> Input parameters
        real(dp), intent(in), dimension(:)   :: theta
        !> Output derived parameters
        real(dp), intent(out),  dimension(:) :: phi

        real(dp) :: loglikelihood

        real(dp), dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        real(dp), dimension(size(theta)) :: mu    ! Mean


        ! Initialise the mean and standard deviation
        mu    = 5d-1   ! mean in the center
        sigma = 1d-1  ! all sigma set relatively small

        ! Gaussian normalisation
        loglikelihood = - sum( log( sigma ) + logTwoPi/2d0 ) 

        ! theta dependence
        loglikelihood = loglikelihood - sum( ( ( theta - mu ) / sigma ) ** 2d0 ) / 2d0

        ! The radius
        phi(1) = sqrt(sum((theta-mu)**2))
        phi(2) = log(phi(1)**size(theta) * Vn(size(theta)))


    end function loglikelihood

    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings

    end subroutine setup_loglikelihood

end module loglikelihood_module
