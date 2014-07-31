module example_likelihoods
    use model_module , only: model

    contains

    !> Basic Gaussian likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, and all sigmas at 0.01

    function gaussian_loglikelihood(M,theta)
        class(model),intent(in)                    :: M
        double precision, intent(in), dimension(:) :: theta


        double precision :: gaussian_loglikelihood

        double precision, dimension(M%nDims) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(M%nDims) :: mu    ! Mean

        double precision, parameter :: TwoPi = 8*atan(1d0) ! 2\pi in double precision
        
        ! Initialise the mean and standard deviation
        mu    = 0.5   ! mean in the center
        sigma = 0.01  ! all sigma set relatively small

        ! Gaussian normalisation
        gaussian_loglikelihood = - M%nDims / 2d0 * log( TwoPi ) - sum( log( sigma ) ) 

        ! theta dependence
        gaussian_loglikelihood = gaussian_loglikelihood - sum( ( ( theta - mu ) / sigma ) ** 2d0 ) / 2d0

    end function gaussian_loglikelihood




end module example_likelihoods
