module example_likelihoods
    use model_module , only: model

    contains

    function gaussian_loglikelihood(M,theta)
        class(model),intent(in)                    :: M
        double precision, intent(in), dimension(:) :: theta

        double precision :: gaussian_loglikelihood

        double precision, dimension(M%nDims) :: sigma  
        double precision, dimension(M%nDims) :: mu

        double precision, parameter :: TwoPi = 6.2831853d0
        
        sigma = 0.01
        mu    = 0.5

        if (size(theta) .ne. M%nDims) then
            write(*,*) 'error in log likelihood, size(theta) != nDims'
            return
        end if

        ! Gaussian normalisation
        gaussian_loglikelihood = - M%nDims / 2d0 * log( TwoPi ) - sum( log( sigma ) ) 

        ! theta dependence
        gaussian_loglikelihood = gaussian_loglikelihood &
            - sum( ( ( theta - mu ) / sigma ) ** 2d0 ) / 2d0

    end function gaussian_loglikelihood




end module example_likelihoods
