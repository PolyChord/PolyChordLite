module example_likelihoods
    use model_module , only: model,logzero

    contains

    !> Basic Gaussian likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, and all sigmas at 0.01

    function gaussian_loglikelihood(M,theta,feedback)
        implicit none
        class(model),     intent(in)               :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: gaussian_loglikelihood

        double precision, dimension(M%nDims) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(M%nDims) :: mu    ! Mean

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision

        
        ! Initialise the mean and standard deviation
        mu    = 5d-1   ! mean in the center
        sigma = 1d-2  ! all sigma set relatively small

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : Gaussian" )')
            end if
            if(feedback>=2) then
                write(*,'( "     mean: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') mu
                write(*,'( "     sigma: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') sigma
            end if
            gaussian_loglikelihood = logzero
            return
        end if


        ! Gaussian normalisation
        gaussian_loglikelihood = - M%nDims/2d0 * log( TwoPi ) - sum( log( sigma ) ) 

        ! theta dependence
        gaussian_loglikelihood = gaussian_loglikelihood - sum( ( ( theta - mu ) / sigma ) ** 2d0 ) / 2d0


    end function gaussian_loglikelihood





    !> Pyramidal likelihood centered on 0.5.
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function pyramidal_loglikelihood(M,theta,feedback)
        implicit none
        class(model),intent(in)                    :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: pyramidal_loglikelihood

        double precision, dimension(M%nDims) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(M%nDims) :: center    ! Mean
        
        center    = 5d-1   ! mean in the center
        sigma = 1d-2 

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : Pyramidal" )')
            end if
            if(feedback>=2) then
                write(*,'( "     center: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') center
            end if
            return
        end if



        ! normalisation
        pyramidal_loglikelihood =   -(M%nDims)*log(2d0) - log(gamma(1d0+M%nDims/2d0)) -sum(log(sigma))

        ! theta dependence
        pyramidal_loglikelihood = pyramidal_loglikelihood - maxval(abs(theta-center)/sigma)**2

    end function pyramidal_loglikelihood



end module example_likelihoods
