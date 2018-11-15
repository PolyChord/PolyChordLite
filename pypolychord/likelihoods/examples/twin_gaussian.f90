module loglikelihood_module
    use utils_module, only: dp

    contains

    !> Twin gaussian peaks likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, apart from the x and y directions, where
    !! they are separated by 20sqrt(2) sigma widths and all sigmas at 0.01

    function loglikelihood(theta,phi)
        use utils_module, only: logaddexp,logTwoPi
        implicit none
        real(dp), intent(in),  dimension(:) :: theta         !> Input parameters
        real(dp), intent(out), dimension(:) :: phi           !> Output derived parameters
        real(dp)                            :: loglikelihood ! loglikelihood value to output

        real(dp) :: loglikelihood1
        real(dp) :: loglikelihood2

        real(dp), dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        real(dp), dimension(size(theta)) :: mu1   ! Mean
        real(dp), dimension(size(theta)) :: mu2   ! Mean


        ! Initialise the mean and standard deviation
        sigma  = 1d-1 ! all sigma set relatively small
        mu1    = 0d0  ! mean in the center
        mu1(1) = -5d-1
        mu1(2) = -5d-1
        mu2    = 0d0  ! mean in the center
        mu2(1) = +5d-1
        mu2(2) = +5d-1

        ! Gaussian normalisation
        loglikelihood1 = - sum( log( sigma ) + logTwoPi/2d0 ) 
        loglikelihood2 = loglikelihood1

        ! theta dependence
        loglikelihood1 = loglikelihood1 - sum( ( ( theta - mu1 ) / sigma ) ** 2d0 ) / 2d0
        loglikelihood2 = loglikelihood2 - sum( ( ( theta - mu2 ) / sigma ) ** 2d0 ) / 2d0

        loglikelihood = logaddexp(loglikelihood1,loglikelihood2) - log(2d0)

        ! One derived parameter, which is +1 if its in the right hand cluster,
        ! or -1 if its in the left hand cluster
        if(theta(1) >0.5d0) then
            phi(1)=1d0
        else
            phi(1)=-1d0
        end if

    end function loglikelihood

    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings

    end subroutine setup_loglikelihood

end module loglikelihood_module
