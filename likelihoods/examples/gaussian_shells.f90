module loglikelihood_module
    use utils_module, only: dp

    real(dp), parameter :: sigma  = 0.1d0   ! all sigma set relatively small
    real(dp), parameter :: radius = 2d0

    real(dp) :: A = -huge(1d0)

    contains

    function loglikelihood(theta,phi)
        use utils_module, only: logTwoPi,Vn,loggamma,logincexp,Hypergeometric1F1
        implicit none
        !> Input parameters
        real(dp), intent(in), dimension(:)   :: theta
        !> Output derived parameters
        real(dp), intent(out),  dimension(:) :: phi

        real(dp) :: loglikelihood
        real(dp) :: loglikelihood_temp

        real(dp), dimension(size(theta)) :: mu    ! Mean

        real(dp),parameter :: logsqrttwopi = log(sqrt(8d0*atan(1d0)))
        real(dp),parameter :: logpi = log(4d0*atan(1d0))

        integer :: ndims
        integer :: i

        real(dp) :: r0,sigma0,logf0

        ndims  = size(theta)

        if(A==-huge(1d0)) then
            r0 = (radius + sqrt(radius**2 + 4*(ndims-1)*sigma**2))/2
            logf0 = - (radius - r0)**2/2/sigma**2 + (ndims-1)*log(r0) + log(ndims+0d0) + ndims/2d0*logpi -loggamma(1+ndims/2d0) 
            sigma0 = sigma*sqrt((1+radius/sqrt(radius**2 + 4*(ndims-1)*sigma**2))/2d0)
            A = logf0 + logsqrttwopi + log(sigma0)
        end if


        mu=0
        ! Initialise the mean and standard deviation
        mu(1)  = -3.5   ! mean in the center

        ! Gaussian normalisation
        loglikelihood_temp = -A - ( (sqrt( sum( (mu-theta)**2 ) ) -radius)**2d0 /(2d0*sigma*sigma) )

        loglikelihood = loglikelihood_temp

        mu(1)  = +3.5   ! mean in the center

        loglikelihood_temp = -A - ( (sqrt( sum( (mu-theta)**2 ) ) -radius)**2d0 /(2d0*sigma*sigma) )

        call logincexp(loglikelihood,loglikelihood_temp)
        loglikelihood = loglikelihood - log(2d0)

    end function loglikelihood

    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings

    end subroutine setup_loglikelihood


end module loglikelihood_module
