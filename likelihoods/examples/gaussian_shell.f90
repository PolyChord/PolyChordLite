module loglikelihood_module
    use utils_module, only: dp
    use utils_module, only: logzero

    real(dp), parameter :: sigma  = 0.1d0   ! all sigma set relatively small
    real(dp), parameter :: radius = 2d0

    real(dp) :: A = logzero

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

        if(A==logzero) then
            r0 = (radius + sqrt(radius**2 + 4*(ndims-1)*sigma**2))/2
            logf0 = - (radius - r0)**2/2/sigma**2 + (ndims-1)*log(r0) + log(ndims+0d0) + ndims/2d0*logpi -loggamma(1+ndims/2d0) 
            sigma0 = sigma*sqrt((1+radius/sqrt(radius**2 + 4*(ndims-1)*sigma**2))/2d0)
            A = logf0 + logsqrttwopi + log(sigma0)
        end if


        mu=0
        ! Gaussian normalisation
        loglikelihood = -A - ( (sqrt( sum( (mu-theta)**2 ) ) -radius)**2d0 /(2d0*sigma*sigma) )

        phi(1) = sqrt(sum(theta**2)) ! Radius
        do i=2,size(phi)
            phi(i) = acos(theta(i-1)/sqrt(sum(theta(i-1:)**2)))
        end do

    end function loglikelihood

    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings

    end subroutine setup_loglikelihood


end module loglikelihood_module
