module loglikelihood_module
    use utils_module, only: logzero

    double precision, parameter :: sigma  = 0.1d0   ! all sigma set relatively small
    double precision, parameter :: radius = 2d0

    double precision :: A = logzero

    contains

    function loglikelihood(theta,phi)
        use utils_module, only: logTwoPi,Vn,loggamma,logincexp,Hypergeometric1F1
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi

        double precision :: loglikelihood
        double precision :: loglikelihood_temp

        double precision, dimension(size(theta)) :: mu    ! Mean

        double precision,parameter :: logsqrttwopi = log(sqrt(8d0*atan(1d0)))
        double precision,parameter :: logpi = log(4d0*atan(1d0))

        integer :: ndims
        integer :: i

        double precision :: r0,sigma0,logf0

        ndims  = size(theta)

        if(A==logzero) then
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

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

    end subroutine setup_loglikelihood


end module loglikelihood_module
