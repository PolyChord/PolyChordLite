module loglikelihood_module

    contains

    !> Basic Gaussian likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, and all sigmas at 0.01
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

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

    end subroutine setup_loglikelihood

end module loglikelihood_module
