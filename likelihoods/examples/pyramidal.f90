module loglikelihood_module

    contains

    !> Pyramidal likelihood centered on 0.5.
    !!
    !! This is effectively a gaussian loglikelihood with an L_\infty norm
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! Note that there is a dimensional component to the argument of the exponential function
    !!
    !! This is chosen so as to ensure that the entropy H(N) ~ O(N)
    !!
    function loglikelihood(theta,phi)
        use utils_module, only: loggamma
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi

        integer,save :: nDims=0
        double precision, save :: factor=1d0

        double precision :: loglikelihood

        double precision, dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(size(theta)) :: mu    ! Mean

        double precision, parameter :: logSqrtTwoPi = log(sqrt(8d0*atan(1d0)))
        double precision, parameter :: PiOverTwo = 2d0*atan(1d0)



        mu= 5d-1   ! mean in the center
        sigma = 1d-1 
        if(size(theta)/=nDims) then
            nDims = size(theta)
            factor = exp(-2d0/nDims * loggamma(1d0 + nDims/2d0)) * PiOverTwo
        end if

        ! normalisation
        loglikelihood =   - sum(logSqrtTwoPi+log(sigma))

        ! theta dependence
        loglikelihood = loglikelihood - maxval(abs(theta-mu)/sigma)**2  / factor


    end function loglikelihood

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

    end subroutine setup_loglikelihood

end module loglikelihood_module
