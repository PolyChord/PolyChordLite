module loglikelihood_module

    double precision :: normalisation

    ! Parameters of the rosenbrock function
    double precision, parameter :: a = 1d0
    double precision, parameter :: b = 1d2
    double precision, parameter :: pi = 4d0*atan(1d0) ! \pi in double precision

    ! number of dimensions
    integer :: nDims

    contains

    !> Upside down [Rosenbrock function](http://en.wikipedia.org/wiki/Rosenbrock_function).
    !!
    !! \f[ -\log\mathcal{L}(\theta) = \sum_{i=1}^{N-1}   (a-\theta_i)^2+ b (\theta_{i+1} -\theta_i^2 )^2 \f]
    !!
    !! This is the industry standard 'Banana'. It is useful for testing whether the
    !! algorithm is capable of navigating a curving degenerate space and able to
    !! find the global maximum. As is conventional, we choose \f$a=1\f$ and
    !! \f$b=100\f$
    !!
    !! It is only valid in nDims>2.
    !!
    !! To be precise, we choose our log likelihood as the negative of the 'true'
    !! Rosenbrock function, as our algorithm is a maximiser.
    !!
    !! In addition, we add an offset to normalise the loglikelihood so that in
    !! 2D it should integrate to 1. (There is no analytic formula for ND).
    !! 
    !! The global maximum is atop a long, narrow, parabolic shaped ridge.
    !!
    !! dimension | extrema
    !! ----------|------
    !! \f$ 2 \f$ | One maximum at \f$(1,1)\f$
    !! \f$ 3 \f$ | One maximum at \f$(1,1,1)\f$ 
    !! \f$4-7\f$ | One global maximum at \f$(1,\ldots,1)\f$, one local near \f$(-1,1,\ldots,1)\f$
    !!
    !! \f$(1,\ldots,1)\f$ is always the
    !! global maximum, but it is unclear analytically how many (if any) local maxima there
    !! are elsewhere in dimensions higher than 7.
    !!
    function loglikelihood(theta,phi)
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision                            :: loglikelihood ! loglikelihood value to output


        ! Normalisation for 2D
        loglikelihood = normalisation

        ! Sum expressed with fortran intrinsics
        loglikelihood =  loglikelihood  - sum( (a-theta(1:nDims-1))**2 + b*(theta(2:nDims) - theta(1:nDims-1)**2)**2 )

    end function loglikelihood

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

        ! Get the number of dimensions from the settings variable
        nDims = settings%nDims

        ! 
        normalisation = -0.5d0 * log( pi**nDims / det(nDims) )

    end subroutine setup_loglikelihood

    function det(n)
        integer,intent(in) :: n
        double precision   :: det

        det = abs(-2d0*b*recur(n-1) - 16*b*b*recur(n-2))

    end function

    recursive function recur(n) result(ans)
        integer,intent(in) :: n
        double precision   :: ans

        if(n<=0) then
            ans = 0d0
            return
        else if(n==1) then
            ans = 1d0
            return
        else
            ans = (-2d0-10d0*b) * recur(n-1) - 16d0*b*b * recur(n-2)
        end if

    end function recur

end module loglikelihood_module
