module loglikelihood_module

    contains
    !> Upside down [Rastrigin function](http://en.wikipedia.org/wiki/Rastrigin_function).
    !!
    !! \f[ -\log L(\theta) = f(\mathbf{x}) = A n + \sum_{i=1}^n \left[\theta_i^2 - A\cos(2 \pi \theta_i)\right] \f]
    !!
    !! where  \f$ A=10\f$  and \f$ \theta_i\in[-5.12,5.12] \f$ . It has a global minimum 
    !! at \f$ \theta = \mathbf{0}\f$  where \f$ f(\theta)=0\f$ .
    !!
    !! In mathematical optimization, the Rastrigin function is a
    !! non-convex function used as a performance test problem for optimization
    !! algorithms. It is a typical example of non-linear multimodal function. It
    !! was first proposed by Rastrigin as a 2-dimensional function and has been generalized
    !! by Mühlenbein et al. Finding the minimum of this function is a
    !! fairly difficult problem due to its large search space and its large number
    !! of local minima.
    !!
    function loglikelihood(theta,phi)
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi

        double precision, parameter :: A=10d0
        double precision, parameter :: TwoPi = 8d0*atan(1d0)

        ! The return value
        double precision :: loglikelihood

        loglikelihood =  - sum( log(4991.21750d0) + theta**2 - A*cos(TwoPi*theta) )

    end function loglikelihood

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

    end subroutine setup_loglikelihood

end module loglikelihood_module
