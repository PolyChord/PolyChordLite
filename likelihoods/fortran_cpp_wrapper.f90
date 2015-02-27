module loglikelihood_module

    interface 
        function cpp_loglikelihood(theta,phi) bind(c)
            use iso_c_binding
            implicit none
            real (c_double),dimension(*),intent(in)  :: theta
            real (c_double),dimension(*),intent(out) :: phi
            real (c_double)                          :: cpp_loglikelihood
        end function cpp_loglikelihood
    end interface

    interface 
        subroutine cpp_loglikelihood_setup() bind(c)
            use iso_c_binding
            implicit none
        end subroutine cpp_loglikelihood_setup
    end interface

    contains

    function loglikelihood(theta,phi)
        use iso_c_binding
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi

        double precision :: loglikelihood

        real (c_double),dimension(size(theta)) :: c_theta
        real (c_double),dimension(size(phi))   :: c_phi

        real (c_double) :: c_loglike

        ! convert inputs to c
        c_theta   = theta

        ! Call the c loglikelihood function
        c_loglike = cpp_loglikelihood(c_theta,c_phi)

        ! Convert outputs to fortran
        loglikelihood = c_loglike
        phi = c_phi

    end function loglikelihood

    subroutine setup_loglikelihood
        implicit none

        call cpp_loglikelihood_setup

    end subroutine setup_loglikelihood

end module loglikelihood_module
