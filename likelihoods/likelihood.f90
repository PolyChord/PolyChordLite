module loglikelihood_module

    contains

    function loglikelihood(theta,phi,context)
        use utils_module, only: TwoPi
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        double precision :: loglikelihood
        
        loglikelihood = maxval(-abs(theta-0.5d0))

        phi=0d0

    end function loglikelihood

    subroutine setup_loglikelihood
        implicit none

    end subroutine setup_loglikelihood

end module loglikelihood_module
