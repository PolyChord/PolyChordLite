module loglikelihood_module

    contains

    function loglikelihood(theta,phi)
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision                            :: loglikelihood ! loglikelihood value to output




    end function loglikelihood

    subroutine setup_loglikelihood
        implicit none

    end subroutine setup_loglikelihood

end module loglikelihood_module
