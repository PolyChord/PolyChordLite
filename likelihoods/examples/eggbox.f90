module loglikelihood_module
    use utils_module, only: dp

    contains

    !> Eggbox likelihood
    !!
    !! \f[ -\log L(x, y) = (2+ \prod_i^n \cos(\theta_i/2) )^5 \f]
    !!
    function loglikelihood(theta,phi)
        implicit none
        real(dp), intent(in),  dimension(:) :: theta         !> Input parameters
        real(dp), intent(out), dimension(:) :: phi           !> Output derived parameters
        real(dp)                            :: loglikelihood ! loglikelihood value to output

        ! No normalisation implemented yet
        loglikelihood = 0

        ! calculate
        loglikelihood =  loglikelihood  -  (2 + product(cos(theta / 2d0) ) )**5

    end function loglikelihood

    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings

    end subroutine setup_loglikelihood

end module loglikelihood_module
