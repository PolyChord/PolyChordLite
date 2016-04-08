module loglikelihood_module
    use utils_module, only: dp

    interface 
        function cpp_loglikelihood(theta,nDims,phi,nDerived) bind(c)
            use iso_c_binding
            implicit none
            real (c_double),dimension(nDims),intent(in)     :: theta
            integer (c_int),intent(in)                      :: nDims
            real (c_double),dimension(nDerived),intent(out) :: phi
            integer (c_int),intent(in)                      :: nDerived
            real (c_double)                                 :: cpp_loglikelihood
        end function cpp_loglikelihood
    end interface

    interface 
        subroutine cpp_prior(cube,theta,nDims) bind(c)
            use iso_c_binding
            implicit none
            integer (c_int),intent(in)                      :: nDims
            real (c_double),dimension(nDims),intent(in)     :: cube
            real (c_double),dimension(nDims),intent(out)    :: theta
        end subroutine cpp_prior
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
        real(dp), intent(in), dimension(:)   :: theta
        !> Output derived parameters
        real(dp), intent(out),  dimension(:) :: phi

        real(dp) :: loglikelihood

        real (c_double),dimension(size(theta)) :: c_theta
        integer (c_int)                        :: nDims
        real (c_double),dimension(size(phi))   :: c_phi
        integer (c_int)                        :: nDerived

        real (c_double) :: c_loglike

        ! convert inputs to c
        c_theta   = theta

        ! Get the sizes of the arrays
        nDims     = size(theta)
        nDerived  = size(phi)

        ! Call the c loglikelihood function
        c_loglike = cpp_loglikelihood(c_theta,nDims,c_phi,nDerived)

        ! Convert outputs to fortran
        loglikelihood = c_loglike
        phi = c_phi

    end function loglikelihood


    function prior(cube) result(theta)
        use iso_c_binding
        implicit none
        !> Input parameters
        real(dp), intent(in), dimension(:)   :: cube
        !> Output derived parameters
        real(dp),  dimension(size(cube)) :: theta

        real(dp) :: loglikelihood

        real (c_double),dimension(size(cube))  :: c_cube
        real (c_double),dimension(size(theta)) :: c_theta
        integer (c_int)                        :: nDims

        ! convert inputs to c
        c_cube   = cube

        ! Get the sizes of the arrays
        nDims     = size(cube)

        ! Call the c loglikelihood function
        call cpp_prior(c_cube,c_theta,nDims)

        theta = c_theta

    end function prior

    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings

        call cpp_loglikelihood_setup

    end subroutine setup_loglikelihood

end module loglikelihood_module
