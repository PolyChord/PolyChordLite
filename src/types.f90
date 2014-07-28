module settings_module
    implicit none

    Type, public :: program_settings

        integer :: nlive

    end type program_settings


end module settings_module 


module model_module
    implicit none

    !> Type to encode all of the information about the priors.
    type :: model

        
        integer :: nDims       !> Dimensionality of the space

        integer :: nDerived    !> Number of derived parameters

        integer :: nTotal      !> nDerived + nDims + 1

        contains
        procedure :: loglikelihood => gaussian_loglikelihood

    end type model

    contains

    function gaussian_loglikelihood(this,theta)
        class(model),intent(in)                             :: this
        double precision, intent(in), dimension(this%nDims) :: theta

        double precision :: gaussian_loglikelihood

        double precision, dimension(this%nDims) :: sigma  
        double precision, dimension(this%nDims) :: mu

        double precision :: TwoPi 

        integer :: i


        TwoPi = 6.2831853d0
        sigma = 0.01
        mu    = 0.5

        gaussian_loglikelihood = - this%nDims / 2d0 * log( TwoPi )

        do i = 1, this%nDims
            gaussian_loglikelihood = gaussian_loglikelihood - log( sigma(i) )
        enddo

        gaussian_loglikelihood = gaussian_loglikelihood &
            - sum( ( ( theta( 1:this%nDims ) - mu(1:this%nDims ) ) / sigma( 1:this%nDims ) ) ** 2d0 ) / 2d0

    end function gaussian_loglikelihood

end module model_module 
