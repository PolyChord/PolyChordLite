module settings_module
    implicit none

    Type, public :: program_settings

        integer :: nlive

    end type program_settings


end module settings_module 


module model_module
    implicit none

    !> Type to encode all of the information about the priors.
    Type, public :: model_details

        
        integer :: nDims       !> Dimensionality of the space

        integer :: num_derived !> Number of derived parameters

    end type model_details


end module model_module 
