module settings_module
    use model_module, only: model
    implicit none

    Type, public :: program_settings

        integer :: nlive

        procedure(samp), nopass, pointer :: sampler

    end type program_settings

    interface
        subroutine samp(new_point, live_data, likelihood_bound, M) 
            import :: model
            implicit none
            double precision, intent(out),    dimension(:)   :: new_point
            double precision, intent(in), dimension(:,:) :: live_data
            double precision, intent(in) :: likelihood_bound
            type(model),            intent(in) :: M
        end subroutine samp
    end interface


end module settings_module 






module model_module
    implicit none

    !> Type to encode all of the information about the priors.
    type :: model


        integer :: nDims       !> Dimensionality of the space
        integer :: nDerived    !> Number of derived parameters
        integer :: nTotal      !> 2*ndims + nDerived + 1

        !> Indices for the sections of a live_points array
        integer :: h0,h1       !> hypercube indices
        integer :: p0,p1       !> physical indices
        integer :: d0,d1       !> derived indices
        integer :: l0          !> likelihood index

        procedure(loglike), pass(M), pointer :: loglikelihood 

    end type model

    interface
        function loglike(M,theta)
            import :: model
            class(model),     intent(in)               :: M
            double precision, intent(in), dimension(:) :: theta

            double precision :: loglike
        end function
    end interface


    contains

    subroutine hypercube_to_physical(M, live_data)
        type(model),     intent(in)                   :: M
        double precision, intent(inout) , dimension(:) :: live_data

        double precision, dimension(M%nDims) :: hypercube_coords
        double precision, dimension(M%nDims) :: physical_coords

        ! copy the hypercube coordinates to a temporary variable
        hypercube_coords = live_data(M%h0:M%h1)

        ! Transform to physical coordinates
        physical_coords = hypercube_coords

        ! copy the physical coordinates back to live_data
        live_data(M%p0:M%p1) = physical_coords

    end subroutine hypercube_to_physical



    subroutine calculate_derived_parameters(M, live_data)
        type(model),      intent(in)                   :: M
        double precision, intent(inout) , dimension(:) :: live_data

        double precision, dimension(M%nDims)    :: physical_coords
        double precision                        :: loglike
        double precision, dimension(M%nDerived) :: derived_parameters

        ! Copy the physical coordinates and loglike to temporary variables
        physical_coords    = live_data(M%p0:M%p1)
        loglike            = live_data(M%l0) 

        derived_parameters = 0.0

        ! transfer the derived parameter back to live_data
        live_data(M%d0:M%d1) = derived_parameters

    end subroutine calculate_derived_parameters



    subroutine initialise_model(M)
        type(model), intent(inout) :: M

        ! Total number of parameters
        M%nTotal = 2*M%nDims+M%nDerived+1

        ! Hypercube parameter indices
        M%h0=1
        M%h1=M%nDims

        ! Physical parameter indices
        M%p0=M%nDims+1
        M%p1=2*M%nDims

        ! Derived parameter indices
        M%d0=2*M%nDims+1
        M%d1=2*M%nDims+M%nDerived

        ! Loglikelihood index
        M%l0=M%nTotal


    end subroutine initialise_model



end module model_module 
