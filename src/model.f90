module model_module
    implicit none
    double precision, parameter :: logzero = -1d20

    !> Type to encode all of the information about the priors.
    type :: model
        !> Dimensionality of the space
        integer :: nDims       
        !> Number of derived parameters
        integer :: nDerived    
        !> 2*ndims + nDerived + 1
        integer :: nTotal      

        ! Indices for the sections of a live_points array
        !> @todo make array indices read only
        !> hypercube indices
        integer :: h0,h1       
        !> physical indices
        integer :: p0,p1       
        !> derived indices
        integer :: d0,d1       
        !> likelihood index
        integer :: l0          

        ! Prior details:
        ! Separable priors:
        !> Uniform prior indices
        integer                                        :: uniform_num   = 0
        double precision, allocatable, dimension(:,:)  :: uniform_params
        integer                                        :: uniform_index = 1

        !> log-uniform indices
        integer                                        :: log_uniform_num   = 0   
        double precision, allocatable, dimension(:,:)  :: log_uniform_params
        integer                                        :: log_uniform_index = 1


        procedure(loglike), pass(M), pointer :: loglikelihood 



    end type model


    interface
        function loglike(M,theta,feedback)
            import :: model
            class(model),     intent(in)                :: M
            double precision, intent(in),  dimension(:) :: theta
            integer,          intent(in),  optional     :: feedback

            double precision :: loglike
        end function
    end interface


    contains


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


    subroutine calculate_point(M, live_data)
        type(model),     intent(in)                   :: M
        double precision, intent(inout) , dimension(:) :: live_data

        if ( any(live_data(M%h0:M%h1)<0d0) .or. any(live_data(M%h0:M%h1)>1d0) )  then
            live_data(M%p0:M%p1) = 0
            live_data(M%l0) = logzero
        else
            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, live_data(:) )

            ! Calculate the likelihood and store it in the last index
            live_data(M%l0) = M%loglikelihood( live_data(M%p0:M%p1))

            ! Calculate the derived parameters
            call calculate_derived_parameters( M, live_data(:) )
        end if

    end subroutine calculate_point



    subroutine hypercube_to_physical(M, live_data)
        type(model),     intent(in)                   :: M
        double precision, intent(inout) , dimension(:) :: live_data

        double precision, dimension(M%nDims) :: hypercube_coords
        double precision, dimension(M%nDims) :: physical_coords

        ! copy the hypercube coordinates to a temporary variable
        hypercube_coords = live_data(M%h0:M%h1)

        ! Transform to physical coordinates

        ! Uniform priors:
        ! This is a fairly simple transformation, each parameter is transformed as
        ! hypercube_coord -> min + hypercube_coord * (max-min)
        ! Param 1 is the minimum bound, Param 2 is the maximum
        if(M%uniform_num>0) &
            physical_coords(M%uniform_index:) = &
            M%uniform_params(:,1) + (M%uniform_params(:,2) - M%uniform_params(:,1) ) * hypercube_coords(M%uniform_index:)

        ! log-uniform priors:
        ! hypercube_coord -> min * (max/min)**hypercube_coord
        ! Param 1 is the minimum bound, Param 2 is the maximum
        if(M%log_uniform_num>0) &
            physical_coords(M%log_uniform_index:) = &
            M%log_uniform_params(:,1) * (M%log_uniform_params(:,2)/M%log_uniform_params(:,1) ) ** hypercube_coords(M%log_uniform_index:)



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

        derived_parameters = live_data(M%d0:M%d1)

        ! accumulate the number of likelihood calls that we've made
        derived_parameters(1) = derived_parameters(1) + 1

        ! transfer the derived parameter back to live_data
        live_data(M%d0:M%d1) = derived_parameters

    end subroutine calculate_derived_parameters

    function prior_log_volume(M) result(log_volume)
        implicit none
        type(model),      intent(in)                   :: M

        double precision :: log_volume

        log_volume = 0

        ! Uniform contribution
        if (M%uniform_num >0 ) &
            log_volume = log_volume + sum(log(M%uniform_params(:,2)-M%uniform_params(:,1) ))
        if (M%log_uniform_num >0 ) &
            log_volume = log_volume + sum( log(log(M%log_uniform_params(:,2)/M%log_uniform_params(:,1))) )
        
    end function prior_log_volume




end module model_module 
