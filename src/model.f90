!include 'mkl_vml.f90'   ! include the vector mathematical library, in order to
                        ! use the vector functions for the priors


module model_module
    use utils_module, only: logzero
    implicit none

    !> Type to encode all of the information about the priors.
    type :: model
        !> Dimensionality of the space
        integer :: nDims       
        !> Number of derived parameters
        integer :: nDerived    
        !> 2*ndims + nDerived + 1
        integer :: nTotal      

        ! Indices for the sections of a live_points array

        !> hypercube indices
        !!
        !! These are the coordinates in the unit hypercube.
        !! h0 indicates the starting index, h1 the finishing index.
        integer :: h0,h1       

        !> algorithm indices
        !!
        !! These are any variables attached to each live point, e.g. the number
        !! of likelihood calls, or the estimated slice size.
        !! a0 indicates the starting index, a1 the finishing index.
        integer :: a0,a1       

        !> physical indices
        !!
        !! These are the coordinates in the physical space 
        !! These are linked to the hypercube coordinates via the prior.
        !! p0 indicates the starting index, p1 the finishing index
        integer :: p0,p1       

        !> derived indices
        !!
        !! These are the indices of the derived parameters to be calculated
        !! during the run.
        !! d0 indicates the starting index, d1 the finishing index
        integer :: d0,d1       

        !> likelihood index
        !!
        !! This is the likelihood attached to each live point
        integer :: l0          

        ! Prior details:
        !
        ! Separable priors:
        !
        ! Uniform priors:
        !> The number of parameters with uniform priors
        integer                                        :: uniform_num   = 0
        !> The maximum and minimum of the uniform prior
        !!
        !! This is a two dimensional array of shape (uniform_num,2).
        !! The first row is the minimum, the second the maximum
        double precision, allocatable, dimension(:,:)  :: uniform_params
        !> The starting index of the coordinates of the uniform parameters (in
        !! both physical and hypercube coordinates)
        integer                                        :: uniform_index = 1

        ! Log-uniform
        !> The number of parameters with log-uniform priors
        integer                                        :: log_uniform_num   = 0   
        !> The maximum and minimum of the log-uniform priors
        !!
        !! This is a two dimensional array of shape (uniform_num,2).
        !! The first row is the minimum, the second the maximum
        double precision, allocatable, dimension(:,:)  :: log_uniform_params
        !> The starting index of the coordinates of the log-uniform parameters (in
        !! both physical and hypercube coordinates)
        integer                                        :: log_uniform_index = 1

        ! Gaussian 
        !> The number of parameters with Gaussian priors
        integer                                        :: gaussian_num   = 0   
        !> The mean and standard deviation of the Gaussian priors
        !!
        !! This is a two dimensional array of shape (uniform_num,2).
        !! The first row is the mean, the second the standard deviation
        double precision, allocatable, dimension(:,:)  :: gaussian_params
        !> The starting index of the coordinates of the Gaussian parameters (in
        !! both physical and hypercube coordinates)
        integer                                        :: gaussian_index = 1


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
        if(M%uniform_num>0) then
            physical_coords(M%uniform_index:) = M%uniform_params(:,1) &
                + (M%uniform_params(:,2) - M%uniform_params(:,1) ) * hypercube_coords(M%uniform_index:)
        end if

        ! log-uniform priors:
        ! hypercube_coord -> min * (max/min)**hypercube_coord
        ! Param 1 is the minimum bound, Param 2 is the maximum
        if(M%log_uniform_num>0) then
            physical_coords(M%log_uniform_index:) = M%log_uniform_params(:,1) &
                * (M%log_uniform_params(:,2)/M%log_uniform_params(:,1) ) ** hypercube_coords(M%log_uniform_index:)
        end if

        ! Gaussian priors:
        ! We use the intel function vcdfnorminv to transform from [0,1]->[-inf,inf] via a variant of the error function
        ! We then scale by the gaussian standard deviations in gaussian_params(:,2) and shift by the mean in gaussian_params(:,1).
        if(M%gaussian_num>0) then

            ! Transform via the intel function cdfnorm inverse (v=vector, d=double)
            ! https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-67369FA5-ABFD-4B5D-82D4-E6A5E4AB565B.htm#GUID-67369FA5-ABFD-4B5D-82D4-E6A5E4AB565B
            call vdcdfnorminv(M%gaussian_num,hypercube_coords(M%gaussian_index:),physical_coords(M%gaussian_index:) )

            ! Scale by the standard deviation and shift by the mean
            physical_coords(M%gaussian_index:) = M%gaussian_params(:,1) + M%gaussian_params(:,2) * physical_coords(M%gaussian_index:)

        end if



        ! copy the physical coordinates back to live_data
        live_data(M%p0:M%p1) = physical_coords

    end subroutine hypercube_to_physical



    subroutine calculate_derived_parameters(M, live_data)
        implicit none
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
        use utils_module, only: TwoPi
        implicit none
        type(model),      intent(in)                   :: M

        double precision :: log_volume

        log_volume = 0

        ! Uniform contribution
        if (M%uniform_num >0 ) &
            log_volume = log_volume + sum(log(M%uniform_params(:,2)-M%uniform_params(:,1) ))
        if (M%log_uniform_num >0 ) &
            log_volume = log_volume + sum( log(log(M%log_uniform_params(:,2)/M%log_uniform_params(:,1))) )
        if (M%gaussian_num >0 ) &
            log_volume = log_volume + sum( 0.5d0*log(TwoPi) + log(M%gaussian_params(:,2)) )
        
    end function prior_log_volume




end module model_module 
