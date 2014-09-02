module model_module
    use utils_module, only: logzero
    implicit none

    !> Type to encode all of the information about the priors.
    type :: model
        !> Dimensionality of the space
        integer :: nDims = 1
        !> Number of derived parameters
        integer :: nDerived = 0   
        !> 2*ndims + nDerived + 1
        integer :: nTotal      

        ! Indices for the sections of a live_points array

        !> hypercube indices
        !!
        !! These are the coordinates in the unit hypercube.
        !! h0 indicates the starting index, h1 the finishing index.
        integer :: h0,h1       

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

        ! Algorithm indices
        !> The number of likelihood evaluations required to calculate this point
        integer :: nlike
        !> The last chord length used in calculating this point
        integer :: last_chord
#ifdef MPI
        !> Pointer to any daughter points
        integer :: daughter
#endif

        !> likelihood index
        !!
        !! This is the likelihood evaluated at the position of the live point
        integer :: l0          
        !> This is the likelihood contour which generated the live point
        integer :: l1          

        !> Pointer to any additional files that need to be stored in between
        !! evaluations (only really important for C likelihoods)
        integer :: context


        !> Grades of parameters
        !!
        !! In the case of 'fast-slow' parameters, these indicate the 'grade' of
        !! parameter. grade=1 is slowest, 
        integer, dimension(:), allocatable :: grade


        !==========================================================================
        ! Prior details:
        !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Separable priors:
        !
        ! These are priors which are separable in each parameter {x_i}:
        ! pi(x_1,x_2,...) = pi_1(x_1) pi_2(x_2) ...
        !
        ! These may be subdivided into categories of prior. Currently supported
        ! are:
        ! 1) Uniform priors: 
        !     - parameters  a, b 
        !                       pi(x) = 1/(b-a)  if a<x<b
        !                               0        otherwise
        ! 2) Log-uniform priors
        !     - parameters  a, b 
        !                       pi(x) = log(a/b) * 1/x      if a<x<b
        !                               0                   otherwise
        ! 3) Gaussian priors
        !     - parameters  mu, sigma 
        !      pi(x) = 1/(sqrt(2pi)sigma) exp( -(x-mu)^2/(2*sigma^2) )
        !
        ! In general there will be <prior type>_num of each of these, and their
        ! parameters are stored in the 2D array <prior type>_params, which has
        ! shape (<prior type>_num, number of parameters)
        !
        ! These are stored sequentially in the hypercube array and physical
        ! parameters array, in between the indices: 
        !      <prior type>_index : <prior type>_index +<prior type>_num -1
        !
        ! If you wish to add a new separable prior type, follow the examples
        ! above. You need to:
        !   1) add <my type>_num, <my type>_params, <my type>_index as below
        !   2) add code to allocate the array in the subroutine initialise_model
        !   3) work out the transformation from the unit hypercube to the
        !      physical parameters. To do this, you should calculate the inverse
        !      function F^-1(x) of the cumulative distribution function of the
        !      prior pi(x):
        !       F(x) = int_{-inf}^x  pi(t) dt
        !   4) Implement this transformation in the subroutine hypercube_to_physical
        !   5) initialise the prior in main.f90 like so:
        !      M%<my type>_num = <number of paramers with this type>
        !      call initialise_model(M)
        !      M%<my_type>_params = <parameter choices>
        !
        !--------------------------------------------------------------------------
        ! Uniform priors:
        !
        !> The number of parameters with uniform priors
        integer                                        :: uniform_num   = 0
        !> The maximum and minimum of the uniform prior
        !!
        !! This is a two dimensional array of shape (uniform_num,2).
        !! The first row is the minimum, the second the maximum
        double precision, allocatable, dimension(:,:)  :: uniform_params
        !> The starting index of the coordinates of the uniform parameters (in
        !! both physical and hypercube coordinates)
        integer                                        :: uniform_index = 0
        !
        !--------------------------------------------------------------------------
        ! Log-uniform
        !
        !> The number of parameters with log-uniform priors
        integer                                        :: log_uniform_num   = 0   
        !> The maximum and minimum of the log-uniform priors
        !!
        !! This is a two dimensional array of shape (sorted_uniform_num,2).
        !! The first row is the minimum, the second the maximum
        double precision, allocatable, dimension(:,:)  :: log_uniform_params
        !> The starting index of the coordinates of the log-uniform parameters (in
        !! both physical and hypercube coordinates)
        integer                                        :: log_uniform_index = 0
        !
        !--------------------------------------------------------------------------
        ! Gaussian 
        !
        !> The number of parameters with Gaussian priors
        integer                                        :: gaussian_num   = 0   
        !> The mean and standard deviation of the Gaussian priors
        !!
        !! This is a two dimensional array of shape (gaussian_num,2).
        !! The first row is the mean, the second the standard deviation
        double precision, allocatable, dimension(:,:)  :: gaussian_params
        !> The starting index of the coordinates of the Gaussian parameters (in
        !! both physical and hypercube coordinates)
        integer                                        :: gaussian_index = 0
        !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Non separable priors
        !
        ! These priors cannot be separated out, and must be transformed together
        ! using joint conditionals.
        !
        ! These may be sub-divided into categories of multivariate prior.
        ! Currently supported are:
        !
        ! 1) Sorted uniform priors
        !     - parameters  x_min,x_max,N
        !     
        !    These are a set of variables {x_i} distributed uniformly in
        !    [x_min,x_max] such that x_1<x_2<...<x_N
        !    We transform using the distributions which are marginalised over
        !    the greater parameters and conditioned on the lower parameters:
        !
        !    pi(x_1)                 = N/(x_N-x_min) * x_1^(N-1)     xmin<x<xmax
        !                              0                             otherwise
        !    pi(x_{i+1}|x_i,...,x_1) = 1/(x_N-x_i)   * x_1^(N-i-1)   x_i <x<xmax
        !                              0                             otherwise
        !--------------------------------------------------------------------------
        ! Sorted uniform priors
        !
        !> The number of sets of sorted uniform priors.
        integer                                        :: sorted_uniform_num   = 0   
        !
        !> The maximum and minimum of the sorted uniform priors
        !! This is a two dimensional array of shape (sorted_uniform_num,3).
        !!
        !! The first row is the minimum, the second the maximum, and the third
        !! is the number of parameters in the set.
        !! columns indicate which set of parameters these belong to.
        double precision, allocatable, dimension(:,:)  :: sorted_uniform_params
        !> The starting index of the coordinates of the sorted uniform parameters (in
        !! both physical and hypercube coordinates)
        integer, allocatable, dimension(:)             :: sorted_uniform_index
        !
        !==========================================================================


    end type model


    contains


    subroutine allocate_live_indices(M)
        implicit none
        type(model), intent(inout) :: M


        ! Hypercube parameter indices
        M%h0=1
        M%h1=M%h0+M%nDims-1

        ! Physical parameter indices
        M%p0=M%h1+1
        M%p1=M%p0+M%nDims-1 

        ! Derived parameter indices
        M%d0=M%p1+1
        M%d1=M%d0+M%nDerived-1  

        ! Algorithm indices
        M%nlike=M%d1+1
        M%last_chord=M%nlike+1
#ifdef MPI
        M%daughter=M%last_chord+1
        ! Loglikelihood indices
        M%l0=M%daughter+1
#else
        M%l0=M%last_chord+1
#endif

        M%l1=M%l0+1

        ! Total number of parameters
        M%nTotal = M%l1

        ! grades
        if(.not. allocated(M%grade)) allocate(M%grade(M%nDims))

    end subroutine allocate_live_indices

    subroutine allocate_prior_arrays(M)
        implicit none
        type(model), intent(inout) :: M

        ! Allocate the prior arrays if they're required, and not already
        ! allocated
        if(M%uniform_num>0 .and. .not. allocated(M%uniform_params) ) then
            allocate( M%uniform_params(M%uniform_num,2) )
            M%uniform_params(:,1) = 0d0
            M%uniform_params(:,2) = 1d0
        end if
        if(M%log_uniform_num>0 .and. .not. allocated(M%log_uniform_params) ) then
            allocate( M%log_uniform_params(M%log_uniform_num,2) )
            M%log_uniform_params(:,1) = exp(0d0)
            M%log_uniform_params(:,2) = exp(1d0)
        end if
        if(M%gaussian_num>0 .and. .not. allocated(M%gaussian_params) ) then
            allocate( M%gaussian_params(M%gaussian_num,2) )
            M%log_uniform_params(:,1) = 0d0
            M%log_uniform_params(:,2) = 1d0
        end if

        if(M%sorted_uniform_num>0) then
            if (.not. allocated(M%sorted_uniform_params) ) allocate( M%sorted_uniform_params(M%sorted_uniform_num,3) )
            if (.not. allocated(M%sorted_uniform_index)  ) allocate( M%sorted_uniform_index(M%sorted_uniform_num) )
            M%sorted_uniform_params(:,1) = 0d0
            M%sorted_uniform_params(:,2) = 1d0
            M%sorted_uniform_params(:,3) = nint(0d0)
        end if
    end subroutine allocate_prior_arrays
        
        
    subroutine set_up_prior_indices(M)
        implicit none
        type(model), intent(inout) :: M
        integer :: i

        ! Set up the indices
        M%uniform_index        = 1
        M%log_uniform_index    = M%uniform_index     + M%uniform_num
        M%gaussian_index       = M%log_uniform_index + M%log_uniform_num

        do i=1,M%sorted_uniform_num
            M%sorted_uniform_index(i) = M%gaussian_index  + sum(nint(M%sorted_uniform_params(:i,3)))
        end do

    end subroutine set_up_prior_indices



    subroutine calculate_point(loglikelihood, M, live_data)
        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        type(model),     intent(in)                   :: M
        double precision, intent(inout) , dimension(:) :: live_data

        if ( any(live_data(M%h0:M%h1)<0d0) .or. any(live_data(M%h0:M%h1)>1d0) )  then
            live_data(M%p0:M%p1) = 0
            live_data(M%l0) = logzero
        else
            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, live_data(:) )

            ! Calculate the likelihood and store it in the last index
            live_data(M%l0) = loglikelihood( live_data(M%p0:M%p1), live_data(M%d0:M%d1),M%context)

            ! accumulate the number of likelihood calls that we've made
            live_data(M%nlike) = live_data(M%nlike)+1
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

        ! Sorted uniform priors:
        ! 



        ! copy the physical coordinates back to live_data
        live_data(M%p0:M%p1) = physical_coords

    end subroutine hypercube_to_physical




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


    function gradloglike(loglikelihood,M,theta,loglike,delta)
        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface
        type(model),      intent(in)                     :: M
        double precision, intent(in), dimension(M%nDims) :: theta
        double precision, intent(in)                     :: loglike
        double precision, intent(in)                     :: delta

        double precision, dimension(M%nDims)             :: gradloglike


        double precision, dimension(M%nDerived)               :: derived
        integer :: context

        double precision, dimension(M%nDims) :: delta_vec

        integer i

        do i=1,M%nDims
            delta_vec = 0
            delta_vec(i) = delta
            gradloglike(i) = ( loglikelihood(theta+delta_vec,derived,context) - loglike )/ delta
        end do

    end function gradloglike




end module model_module 
