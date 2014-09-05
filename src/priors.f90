module priors_module
    use utils_module, only: stdout_unit
    implicit none

    integer, parameter :: uniform_type        = 1
    integer, parameter :: gaussian_type       = 2
    integer, parameter :: log_uniform_type    = 3
    integer, parameter :: sorted_uniform_type = 4

    type prior
        integer :: npars = 0
        integer :: prior_type = 1
        integer, dimension(:), allocatable :: hypercube_indices
        integer, dimension(:), allocatable :: physical_indices

        double precision, dimension(:), allocatable :: parameters

    end type prior



    contains

    ! Prior transformations

    function hypercube_to_physical(hypercube_coords,priors) result(physical_coords)
        implicit none
        type(prior), dimension(:), intent(in) :: priors
        double precision, intent(in), dimension(:) :: hypercube_coords

        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: i

        physical_coords=0d0

        do i=1,size(priors)
            select case(priors(i)%prior_type)
            case(uniform_type)
                physical_coords(priors(i)%physical_indices)= uniform_htp    (hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            case(gaussian_type)
                physical_coords(priors(i)%physical_indices)= gaussian_htp   (hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            case(log_uniform_type)
                physical_coords(priors(i)%physical_indices)= log_uniform_htp(hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            end select
        end do

    end function hypercube_to_physical


    function prior_log_volume(priors) result(log_volume)
        use utils_module, only: TwoPi
        implicit none
        type(prior), dimension(:), intent(in) :: priors

        double precision :: log_volume
        integer :: i

        log_volume = 0

        do i=1,size(priors)
            select case(priors(i)%prior_type)
            case(uniform_type)
                log_volume = log_volume + sum( log(priors(i)%parameters(priors(i)%npars+1:)- priors(i)%parameters(:priors(i)%npars) )) 
            case(gaussian_type)
                log_volume = log_volume + sum( 0.5d0*log(TwoPi) + log(priors(i)%parameters(priors(i)%npars+1:)) )
            case(log_uniform_type)
                log_volume = log_volume + sum( log(log( priors(i)%parameters(priors(i)%npars+1:)/priors(i)%parameters(:priors(i)%npars) )) ) 
            end select
        end do

        ! Uniform contribution
        !if (M%uniform_num >0 ) &
        !    log_volume = log_volume + sum(log(M%uniform_params(:,2)-M%uniform_params(:,1) ))
        !if (M%log_uniform_num >0 ) &
        !    log_volume = log_volume + sum( log(log(M%log_uniform_params(:,2)/M%log_uniform_params(:,1))) )
        !if (M%gaussian_num >0 ) &
        !    log_volume = log_volume + sum( 0.5d0*log(TwoPi) + log(M%gaussian_params(:,2)) )
        
    end function prior_log_volume

    !> Routine to initialise a set of seperable uniform parameters
    !!
    !! The user should provide:
    !!
    !! minimums & maximums: These are two arrays indicating the lower and upper
    !                       bound of the uniform prior. 
    !! hypercube_indices:   a set of indices indicating where these parameters
    !!                      are to lie in the vector of the unit hypercube
    !! physical_indices:    a set of indices indicating where these parameters
    !!                      are to lie in the vector of physical parameters
    !!
    subroutine initialise_uniform(uniform_prior,hypercube_indices,physical_indices,minimums,maximums)
        implicit none
        !> The prior to be initialised
        type(prior),intent(inout) :: uniform_prior
        !> The hypercube indices
        integer, dimension(:), intent(in) :: hypercube_indices
        !> The physical indices
        integer, dimension(:), intent(in) :: physical_indices
        !> The minimums of the uniform transformation
        double precision, dimension(:), intent(in) :: minimums
        !> The maximums of the uniform transformation
        double precision, dimension(:), intent(in) :: maximums

        integer :: npars

        ! Get the number of parameters
        npars = size(hypercube_indices)

        ! Allocate arrays and pass in values
        allocate(uniform_prior%hypercube_indices(npars))
        allocate(uniform_prior%physical_indices(npars))
        allocate(uniform_prior%parameters(2*npars))

        uniform_prior%hypercube_indices = hypercube_indices
        uniform_prior%physical_indices = physical_indices
        uniform_prior%parameters(:npars) = minimums
        uniform_prior%parameters(npars+1:) = maximums

        ! Initialise the number of parameters
        uniform_prior%npars = npars

        ! Record its type as uniform
        uniform_prior%prior_type = uniform_type

    end subroutine initialise_uniform

    !> Uniform transformation
    function uniform_htp(hypercube_coords,parameters) result(physical_coords)
        implicit none
        !> The hypercube coordinates to be transformed
        double precision, intent(in), dimension(:) :: hypercube_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(hypercube_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: npars 

        npars=size(hypercube_coords)

        ! This is a fairly simple transformation, each parameter is transformed as
        ! hypercube_coord -> min + hypercube_coord * (max-min)
        ! Lower half of the parameters array are the minimums
        ! Upper half of the parameters array are the maximums

        physical_coords = parameters(:npars) + (parameters(npars+1:) - parameters(:npars) ) * hypercube_coords

    end function uniform_htp



    !> Routine to initialise a set of seperable gaussian parameters
    !!
    !! The user should provide:
    !!
    !! means & stdevs:      These are two arrays indicating the means and
    !!                      standard deviations of the gaussian priors
    !! hypercube_indices:   a set of indices indicating where these parameters
    !!                      are to lie in the vector of the unit hypercube
    !! physical_indices:    a set of indices indicating where these parameters
    !!                      are to lie in the vector of physical parameters
    !!
    subroutine initialise_gaussian(gaussian_prior,hypercube_indices,physical_indices,means,stdevs)
        implicit none
        !> The prior to be initialised
        type(prior),intent(inout) :: gaussian_prior
        !> The hypercube indices
        integer, dimension(:), intent(in) :: hypercube_indices
        !> The physical indices
        integer, dimension(:), intent(in) :: physical_indices
        !> The means of the gaussian transformation
        double precision, dimension(:), intent(in) :: means
        !> The standard deviations of the gaussian transformation
        double precision, dimension(:), intent(in) :: stdevs

        integer :: npars

        ! Get the number of parameters
        npars = size(hypercube_indices)

        ! Allocate arrays and pass in values
        allocate(gaussian_prior%hypercube_indices(npars))
        allocate(gaussian_prior%physical_indices(npars))
        allocate(gaussian_prior%parameters(2*npars))

        gaussian_prior%hypercube_indices = hypercube_indices
        gaussian_prior%physical_indices = physical_indices
        gaussian_prior%parameters(:npars) = means
        gaussian_prior%parameters(npars+1:) = stdevs

        ! Initialise the number of parameters
        gaussian_prior%npars = npars

        ! Record its type as uniform
        gaussian_prior%prior_type = gaussian_type

    end subroutine initialise_gaussian

    !> Gaussian transformation
    function gaussian_htp(hypercube_coords,parameters) result(physical_coords)
        implicit none
        !> The hypercube coordinates to be transformed
        double precision, intent(in), dimension(:) :: hypercube_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(hypercube_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: npars 

        npars=size(hypercube_coords)

        ! Transform via the intel function cdfnorm inverse 
        ! (v=vector, d=double, cdf=cumulative distribution function, norm=normal, inv=inverse)
        ! https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-67369FA5-ABFD-4B5D-82D4-E6A5E4AB565B.htm#GUID-67369FA5-ABFD-4B5D-82D4-E6A5E4AB565B
        ! Lower half of the parameters array are the means
        ! Upper half of the parameters array are the stdevs
        call vdcdfnorminv(npars,hypercube_coords,physical_coords)

        ! Scale by the standard deviation and shift by the mean
        physical_coords = parameters(:npars) + parameters(npars+1:) * physical_coords

    end function gaussian_htp


    




    !> Routine to initialise a set of seperable log-uniform parameters
    !!
    !! The user should provide:
    !!
    !! minimums & maximums: These are two arrays indicating the lower and upper
    !                       bound of the log uniform prior. 
    !! hypercube_indices:   a set of indices indicating where these parameters
    !!                      are to lie in the vector of the unit hypercube
    !! physical_indices:    a set of indices indicating where these parameters
    !!                      are to lie in the vector of physical parameters
    !!
    subroutine initialise_log_uniform(log_uniform_prior,hypercube_indices,physical_indices,minimums,maximums)
        implicit none
        !> The prior to be initialised
        type(prior),intent(inout) :: log_uniform_prior
        !> The hypercube indices
        integer, dimension(:), intent(in) :: hypercube_indices
        !> The physical indices
        integer, dimension(:), intent(in) :: physical_indices
        !> The minimums of the uniform transformation
        double precision, dimension(:), intent(in) :: minimums
        !> The maximums of the uniform transformation
        double precision, dimension(:), intent(in) :: maximums

        integer :: npars

        ! Get the number of parameters
        npars = size(hypercube_indices)

        ! Allocate arrays and pass in values
        allocate(log_uniform_prior%hypercube_indices(npars))
        allocate(log_uniform_prior%physical_indices(npars))
        allocate(log_uniform_prior%parameters(2*npars))

        log_uniform_prior%hypercube_indices = hypercube_indices
        log_uniform_prior%physical_indices = physical_indices
        log_uniform_prior%parameters(:npars) = minimums
        log_uniform_prior%parameters(npars+1:) = maximums

        ! Initialise the number of parameters
        log_uniform_prior%npars = npars

        ! Record its type as uniform
        log_uniform_prior%prior_type = log_uniform_type

    end subroutine initialise_log_uniform


    !> Log-Uniform transformation
    function log_uniform_htp(hypercube_coords,parameters) result(physical_coords)
        implicit none
        !> The hypercube coordinates to be transformed
        double precision, intent(in), dimension(:) :: hypercube_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(hypercube_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: npars 

        npars=size(hypercube_coords)

        ! hypercube_coord -> min * (max/min)**hypercube_coord
        ! Lower half of the parameters array are the minimums
        ! Upper half of the parameters array are the maximums
        physical_coords = parameters(:npars) * (parameters(npars+1:)/parameters(:npars)) ** hypercube_coords

    end function log_uniform_htp


end module priors_module
