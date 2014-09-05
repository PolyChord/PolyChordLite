module priors_module
    use utils_module, only: stdout_unit
    implicit none

    integer, parameter :: uniform_type = 1

    type prior
        integer :: npars = 0
        integer :: prior_type = 1
        integer, dimension(:), allocatable :: hypercube_indices
        integer, dimension(:), allocatable :: physical_indices

        double precision, dimension(:), allocatable :: parameters

    end type prior



    contains

    subroutine initialise_uniform(uniform_prior,hypercube_indices,physical_indices,minimums,maximums)
        implicit none
        !> The prior to be initialised
        type(prior), intent(inout) :: uniform_prior
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

        ! Check that all inputs have the same size
        if(npars /= size(physical_indices)) then
            write(stdout_unit, '("Error in initialise_uniform: size(physical_indices)=",I4,"/= size(hypercube_indices)=",I4)') size(physical_indices),npars
            call exit(1)
        else if(npars /= size(minimums)) then
            write(stdout_unit, '("Error in initialise_uniform: size(minimums)=",I4,"/= size(hypercube_indices)=",I4)') size(minimums),npars
            call exit(1)
        else if(npars /= size(maximums)) then
            write(stdout_unit, '("Error in initialise_uniform: size(maximums)=",I4,"/= size(hypercube_indices)=",I4)') size(maximums),npars
            call exit(1)
        end if


        ! Allocate arrays and pass in values
        allocate(uniform_prior%hypercube_indices(npars))
        uniform_prior%hypercube_indices = hypercube_indices

        allocate(uniform_prior%physical_indices(npars))
        uniform_prior%physical_indices = physical_indices

        allocate(uniform_prior%parameters(2*npars))
        uniform_prior%parameters(:npars) = minimums
        uniform_prior%parameters(npars+1:) = maximums

        ! Initialise the rest
        uniform_prior%npars = npars
        uniform_prior%prior_type = uniform_type

    end subroutine initialise_uniform

    ! Prior transformations

    function hypercube_to_physical(hypercube_coords,priors) result(physical_coords)
        implicit none
        type(prior), dimension(:), intent(in) :: priors
        double precision, intent(in), dimension(:) :: hypercube_coords

        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: i

        do i=1,size(priors)
            select case(priors(i)%prior_type)
            case(uniform_type)
                physical_coords= uniform_htp(hypercube_coords,priors(i)%parameters)
            end select
        end do

    end function hypercube_to_physical












    !> Uniform transformation
    function uniform_htp(hypercube_coords,parameters) result(physical_coords)
        implicit none
        !> The hypercube coordinates to be transformed
        double precision, intent(in), dimension(:) :: hypercube_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(hypercube_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: nDims 

        nDims=size(hypercube_coords)

        ! This is a fairly simple transformation, each parameter is transformed as
        ! hypercube_coord -> min + hypercube_coord * (max-min)
        ! Param 1 is the minimum bound, Param 2 is the maximum

        physical_coords = parameters(1:nDims) + (parameters(nDims+1:2*nDims) - parameters(1:nDims) ) * hypercube_coords

    end function uniform_htp



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
                log_volume = log_volume + sum( log(priors(i)%parameters(:priors(i)%npars)- priors(i)%parameters(priors(i)%npars+1:) )) 
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



end module priors_module
