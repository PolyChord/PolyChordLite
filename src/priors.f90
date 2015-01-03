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

    ! prior transformations

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
                physical_coords(priors(i)%physical_indices)= uniform_htp&
                    (hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            case(gaussian_type)
                physical_coords(priors(i)%physical_indices)= gaussian_htp&
                    (hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            case(log_uniform_type)
                physical_coords(priors(i)%physical_indices)= log_uniform_htp&
                    (hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            case(sorted_uniform_type)
                physical_coords(priors(i)%physical_indices)= sorted_uniform_htp&
                    (hypercube_coords(priors(i)%hypercube_indices),priors(i)%parameters)
            end select
        end do

    end function hypercube_to_physical

    function physical_to_hypercube(physical_coords,priors) result(hypercube_coords)
        implicit none
        type(prior), dimension(:), intent(in) :: priors
        double precision, intent(in), dimension(:) :: physical_coords

        double precision, dimension(size(physical_coords)) :: hypercube_coords

        integer :: i

        hypercube_coords=0d0

        do i=1,size(priors)
            select case(priors(i)%prior_type)
            case(uniform_type)
                hypercube_coords(priors(i)%hypercube_indices)= uniform_pth&
                    (physical_coords(priors(i)%physical_indices),priors(i)%parameters)
            case(gaussian_type)
                hypercube_coords(priors(i)%hypercube_indices)= gaussian_pth&
                    (physical_coords(priors(i)%physical_indices),priors(i)%parameters)
            case(log_uniform_type)
                hypercube_coords(priors(i)%hypercube_indices)= log_uniform_pth&
                    (physical_coords(priors(i)%physical_indices),priors(i)%parameters)
            case(sorted_uniform_type)
                hypercube_coords(priors(i)%hypercube_indices)= sorted_uniform_pth&
                    (physical_coords(priors(i)%physical_indices),priors(i)%parameters)
            end select
        end do

    end function physical_to_hypercube

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
            case(sorted_uniform_type)
                log_volume = log_volume + log(priors(i)%parameters(2)- priors(i)%parameters(1) ) - log(gamma(1d0+priors(i)%npars)) 
            end select
        end do

    end function prior_log_volume





    !============= Uniform (Separable) ================================================


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

    function uniform_pth(physical_coords,parameters) result(hypercube_coords)
        implicit none
        !> The physical coordinates to be transformed
        double precision, intent(in), dimension(:) :: physical_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(physical_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(physical_coords)) :: hypercube_coords

        integer :: npars 

        npars=size(physical_coords)


        hypercube_coords = (physical_coords - parameters(:npars)) / (parameters(npars+1:) - parameters(:npars) )

    end function uniform_pth


    !============= Gaussian (Separable) ===============================================

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
        use utils_module, only: inv_normal_cdf
        implicit none
        !> The hypercube coordinates to be transformed
        double precision, intent(in), dimension(:) :: hypercube_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(hypercube_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer :: npars 

        npars=size(hypercube_coords)

        ! Transform via the inverse normal cumulative distribution function
        physical_coords = inv_normal_cdf(hypercube_coords)

        ! Scale by the standard deviation and shift by the mean
        physical_coords = parameters(:npars) + parameters(npars+1:) * physical_coords

    end function gaussian_htp

    !> Gaussian transformation
    function gaussian_pth(physical_coords,parameters) result(hypercube_coords)
        use utils_module, only: normal_cdf
        implicit none
        !> The physical coordinates to be transformed
        double precision, intent(in), dimension(:) :: physical_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(physical_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(physical_coords)) :: hypercube_coords

        integer :: npars 

        npars=size(physical_coords)

        ! Scale by the standard deviation and shift by the mean
        ! Lower half of the parameters array are the means
        ! Upper half of the parameters array are the stdevs
        
        hypercube_coords = ( physical_coords - parameters(:npars) )/parameters(npars+1:)

        ! Transform via the normal cumulative distribution function
        hypercube_coords = normal_cdf(hypercube_coords)

    end function gaussian_pth

    



    !============= Log-Uniform (Separable) ============================================

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

    function log_uniform_pth(physical_coords,parameters) result(hypercube_coords)
        implicit none
        !> The physical coordinates to be transformed
        double precision, intent(in), dimension(:) :: physical_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2*size(physical_coords)) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(physical_coords)) :: hypercube_coords

        integer :: npars 

        npars=size(physical_coords)
                         
        ! hypercube_coord -> min * (max/min)**hypercube_coord
        ! Lower half of the parameters array are the minimums
        ! Upper half of the parameters array are the maximums
        hypercube_coords = log(physical_coords/parameters(:npars)) / log(parameters(npars+1:)/parameters(:npars))

    end function log_uniform_pth



    !============= Sorted Uniform =====================================================

    !> Routine to initialise a set of sorted uniform parameters
    !!
    !! These are a set of npars parameters \f$\{x_1,\ldots,x_N\}\f$ which are uniformly distributed
    !! between a maximum and a minimum, but are sorted so that \f$x_1<x_2<\ldots<x_N\f$.
    !!
    !! The user should provide:
    !!
    !! minimum & maximum: The minimum and maximum bounds of the sorted uniform
    !!                    prior
    !! hypercube_indices: a set of indices indicating where these parameters
    !!                    are to lie in the vector of the unit hypercube
    !! physical_indices:  a set of indices indicating where these parameters
    !!                    are to lie in the vector of physical parameters
    !!
    subroutine initialise_sorted_uniform(sorted_uniform_prior,hypercube_indices,physical_indices,minimum,maximum)
        implicit none
        !> The prior to be initialised
        type(prior),intent(inout) :: sorted_uniform_prior
        !> The hypercube indices
        integer, dimension(:), intent(in) :: hypercube_indices
        !> The physical indices
        integer, dimension(:), intent(in) :: physical_indices
        !> The minimums of the uniform transformation
        double precision, intent(in) :: minimum
        !> The maximums of the uniform transformation
        double precision, intent(in) :: maximum

        integer :: npars

        ! Get the number of parameters
        npars = size(hypercube_indices)

        ! Allocate arrays and pass in values
        allocate(sorted_uniform_prior%hypercube_indices(npars))
        allocate(sorted_uniform_prior%physical_indices(npars))
        allocate(sorted_uniform_prior%parameters(2))

        sorted_uniform_prior%hypercube_indices = hypercube_indices
        sorted_uniform_prior%physical_indices = physical_indices
        sorted_uniform_prior%parameters(1) = minimum
        sorted_uniform_prior%parameters(2) = maximum

        ! Initialise the number of parameters
        sorted_uniform_prior%npars = npars

        ! Record its type as uniform
        sorted_uniform_prior%prior_type = sorted_uniform_type

    end subroutine initialise_sorted_uniform

    !> This transforms the unit hypercube to a "forced identifiablity" prior.
    !! This means that the \f$(\theta_1,\theta_2,\ldots,\theta_n)\f$ variables are uniformly distributed in the
    !! physical prior space between \f$\theta_\mathrm{min}\f$ and \f$\theta_\mathrm{max}\f$, but have been sorted so
    !! that \f$(\theta_1<\theta_2<\ldots<\theta_n)\f$. This amounts to choosing a non-separable prior such that:
    !!
    !! \f[ \pi_n(\theta_n)            
    !!     =  n  \frac{(\theta_n-\theta_\mathrm{min})^{n-1}}{(\theta_\mathrm{max}-\theta_\mathrm{min})^n}
    !!     \qquad  \theta_\mathrm{min}<\theta_n<\theta_\mathrm{max} \f]
    !! \f[ \pi_{n-1}(\theta_{n-1}|\theta_n)            
    !!     =(n-1)\frac{(\theta_{n-1}-\theta_\mathrm{min})^{n-2}}{(\theta_\mathrm{max}-\theta_\mathrm{min})^{n-1}}
    !!     \qquad  \theta_\mathrm{min}<\theta_{n-1}<\theta_n \f]
    !! \f[ \pi_{n-2}(\theta_{n-2}|\theta_n,\theta{n-1})            
    !!     =(n-2)\frac{(\theta_{n-2}-\theta_\mathrm{min})^{n-3}}{(\theta_\mathrm{max}-\theta_\mathrm{min})^{n-2}}
    !!     \qquad  \theta_\mathrm{min}<\theta_{n-2}<\theta_{n-1} \f]
    !! \f[...\f]
    !! \f[ \pi_1(\theta_1|\theta_n,\ldots,\theta_2)            
    !!     =\frac{1}{(\theta_\mathrm{max}-\theta_\mathrm{min})}
    !!     \qquad  \theta_\mathrm{min}<\theta_1<\theta_2 \f]
    !!
    !! The first of these is the probability density for the largest of n points
    !! in \f$[\theta_\mathrm{min},\theta_\mathrm{max}]\f$.
    !!
    !! For the next highest point it is the smallest of \f$n-1\f$ points in the
    !! range \f$[\theta_\mathrm{min},\theta_n]\f$. 
    !!
    !! To perform this transformation, it is cleanest to do this in two
    !! steps. First perform the above transformations using the inverse of the cumulative
    !! distribution function:
    !! \f[ CDF^{-1}(x) = x^{1/n} \f]
    !!
    !! Then transform the unit hypercube into the physical space with a linear
    !! rescaling

    function sorted_uniform_htp(hypercube_coords,parameters) result(physical_coords)
        implicit none

        !> The hypercube coordinates to be transformed
        double precision, intent(in), dimension(:) :: hypercube_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(hypercube_coords)) :: physical_coords

        integer n_prior ! the dimension

        ! Get the size of the array
        n_prior = size(hypercube_coords)

        ! Transform the largest index to the largest of n_prior variables in [0,1]
        physical_coords(n_prior) = hypercube_coords(n_prior)**(1d0/n_prior)

        ! Then for the remaining variables, transform them to the largest of the
        ! remaining variables, and rescale so that the variable one larger is
        ! the maximum
        do n_prior=n_prior-1,1,-1
            physical_coords(n_prior) = hypercube_coords(n_prior)**(1d0/n_prior)*physical_coords(n_prior+1)
        end do

        ! Rescale using the parameters
        physical_coords = parameters(1) + physical_coords*(parameters(2)-parameters(1))

    end function sorted_uniform_htp


    function sorted_uniform_pth(physical_coords,parameters) result(hypercube_coords)
        implicit none

        !> The physical coordinates to be transformed
        double precision, intent(in), dimension(:) :: physical_coords

        !> The parameters of the transformation
        double precision, intent(in), dimension(2) :: parameters

        !> The transformed coordinates
        double precision, dimension(size(physical_coords)) :: hypercube_coords

        integer n_prior ! the dimension
        integer i_prior ! the dimension

        ! Get the size of the array
        n_prior = size(physical_coords)


        ! Rescale back to [0,1]
        hypercube_coords =  (physical_coords - parameters(1)) / (parameters(2)-parameters(1))  

        ! Undo the trasformation piece by piece
        do i_prior = 1,n_prior-1
            hypercube_coords(i_prior) = ( hypercube_coords(i_prior)/hypercube_coords(i_prior+1) )**i_prior
        end do
        hypercube_coords(n_prior) = hypercube_coords(n_prior)**n_prior


    end function sorted_uniform_pth







end module priors_module
