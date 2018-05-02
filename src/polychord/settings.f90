!> This module encodes the type 'program_settings' which contains all of the
!! details required to perform a nested sampling run.
module settings_module
    use utils_module, only: dp
    use utils_module,   only: STR_LENGTH
    implicit none

    integer, parameter :: live_type    = 1
    integer, parameter :: blank_type   = 0
    integer, parameter :: phantom_type =-1

    !> Type to contain all of the parameters involved in a nested sampling run
    Type :: program_settings

        integer :: nlive = 500                  !> The number of live points
        integer :: num_repeats = -1             !> The number of slow chords to draw
        logical :: do_clustering = .false.      !> Whether to do clustering or not

        integer :: feedback = 1                 !> The degree of feedback to provide

        real(dp) :: precision_criterion = 1d-3  !> The stopping criterion
        real(dp) :: logzero = -1d30             !> The threshold for 'log(0) values'

        !> The maximum number of dead points/samples
        !!
        !! Set equal to -1 for no maximum number
        integer :: max_ndead = -1

        !> What factor should we bulk up the posterior points by (using
        !! inter-chain points)
        !!
        !! set to <=0 to use all of them
        real(dp) :: boost_posterior = 0d0

        logical :: posteriors = .false.         !> Whether to calculate weighted posteriors
        logical :: equals     = .true.          !> Whether to calculate equally weighted posteriors
        logical :: cluster_posteriors = .false. !> Whether to calculate clustered posteriors


        
        logical :: write_resume = .false.     !> Whether or not to write resume files
        logical :: write_paramnames = .false. !> Whether or not to write paramnames file
        logical :: read_resume = .false.      !> Whether or not to resume from file
        logical :: write_stats = .true.       !> Whether or not to write stats file
        logical :: write_live = .true.        !> Whether or not to write phys_live points
        logical :: write_dead = .true.        !> Whether or not to write dead points

        integer, allocatable,dimension(:) :: grade_dims          !> The number of parameters in each grade
        real(dp), allocatable,dimension(:) :: grade_frac !> The fraction of time spent in each grade


        real(dp), allocatable,dimension(:) :: seed_point
        logical :: generate_from_seed = .false.
        integer :: nprior = -1
        integer :: nprior_repeat = -1
        logical :: write_prior = .false.

        integer :: seed = -1


        character(STR_LENGTH) :: base_dir='chains' !> The directory to put outputs in
        character(STR_LENGTH) :: file_root='test'  !> The file root for outputs


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

        !> Birth index
        !!
        !! This is the moment that the point was born
        integer :: b0          

        !> likelihood index
        !!
        !! This is the likelihood evaluated at the position of the live point
        integer :: l0          


        ! Posterior indices

        ! Weight index
        integer :: pos_w
        ! Cumulative weight index
        integer :: pos_Z
        ! Loglikelihood index
        integer :: pos_l
        ! Volume index
        integer :: pos_X
        ! physical parameter indices
        integer :: pos_p0
        integer :: pos_p1
        ! derived parameter indices
        integer :: pos_d0
        integer :: pos_d1

        ! Number of posterior parameters
        integer :: nposterior

        ! Final posterior indices
        integer :: p_w
        integer :: p_2l
        integer :: p_p0
        integer :: p_p1
        integer :: p_d0
        integer :: p_d1
        integer :: np



        ! Sub clustering dimensions
        integer,dimension(:),allocatable :: sub_clustering_dimensions

        integer,  dimension(:), allocatable :: nlives     !> The number of live points per contour
        real(dp), dimension(:), allocatable :: loglikes   !> The contours for nlive_list
        
        real(dp) :: compression_factor = exp(-1d0)

    end type program_settings

    contains

    !> This subroutine initialises all of the default settings
    !!
    !! It deals with all the fiddly stuff that's important, such as where to
    !! store things
    !!
    subroutine initialise_settings(settings)
        use abort_module, only: halt_program
        use utils_module, only: sort_doubles
        implicit none
        type(program_settings), intent(inout) :: settings


        ! Hypercube parameter indices
        settings%h0=1
        settings%h1=settings%h0+settings%nDims-1

        ! Physical parameter indices
        settings%p0=settings%h1+1
        settings%p1=settings%p0+settings%nDims-1 

        ! Derived parameter indices
        settings%d0=settings%p1+1
        settings%d1=settings%d0+settings%nDerived-1  

        ! birth index
        settings%b0=settings%d1+1

        ! Loglikelihood index
        settings%l0=settings%b0+1

        ! Total number of parameters
        settings%nTotal = settings%l0



        ! Posterior indices

        ! Volume index
        settings%pos_X = 1
        ! Loglikelihood index
        settings%pos_l = settings%pos_X+1
        ! Weight index
        settings%pos_w = settings%pos_l+1
        ! Cumulative weight index
        settings%pos_Z = settings%pos_w+1
        ! physical parameter indices
        settings%pos_p0= settings%pos_Z+1
        settings%pos_p1= settings%pos_Z+settings%nDims
        ! derived parameter indices
        settings%pos_d0= settings%pos_p1+1
        settings%pos_d1= settings%pos_p1+settings%nDerived

        ! Number of posterior parameters
        settings%nposterior = settings%pos_d1


        settings%p_w = 1
        settings%p_2l = settings%p_w+1
        settings%p_p0 = settings%p_2l+1
        settings%p_p1 = settings%p_2l+settings%nDims
        settings%p_d0 = settings%p_p1+1
        settings%p_d1 = settings%p_p1+settings%nDerived

        settings%np   = settings%p_d1

        if( settings%num_repeats < 1 ) call halt_program("You need to set num_repeats. Suggestion: 5*nDims")

        if(.not. allocated(settings%grade_dims)) then
            allocate(settings%grade_dims(1))
            settings%grade_dims = [settings%nDims]
        end if

        if(.not. allocated(settings%grade_frac)) then
            allocate(settings%grade_frac(1))
            settings%grade_frac = [1.0d0]
        end if

        if (.not. allocated(settings%nlives)) then
            allocate(settings%nlives(1),settings%loglikes(1))
            settings%nlives(1) = settings%nlive
            settings%loglikes(1) = settings%logzero
        end if

        settings%nlives = settings%nlives(sort_doubles(settings%loglikes))
        settings%loglikes = settings%loglikes(sort_doubles(settings%loglikes))



    end subroutine initialise_settings


end module settings_module 
