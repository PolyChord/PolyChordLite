!> This module encodes the type 'program_settings' which contains all of the
!! details required to perform a nested sampling run.
module settings_module
    use priors_module,   only: prior
    use utils_module,   only: STR_LENGTH
    use grades_module,  only: parameter_grades
    implicit none

    integer, parameter :: live_type    = 1
    integer, parameter :: blank_type   = 0
    integer, parameter :: phantom_type =-1

    !> Type to contain all of the parameters involved in a nested sampling run
    Type :: program_settings

        !> The directory to put outputs in
        character(STR_LENGTH) :: base_dir='chains'
        !> The file root for outputs
        character(STR_LENGTH) :: file_root='test'

        !> The number of live points
        integer :: nlive =500

        !> The number of live points
        integer :: nstack =500*10

        !> The degree of feedback to provide
        integer :: feedback = 1

        !> The degree of precision in the final answer
        double precision :: precision_criterion = 1d-3

        !> The maximum number of dead points/samples
        !!
        !! Set equal to -1 for no maximum number
        integer :: max_ndead = -1

        !> The maximum number of posterior points
        !!
        !! This is for memory allocation purposes, it won't necessarily have
        !! this many points if they're not 'good enough'
        integer :: sigma_posterior = 5

        !> What to thin the posterior points (i.e. probability of keeping
        !! phantom points)
        double precision :: thin_posterior = 0

        !> Whether or not to calculate the posterior
        logical :: calculate_posterior = .false.

        !> Whether or not to write resume files
        logical :: write_resume = .false.

        !> How often to update the resume file
        integer :: update_resume = 500

        !> How many baby points are generated in each iteration
        integer :: num_babies

        !> The number of chords to draw
        integer :: num_repeats

        !> Whether or not to resume from file
        logical :: read_resume = .false.

        !> Whether or not to write phys_live points
        logical :: write_live = .false.

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

        !> likelihood index
        !!
        !! This is the likelihood evaluated at the position of the live point
        integer :: l0          
        !> This is the likelihood contour which generated the live point
        integer :: l1          


        ! Posterior indices

        ! Weight index
        integer :: pos_w
        ! Cumulative weight index
        integer :: pos_Z
        ! Loglikelihood index
        integer :: pos_l
        ! physical parameter indices
        integer :: pos_p0
        integer :: pos_p1
        ! derived parameter indices
        integer :: pos_d0
        integer :: pos_d1

        ! Number of posterior parameters
        integer :: nposterior



        !> Pointer to any additional files that need to be stored in between
        !! evaluations (only really important for C likelihoods)
        integer :: context

        !> Grades of parameters
        logical :: do_grades = .false.
        type(parameter_grades) :: grades

        !> Whether to time likelihood calls for grades
        logical :: do_timing = .false.
        !> How long to wait in between the printing of times
        integer :: print_timing = 100

        !> Whether to do clustering or not
        logical :: do_clustering = .false.

        !> Clustering parameters
        integer :: SNN_k = 20

        !> Memory storage considerations; The maximum number of active clusters
        integer :: ncluster = 30

        !> To generate from a seed point instead
        logical :: generate_from_seed = .false.
        integer :: ngenerate = 500*10
        integer :: generate_burn_in = 0
        double precision, allocatable,dimension(:) :: seed_point

    end type program_settings

    contains

    !> This subroutine initialises all of the default settings
    !!
    !! It deals with all the fiddly stuff that's important, such as where to
    !! store things
    !!
    subroutine initialise_settings(settings)
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

        ! Algorithm indices
        settings%nlike=settings%d1+1
        settings%last_chord=settings%nlike+1

        ! Loglikelihood indices
        settings%l0=settings%last_chord+1
        settings%l1=settings%l0+1

        ! Total number of parameters
        settings%nTotal = settings%l1



        ! Number of babies to generate,
        settings%num_babies           = settings%nDims*settings%num_repeats
        ! Space in stack for excess phantoms
        settings%nstack               = settings%nlive*settings%num_babies*2


        ! Posterior indices

        ! Weight index
        settings%pos_w = 1
        ! Cumulative weight index
        settings%pos_Z = settings%pos_w+1
        ! Loglikelihood index
        settings%pos_l = settings%pos_Z+1
        ! physical parameter indices
        settings%pos_p0= settings%pos_l+1
        settings%pos_p1= settings%pos_l+settings%nDims
        ! derived parameter indices
        settings%pos_d0= settings%pos_p1+1
        settings%pos_d1= settings%pos_p1+settings%nDerived

        ! Number of posterior parameters
        settings%nposterior = settings%pos_d1


    end subroutine initialise_settings


end module settings_module 
