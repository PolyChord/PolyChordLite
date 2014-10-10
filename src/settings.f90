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

    ! Samplers
    integer, parameter :: sampler_covariance=0
    integer, parameter :: sampler_graded_covariance=1
    integer, parameter :: sampler_adaptive_parallel=2

    !> Type to contain all of the parameters involved in a nested sampling run
    Type :: program_settings

        !> The file root for outputs
        character(STR_LENGTH) :: file_root='chains/test'

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
        integer :: nmax_posterior = 100000

        !> Whether or not to calculate the posterior
        logical :: calculate_posterior = .true.

        !> Whether or not to write resume files
        logical :: write_resume = .true.

        !> How often to update the resume file
        integer :: update_resume = 10000

        !> The number of chords to draw
        integer :: chain_length

        !> Whether or not to resume from file
        logical :: read_resume = .true.

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
        !> Whether or not this is a 'true' live point
        !!
        !! index |    name      | meaning
        !!-------|--------------|--------------------------------------------
        !!   -1  | phantom_type | 'blank' live point (blank slot to be filled)
        !!    0  | blank_type   | 'fake' live point (not independent)
        !!    1  | live_type    | 'true' live point
        !!
        integer :: point_type

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
        type(parameter_grades) :: grades

        !> Which sampling algorithm to use
        integer :: sampler = sampler_covariance

    end type program_settings

    contains

    subroutine allocate_indices(settings)
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
        settings%point_type=settings%last_chord+1

        ! Loglikelihood indices
        settings%l0=settings%point_type+1
        settings%l1=settings%l0+1

        ! Total number of parameters
        settings%nTotal = settings%l1

    end subroutine allocate_indices


end module settings_module 
