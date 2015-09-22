!> This module encodes the type 'program_settings' which contains all of the
!! details required to perform a nested sampling run.
module settings_module
    use utils_module,   only: STR_LENGTH
    implicit none

    integer, parameter :: live_type    = 1
    integer, parameter :: blank_type   = 0
    integer, parameter :: phantom_type =-1

    !> Type to contain all of the parameters involved in a nested sampling run
    Type :: program_settings

        integer :: nlive =500 !> The number of live points
        integer :: num_repeats !> The number of slow chords to draw
        logical :: do_clustering = .false.  !> Whether to do clustering or not

        integer :: feedback = 1 !> The degree of feedback to provide

        double precision :: precision_criterion = 1d-3 !> The stopping criterion

        !> The maximum number of dead points/samples
        !!
        !! Set equal to -1 for no maximum number
        integer :: max_ndead = -1

        !> What factor should we bulk up the posterior points by (using
        !! inter-chain points)
        !!
        !! set to <=0 to use all of them
        double precision :: boost_posterior = 5d0

        logical :: posteriors = .false.        !> Whether to calculate weighted posteriors
        logical :: equals     = .true.         !> Whether to calculate equally weighted posteriors
        logical :: cluster_posteriors = .false.!> Whether to calculate clustered posteriors


        
        logical :: write_resume = .false. !> Whether or not to write resume files
        logical :: write_paramnames = .false. !> Whether or not to write paramnames file
        logical :: read_resume = .false.  !> Whether or not to resume from file
        logical :: write_stats = .true.   !> Whether or not to write stats file
        logical :: write_live = .false.   !> Whether or not to write phys_live points

        integer :: update_files = 500 !> How often to update the resume file

        integer, allocatable,dimension(:) :: grade_dims          !> The number of parameters in each grade
        double precision, allocatable,dimension(:) :: grade_frac !> The fraction of time spent in each grade




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

        ! Loglikelihood index
        settings%l0=settings%d1+1

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

    end subroutine initialise_settings


end module settings_module 
