!> This allows for a simple C interface... 


module interfaces_module


    contains

subroutine simple_interface(loglikelihood, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, &
    boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, &
    write_dead, update_files, nDims, nDerived, grade_dims, grade_frac)

    use ini_module,               only: read_params,initialise_program,default_params
    use params_module,            only: add_parameter,param_type
    use priors_module
    use settings_module,          only: program_settings,initialise_settings
    use random_module,            only: initialise_random
    use nested_sampling_module,   only: NestedSampling
    use utils_module,             only: STR_LENGTH
    use abort_module,             only: halt_program
#ifdef MPI
    use mpi_module,               only: initialise_mpi, finalise_mpi
    use mpi,                      only: MPI_COMM_WORLD
#endif

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    interface
        function loglikelihood(theta,phi)
            double precision, intent(in),  dimension(:) :: theta
            double precision, intent(out),  dimension(:) :: phi
            double precision :: loglikelihood
        end function
    end interface

    integer, intent(in)          :: nlive
    integer, intent(in)          :: num_repeats
    logical, intent(in)          :: do_clustering
    integer, intent(in)          :: feedback
    double precision, intent(in) :: precision_criterion
    integer, intent(in)          :: max_ndead
    double precision, intent(in) :: boost_posterior
    logical, intent(in)          :: posteriors
    logical, intent(in)          :: equals
    logical, intent(in)          :: cluster_posteriors
    logical, intent(in)          :: write_resume
    logical, intent(in)          :: write_paramnames
    logical, intent(in)          :: read_resume
    logical, intent(in)          :: write_stats
    logical, intent(in)          :: write_live
    logical, intent(in)          :: write_dead
    integer, intent(in)          :: update_files
    integer, intent(in)          :: nDims
    integer, intent(in)          :: nDerived

    integer, intent(in), dimension(:) :: grade_dims          !> The number of parameters in each grade
    double precision, intent(in), dimension(:) :: grade_frac !> The fraction of time spent in each grade




    double precision, dimension(5) :: output_info

    type(program_settings)    :: settings  ! The program settings 
    type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

    character(len=STR_LENGTH)                 :: input_file     ! input file
    type(param_type),dimension(:),allocatable :: params         ! Parameter array
    type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

    ! Temporary variables for initialising loglikelihoods
    double precision :: loglike

    ! Basic initialisation
#ifdef MPI
    call initialise_mpi()
#endif
    call initialise_random()

    settings%nlive               = nlive                
    settings%num_repeats         = num_repeats          
    settings%do_clustering       = do_clustering        
    settings%feedback            = feedback             
    settings%precision_criterion = precision_criterion  
    settings%max_ndead           = max_ndead            
    settings%boost_posterior     = boost_posterior      
    settings%posteriors          = posteriors           
    settings%equals              = equals               
    settings%cluster_posteriors  = cluster_posteriors   
    settings%write_resume        = write_resume         
    settings%write_paramnames    = write_paramnames     
    settings%read_resume         = read_resume          
    settings%write_stats         = write_stats          
    settings%write_live          = write_live           
    settings%write_dead          = write_dead           
    settings%update_files        = update_files         
    settings%nDims               = nDims
    settings%nDerived            = nDerived

    allocate(settings%grade_dims(size(grade_dims)))
    settings%grade_dims          = grade_dims           
    allocate(settings%grade_frac(size(grade_frac)))
    settings%grade_frac          = grade_frac           

    params = default_params(settings%nDims,'theta')
    derived_params = default_params(settings%nDerived,'phi')

    call initialise_program(settings,priors,params,derived_params)

#ifdef MPI
    output_info = NestedSampling(loglikelihood,priors,settings,MPI_COMM_WORLD) 
#else
    output_info = NestedSampling(loglikelihood,priors,settings,0) 
#endif

#ifdef MPI
    call finalise_mpi
#endif

end subroutine simple_interface

end module interfaces_module
