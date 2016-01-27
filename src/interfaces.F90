!> This allows for a simple C interface... 


module interfaces_module


    contains

subroutine simple_interface(loglikelihood_wrapper, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, &
    boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, &
    write_dead, update_files, nDims, nDerived, base_dir, file_root, grade_dims, grade_frac)

    use ini_module,               only: read_params,initialise_program,default_params
    use params_module,            only: add_parameter,param_type
    use priors_module
    use settings_module,          only: program_settings,initialise_settings
    use random_module,            only: initialise_random
    use nested_sampling_module,   only: NestedSampling
    use abort_module,             only: halt_program

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    interface
        function loglikelihood_wrapper(theta,phi,nDims,nDerived)
            integer, intent(in) :: nDims
            integer, intent(in) :: nDerived
            double precision, intent(in),  dimension(nDims) :: theta
            double precision, intent(out),  dimension(nDerived) :: phi
            double precision :: loglikelihood_wrapper
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

    character(*), intent(in)     :: base_dir
    character(*), intent(in)     :: file_root

    integer, intent(in), dimension(:) :: grade_dims
    double precision, intent(in), dimension(:) :: grade_frac





    double precision, dimension(5) :: output_info

    type(program_settings)    :: settings  ! The program settings 
    type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

    type(param_type),dimension(:),allocatable :: params         ! Parameter array
    type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

    ! Temporary variables for initialising loglikelihoods
    double precision :: loglike

    double precision, dimension(nDims)    :: theta0
    double precision, dimension(nDerived) :: phi0

    ! Basic initialisation
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

    settings%base_dir            = base_dir
    settings%file_root           = file_root

    allocate(settings%grade_dims(size(grade_dims)))
    settings%grade_dims          = grade_dims           
    allocate(settings%grade_frac(size(grade_frac)))
    settings%grade_frac          = grade_frac           

    params = default_params(settings%nDims,'theta')
    derived_params = default_params(settings%nDerived,'phi')

    call initialise_program(settings,priors,params,derived_params)

    output_info = NestedSampling(loglikelihood,priors,settings,0) 

    contains
        function loglikelihood(theta,phi)
            double precision, intent(in),  dimension(:) :: theta
            double precision, intent(out),  dimension(:) :: phi
            double precision :: loglikelihood
            loglikelihood = loglikelihood_wrapper(theta,phi,size(theta),size(phi))
        end function

end subroutine simple_interface

end module interfaces_module
