!> This allows for a simple C interface... 


module interfaces_module
    use utils_module, only: dp


contains

    subroutine simple_interface(loglikelihood_input, prior_input, nlive, num_repeats, do_clustering, feedback, precision_criterion, &
            max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, &
            write_stats, write_live, write_dead, update_files, nDims, nDerived, base_dir, file_root, grade_dims, grade_frac)

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
            function loglikelihood_input(theta,nDims,phi,nDerived)
                import :: dp
                integer, intent(in) :: nDims
                integer, intent(in) :: nDerived
                real(dp), intent(in),  dimension(nDims) :: theta
                real(dp), intent(out),  dimension(nDerived) :: phi
                real(dp) :: loglikelihood_input
            end function loglikelihood_input
        end interface
        interface
            function prior_input(cube,nDims) result(theta)
                import :: dp
                integer, intent(in) :: nDims
                real(dp), intent(in),  dimension(nDims) :: cube
                real(dp),  dimension(nDims) :: theta
            end function prior_input
        end interface

        integer, intent(in)          :: nlive
        integer, intent(in)          :: num_repeats
        logical, intent(in)          :: do_clustering
        integer, intent(in)          :: feedback
        real(dp), intent(in) :: precision_criterion
        integer, intent(in)          :: max_ndead
        real(dp), intent(in) :: boost_posterior
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
        real(dp), intent(in), dimension(:) :: grade_frac

        real(dp), dimension(4) :: output_info

        type(program_settings)    :: settings  ! The program settings 
        type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

        type(param_type),dimension(:),allocatable :: params         ! Parameter array
        type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array


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

        output_info = NestedSampling(loglikelihood_wrapper,prior_wrapper,settings,0) 

    contains
        function loglikelihood_wrapper(theta,phi)
            real(dp), intent(in),  dimension(:) :: theta
            real(dp), intent(out),  dimension(:) :: phi
            real(dp) :: loglikelihood_wrapper
            loglikelihood_wrapper = loglikelihood_input(theta,size(theta),phi,size(phi))
        end function loglikelihood_wrapper
        function prior_wrapper(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = prior_input(cube,size(cube))
        end function prior_wrapper

    end subroutine simple_interface

end module interfaces_module
