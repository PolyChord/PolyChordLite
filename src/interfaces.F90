!> This allows for a simple C interface... 


module interfaces_module
    use utils_module, only: dp

    implicit none
    interface run_polychord
        module procedure run_polychord_full, run_polychord_no_prior, run_polychord_no_setup, run_polychord_no_prior_no_setup
    end interface run_polychord

contains


    subroutine run_polychord_full(loglikelihood, prior, setup_loglikelihood, settings_in)
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
#ifdef MPI
        use mpi_module,               only: initialise_mpi, finalise_mpi
        use mpi,                      only: MPI_COMM_WORLD
#endif
        implicit none

        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface
        interface
            function prior(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function prior
        end interface
        interface
            subroutine setup_loglikelihood(settings,mpi_communicator)
                import :: program_settings
                type(program_settings), intent(in) :: settings
                integer,intent(in) :: mpi_communicator
            end subroutine setup_loglikelihood
        end interface

        type(program_settings),intent(in)    :: settings_in
        type(program_settings)               :: settings 

        real(dp), dimension(4) :: output_info

#ifdef MPI
        call initialise_mpi
#endif

        call initialise_random()
        settings = settings_in
        call initialise_settings(settings)   

#ifdef MPI
        call setup_loglikelihood(settings,MPI_COMM_WORLD)
        output_info = NestedSampling(loglikelihood,prior,settings,MPI_COMM_WORLD) 
        call finalise_mpi
#else
        call setup_loglikelihood(settings,0)
        output_info = NestedSampling(loglikelihood,prior,settings,0) 
#endif

    end subroutine run_polychord_full


    !===================== INTERFACE ===============================================
    subroutine run_polychord_no_setup(loglikelihood, prior, settings)
        use settings_module,          only: program_settings
        implicit none
        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface
        interface
            function prior(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function prior
        end interface
        type(program_settings),intent(in)    :: settings  ! The program settings 

        call run_polychord(loglikelihood,prior,setup_loglikelihood,settings)
    contains

        subroutine setup_loglikelihood(settings,mpi_communicator)
            implicit none
            type(program_settings), intent(in) :: settings
            integer,intent(in) :: mpi_communicator
        end subroutine setup_loglikelihood

    end subroutine run_polychord_no_setup

    !===================== INTERFACE ===============================================
    subroutine run_polychord_no_prior(loglikelihood, setup_loglikelihood, settings)
        use settings_module,          only: program_settings
        implicit none
        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface
        interface
            subroutine setup_loglikelihood(settings,mpi_communicator)
                import :: program_settings
                type(program_settings), intent(in) :: settings
                integer,intent(in) :: mpi_communicator
            end subroutine setup_loglikelihood
        end interface
        type(program_settings),intent(in)    :: settings  ! The program settings 

        call run_polychord(loglikelihood,prior,setup_loglikelihood,settings)
    contains

        function prior(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = cube
        end function prior

    end subroutine run_polychord_no_prior

    !===================== INTERFACE ===============================================
    subroutine run_polychord_no_prior_no_setup(loglikelihood, settings)
        use settings_module,          only: program_settings
        implicit none
        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface
        type(program_settings),intent(in)    :: settings  ! The program settings 

        call run_polychord(loglikelihood,prior,setup_loglikelihood,settings)
    contains

        function prior(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = cube
        end function prior
        subroutine setup_loglikelihood(settings,mpi_communicator)
            implicit none
            type(program_settings), intent(in) :: settings
            integer,intent(in) :: mpi_communicator
        end subroutine setup_loglikelihood

    end subroutine run_polychord_no_prior_no_setup


    subroutine simple_c_interface(c_loglikelihood_ptr, nlive, num_repeats, do_clustering, feedback, &
        precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, &
            write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived, &
            base_dir, file_root) &
            bind(c,name='polychord_simple_c_interface')

        use iso_c_binding
        use utils_module,             only: STR_LENGTH, convert_c_string
        use ini_module,               only: default_params
        use params_module,            only: param_type
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
        use read_write_module,        only: write_paramnames_file

        ! ~~~~~~~ Local Variable Declaration ~~~~~~~
        implicit none

        interface
            function c_loglikelihood(theta,nDims,phi,nDerived)
                use iso_c_binding
                integer(c_int), intent(in), value :: nDims
                integer(c_int), intent(in), value :: nDerived
                real(c_double), intent(in),  dimension(nDims) :: theta
                real(c_double), intent(out),  dimension(nDerived) :: phi
                real(c_double) :: c_loglikelihood
            end function c_loglikelihood
        end interface

        type(c_funptr), intent(in), value   :: c_loglikelihood_ptr
        integer(c_int), intent(in), value   :: nlive
        integer(c_int), intent(in), value   :: num_repeats
        logical(c_bool), intent(in), value  :: do_clustering
        integer(c_int), intent(in), value   :: feedback
        real(c_double), intent(in), value   :: precision_criterion
        integer(c_int), intent(in), value   :: max_ndead
        real(c_double), intent(in), value   :: boost_posterior
        logical(c_bool), intent(in), value  :: posteriors
        logical(c_bool), intent(in), value  :: equals
        logical(c_bool), intent(in), value  :: cluster_posteriors
        logical(c_bool), intent(in), value  :: write_resume
        logical(c_bool), intent(in), value  :: write_paramnames
        logical(c_bool), intent(in), value  :: read_resume
        logical(c_bool), intent(in), value  :: write_stats
        logical(c_bool), intent(in), value  :: write_live
        logical(c_bool), intent(in), value  :: write_dead
        integer(c_int), intent(in), value   :: update_files
        integer(c_int), intent(in), value   :: nDims
        integer(c_int), intent(in), value   :: nDerived
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: base_dir
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: file_root

        type(program_settings)    :: settings  ! The program settings 

        type(param_type),dimension(:),allocatable :: params         ! Parameter array
        type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

        procedure(c_loglikelihood), pointer :: f_loglikelihood

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

        settings%base_dir            = convert_c_string(base_dir)
        settings%file_root           = convert_c_string(file_root) 

        if(settings%write_paramnames) then
            params = default_params(settings%nDims,'theta','\theta')
            derived_params = default_params(settings%nDerived,'phi','\phi')
            call write_paramnames_file(settings,params,derived_params)
        end if

        call c_f_procpointer(c_loglikelihood_ptr, f_loglikelihood)

        call run_polychord(loglikelihood,settings) 

    contains
        function loglikelihood(theta,phi)
            implicit none
            real(dp), intent(in),  dimension(:) :: theta
            real(dp), intent(out),  dimension(:) :: phi
            real(dp) :: loglikelihood

            real (c_double),dimension(size(theta)) :: c_theta
            integer (c_int)                        :: c_nDims
            real (c_double),dimension(size(phi))   :: c_phi
            integer (c_int)                        :: c_nDerived
            real (c_double)                        :: c_loglike

            c_nDims = size(theta)
            c_nDerived = size(phi)
            c_theta = theta
            c_phi = phi
            c_loglike = f_loglikelihood(c_theta,c_nDims,c_phi,c_nDerived)
            loglikelihood = c_loglike

        end function loglikelihood


    end subroutine simple_c_interface


    subroutine simple_interface(loglikelihood_input, prior_input, nlive, num_repeats, do_clustering, feedback, &
        precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, &
            write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived, &
            base_dir, file_root, grade_dims, grade_frac)

        use ini_module,               only: default_params
        use params_module,            only: param_type
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
        use read_write_module,        only: write_paramnames_file

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

        type(program_settings)    :: settings  ! The program settings 

        type(param_type),dimension(:),allocatable :: params         ! Parameter array
        type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array


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

        if(settings%write_paramnames) then
            params = default_params(settings%nDims,'theta','\theta')
            derived_params = default_params(settings%nDerived,'phi','\phi')
            call write_paramnames_file(settings,params,derived_params)
        end if

        call run_polychord(loglikelihood_wrapper,prior_wrapper,settings) 

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
