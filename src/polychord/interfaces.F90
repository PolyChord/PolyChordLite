!> This allows for a simple C interface... 


module interfaces_module
    use utils_module, only: dp

    implicit none
    interface run_polychord
        module procedure run_polychord_full, run_polychord_no_prior
    end interface run_polychord

contains


    subroutine run_polychord_full(loglikelihood, prior, settings_in)
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
        output_info = NestedSampling(loglikelihood,prior,settings,MPI_COMM_WORLD) 
        call finalise_mpi
#else
        output_info = NestedSampling(loglikelihood,prior,settings,0) 
#endif

    end subroutine run_polychord_full


    !===================== INTERFACE ===============================================
    subroutine run_polychord_no_prior(loglikelihood, settings)
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
        call run_polychord(loglikelihood,prior,settings)
    contains
        function prior(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = cube
        end function prior
    end subroutine run_polychord_no_prior










    subroutine polychord_c_interface(c_loglikelihood_ptr, c_prior_ptr, nlive, num_repeats, do_clustering, feedback, &
        precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, &
            write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived, &
            base_dir, file_root) &
            bind(c,name='polychord_c_interface')

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
            function c_loglikelihood(theta,nDims,phi,nDerived) bind(c)
                use iso_c_binding
                integer(c_int), intent(in), value :: nDims
                integer(c_int), intent(in), value :: nDerived
                real(c_double), intent(in),  dimension(nDims) :: theta
                real(c_double), intent(out),  dimension(nDerived) :: phi
                real(c_double) :: c_loglikelihood
            end function c_loglikelihood
        end interface
        interface
            subroutine c_prior(cube,theta,nDims) bind(c)
                use iso_c_binding
                integer(c_int), intent(in), value :: nDims
                real(c_double), intent(in),  dimension(nDims) :: cube
                real(c_double), intent(out), dimension(nDims) :: theta
            end subroutine c_prior
        end interface

        type(c_funptr), intent(in), value   :: c_loglikelihood_ptr
        type(c_funptr), intent(in), value   :: c_prior_ptr
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

        procedure(c_loglikelihood), pointer :: f_loglikelihood_ptr
        procedure(c_prior), pointer         :: f_prior_ptr

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

        call c_f_procpointer(c_loglikelihood_ptr, f_loglikelihood_ptr)
        call c_f_procpointer(c_prior_ptr, f_prior_ptr)

        call run_polychord(loglikelihood,prior,settings) 

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
            c_loglike = f_loglikelihood_ptr(c_theta,c_nDims,c_phi,c_nDerived)
            phi = c_phi
            loglikelihood = c_loglike

        end function loglikelihood

        function prior(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta

            integer (c_int)                       :: c_nDims
            real (c_double),dimension(size(cube)) :: c_cube
            real (c_double),dimension(size(cube)) :: c_theta

            c_nDims = size(cube)
            c_cube = cube
            call f_prior_ptr(c_cube,c_theta,c_nDims)
            theta = c_theta

        end function prior

    end subroutine polychord_c_interface



end module interfaces_module
