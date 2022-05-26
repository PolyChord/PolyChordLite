!> This allows for a simple C interface... 
module interfaces_module
    use utils_module, only: dp
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    interface run_polychord
        module procedure run_polychord_full, run_polychord_no_cluster, run_polychord_no_prior, run_polychord_no_dumper, run_polychord_no_prior_no_dumper, run_polychord_ini
    end interface run_polychord

contains


    subroutine run_polychord_full(loglikelihood, prior_transform, dumper, cluster, settings_in, mpi_communicator)
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
#ifdef MPI
        use mpi_module,               only: initialise_mpi, finalise_mpi
#endif
        implicit none
#ifdef MPI
        include 'mpif.h'
#endif

        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface
        interface
            function prior_transform(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function prior_transform
        end interface
        interface
            subroutine dumper(live, dead, logweights, logZ, logZerr)
                import :: dp
                real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
                real(dp), intent(in) :: logZ, logZerr
            end subroutine dumper
        end interface
        interface
            function cluster(points) result(cluster_list)
                import :: dp
                real(dp), intent(in), dimension(:,:) :: points
                integer, dimension(size(points,2)) :: cluster_list
            end function
        end interface

        type(program_settings),intent(in)    :: settings_in
        type(program_settings)               :: settings 

        !> MPI handle
        integer, intent(in), optional :: mpi_communicator
        integer :: comm

        real(dp), dimension(4) :: output_info

#ifdef MPI
        if (present(mpi_communicator)) then
            comm = mpi_communicator
        else 
            comm = MPI_COMM_WORLD
        end if
        call initialise_mpi(settings_in%feedback)
#else 
        comm =0
#endif
        if (settings_in%seed >= 0) then
            call initialise_random(comm,settings_in%seed)
        else
            call initialise_random(comm)
        end if
        settings = settings_in
        call initialise_settings(settings)   
#ifdef MPI
        output_info = NestedSampling(loglikelihood,prior_transform,dumper,cluster,settings,comm) 
        call finalise_mpi
#else
        output_info = NestedSampling(loglikelihood,prior_transform,dumper,cluster,settings,0) 
#endif

    end subroutine run_polychord_full


    !===================== INTERFACE ===============================================
    subroutine run_polychord_no_cluster(loglikelihood,prior_transform,dumper,settings,mpi_communicator)
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
            function prior_transform(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function prior_transform
        end interface
        interface
            subroutine dumper(live, dead, logweights, logZ, logZerr)
                import :: dp
                real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
                real(dp), intent(in) :: logZ, logZerr
            end subroutine dumper
        end interface
        type(program_settings),intent(in)    :: settings  ! The program settings 
        integer, intent(in), optional :: mpi_communicator
        integer :: comm
#ifdef MPI
        if (present(mpi_communicator)) then
            comm = mpi_communicator
        else 
            comm = MPI_COMM_WORLD
        end if
#else
            comm = 0
#endif
        call run_polychord(loglikelihood,prior_transform,dumper,cluster,settings,comm)
    contains
        function cluster(points) result(cluster_list)
            real(dp), intent(in), dimension(:,:) :: points
            integer, dimension(size(points,2)) :: cluster_list
            cluster_list = 0
        end function
    end subroutine run_polychord_no_cluster

    subroutine run_polychord_no_prior(loglikelihood,dumper,settings,mpi_communicator)
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
            subroutine dumper(live, dead, logweights, logZ, logZerr)
                import :: dp
                real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
                real(dp), intent(in) :: logZ, logZerr
            end subroutine dumper
        end interface
        type(program_settings),intent(in)    :: settings  ! The program settings 
        integer, intent(in), optional :: mpi_communicator
        integer :: comm
#ifdef MPI
        if (present(mpi_communicator)) then
            comm = mpi_communicator
        else 
            comm = MPI_COMM_WORLD
        end if
#else
            comm = 0
#endif
        call run_polychord(loglikelihood,prior_transform,dumper,settings,comm)
    contains
        function prior_transform(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = cube
        end function prior_transform
    end subroutine run_polychord_no_prior

    subroutine run_polychord_no_dumper(loglikelihood,prior_transform,settings,mpi_communicator)
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
            function prior_transform(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function prior_transform
        end interface
        type(program_settings),intent(in)    :: settings  ! The program settings 
        integer, intent(in), optional :: mpi_communicator
        integer :: comm
#ifdef MPI
        if (present(mpi_communicator)) then
            comm = mpi_communicator
        else 
            comm = MPI_COMM_WORLD
        end if
#else
            comm = 0
#endif
        call run_polychord(loglikelihood,prior_transform,dumper,settings,comm)
    contains
        subroutine dumper(live, dead, logweights, logZ, logZerr)
            real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
            real(dp), intent(in) :: logZ, logZerr
        end subroutine dumper
    end subroutine run_polychord_no_dumper

    subroutine run_polychord_no_prior_no_dumper(loglikelihood, settings,mpi_communicator)
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
        !> MPI handle
        integer, intent(in), optional :: mpi_communicator
        integer :: comm
#ifdef MPI
        if (present(mpi_communicator)) then
            comm = mpi_communicator
        else 
            comm = MPI_COMM_WORLD
        end if
#else
            comm = 0
#endif
        call run_polychord(loglikelihood,prior_transform,settings, comm)
    contains
        function prior_transform(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = cube
        end function prior_transform
        subroutine dumper(live, dead, logweights, logZ, logZerr)
            real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
            real(dp), intent(in) :: logZ, logZerr
        end subroutine dumper
    end subroutine run_polychord_no_prior_no_dumper


    subroutine run_polychord_ini(loglikelihood, setup_loglikelihood, input_file,mpi_communicator)
        use ini_module,               only: read_params
        use read_write_module,        only: write_paramnames_file
        use params_module,            only: add_parameter,param_type
        use priors_module
        use settings_module,          only: program_settings
        use utils_module,             only: STR_LENGTH
        implicit none
        type(program_settings)                    :: settings  ! The program settings 
        type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

        character(len=STR_LENGTH), intent(in)     :: input_file     ! input file
        type(param_type),dimension(:),allocatable :: params         ! Parameter array
        type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array
        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),  dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface
        interface
            subroutine setup_loglikelihood(settings)
                import :: program_settings
                type(program_settings), intent(in) :: settings
            end subroutine setup_loglikelihood
        end interface
        integer, intent(in), optional :: mpi_communicator
        integer :: comm
#ifdef MPI
        if (present(mpi_communicator)) then
            comm = mpi_communicator
        else 
            comm = MPI_COMM_WORLD
        end if
#else
            comm = 0
#endif
        call read_params(trim(input_file),settings,params,derived_params)
        if(settings%write_paramnames) call write_paramnames_file(settings,params,derived_params)
        call create_priors(priors,params,settings)
        call setup_loglikelihood(settings)
        call run_polychord(loglikelihood,prior_wrapper,settings,comm) 

        contains
            function prior_wrapper(cube) result(theta)
                implicit none
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
                theta = hypercube_to_physical(cube,priors)
            end function prior_wrapper
    end subroutine run_polychord_ini


    subroutine polychord_c_interface(&
            c_loglikelihood_ptr,&
            c_prior_ptr,&
            c_dumper_ptr,&
            c_cluster_ptr,&
            nlive,&
            num_repeats,&
            nprior,&
            nfail,&
            do_clustering,&
            feedback, &
            precision_criterion,&
            logzero,&
            max_ndead,&
            boost_posterior,&
            posteriors,&
            equals,&
            cluster_posteriors,&
            write_resume, &
            write_paramnames,&
            read_resume,&
            write_stats,&
            write_live,&
            write_dead,&
            write_prior,&
            maximise,&
            compression_factor,&
            synchronous,&
            nDims,&
            nDerived, &
            base_dir,&
            file_root,&
            nGrade,&
            grade_frac,&
            grade_dims,&
            n_nlives,&
            loglikes,&
            nlives,&
            seed,&
            comm) &
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
        interface
            subroutine c_dumper(ndead, nlive, npars, live, dead, logweights, logZ, logZerr) bind(c)
                use iso_c_binding
                integer(c_int), intent(in), value :: ndead, nlive, npars
                real(c_double), intent(in), dimension(npars,nlive) :: live, dead(npars,ndead), logweights(ndead)
                real(c_double), intent(in), value :: logZ, logZerr
            end subroutine c_dumper
        end interface
        interface
            subroutine c_cluster(points,cluster_list,nDims,nPoints) bind(c)
                use iso_c_binding
                integer(c_int), intent(in), value :: nDims, nPoints
                real(c_double), intent(in),  dimension(nDims,nPoints) :: points
                integer(c_int), intent(out), dimension(nPoints) :: cluster_list
            end subroutine c_cluster
        end interface

        type(c_funptr), intent(in), value   :: c_loglikelihood_ptr
        type(c_funptr), intent(in), value   :: c_prior_ptr
        type(c_funptr), intent(in), value   :: c_dumper_ptr
        type(c_funptr), intent(in), value   :: c_cluster_ptr
        integer(c_int), intent(in), value   :: nlive
        integer(c_int), intent(in), value   :: num_repeats
        integer(c_int), intent(in), value   :: nprior
        integer(c_int), intent(in), value   :: nfail
        logical(c_bool), intent(in), value  :: do_clustering
        integer(c_int), intent(in), value   :: feedback
        real(c_double), intent(in), value   :: precision_criterion
        real(c_double), intent(in), value   :: logzero
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
        logical(c_bool), intent(in), value  :: write_prior
        logical(c_bool), intent(in), value  :: maximise
        real(c_double),  intent(in), value  :: compression_factor
        logical(c_bool), intent(in), value  :: synchronous
        integer(c_int), intent(in), value   :: nDims
        integer(c_int), intent(in), value   :: nDerived
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: base_dir
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: file_root

        integer(c_int), intent(in), value             :: nGrade
        real(c_double), intent(in), dimension(nGrade) :: grade_frac
        integer(c_int), intent(in), dimension(nGrade) :: grade_dims

        integer(c_int), intent(in), value               :: n_nlives
        real(c_double), intent(in), dimension(n_nlives) :: loglikes
        integer(c_int), intent(in), dimension(n_nlives) :: nlives
        integer(c_int), intent(in), value               :: seed

        type(program_settings)    :: settings  ! The program settings 

        type(param_type),dimension(:),allocatable :: params         ! Parameter array
        type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

        procedure(c_loglikelihood), pointer :: f_loglikelihood_ptr
        procedure(c_prior), pointer         :: f_prior_ptr
        procedure(c_cluster), pointer       :: f_cluster_ptr
        procedure(c_dumper), pointer        :: f_dumper_ptr

        integer, intent(in) :: comm

        settings%nlive               = nlive                
        settings%num_repeats         = num_repeats          
        settings%nprior              = nprior
        settings%nfail               = nfail
        settings%do_clustering       = do_clustering        
        settings%feedback            = feedback             
        settings%precision_criterion = precision_criterion  
        settings%logzero             = logzero  
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
        settings%write_prior         = write_prior
        settings%maximise            = maximise
        settings%compression_factor  = compression_factor         
        settings%synchronous         = synchronous         
        settings%nDims               = nDims
        settings%nDerived            = nDerived

        settings%base_dir            = convert_c_string(base_dir)
        settings%file_root           = convert_c_string(file_root) 

        settings%seed                = seed

        allocate(settings%grade_frac(nGrade),settings%grade_dims(nGrade))
        settings%grade_frac = grade_frac
        settings%grade_dims = grade_dims

        allocate(settings%loglikes(n_nlives),settings%nlives(n_nlives))
        settings%loglikes = loglikes
        settings%nlives = nlives

        if(settings%write_paramnames) then
            params = default_params(settings%nDims,'theta','\theta')
            derived_params = default_params(settings%nDerived,'phi','\phi')
            call write_paramnames_file(settings,params,derived_params)
        end if

        call c_f_procpointer(c_loglikelihood_ptr, f_loglikelihood_ptr)
        call c_f_procpointer(c_prior_ptr, f_prior_ptr)
        call c_f_procpointer(c_cluster_ptr, f_cluster_ptr)
        call c_f_procpointer(c_dumper_ptr, f_dumper_ptr)

        call run_polychord(loglikelihood,prior_transform,dumper,cluster,settings,comm)

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

        function prior_transform(cube) result(theta)
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

        end function prior_transform

        subroutine dumper(live, dead, logweights, logZ, logZerr)
            implicit none
            real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
            real(dp), intent(in) :: logZ, logZerr

            integer(c_int) :: c_ndead, c_nlive, c_npars
            real(c_double) :: c_live(size(live,1),size(live,2)), c_dead(size(dead,1),size(dead,2)), c_logweights(size(logweights))
            real(c_double) :: c_logZ, c_logZerr

            c_npars = size(live,1)
            c_nlive = size(live,2)
            c_ndead = size(dead,2)
            c_live = live
            c_dead = dead
            c_logweights = logweights
            c_logZ = logZ
            c_logZerr = logZerr
            call f_dumper_ptr(c_ndead, c_nlive, c_npars, c_live, c_dead, c_logweights, c_logZ, c_logZerr)
        end subroutine dumper

        function cluster(points) result(cluster_list)
            implicit none
            real(dp), intent(in), dimension(:,:) :: points
            integer, dimension(size(points,2)) :: cluster_list

            integer(c_int) :: c_npoints, c_ndims
            real(c_double),  dimension(size(points,1),size(points,2)) :: c_points
            integer(c_int), dimension(size(points,2)) :: c_cluster_list
            c_ndims = size(points,1)
            c_npoints = size(c_cluster_list)
            c_points = points
            call f_cluster_ptr(c_points,c_cluster_list,c_ndims,c_npoints)
            cluster_list = c_cluster_list

        end function cluster


    end subroutine polychord_c_interface


    subroutine polychord_c_interface_ini(c_loglikelihood_ptr, c_setup_loglikelihood_ptr, input_file_c, comm)&
            bind(c,name='polychord_c_interface_ini')

        use iso_c_binding
        use utils_module,             only: STR_LENGTH, convert_c_string

        ! ~~~~~~~ Local Variable Declaration ~~~~~~~
        implicit none
        integer, intent(in) :: comm

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
            subroutine c_setup_loglikelihood() bind(c)
            end subroutine c_setup_loglikelihood
        end interface
        character(len=STR_LENGTH)     :: input_file     ! input file

        type(c_funptr), intent(in), value   :: c_loglikelihood_ptr
        type(c_funptr), intent(in), value   :: c_setup_loglikelihood_ptr
        character(len=1,kind=c_char), intent(in), dimension(STR_LENGTH) :: input_file_c

        procedure(c_loglikelihood), pointer :: f_loglikelihood_ptr
        procedure(c_setup_loglikelihood), pointer :: f_setup_loglikelihood_ptr

        input_file = convert_c_string(input_file_c)

        call c_f_procpointer(c_loglikelihood_ptr, f_loglikelihood_ptr)
        call c_f_procpointer(c_setup_loglikelihood_ptr, f_setup_loglikelihood_ptr)

        call run_polychord(loglikelihood, setup_loglikelihood, input_file, comm) 

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

        subroutine setup_loglikelihood(settings)
            use settings_module,          only: program_settings
            implicit none
            type(program_settings), intent(in)    :: settings
            call f_setup_loglikelihood_ptr()
        end subroutine setup_loglikelihood

    end subroutine polychord_c_interface_ini


end module interfaces_module
