!> This is a log-likelihood wrapper
!!
!! we define this, so we can easily switch in and out differing likelihoods
function loglikelihood(theta,phi,context)
    use example_likelihoods
    implicit none
    double precision, intent(in),  dimension(:) :: theta
    double precision, intent(out),  dimension(:) :: phi
    integer,          intent(in)                 :: context
    double precision :: loglikelihood
    !       These can be written into src/example_likelihoods.f90, and should
    !       have the interface:
    !
    !interface
    !    function loglikelihood(theta,phi,context)
    !        double precision, intent(in),  dimension(:) :: theta
    !        double precision, intent(out),  dimension(:) :: phi
    !        integer,          intent(in)                 :: context
    !        double precision :: loglikelihood
    !    end function
    !end interface
    !
    ! where theta are the input params, phi are any derived params and 
    ! context is an integer which can be useful as a C pointer
    !
    !
    ! Possible example likelihoods already written are:
    !       - gaussian_loglikelihood
    !       - gaussian_shell
    !       - rosenbrock_loglikelihood
    !       - himmelblau_loglikelihood
    !       - rastrigin_loglikelihood
    !       - eggbox_loglikelihood
    !       - gaussian_loglikelihood_corr
    !       - gaussian_loglikelihood_cluster
    !       - twin_gaussian_loglikelihood 
    !

    loglikelihood = gaussian_loglikelihood(theta,phi,context)

end function



!> This is the main driving routine of the nested sampling algorithm
program main

    ! ~~~~~~~ Loaded Modules ~~~~~~~

    use ini_module,               only: read_params,initialise_program
    use params_module,            only: add_parameter,param_type
    use priors_module
    use settings_module,          only: program_settings,initialise_settings
    use random_module,            only: initialise_random
    use nested_sampling_module,   only: NestedSampling
    use utils_module,             only: STR_LENGTH
    use abort_module,             only: halt_program
#ifdef MPI
    use mpi_module
#endif

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    ! Output of the program
    ! 1) log(evidence)
    ! 2) error(log(evidence))
    ! 3) ndead
    ! 4) number of likelihood calls
    ! 5) log(evidence) + log(prior volume)
    double precision, dimension(5)            :: output_info

    type(program_settings)                    :: settings  ! The program settings 
    type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

    character(len=STR_LENGTH)                 :: input_file     ! input file
    type(param_type),dimension(:),allocatable :: params         ! Parameter array
    type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

    ! Interface for the loglikelihood function
    interface
        function loglikelihood(theta,phi,context)
            double precision, intent(in),  dimension(:) :: theta
            double precision, intent(out),  dimension(:) :: phi
            integer,          intent(in)                 :: context
            double precision :: loglikelihood
        end function
    end interface



    ! ======= (1) Initialisation =======
    ! We need to initialise:
    ! a) mpi threads
    ! b) random number generator
    ! c) model
    ! d) program settings
    ! e) likelihoods


    ! ------- (1a) Initialise MPI threads -------------------
#ifdef MPI
    call initialise_mpi
#endif

    ! ------- (1b) Initialise random number generator -------
    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
    call initialise_random()



    ! ------- (1c) Define the parameters of the loglikelihood and the system settings -------
    ! The most convenient way to do this is via a .ini file. See examples
    ! provided for more details. The .ini file should be passed as a command
    ! line argument to PolyChord
    !
    ! Failing this, you can manually set up the parameters as detailed below
    if(iargc()==1) then

        ! ------ (1ci) read from .ini file ------
        call getarg(1,input_file) ! get the input filename from the command line

        ! Call this subroutine to read the input file and set up the settings and priors
        call read_params(trim(input_file),settings,params,derived_params)

    else if (iargc()==0) then

        ! ------ (1ci) manual setup ------
        ! Here we initialise the array params with all of the details we need
        allocate(params(0),derived_params(0))
        ! The argument to add_parameter are:
        ! 1) params:            the parameter array to add to
        ! 2) name:              paramname for getdist processing
        ! 3) latex:             latex name for getdist processing
        ! 4) speed:             The speed of this parameter (lower => slower)
        ! 5) prior_type:        what kind of prior it is
        ! 6) prior_block:       what other parameters are associated with it
        ! 7) prior_parameters:  parameters of the prior
        !                  array   name     latex     speed  prior_type   prior_block prior_parameters
        call add_parameter(params,'param1','\theta_1',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param2','\theta_2',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param3','\theta_3',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param4','\theta_4',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param5','\theta_5',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param6','\theta_6',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param7','\theta_7',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param8','\theta_8',1,     uniform_type,1,          [ 0d0 , 1d0 ])

        call add_parameter(derived_params,'param9','r')

        ! Now initialise the rest of the system settings
        settings%nlive         = 500
        settings%num_repeats   = 8
        settings%do_clustering = .true.

        settings%base_dir      = 'chains'
        settings%file_root     = 'gaussian'

        settings%calculate_posterior = .true.
        
        settings%write_resume  = .true.
        settings%read_resume   = .false.
        settings%write_live    = .true.

        settings%feedback      = 1
        settings%update_resume = settings%nlive

        settings%thin_posterior= 0d0
        allocate(settings%grade_frac(2)) 
        settings%grade_frac=[1d0,1d-1]

    else
        call halt_program('PolyChord should be called with at most one argument, the input file')
    end if

    ! Initialise the program
    call initialise_program(settings,priors,params,derived_params)

    ! ======= (2) Perform Nested Sampling =======
    ! Call the nested sampling algorithm on our chosen likelihood and priors
#ifdef MPI
    output_info = NestedSampling(loglikelihood,priors,settings,MPI_COMM_WORLD) 
#else
    output_info = NestedSampling(loglikelihood,priors,settings,0) 
#endif


    ! ======= (3) De-initialise =======

    ! Finish off all of the threads
#ifdef MPI
    call finalise_mpi
#endif


end program main
