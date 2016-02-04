!> This is the main driving routine of the nested sampling algorithm
program PolyChord

    ! ~~~~~~~ Loaded Modules ~~~~~~~
    use ini_module,               only: read_params
    use params_module,            only: add_parameter,param_type
    use priors_module
    use settings_module,          only: program_settings
    use interfaces_module,        only: run_polychord
    use utils_module,             only: STR_LENGTH
    use abort_module,             only: halt_program
    use loglikelihood_module,     only: loglikelihood, setup_loglikelihood

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    type(program_settings)                    :: settings  ! The program settings 
    type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

    character(len=STR_LENGTH)                 :: input_file     ! input file
    type(param_type),dimension(:),allocatable :: params         ! Parameter array
    type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

    ! -------  Define the parameters of the loglikelihood and the system settings -------
    ! The most convenient way to do this is via a .ini file. See examples
    ! provided for more details. The .ini file should be passed as a command
    ! line argument to PolyChord
    !
    ! Failing this, you can manually set up the parameters as detailed below
    !
    if(iargc()==1) then

        ! ------ read from .ini file ------
        ! Put the command line provided file name into input_file
        call getarg(1,input_file) 

        ! Call this subroutine to read the input file and set up the settings and priors
        call read_params(trim(input_file),settings,params,derived_params)

    else if (iargc()==0) then

        ! ------ manual setup ------
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
        call add_parameter(params,'param1','\theta_{1}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param2','\theta_{2}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param3','\theta_{3}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param4','\theta_{4}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param5','\theta_{5}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param6','\theta_{6}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param7','\theta_{7}',1,     uniform_type,1,          [ 0d0 , 1d0 ])
        call add_parameter(params,'param8','\theta_{8}',1,     uniform_type,1,          [ 0d0 , 1d0 ])

        call add_parameter(derived_params,'radius','r')
        call add_parameter(derived_params,'logVolume','\log X')

        ! Now initialise the rest of the system settings
        settings%nDims         = size(params)
        settings%nDerived      = size(derived_params)
        settings%nlive         = 500
        settings%num_repeats   = 16
        settings%do_clustering = .false.

        settings%base_dir      = 'chains'
        settings%file_root     = 'test'

        settings%write_resume  = .false.
        settings%read_resume   = .false.
        settings%write_live    = .false.
        settings%write_stats   = .false.

        settings%equals        = .false.
        settings%posteriors    = .false.
        settings%cluster_posteriors = .false.

        settings%feedback      = 1
        settings%update_files  = settings%nlive

        settings%boost_posterior= 5d0
        allocate(settings%grade_frac(1)) 
        settings%grade_frac=[1d0]

    else
        call halt_program('PolyChord should be called with at most one argument, the input file')
    end if

    call create_priors(priors,params,settings)

    call run_polychord(loglikelihood,prior_wrapper, setup_loglikelihood, settings) 

contains

    function prior_wrapper(cube) result(theta)
        implicit none
        double precision, intent(in), dimension(:) :: cube
        double precision, dimension(size(cube))    :: theta
        theta = hypercube_to_physical(cube,priors)
    end function prior_wrapper

end program PolyChord
