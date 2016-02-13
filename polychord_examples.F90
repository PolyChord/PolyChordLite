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

    else
        call halt_program('PolyChord should be called with at most one argument, the input file')
    end if

    call create_priors(priors,params,settings)

    call setup_loglikelihood(settings)
    call run_polychord(loglikelihood,prior_wrapper, settings) 

contains

    function prior_wrapper(cube) result(theta)
        implicit none
        double precision, intent(in), dimension(:) :: cube
        double precision, dimension(size(cube))    :: theta
        theta = hypercube_to_physical(cube,priors)
    end function prior_wrapper

end program PolyChord
