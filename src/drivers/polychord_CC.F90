program PolyChord

    ! ~~~~~~~ Loaded Modules ~~~~~~~
    use utils_module,             only: dp
    use settings_module,          only: program_settings
    use interfaces_module,        only: run_polychord
    use loglikelihood_module,     only: loglikelihood,prior,setup_loglikelihood

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none
    type(program_settings)                    :: settings  ! The program settings 

    ! Initialise the system settings
    settings%nDims         = 20
    settings%nDerived      = 1
    settings%nlive         = 500
    settings%num_repeats   = settings%nDims*5
    settings%do_clustering = .false.

    settings%precision_criterion = 1d-3

    settings%base_dir      = 'chains'
    settings%file_root     = 'demo_gaussian'

    settings%write_resume  = .false.
    settings%read_resume   = .false.
    settings%write_live    = .true.
    settings%write_dead    = .false.
    settings%write_stats   = .false.

    settings%equals        = .false.
    settings%posteriors    = .false.
    settings%cluster_posteriors = .false.

    settings%feedback      = 1
    settings%update_files  = settings%nlive

    settings%boost_posterior= 5.0_dp

    call setup_loglikelihood(settings)
    call run_polychord(loglikelihood, prior, settings) 


end program PolyChord
