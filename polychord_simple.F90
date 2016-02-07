program demo_gaussian

    ! ~~~~~~~ Loaded Modules ~~~~~~~
    use utils_module,             only: dp
    use settings_module,          only: program_settings
    use interfaces_module,        only: run_polychord

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
    settings%write_stats   = .false.

    settings%equals        = .false.
    settings%posteriors    = .false.
    settings%cluster_posteriors = .false.

    settings%feedback      = 1
    settings%update_files  = settings%nlive

    settings%boost_posterior= 5.0_dp

    call run_polychord(loglikelihood, prior, settings) 


contains

    function loglikelihood(theta, phi)
        implicit none
        real(dp), intent(in), dimension(:)  :: theta
        real(dp), intent(out), dimension(:) :: phi
        real(dp) :: loglikelihood

        real(dp), parameter :: mu = 0.0_dp
        real(dp), parameter :: sigma = 1.0_dp
        real(dp), parameter :: pi = atan(1.0_dp)*4

        real(dp) :: radius_squared

        integer :: nDims

        nDims = size(theta)

        radius_squared = sum((theta-mu)*(theta-mu))

        loglikelihood = -log(2*pi*sigma*sigma)*nDims/2.0
        loglikelihood = loglikelihood - radius_squared/2/sigma/sigma

        phi(1) = radius_squared

    end function loglikelihood

    function prior(cube) result(theta)
        implicit none
        real(dp), intent(in), dimension(:) :: cube
        real(dp), dimension(size(cube))    :: theta

        real(dp), parameter :: xmin = -10.0_dp
        real(dp), parameter :: xmax = +10.0_dp

        theta = xmin + cube*(xmax-xmin)

    end function prior


end program demo_gaussian
