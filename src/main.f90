!> This is the main driving routine of the nested sampling algorithm
program main

    ! ~~~~~~~ Loaded Modules ~~~~~~~

    use priors_module
    use settings_module
    use random_module,          only: initialise_random, deinitialise_random
    use example_likelihoods
    use feedback_module
    use grades_module,          only: allocate_grades,calc_num_babies
    use nested_sampling_module,   only: NestedSampling
    use mpi_module

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    ! Output of the program
    ! 1) log(evidence)
    ! 2) error(log(evidence))
    ! 3) ndead
    ! 4) number of likelihood calls
    double precision, dimension(5) :: output_info

    type(program_settings)    :: settings  ! The program settings 
    type(prior), dimension(1) :: priors

    pointer loglikelihood
    double precision :: loglike

    double precision, allocatable, dimension(:) :: theta
    double precision, allocatable, dimension(:) :: phi

    double precision, allocatable, dimension(:) :: minimums 
    double precision, allocatable, dimension(:) :: maximums
    integer, allocatable, dimension(:) :: hypercube_indices
    integer, allocatable, dimension(:) :: physical_indices
    integer :: i



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


    ! ------- (1a) Initialise MPI threads -------------------
    call MPI_INIT(mpierror)

    ! ------- (1b) Initialise random number generator -------
    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
    call initialise_random()


    ! ------- (1c) Initialise the model -------
    ! (i) Choose the loglikelihood
    !       Possible example likelihoods:
    !       - gaussian_loglikelihood
    !       - gaussian_shell
    !       - rosenbrock_loglikelihood
    !       - himmelblau_loglikelihood
    !       - rastrigin_loglikelihood
    !       - eggbox_loglikelihood
    !       - gaussian_loglikelihood_corr
    !       - gaussian_loglikelihood_cluster
    !       - twin_gaussian_loglikelihood 
    loglikelihood => rastrigin_loglikelihood!gaussian_loglikelihood_corr!rastrigin_loglikelihood

    ! (ii) Set the dimensionality
    settings%nDims= 2                ! Dimensionality of the space
    settings%nDerived = 0             ! Assign the number of derived parameters

    ! (iii) Assign the priors
    call allocate_indices(settings)

    ! (v) Set up priors
    allocate(minimums(settings%nDims))
    allocate(maximums(settings%nDims))
    allocate(physical_indices(settings%nDims))
    allocate(hypercube_indices(settings%nDims))

    !minimums=0.5-1d-2*5
    !maximums=0.5+1d-2*5
    !minimums=0.5-1d-2*20
    !maximums=0.5+1d-2*20
    minimums = -5
    maximums =  5

    do i=1,settings%nDims
        physical_indices(i)  = i
        hypercube_indices(i) = i
    end do

    call initialise_uniform(priors(1),hypercube_indices,physical_indices,minimums,maximums)




    ! ------- (1d) Initialise the program settings -------
    settings%nlive                = 5000!25*settings%nDims        !number of live points
    settings%num_repeats          = 5                        !Number of chords to draw

    settings%num_babies           = settings%nDims*settings%num_repeats
    settings%nstack               = settings%nlive*settings%num_babies*2
    settings%file_root            =  'chains/test'           !file root
    settings%feedback             = 1                        !degree of feedback

    ! stopping criteria
    settings%precision_criterion  =  1d-3                    !degree of precision in answer 
    settings%max_ndead            = -1                       !maximum number of samples 
    ! posterior calculation
    settings%nmax_posterior       = 100000                   !max number of posterior points
    settings%calculate_posterior  = .false.                  !calculate the posterior (slows things down at the end of the run)

    ! reading and writing
    settings%read_resume          = .false.                  !whether or not to resume from file
    settings%write_resume         = .false.                  !whether or not to write resume files
    settings%update_resume        = settings%nlive           !How often to update the resume files
    settings%write_live           = .true.                   !write out the physical live points?

    settings%do_clustering = .true.
    settings%SNN_k = 10 !settings%nDims
    settings%SNN_kt = 5

    settings%ncluster = 100
    settings%nclustertot = 400


    ! Initialise the loglikelihood
    allocate(theta(settings%nDims),phi(settings%nDerived))
    loglike = loglikelihood(theta,phi,0)

    ! Sort out the grades
    !call allocate_grades(settings%grades,(/1,1,1,1,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4/) ) 
    !settings%grades%num_repeats(2)= 5
    !settings%grades%num_repeats(4)= 5
    !settings%num_babies = calc_num_babies(settings%grades)
    !settings%nstack               = settings%nlive*settings%num_babies*2
    settings%do_grades=.false.
    settings%do_timing=.false.
    !nest_settings%thin_posterior = max(0d0,&
    !    nest_settings%grades%num_repeats(1)*nest_settings%grades%grade_nDims(1)/&
    !    (sum(nest_settings%grades%num_repeats*nest_settings%grades%grade_nDims)+0d0)&
    !    )


    ! ======= (2) Perform Nested Sampling =======
    ! Call the nested sampling algorithm on our chosen likelihood and priors

    do i=1,1
        output_info = NestedSampling(loglikelihood,priors,settings,MPI_COMM_WORLD) 
        write(*,'(2E17.8)') output_info(5),output_info(2)
    end do



    ! ======= (3) De-initialise =======
    ! De-initialise the random number generator 
    call deinitialise_random()

    call MPI_FINALIZE(mpierror)


end program main
