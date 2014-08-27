!> This is the main driving routine of the nested sampling algorithm
program main

    ! ~~~~~~~ Loaded Modules ~~~~~~~

    use model_module,           only: model, allocate_live_indices, allocate_prior_arrays, set_up_prior_indices
    
    use settings_module,        only: program_settings
    use random_module,          only: initialise_random, deinitialise_random

    use chordal_module,      only: ChordalSampling
    use evidence_module
    use example_likelihoods
    use feedback_module
#ifdef MPI
    use mpi_module
    use chordal_module,                  only: ChordalSampling
    use nested_sampling_parallel_module, only: NestedSamplingP
    use nested_sampling_linear_module,   only: NestedSamplingL
#else
    use nested_sampling_linear_module, only: NestedSamplingL
#endif
    use utils_module, only: TwoPi

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    type(program_settings) :: settings  ! The program settings 
    type(model)            :: M         ! The model details




    ! ======= (1) Initialisation =======
    ! We need to initialise:
    ! a) mpi threads
    ! b) random number generator
    ! c) model
    ! d) program settings


    ! ------- (1a) Initialise MPI threads -------------------
#ifdef MPI
    call mpi_initialise()
#endif

    ! ------- (1b) Initialise random number generator -------
    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
    call initialise_random()


    ! ------- (1c) Initialise the model -------
    ! (i) Assign the likelihood function
    M%loglikelihood => gaussian_loglikelihood
    !M%loglikelihood => himmelblau_loglikelihood
    !M%loglikelihood => rastrigin_loglikelihood
    !M%loglikelihood => rosenbrock_loglikelihood
    !M%loglikelihood => eggbox_loglikelihood

    M%calc_derived => zero_derived
    

    ! (ii) Set the dimensionality
    M%nDims=20                 ! Dimensionality of the space
    M%nDerived = 0             ! Assign the number of derived parameters
    ! There are two derived parameters:
    ! 1) the number of likelihood evaluations required for the calculation of the
    ! new point
    ! 2) the smallest plausible chord length
        

    ! (iii) Assign the priors
    !       
    M%uniform_num = M%nDims

    call allocate_live_indices(M)

    ! (iv) Assign the grades

    !M%grade=1
    M%grade(:4) = 1
    M%grade(5:6) = 2
    M%grade(7:) = 3

    ! (v) Set up prior arrays
    call allocate_prior_arrays(M)

    !       - settings of priors
    M%uniform_params(:,1) = 0.45
    M%uniform_params(:,2) = 0.55
    !M%uniform_params(:,1) = -5d0
    !M%uniform_params(:,2) =  5d0

    call set_up_prior_indices(M)




    ! ------- (1d) Initialise the program settings -------
    settings%file_root            =  'chains/test'           !file root
    settings%nlive                =  100*M%nDims             !number of live points
    settings%sampler              => ChordalSampling         !Sampler choice
    settings%evidence_calculator  => KeetonEvidence          !evidence calculator
    settings%feedback             =  1                       !degree of feedback
    settings%precision_criterion  =  1d-1                    !degree of precision in answer
    settings%max_ndead            =  -1                      !maximum number of samples
    settings%nmax_posterior       = 100000                   !max number of posterior points
    settings%minimum_weight       = 1d-50                    !minimum weight of the posterior points
    settings%calculate_posterior  = .true.                   !calculate the posterior (slows things down at the end of the run)
    settings%write_resume         = .true.                   !whether or not to write resume files
    settings%update_resume        = settings%nlive*100       !How often to update the resume files

    allocate(settings%num_chords(3))
    settings%num_chords(1) = count(M%grade==1)
    settings%num_chords(2) = count(M%grade==2)
    settings%num_chords(3) = count(M%grade==3)
    !allocate(settings%num_chords(1))
    !settings%num_chords=20


    ! ======= (2) Perform Nested Sampling =======
    ! Call the nested sampling algorithm on our chosen likelihood and priors
#ifdef MPI
    if (mpi_size()>1) then
        call NestedSamplingP(M,settings)
    else
        call NestedSamplingL(M,settings) 
    end if
#else
    call NestedSamplingL(M,settings) 
#endif 



    ! ======= (3) De-initialise =======
    ! De-initialise the random number generator 
    call deinitialise_random()

#ifdef MPI
    call mpi_finalise()
#endif


end program main
