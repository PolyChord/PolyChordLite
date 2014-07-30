program galileo

    use nested_sampling_module, only: NestedSampling
    use model_module,           only: model, initialise_model
    use settings_module,        only: program_settings
    use random_module,          only: initialise_random, deinitialise_random

    use test_sampler_module
!    use evidence
    use example_likelihoods

    implicit none

    type(program_settings) :: settings
    type(model)            :: M

    integer, parameter :: nDims = 2

    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
    call initialise_random()

    !* Initialise the model
    M%loglikelihood => gaussian_loglikelihood  ! Assign the likelihood 
    M%nDims = nDims  ! Assign the dimensionality
                     ! Assign the priors
    M%nDerived = 1   ! Assign the number of derived parameters
    
    call initialise_model(M)   ! Configure the rest


    !* Initialise the program settings
    settings%nlive=1024  ! number of live points

    settings%sampler              => SphericalCenterSampling !Sampler choice
!    settings%evidence_calculator  => KeetonCalculator        !evidence calculator


    ! Call the nested sampling algorithm on our chosen likelihood and priors
    call NestedSampling(M,settings)

    ! De-initialise the random number generator 
    call deinitialise_random

end program galileo
