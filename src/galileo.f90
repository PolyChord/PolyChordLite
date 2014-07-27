program galileo

    use nested_sampling_module
    !use example_likelihoods
    use model_module
    use settings_module
    use random_module

    implicit none

    type(program_settings) settings
    type(model_details) priors

    double precision, dimension(20) :: nhat


    call initialise_random(1)
    
    call random_direction(nhat,20)

    call deinitialise_random


    ! Initialise the program settings

    ! Initialise the priors


    ! Call the nested sampling algorithm on our chosen likelihood and priors
    !call NestedSampling(gaussian_loglikelihood,priors,settings)


end program galileo
