program galileo

    use nested_sampling_module
    use likelihoods_module
    use model_module
    use settings_module
    use random_module

    implicit none

    type(program_settings) :: settings
    type(model)            :: M

    integer, parameter :: nDims = 2

    double precision, dimension(nDims) :: coords

    !* Initialise the program settings
    settings%nlive=1024
    !* Initialise the model
    M%nDims = nDims
    M%nDerived = nDims+1
    M%nTotal = nDims+1

    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
    call initialise_random()
    
    ! Call the nested sampling algorithm on our chosen likelihood and priors
    call NestedSampling(M,settings)


    ! De-initialise the random number generator 
    call deinitialise_random

end program galileo
