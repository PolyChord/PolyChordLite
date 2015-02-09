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
        loglikelihood = rastrigin_loglikelihood(theta,phi,context)
end function



!> This is the main driving routine of the nested sampling algorithm
program main

    ! ~~~~~~~ Loaded Modules ~~~~~~~

    use priors_module
    use settings_module,          only: program_settings,initialise_settings
    use random_module,            only: initialise_random
    use nested_sampling_module,   only: NestedSampling
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
    double precision, dimension(5) :: output_info

    type(program_settings)    :: settings  ! The program settings 
    type(prior), dimension(1) :: priors    ! The details of the priors

    !pointer loglikelihood                  ! Pointer to a loglikelihood function
    ! (this just makes assigning them easier, one can just pass your chosen function to the nested sampling algorithm)



    ! Temporary variables for initialising priors
    double precision, allocatable, dimension(:) :: minimums 
    double precision, allocatable, dimension(:) :: maximums
    integer, allocatable, dimension(:) :: hypercube_indices
    integer, allocatable, dimension(:) :: physical_indices

    ! Iterator
    integer :: i

    ! Temporary variables for initialising loglikelihoods
    double precision :: loglike
    double precision, allocatable, dimension(:) :: theta
    double precision, allocatable, dimension(:) :: phi



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
    call MPI_INIT(mpierror)
#endif

    ! ------- (1b) Initialise random number generator -------
    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
    call initialise_random()


    ! ------- (1c) Initialise the model -------
    ! (i) Choose the loglikelihood at the top of this file

    ! (ii) Set the dimensionality
    settings%nDims= 2                 ! Dimensionality of the space
    settings%nDerived = 0             ! Assign the number of derived parameters

    ! (iii) Set up priors

    ! If just using uniform priors, then all you need to alter are the
    ! 'minimums' and 'maximums' array. If they're more complicated, then these
    ! are initialised in a very similar way, but using the initialise_<prior choice> 
    ! function in src/priors.f90
    !
    ! Currently have the possibility of
    ! 1) Separable Uniform priors
    ! 2) Separable Gaussian priors
    ! 3) Separable Log-Uniform priors
    ! 4) Sorted Uniform Priors

    allocate(                                &
        minimums(settings%nDims),            &
        maximums(settings%nDims),            &
        physical_indices(settings%nDims),    &
        hypercube_indices(settings%nDims)    &
        )

    minimums=-2
    maximums= 2
    hypercube_indices  = [ (i,i=1,settings%nDims) ]
    physical_indices  = hypercube_indices

    call initialise_uniform(priors(1),hypercube_indices,physical_indices,minimums,maximums)




    ! ------- (1d) Initialise the program settings -------
    settings%nlive                = 500                      !number of live points
    settings%num_repeats          = 1                        !Number of chords to draw (this is multiplied by nDims)

    settings%do_clustering        = .true.                   !whether or not to do clustering
    settings%ncluster             = 30                       !maximum number of clusters + 1 (memory allocation purposes)

    settings%feedback             = 1                        !degree of feedback

    settings%calculate_posterior  = .true.                   !calculate the posterior (slows things down)
    settings%thin_posterior       = 1                        !Factor by which the posterior file should be thinned
                                                             ! 0 uses just live points,
                                                             ! 1 uses all inter-chain points

    ! reading and writing
    settings%file_root            = 'rastrigin'              !file root
    settings%base_dir             = 'chains'                 !directory to put chains in
    settings%read_resume          = .false.                  !whether or not to resume from file
    settings%write_resume         = .true.                   !whether or not to write resume files
    settings%update_resume        = settings%nlive           !How often to update the resume files
    settings%write_live           = .true.                   !write out the physical live points?


    ! Calculate all of the rest of the settings
    call initialise_settings(settings)   



    ! ------- (1e) Initialise loglikelihood -----------------
    ! This is only needed for a few things (e.g. generating a random correlated gaussian)
    allocate(theta(settings%nDims),phi(settings%nDerived))
    theta   = 0d0
    loglike = loglikelihood(theta,phi,0)






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
    call MPI_FINALIZE(mpierror)
#endif


end program main
