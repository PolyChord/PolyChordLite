!> This is the main driving routine of the nested sampling algorithm
program main

    ! ~~~~~~~ Loaded Modules ~~~~~~~

    use priors_module
    use settings_module,        only: program_settings,allocate_indices
    use random_module,          only: initialise_random, deinitialise_random

    use chordal_module,         only: SliceSampling, SliceSampling_Graded, &
                                      HitAndRun, Adaptive_Parallel, &
                                      no_processing, get_live_coordinates
    use evidence_module,        only: KeetonEvidence
    use example_likelihoods
    use feedback_module
#ifdef MPI
    use mpi_module
#endif
    use nested_sampling_linear_module,   only: NestedSamplingL

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    type(program_settings)    :: settings  ! The program settings 
    type(prior), dimension(1) :: priors

    double precision, dimension(5) :: output_info

    pointer loglikelihood

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
#ifdef MPI
    call mpi_initialise()
#endif

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
    loglikelihood => planck_loglikelihood

    ! (ii) Set the dimensionality
    settings%nDims=20                 ! Dimensionality of the space
    settings%nDerived = 0             ! Assign the number of derived parameters

    ! (iii) Assign the priors
    call allocate_indices(settings)

    ! (v) Set up priors
    allocate(minimums(settings%nDims))
    allocate(maximums(settings%nDims))
    allocate(physical_indices(settings%nDims))
    allocate(hypercube_indices(settings%nDims))

    ! Standard minimums and maximums from a normal run
    !minimums( 1) = 0.019  !omegabh2  
    !minimums( 2) = 0.095  !omegach2  
    !minimums( 3) = 1.03   !theta     
    !minimums( 4) = 0.01   !tau       
    !minimums( 5) = 2.5    !logA        
    !minimums( 6) = 0.885  !ns          
    !minimums( 7) = 0      !aps100    
    !minimums( 8) = 0      !aps143    
    !minimums( 9) = 0      !aps217    
    !minimums(10) = 0      !acib143   
    !minimums(11) = 0      !acib217   
    !minimums(12) = 0      !asz143    
    !minimums(13) = 0.0    !psr       
    !minimums(14) = 0.0    !cibr      
    !minimums(15) = -2     !ncib      
    !minimums(16) = 0.98   !cal0      
    !minimums(17) = 0.95   !cal2      
    !minimums(18) = 0      !xi        
    !minimums(19) = 0      !aksz      
    !minimums(20) = -20    !bm_1_1    
     
    !maximums( 1) = 0.025 !omegabh2  
    !maximums( 2) = 0.145 !omegach2  
    !maximums( 3) = 1.05  !theta     
    !maximums( 4) = 0.4   !tau       
    !maximums( 5) = 3.7   !logA        
    !maximums( 6) = 1.04  !ns          
    !maximums( 7) = 360   !aps100    
    !maximums( 8) = 270   !aps143    
    !maximums( 9) = 450   !aps217    
    !maximums(10) = 20    !acib143   
    !maximums(11) = 80    !acib217   
    !maximums(12) = 10    !asz143    
    !maximums(13) = 1.0   !psr       
    !maximums(14) = 1.0   !cibr      
    !maximums(15) = 2     !ncib      
    !maximums(16) = 1.02  !cal0      
    !maximums(17) = 1.05  !cal2      
    !maximums(18) = 1     !xi        
    !maximums(19) = 10    !aksz      
    !maximums(20) = 20    !bm_1_1    

    ! extended minimums and maximums for accurate evidence calculation
    minimums( 1) = 0.019  !omegabh2  
    minimums( 2) = 0.095  !omegach2  
    minimums( 3) = 1.03   !theta     
    minimums( 4) = -0.1   !tau       
    minimums( 5) = 2.5    !logA        
    minimums( 6) = 0.885  !ns          
    minimums( 7) = 0      !aps100    
    minimums( 8) = 0      !aps143    
    minimums( 9) = 0      !aps217    
    minimums(10) = -10    !acib143   
    minimums(11) = 0      !acib217   
    minimums(12) = -4     !asz143    
    minimums(13) = 0.0    !psr       
    minimums(14) = -0.2   !cibr      
    minimums(15) = -2     !ncib      
    minimums(16) = 0.98   !cal0      
    minimums(17) = 0.95   !cal2      
    minimums(18) = -1     !xi        
    minimums(19) = -5     !aksz      
    minimums(20) = -20    !bm_1_1    

    maximums( 1) = 0.025 !omegabh2  
    maximums( 2) = 0.145 !omegach2  
    maximums( 3) = 1.05  !theta     
    maximums( 4) = 0.4   !tau       
    maximums( 5) = 3.7   !logA        
    maximums( 6) = 1.04  !ns          
    maximums( 7) = 360   !aps100    
    maximums( 8) = 270   !aps143    
    maximums( 9) = 450   !aps217    
    maximums(10) = 25    !acib143   
    maximums(11) = 80    !acib217   
    maximums(12) = 12    !asz143    
    maximums(13) = 2.0   !psr       
    maximums(14) = 1.5   !cibr      
    maximums(15) = 2     !ncib      
    maximums(16) = 1.02  !cal0      
    maximums(17) = 1.05  !cal2      
    maximums(18) = 2     !xi        
    maximums(19) = 20    !aksz      
    maximums(20) = 20    !bm_1_1    

    do i=1,settings%nDims
        physical_indices(i)  = i
        hypercube_indices(i) = i
    end do

    call initialise_uniform(priors(1),hypercube_indices,physical_indices,minimums,maximums)




    ! ------- (1d) Initialise the program settings -------
    settings%nlive                = 500                      !number of live points
    settings%chain_length         = 1                        !Number of chords to draw (after each randomisation)
    settings%num_reflections      = 1                        !Number of randomisations to choose, 4 seems fine in most cases

    settings%read_resume          = .false.                  !whether or not to resume from file


    settings%nstack               =  settings%nlive*10       !number of points in the 'stack'
    settings%file_root            =  'chains/test'           !file root
    settings%sampler              => SliceSampling_Graded    !Sampler choice
    settings%get_nhat             => Adaptive_Parallel       !Direction choice
    settings%process_live_points  => get_live_coordinates    !no processing of live points needed

    settings%evidence_calculator  => KeetonEvidence          !evidence calculator
    settings%feedback             =  1                       !degree of feedback
    settings%precision_criterion  =  1d-2                    !degree of precision in answer
    settings%max_ndead            =  -1                      !maximum number of samples
    settings%nmax_posterior       = 100000                   !max number of posterior points
    settings%minimum_weight       = 1d-5                     !minimum weight of the posterior points
    settings%calculate_posterior  = .false.                  !calculate the posterior (slows things down at the end of the run)
    settings%write_resume         = .false.                  !whether or not to write resume files
    settings%update_resume        = settings%nlive           !How often to update the resume files
    settings%save_all             = .false.                  !Save all the dead points?

    settings%grade(:4)=1
    settings%grade(5:6)=2
    settings%grade(7:)=4

    allocate(settings%chain_lengths(maxval(settings%grade)))
    settings%chain_lengths(1) = 4
    settings%chain_lengths(2) = 10
    settings%chain_lengths(3) = 0
    settings%chain_lengths(4) = 10

    ! ======= (2) Perform Nested Sampling =======
    ! Call the nested sampling algorithm on our chosen likelihood and priors

#ifdef MPI
    if (mpi_size()>1) then
!        output_info = NestedSamplingP(loglikelihood,priors,settings)
    else
        settings%nstack=settings%nlive
        output_info = NestedSamplingL(loglikelihood,priors,settings) 
    end if
#else
    settings%nstack=settings%nlive
    output_info = NestedSamplingL(loglikelihood,priors,settings) 
#endif 



    ! ======= (3) De-initialise =======
    ! De-initialise the random number generator 
    call deinitialise_random()

#ifdef MPI
    call mpi_finalise()
#endif


end program main
