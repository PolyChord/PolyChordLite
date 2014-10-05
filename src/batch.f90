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
    use nested_sampling_parallel_module, only: NestedSamplingP
#endif
    use nested_sampling_linear_module,   only: NestedSamplingL

    use utils_module, only: STR_LENGTH

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    ! The number of samples to draw
    integer, parameter :: num_samples = 100000  
    ! How often to update from all the MPI cores
    integer, parameter :: update = 10
    ! The name of the file
    character(STR_LENGTH) :: out_root='corr8_8_1.dat'

    ! Output of the program
    ! 1) log(evidence)
    ! 2) error(log(evidence))
    ! 3) ndead
    ! 4) number of likelihood calls
    double precision, dimension(5,update) :: output_info_local
    double precision, dimension(:,:),allocatable :: output_info

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
    integer :: i,j
    integer :: i_err

    integer :: myrank
    integer :: nprocs


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
    loglikelihood => gaussian_loglikelihood_corr

    ! (ii) Set the dimensionality
    settings%nDims= 8                 ! Dimensionality of the space
    settings%nDerived = 0             ! Assign the number of derived parameters

    ! (iii) Assign the priors
    call allocate_indices(settings)

    ! (v) Set up priors
    allocate(minimums(settings%nDims))
    allocate(maximums(settings%nDims))
    allocate(physical_indices(settings%nDims))
    allocate(hypercube_indices(settings%nDims))

    minimums=0.5-1d-2*5
    maximums=0.5+1d-2*5
    !minimums=-5
    !maximums=5

    do i=1,settings%nDims
        physical_indices(i)  = i
        hypercube_indices(i) = i
    end do

    call initialise_uniform(priors(1),hypercube_indices,physical_indices,minimums,maximums)




    ! ------- (1d) Initialise the program settings -------
    settings%nlive                = 25*settings%nDims        !number of live points
    settings%chain_length         = settings%nDims           !Number of chords to draw

    settings%nstack               =  settings%nlive*10       !number of points in the 'stack'
    settings%file_root            =  'chains/test'           !file root
    settings%sampler              => SliceSampling           !Sampler choice
    settings%get_nhat             => Adaptive_Parallel       !Direction choice
    settings%process_live_points  => get_live_coordinates    !no processing of live points needed

    settings%evidence_calculator  => KeetonEvidence          !evidence calculator
    settings%feedback             =  -1                      !degree of feedback

    ! stopping criteria
    settings%precision_criterion  =  1d-1                    !degree of precision in answer
    settings%max_ndead            =  100000                  !maximum number of samples

    ! posterior calculation
    settings%nmax_posterior       = 100000                   !max number of posterior points
    settings%minimum_weight       = 1d-6                     !minimum weight of the posterior points
    settings%calculate_posterior  = .false.                  !calculate the posterior (slows things down at the end of the run)

    ! reading and writing
    settings%read_resume          = .false.                  !whether or not to resume from file
    settings%write_resume         = .false.                  !whether or not to write resume files
    settings%update_resume        = settings%nlive           !How often to update the resume files
    settings%write_live           = .false.                  !write out the physical live points?
    settings%save_all             = .false.                  !Save all the dead points?

    ! Evidence inference
    settings%infer_evidence       = .false.
    settings%evidence_samples     = 100000 


    ! Initialise the loglikelihood
    allocate(theta(settings%nDims),phi(settings%nDerived))
    loglike = loglikelihood(theta,phi,0)


    ! ======= (2) Perform Nested Sampling =======
    ! Call the nested sampling algorithm on our chosen likelihood and priors

    ! Since we're running each of the nested sampling algorithms in linear mode,
    ! one should set the stack size to the number of live points
    settings%nstack=settings%nlive

    ! Get the processor ID...
    myrank = mpi_rank()
    ! ... and the total number of processors
    nprocs = mpi_size()

    ! Allocate the collecting array on root 
    if(myrank==root) allocate(output_info(5,nprocs*update))

    do i=1,num_samples

        ! Write to the local output_info array the results of a nested sampling run
        output_info_local(:,mod(i-1,update)+1) = NestedSamplingL(loglikelihood,priors,settings) 

        if(mod(i-1,update)+1 == update) then
            ! Gather all of the evidences every update iterations
            call MPI_GATHER(                 &  
                output_info_local,           & ! sending array
                5*update,                    & ! number of elements to be sent
                MPI_DOUBLE_PRECISION,        & ! type of element to be sent
                output_info,                 & ! recieving array
                5*update,                    & ! number of elements to be recieved from each node
                MPI_DOUBLE_PRECISION,        & ! type of element recieved
                root,                        & ! root node address
                MPI_COMM_WORLD,              & ! communication info
                mpierror)                      ! error (from module mpi_module)

            if(myrank==root) then
                write(*,*) i*nprocs
                ! Open the output file for writing, appending to the end
                open(111,file=trim(out_root), action='write', position="append", iostat=i_err)
                ! Write all of the evidences we've recieved 
                ! 5 contains the 'check' i.e. log(evidence*volume)
                ! 2 contains the 'error' i.e. stdev in log(evidence)
                do j=1,update*nprocs
                    if(myrank==root) write(111,'(2E17.8)') output_info(5,j), output_info(2,j)
                end do
                ! Close the output file
                close(111)
            end if
        end if

    end do





    ! ======= (3) De-initialise =======
    ! De-initialise the random number generator 
    call deinitialise_random()

#ifdef MPI
    call mpi_finalise()
#endif


end program main
