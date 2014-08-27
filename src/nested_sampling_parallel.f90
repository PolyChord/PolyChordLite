module nested_sampling_parallel_module
    implicit none

    integer,parameter :: flag_live_waiting = 0
    integer,parameter :: flag_incubator_running = -1
    integer,parameter :: flag_incubator_blank   = -2

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSamplingP(M,settings)
        use mpi_module
        use model_module,      only: model
        use utils_module,      only: logzero,loginf,DBL_FMT,read_resume_unit,stdout_unit
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp
        use read_write_module, only: write_resume_file,write_posterior_file
        use feedback_module

        implicit none
        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings



        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data
        double precision, allocatable, dimension(:,:)        :: live_data_local

        double precision, dimension(M%nTotal,settings%nlive) :: incubating_data
        integer :: incubator_index(1)

        integer :: late_successor_index(1)

        integer :: nprocs
        integer :: myrank
        integer :: nlive_local

        integer :: i_live
        integer :: i_slaves

        integer, parameter :: RUNTAG=0
        integer, parameter :: ENDTAG=1

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        logical :: waiting

        double precision, dimension(M%nDims,2)               :: min_max_array

        double precision, allocatable, dimension(:,:) :: posterior_array
        double precision, dimension(M%nDims+M%nDerived+2) :: posterior_point
        integer :: nposterior
        integer :: insertion_index(1)
        integer :: late_index(1)

        logical :: more_samples_needed

        ! The new-born baby point
        double precision,    dimension(M%nTotal)   :: baby_point
        double precision                           :: baby_likelihood

        ! The recently dead point
        double precision,    dimension(M%nTotal)   :: late_point
        double precision                           :: late_likelihood
        double precision :: late_logweight

        ! Point to seed a new one from
        double precision,    dimension(M%nTotal)   :: seed_point


        ! Evidence info
        double precision, dimension(6)             :: evidence_vec


        logical :: resume=.false.
        ! Means to be calculated
        double precision                           :: mean_likelihood_calls
        integer                                    :: total_likelihood_calls

        integer :: ndead

        double precision :: lognlive 
        double precision :: lognlivep1 
        double precision :: logminimumweight



        nprocs = mpi_size()  ! Get the number of MPI procedures
        myrank = mpi_rank()  ! Get the MPI label of the current processor

        ! Initialise any likelihood details by calling it
        evidence_vec(1) = M%loglikelihood(baby_point(M%p0:M%p1),-1)



        if(myrank==0) then 
            call write_opening_statement(M,settings)

            ! Check to see whether there's a resume file present, and record in the
            ! variable 'resume'
            inquire(file=trim(settings%file_root)//'.resume',exist=resume)

            ! Check if we actually want to resume
            resume = settings%read_resume .and. resume

            if(resume .and. settings%feedback>=0) write(stdout_unit,'("Resuming from previous run")')
        end if


        !======= 1) Initialisation =====================================
        ! (i)   On all nodes generate initial live points in parallel by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables for the master node
        ! (iii) Send out the first nprocs-1 tasks to the slaves

        !~~~ (i) Generate Live Points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(resume) then
            if(myrank==0) then
                ! If there is a resume file present, then load the live points from that
                open(read_resume_unit,file=trim(settings%file_root)//'.resume',action='read')
                ! Read the live data
                read(read_resume_unit,'(<M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
            end if
        else !(not resume)
            if(myrank==0) call write_started_generating(settings%feedback)

            ! Otherwise generate them anew:
            ! Create initial live points on all processors, and then merge them onto
            ! the root with MPI_GATHER


            ! First allocate a local live_data array which is nlive/nprocs on each
            ! of the nprocs nodes
            nlive_local = ceiling(settings%nlive/(nprocs+0d0))
            allocate(live_data_local(M%nTotal,nlive_local))

            ! Generate nlive/nprocs live points on each of the nprocs nodes
            live_data_local = GenerateLivePoints(M,nlive_local)

            ! Gather all of this data onto the root node
            call MPI_GATHER(          &  
                live_data_local,      & ! sending array
                M%nTotal*nlive_local, & ! number of elements to be sent
                MPI_DOUBLE_PRECISION, & ! type of element to be sent
                live_data,            & ! recieving array
                M%nTotal*nlive_local, & ! number of elements to be recieved from each node
                MPI_DOUBLE_PRECISION, & ! type of element recieved
                0,                    & ! root node address
                MPI_COMM_WORLD,       & ! communication info
                mpierror)               ! error (from module mpi_module)

            ! deallocate the now unused local live points array to save memory
            deallocate(live_data_local)

            if(myrank==0) call write_finished_generating(settings%feedback) !Flag to note that we're done generating 
        end if !(resume)



        !~~~ (ii) Initialise all variables on master node ~~~~~~~~~~~~~~
        ! There are several variables used throughout the rest of the
        ! algorithm that need to be initialised here
        !  (a) evidence_vec           | Vector containing the evidence, its error, and any other 
        !                             |  things that need to be accumulated over the run.
        !                             |  we need to initialise its sixth argument.
        !  (b) mean_likelihood_calls  | Mean number of likelihood calls over the past nlive iterations
        !  (c) ndead                  | Number of iterations/number of dead points
        !  (d) min_max_array          | Array of maximums and minimums for each coordinate - allows rescaling
        !  (e) posterior_array        | Array of weighted posterior points
        !  (f) incubating_data        | Points that have been generated from a higher loglikelihood contour, and are
        !                             |  waiting to be 'born' i.e. become live points

        if(myrank==0) then 

            ! (a) 
            if(resume) then
                ! If resuming, get the accumulated stats to calculate the
                ! evidence from the resume file
                read(read_resume_unit,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
            else !(not resume) 
                ! Otherwise compute the average loglikelihood and initialise the evidence vector accordingly
                evidence_vec = logzero
                evidence_vec(6) = logsumexp(live_data(M%l0,:)) - log(settings%nlive+0d0)
            end if !(resume) 

            ! (b) initialise the mean number of likelihood calls to 1
            mean_likelihood_calls = sum(live_data(M%nlike,:))/settings%nlive
            total_likelihood_calls = settings%nlive

            ! (c) get number of dead points
            if(resume) then
                ! If resuming, then get the number of dead points from the resume file
                read(read_resume_unit,'(I)') ndead
            else !(not resume) 
                ! Otherwise no dead points originally
                ndead = 0
            end if !(resume) 

            ! (d) calculate the minimums and maximums of the live data
            min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
            min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

            ! (e) Posterior array

            allocate(posterior_array(M%nDims+2,settings%nmax_posterior))
            nposterior = 0
            ! set all of the loglikelihoods and logweights to be zero initially
            posterior_array(1:2,:) = logzero

            ! set the posterior coordinates to be zero initially
            posterior_array(3:,:) = 0d0

            if(resume) then
                ! Read the actual number we've used so far
                read(read_resume_unit,'(I)') nposterior
                !...followed by the posterior array itself
                read(read_resume_unit,'(<M%nDims+M%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)
            end if !(resume) 

            ! Close the resume file if we've openend it
            if(resume) close(read_resume_unit)

            ! Calculate these global variables so we don't need to again
            lognlive   = log(settings%nlive+0d0)
            lognlivep1 = log(settings%nlive+1d0)
            logminimumweight = log(settings%minimum_weight)

            ! (f) Initialise the incubating stack
            incubating_data = 0d0
            incubating_data(M%incubator,:) = flag_incubator_blank


        end if



        !~~~ (iii) Send out first tasks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! We hand over the first nprocs-1 jobs from the master to the slaves
        ! Note that this is slightly less involved than the manner in which we
        ! hand out tasks later, since these tasks are handed out in sequential
        ! order
        if(myrank==0) then
            do i_slaves=1,nprocs-1

                ! Generate a seed point from live_data and incubating_data, and
                ! update the data arrays accordingly
                seed_point = GenerateSeed(M,settings%nlive,live_data,incubating_data)

                ! Send a seed point to the i_slaves th slave
                call MPI_SEND(            &
                    seed_point,           & ! seed point to be sent
                    M%nTotal,             & ! size of this data
                    MPI_DOUBLE_PRECISION, & ! type of this data
                    i_slaves,             & ! send it to the i_slaves point
                    RUNTAG,               & ! tagging information (not important here)
                    MPI_COMM_WORLD,       & ! communication data
                    mpierror              & ! error information (from mpi_module)
                    )
                ! Send the arrays of minimums and maximums
                call MPI_SEND(              &
                    min_max_array,          & ! seed point to be sent
                    M%nDims*2,              & ! size of this data
                    MPI_DOUBLE_PRECISION,   & ! type of this data
                    i_slaves,               & ! send it to the point we just recieved from
                    RUNTAG,                 & ! tagging information (not important here)
                    MPI_COMM_WORLD,         & ! communication data
                    mpierror                & ! error information (from mpi_module)
                    )
            end do


        end if





        !======= 2) Main loop body =====================================
        !
        ! This parallelised by splitting it into two parts: Master and Slaves
        !
        ! The slaves take the job of generating new points from seed points
        ! (within a given likelihood contour)
        !
        ! The master's job is to collate the newly generated points into the
        ! live points array and to calculate evidence

        if(myrank ==0) call write_started_sampling(settings%feedback)

        ! definitely more samples needed than this
        more_samples_needed = .true.

        do while ( more_samples_needed )

            if(myrank == 0) then
                
                !================================================================
                !===================== MASTER NODE ==============================
                !================================================================
                !
                ! (1) Keep track of the lowest loglikelihood contour that hasn't
                !      yet been sent off for sampling (loglikelihood_bound)
                !
                ! (2) Recieve baby point from any slave, along with the contour
                !      it was generated from and insert it into the incubating
                !      stack in order of contours
                !
                ! (3) Send new seed to the now waiting slave. This seed is drawn 
                !      from the set of live points with a likelihood greater
                !      than loglikelihood_bound
                !
                ! (4) Update the live points by birthing any points that are now
                !      ready from the incubating stack



                do i_slaves=1,nprocs-1

                    call MPI_IPROBE(    &  
                        i_slaves,       & !
                        MPI_ANY_TAG,    & !
                        MPI_COMM_WORLD, & !
                        waiting,        & !
                        mpi_status,     & !
                        mpierror        & !
                        )

                    if (waiting) then

                        ! (2) Recieve newly generated baby point from any slave
                        !
                        call MPI_RECV(            &
                            baby_point,           & ! newly generated point to be receieved
                            M%nTotal,             & ! size of this data
                            MPI_DOUBLE_PRECISION, & ! type of this data
                            i_slaves,             & ! recieve it from any slave
                            MPI_ANY_TAG,          & ! tagging information (not important here)
                            MPI_COMM_WORLD,       & ! communication data
                            mpi_status,           & ! status - important (tells you where it's been recieved from )
                            mpierror              & ! error information (from mpi_module)
                            )

                        ! (2) Insert into incubator
                        !
                        ! get the place in the incubator stack for this point
                        incubator_index = nint(baby_point(M%incubator))
                        ! note that this point hasn't launched any new ones
                        baby_point(M%incubator)=flag_live_waiting
                        ! Insert this into the incubator
                        incubating_data(:,incubator_index(1)) = baby_point


                        ! Generate a seed point from live_data and incubating_data, and
                        ! update the data arrays accordingly
                        seed_point = GenerateSeed(M,settings%nlive,live_data,incubating_data) 

                        ! Send a seed point back to that slave
                        call MPI_SEND(              &
                            seed_point,             & ! seed point to be sent
                            M%nTotal,               & ! size of this data
                            MPI_DOUBLE_PRECISION,   & ! type of this data
                            i_slaves,               & ! send it to the point we just recieved from
                            RUNTAG,                 & ! tagging information (not important here)
                            MPI_COMM_WORLD,         & ! communication data
                            mpierror                & ! error information (from mpi_module)
                            )
                        ! Send the arrays of minimums and maximums
                        call MPI_SEND(              &
                            min_max_array,          & ! seed point to be sent
                            M%nDims*2,              & ! size of this data
                            MPI_DOUBLE_PRECISION,   & ! type of this data
                            i_slaves,               & ! send it to the point we just recieved from
                            RUNTAG,                 & ! tagging information (not important here)
                            MPI_COMM_WORLD,         & ! communication data
                            mpierror                & ! error information (from mpi_module)
                            )
                    end if
                end do


                ! Transfer from incubating stack to live_data if necessary
                do while(.true.)
                    ! Find the point with the lowest likelihood...
                    late_index = minloc(live_data(M%l0,:))
                    ! ...and save it.
                    late_point = live_data(:,late_index(1))
                    ! Get the likelihood contour
                    late_likelihood = late_point(M%l0)
                    ! Calculate the late logweight
                    late_logweight = (ndead-1)*lognlive - ndead*lognlivep1 

                    ! Find the position of the successor to the late point within
                    ! the incubation stack
                    late_successor_index = nint( late_point(M%incubator) )

                    ! Check to see if the late point has a generated index
                    if(late_successor_index(1)==flag_live_waiting) exit
                    if(incubating_data( M%incubator, late_successor_index(1) )<flag_live_waiting ) exit

                    ! Birth the new point from the incubator
                    baby_point = incubating_data(:,late_successor_index(1)) 
                    baby_likelihood  = baby_point(M%l0)
                    ! Delete the point from the incubator
                    incubating_data(M%incubator,late_successor_index(1))=flag_incubator_blank

                    ! (3) Insert the baby point into the set of live points (over the
                    !     old position of the dead points

                    ! Insert the baby point over the late point
                    live_data(:,late_index(1)) = baby_point

                    ! record that we have a new dead point
                    ndead = ndead + 1

                    ! If we've put a limit on the maximum number of iterations, then
                    ! check to see if we've reached this
                    if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

                    ! update the minimum and maximum values of the live points
                    min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
                    min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)



                    ! (4) Calculate the new evidence (and check to see if we're accurate enough)
                    call settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed,evidence_vec)




                    ! (5) Update the set of weighted posteriors
                    if( settings%calculate_posterior .and. late_point(M%l0) + late_logweight - evidence_vec(1) > logminimumweight ) then
                        ! If the late point has a sufficiently large weighting, then we
                        ! should add it to the set of saved posterior points

                        ! calculate a new point for insertion
                        posterior_point(1)  = late_point(M%l0) + late_logweight
                        posterior_point(2)  = late_point(M%l0)
                        posterior_point(3:) = late_point(M%p0:M%d1)

                        if(nposterior<settings%nmax_posterior) then
                            ! If we're still able to use a restricted array,

                            ! Find the closest point in the array which is beneath the minimum weight
                            insertion_index = minloc(posterior_array(1,:nposterior),mask=posterior_array(1,:nposterior)<logminimumweight+evidence_vec(1))

                            if(insertion_index(1)==0) then
                                ! If there are no points to overwrite, then we should
                                ! expand the available storage array
                                nposterior=nposterior+1
                                posterior_array(:,nposterior) = posterior_point
                            else
                                ! Otherwise overwrite the 
                                posterior_array(:,insertion_index(1)) = posterior_point
                            end if

                        else
                            ! Otherwise we have to overwrite the smallest element
                            insertion_index = minloc(posterior_array(1,:nposterior))
                            posterior_array(:,insertion_index(1)) = posterior_point
                        end if

                    end if


                    ! (6) Command line feedback

                    ! update the mean number of likelihood calls
                    mean_likelihood_calls = mean_likelihood_calls + (baby_point(M%nlike) - late_point(M%nlike) ) / (settings%nlive + 0d0)

                    ! update the total number of likelihood calls
                    total_likelihood_calls = total_likelihood_calls + baby_point(M%nlike)


                    ! Feedback to command line every nlive iterations
                    if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                        write(stdout_unit,'("ndead     = ", I20                  )') ndead
                        write(stdout_unit,'("efficiency= ", F20.2                )') mean_likelihood_calls
                        write(stdout_unit,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                        write(stdout_unit,'("")')
                    end if

                    ! (7) Update the resume and posterior files every update_resume iterations, or at program termination
                    if (mod(ndead,settings%update_resume) .eq. 0 .or.  more_samples_needed==.false.)  then
                        if(settings%write_resume) call write_resume_file(settings,M,live_data,evidence_vec,ndead,nposterior,posterior_array) 
                        if(settings%calculate_posterior) call write_posterior_file(settings,M,posterior_array,evidence_vec(1),nposterior)  
                    end if

                end do

            else
                !================================================================
                !===================== SLAVE NODES ==============================
                !================================================================

                ! Listen for a signal from the master
                call MPI_RECV(            &
                    seed_point,           & ! seed point to be recieved
                    M%nTotal,             & ! size of this data
                    MPI_DOUBLE_PRECISION, & ! type of this data
                    0,                    & ! recieve it from the master
                    MPI_ANY_TAG,          & ! recieve any tagging information
                    MPI_COMM_WORLD,       & ! communication data
                    mpi_status,           & ! status (not important here)
                    mpierror              & ! error information (from mpi_module)
                    )

                ! If we receive a kill signal, then exit the loop
                if(mpi_status(MPI_TAG)==ENDTAG) then
                    more_samples_needed=.false.
                    exit
                end if

                call MPI_RECV(            &
                    min_max_array,        & ! seed point to be recieved
                    M%nDims*2,            & ! size of this data
                    MPI_DOUBLE_PRECISION, & ! type of this data
                    0,                    & ! recieve it from the master
                    MPI_ANY_TAG,          & ! recieve any tagging information
                    MPI_COMM_WORLD,       & ! communication data
                    mpi_status,           & ! status (not important here)
                    mpierror              & ! error information (from mpi_module)
                    )

                ! Calculate a new baby point from the seed point
                baby_point = settings%sampler(seed_point, min_max_array, M)

                ! Send the baby point back
                call MPI_SEND(            &
                    baby_point,           & ! baby point to be sent
                    M%nTotal,             & ! size of this data
                    MPI_DOUBLE_PRECISION, & ! type of this data
                    0,                    & ! send it to the master
                    RUNTAG,               & ! tagging information (not important here)
                    MPI_COMM_WORLD,       & ! communication data
                    mpierror              & ! error information (from mpi_module)
                    )

            end if
            
        end do ! End main loop



        if (myrank==0) then
            
            ! If we're done, then clean up by receiving the last piece of
            ! data from each node (and throw it away) and then send a kill signal back to it
            if(more_samples_needed==.false.) then
                do i_live=1,nprocs-1
                    call MPI_RECV(            &
                        baby_point,           & ! newly generated point to be receieved
                        M%nTotal,             & ! size of this data
                        MPI_DOUBLE_PRECISION, & ! type of this data
                        MPI_ANY_SOURCE,       & ! recieve it from any slave
                        MPI_ANY_TAG,          & ! tagging information (not important here)
                        MPI_COMM_WORLD,       & ! communication data
                        mpi_status,           & ! status - important (tells you where it's been recieved from )
                        mpierror              & ! error information (from mpi_module)
                        )
                    call MPI_SEND(              &
                        min_max_array,          & ! seed point to be sent
                        M%nDims*2,              & ! size of this data
                        MPI_DOUBLE_PRECISION,   & ! type of this data
                        mpi_status(MPI_SOURCE), & ! send it to the point we just recieved from
                        ENDTAG,                 & ! tagging information (not important here)
                        MPI_COMM_WORLD,         & ! communication data
                        mpierror                & ! error information (from mpi_module)
                        )
                end do


            end if




            call write_final_results(M,evidence_vec,ndead,total_likelihood_calls,settings%feedback)  
        end if

    end subroutine NestedSamplingP




    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePoints(M,nlive) result(live_data)
        use model_module,    only: model, calculate_point
        use random_module,   only: random_reals
        use utils_module,    only: logzero
        implicit none

        !> The model details (loglikelihood, priors, ndims etc...)
        type(model), intent(in) :: M

        !> The number of points to be generated
        integer, intent(in) :: nlive

        !live_data(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,nlive) :: live_data

        ! Loop variable
        integer i_live

        ! initialise live points at zero
        live_data = 0d0

        do i_live=1,nlive

            ! Generate a random coordinate
            live_data(:,i_live) = random_reals(M%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, live_data(:,i_live) )

        end do

        ! Set the number of likelihood calls for each point to 1
        live_data(M%nlike,:) = 1

        ! Set the initial trial values of the chords as the diagonal of the hypercube
        live_data(M%last_chord,:) = sqrt(M%nDims+0d0)

        ! Initially, none of the points have been calculated yet
        live_data(M%incubator,:) = flag_live_waiting

        ! Set the likelihood contours to logzero for now
        live_data(M%l1,:) = logzero


    end function GenerateLivePoints

    function GenerateSeed(M,nlive,live_data,incubating_data) result(seed_point)
        use model_module,      only: model
        use random_module,     only: random_integer
        implicit none
        type(model),      intent(in) :: M
        integer,          intent(in) :: nlive
        double precision, intent(inout), dimension(M%nTotal,nlive) :: live_data
        double precision, intent(inout), dimension(M%nTotal,nlive) :: incubating_data

        ! Point to seed a new one from
        double precision,    dimension(M%nTotal)   :: seed_point


        integer :: incubator_index(1)
        integer :: live_index(1)

        double precision :: loglikelihood_bound


        ! Find the lowest likelihood point whose contour is waiting to
        ! be generated from
        live_index = minloc(live_data(M%l0,:),mask=nint(live_data(M%incubator,:))==flag_live_waiting) 

        ! Find a place in the incubator stack for the generated point
        ! We search through the incubator stack to find the first index 
        incubator_index = minloc(incubating_data(M%incubator,:),mask=nint(incubating_data(M%incubator,:))==flag_incubator_blank) 

        ! Give this place to the point that generated the contour
        live_data(M%incubator,live_index(1)) = incubator_index(1)

        ! Note at this place in the incubator that we're waiting on a
        ! point to be generated
        incubating_data(M%incubator,incubator_index(1))=flag_incubator_running
             
        ! Select a seed point for the generator
        !  -excluding the points which have likelihoods equal to the
        !   loglikelihood bound
        loglikelihood_bound = live_data(M%l0,live_index(1))
        seed_point(M%l0)=loglikelihood_bound

        do while (seed_point(M%l0)<=loglikelihood_bound )
            ! get a random integer in [1,nlive]
            ! get this point from live_data 
            seed_point = live_data(:,random_integer(nlive))
        end do

        ! Record the likelihood bound which this seed will generate from
        seed_point(M%l1) = loglikelihood_bound

        ! Record the eventual position in the incubator stack
        seed_point(M%incubator) = incubator_index(1)

    end function GenerateSeed




end module nested_sampling_parallel_module
