module nested_sampling_parallel_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSamplingP(M,settings)
        use mpi_module
        use model_module,      only: model
        use utils_module,      only: logzero,loginf,DBL_FMT,read_resume_unit,dbleq
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp
        use random_module,     only: random_integer
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

        double precision, dimension(M%nDims,2)               :: min_max_array

        double precision, allocatable, dimension(:,:) :: posterior_array
        double precision, dimension(M%nDims+2) :: posterior_point
        integer :: nposterior
        integer :: nremove

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

        integer :: nprocs
        integer :: myrank
        integer :: nlive_local

        integer :: i_live
        integer :: nslaves

        integer, parameter :: RUNTAG=0
        integer, parameter :: ENDTAG=1

        integer, dimension(MPI_STATUS_SIZE) :: mpi_status

        double precision :: loglikelihood_bound





        nprocs = mpi_size()  ! Get the number of MPI procedures
        myrank = mpi_rank()  ! Get the MPI label of the current processor

        ! Initialise any likelihood details by calling it
        loglikelihood_bound = M%loglikelihood(baby_point(M%p0:M%p1),-1)



        if(myrank==0) then 
            call write_opening_statement(M,settings)
            inquire(file=trim(settings%file_root)//'.resume',exist=resume)
            if(resume) write(*,*) "Resuming from previous run"
        end if


        !======= 1) Initialisation =====================================
        ! (i)   On all nodes generate initial live points in parallel by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables for the master node
        ! (iii) Send out the first nprocs-1 messages to the slaves



        !~~~ (i) Generate Live Points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(resume) then
            if(myrank==0) then
                ! If there is a resume file present, then load the live points from that
                write(*,*) "Reading live data"
                open(read_resume_unit,file=trim(settings%file_root)//'.resume',action='read')
                read(read_resume_unit,'(<M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
            end if
        else
            if(myrank==0) call write_started_generating(settings%feedback)
            ! Create initial live points on all processors, and then merge them onto
            ! the root with MPI_GATHER

            ! First allocate a local live_data array which is nlive/nprocs on each
            ! of the nprocs nodes
            nlive_local = ceiling(settings%nlive/(nprocs+0d0))
            allocate(live_data_local(M%nTotal,nlive_local))

            ! Generate nlive/nprocs live points on each of the nprocs nodes
            live_data_local = GenerateLivePoints(M,nlive_local)

            ! Sort the live points in order of likelihood, first point lowest, last point highest on every node
            call quick_sort(live_data_local,M%l0)

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

            ! Sort the live points in order of likelihood, first point lowest, last point highest on the root node
            ! (note that this is made easier by the fact that this array is semi-sorted
            if(myrank==0) call quick_sort(live_data,M%l0)

            if(myrank==0) call write_finished_generating(settings%feedback) 
        end if





        !~~~ (ii) Initialise master node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! There are several variables used throughout the rest of the
        ! algorithm that need to be initialised here
        !  (a) evidence_vec           | Vector containing the evidence, its error, and any other 
        !                             |  things that need to be accumulated over the run.
        !                             |  we need to initialise its sixth argument.
        !  (b) mean_likelihood_calls  | Mean number of likelihood calls over the past nlive iterations
        !  (c) ndead                  | Number of iterations/number of dead points
        !  (d) incubating_data        | Points that have been generated from a higher loglikelihood contour, and are
        !                             |  waiting to be 'born' i.e. become live points
        !  (f) min_max_array          | Array of maximums and minimums for each coordinate - allows rescaling
        !  (g) posterior_array        | Array of weighted posterior points

        if(myrank==0) then 

            ! (a) 
            if(resume) then
                ! If resuming, get the accumulated stats to calculate the
                ! evidence from the resume file
                read(read_resume_unit,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
            else
                ! Otherwise compute the average loglikelihood and initialise the evidence vector accordingly
                evidence_vec = logzero
                evidence_vec(6) = logsumexp(live_data(M%l0,:)) - log(settings%nlive+0d0)
            end if

            ! (b) initialise the mean number of likelihood calls to 1
            mean_likelihood_calls = sum(live_data(M%d0,:))/settings%nlive
            total_likelihood_calls = settings%nlive

            ! (c) record the first dead point 
            late_point = live_data(:,1)

            if(resume) then
                ! If resuming, then get the number of dead points from the resume file
                write(*,*) "Reading ndead"
                read(read_resume_unit,'(I)') ndead
            else
                ! Otherwise initialise the number of dead points at 1
                ndead = 1
            end if

            ! Get the likelihood contour
            late_likelihood = late_point(M%l0)

            ! (d) Initialise the incubating stack
            incubating_data = 0d0
            incubating_data(M%l1,:) = loginf

            ! (f) calculate the minimums and maximums of the live data
            min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
            min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

            ! Calculate these global variables so we don't need to again
            lognlive   = log(settings%nlive+0d0)
            lognlivep1 = log(settings%nlive+1d0)

            ! (g) Posterior array
            allocate(posterior_array(M%nDims+2,settings%nmax_posterior))
            nposterior=0
            posterior_array(1:2,:) = logzero
            posterior_array(3:,:) = 0d0
            if(resume) then
                read(read_resume_unit,'(I)') nposterior
                read(read_resume_unit,'(<M%nDims+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)
                close(read_resume_unit)
            end if

        end if




        !~~~ (iii) Send out first tasks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! We hand over the first nprocs-1 jobs from the master to the slaves
        ! Note that this is slightly less involved than the manner in which we
        ! hand out tasks later, since these tasks are handed out in sequential
        ! order
        if(myrank==0) then
            do nslaves=1,nprocs-1

                ! Note the loglikelihood bound that this point will generate from
                loglikelihood_bound = live_data(M%l0,nslaves)

                ! Select a seed point for the generator
                !  -excluding the points which have likelihoods equal to the
                !   loglikelihood bound

                seed_point(M%l0)=loglikelihood_bound
                do while (seed_point(M%l0)<=loglikelihood_bound )
                    ! get a random integer in [1,nlive]
                    ! get this point from live_data 
                    seed_point = live_data(:,random_integer(settings%nlive))
                end do

                ! Record the likelihood bound which this seed will generate from
                seed_point(M%l1) = loglikelihood_bound

                ! Record that this likelihood is now active
                live_data(M%d0+2,nslaves) = 1d0


                ! Send a seed point to the nslaves th slave
                call MPI_SEND(            &
                    seed_point,           & ! seed point to be sent
                    M%nTotal,             & ! size of this data
                    MPI_DOUBLE_PRECISION, & ! type of this data
                    nslaves,              & ! send it to the nslaves point
                    RUNTAG,               & ! tagging information (not important here)
                    MPI_COMM_WORLD,       & ! communication data
                    mpierror              & ! error information (from mpi_module)
                    )
                ! Send the arrays of minimums and maximums
                call MPI_SEND(              &
                    min_max_array,          & ! seed point to be sent
                    M%nDims*2,              & ! size of this data
                    MPI_DOUBLE_PRECISION,   & ! type of this data
                    nslaves,                & ! send it to the point we just recieved from
                    RUNTAG,                 & ! tagging information (not important here)
                    MPI_COMM_WORLD,         & ! communication data
                    mpierror                & ! error information (from mpi_module)
                    )
            end do

        end if





        !======= 2) Main Loop ==========================================
        !
        ! This parallelised by splitting it into two parts: Master and Slaves
        !
        ! The slaves take the job of generating new points from seed points
        !
        ! The master's job is to collate the newly generated points into the
        ! live points array and to calculate evidence



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



                ! (1) Find the lowest likelihood which hasn't been sent off for
                ! sampling
                do i_live=1,settings%nlive

                    if( live_data(M%d0+2,i_live) < 0 ) then
                        loglikelihood_bound= live_data(M%l0,i_live)

                        ! Record that this likelihood is now running
                        live_data(M%d0+2,i_live) = 1d0
                        exit
                    end if

                end do

                if (i_live==settings%nlive) write(*,*) 'oh no!'

                ! (2) Recieve newly generated baby point from any slave
                !
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


                ! add the new baby point to the incubation stack
                call insert_point_into_incubating_data(baby_point,incubating_data,M%l1)

                ! (3) Select a seed point for the generator
                !  -excluding the points which have likelihoods equal to the
                !   loglikelihood bound
                seed_point(M%l0)=loglikelihood_bound
                do while ( seed_point(M%l0)<=loglikelihood_bound )
                    ! get a random number in [1,nlive]
                    ! get this point from live_data 
                    seed_point = live_data(:,random_integer(settings%nlive))
                end do

                ! Record the likelihood bound which this seed will generate from
                seed_point(M%l1) = loglikelihood_bound

                ! Send a seed point back to that slave
                call MPI_SEND(              &
                    seed_point,             & ! seed point to be sent
                    M%nTotal,               & ! size of this data
                    MPI_DOUBLE_PRECISION,   & ! type of this data
                    mpi_status(MPI_SOURCE), & ! send it to the point we just recieved from
                    RUNTAG,                 & ! tagging information (not important here)
                    MPI_COMM_WORLD,         & ! communication data
                    mpierror                & ! error information (from mpi_module)
                    )
                ! Send the arrays of minimums and maximums
                call MPI_SEND(              &
                    min_max_array,          & ! seed point to be sent
                    M%nDims*2,              & ! size of this data
                    MPI_DOUBLE_PRECISION,   & ! type of this data
                    mpi_status(MPI_SOURCE), & ! send it to the point we just recieved from
                    RUNTAG,                 & ! tagging information (not important here)
                    MPI_COMM_WORLD,         & ! communication data
                    mpierror                & ! error information (from mpi_module)
                    )


                ! (4) transfer from incubating stack to live_data if necessary
                !
                do while( dbleq( incubating_data(M%l1,1) , late_likelihood ) )

                    ! birth the new point from the incubator
                    baby_point = incubating_data(:,1)             
                    baby_likelihood  = baby_point(M%l0)
                    ! pop the stack
                    incubating_data(:,:settings%nlive-1) = incubating_data(:,2:)

                    ! Record that this likelihood hasn't been set running yet
                    baby_point(M%d0+2) = -1d0

                    ! Insert the new point
                    call insert_point_into_live_data(baby_point,live_data,M%l0)

                    ! Calculate the new evidence (and check to see if we're accurate enough)
                    call settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed,evidence_vec)

                    ! update the mean number of likelihood calls
                    mean_likelihood_calls = mean_likelihood_calls + (baby_point(M%d0) - late_point(M%d0) ) / (settings%nlive + 0d0)

                    ! update the total number of likelihood calls
                    total_likelihood_calls = total_likelihood_calls + baby_point(M%d0)

                    ! record the next point thats about to die
                    late_point = live_data(:,1)
                    ndead = ndead + 1

                    ! Get the new likelihood contour
                    late_likelihood = late_point(M%l0)

                    ! Calculate the late logweight
                    late_logweight = (ndead-1)*lognlive - ndead*lognlivep1 

                    ! Update the set of weighted posteriors
                    if(settings%calculate_posterior .and. &
                        late_point(M%l0) + late_logweight - evidence_vec(1) > log(settings%minimum_weight) ) then

                        ! First trim off any points that are now under the minimum_weight limit
                        nremove = count( posterior_array(1,1:nposterior)-evidence_vec(1)<=log(settings%minimum_weight) )
                        posterior_array(1:2,nposterior-nremove+1:nposterior) = logzero
                        posterior_array(3:, nposterior-nremove+1:nposterior) = 0d0
                        nposterior = nposterior-nremove


                        ! Now add the new point
                        nposterior=min(nposterior+1,settings%nmax_posterior)
                        posterior_point(1) = late_point(M%l0) + late_logweight
                        posterior_point(2) = late_point(M%l0)
                        posterior_point(3:) = late_point(M%p0:M%p1)
                        call insert_into_posterior(posterior_point,posterior_array(:,1:nposterior))
                    end if


                    ! Feedback to command line every nlive iterations
                    if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                        write(*,'("ndead     = ", I20                  )') ndead
                        write(*,'("efficiency= ", F20.2                )') mean_likelihood_calls
                        write(*,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                        write(*,*)
                    end if

                end do

                ! update the minimum and maximum values of the live points
                min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
                min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)





                ! Halt if we've reached the desired maximum iterations
                if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.


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

                    call write_final_results(M,evidence_vec,ndead,total_likelihood_calls,settings%feedback)
                end if


                ! Update the resume and posterior files every update_resume iterations, or at program termination
                if (mod(ndead,settings%update_resume) .eq. 0 .or.  more_samples_needed==.false.)  then
                    call write_resume_file(settings,M,live_data,evidence_vec,ndead,posterior_array(:,:nposterior),nposterior) 
                    call write_posterior_file(settings,M,posterior_array,evidence_vec(1),nposterior) 
                end if



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



        end do

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
        live_data(M%d0,:) = 1

        ! Set the initial trial values of the chords as the diagonal of the hypercube
        live_data(M%d0+1,:) = sqrt(M%nDims+0d0)

        ! Initially, none of the points have been calculated yet
        live_data(M%d0+2,:) = -1d0

        ! Set the likelihood contours to logzero for now
        live_data(M%l1,:) = logzero


    end function GenerateLivePoints



    !> Sort nlive vectors in data_array(:,nlive) by the second to last element of each
    !! vector
    recursive subroutine quick_sort(data_array,colnum)

        !> an array of nlive  (coord,derived,like) to be sorted by the
        !! likelihood value at the end of each vector
        double precision, dimension(:,:), intent(inout) :: data_array

        !> The column with which to sort by
        integer, intent(in) :: colnum

        integer :: iq

        if( size(data_array, 2 ) > 1) then
            call partition(data_array,iq)
            call quick_sort( data_array(:,:iq-1), colnum )
            call quick_sort( data_array(:,iq:),   colnum )
        endif
        contains

        !> @todo comment this quick_sort partition
        subroutine partition(data_array, marker)

            double precision, intent(inout), dimension(:,:) :: data_array
            integer, intent(out) :: marker
            integer :: i, j
            double precision, dimension(size(data_array,1)) :: temp,x

            x(:) = data_array(:,1)
            i= 0
            j= size( data_array, 2 ) + 1

            do
                j = j-1
                do
                    if (data_array(colnum,j) <= x(colnum)) exit
                    j = j-1
                end do

                i = i+1
                do
                    if (data_array(colnum,i) >= x(colnum)) exit
                    i = i+1
                end do

                if (i < j) then
                    ! exchange data_array(:,i) and data_array(:,j)
                    temp(:)         = data_array(:,i)
                    data_array(:,i) = data_array(:,j)
                    data_array(:,j) = temp(:)

                    elseif (i == j) then
                    marker = i+1
                    return
                else
                    marker = i
                    return
                endif
            end do

        end subroutine partition

    end subroutine quick_sort





    !> Insert the new point into the live data array
    !!
    !! Since the data array is already sorted, one can insert a new point using
    !! binary search insertion algorithm
    subroutine insert_point_into_live_data(baby_point,live_data,loglike_pos)
        !> The point to be inserted by order of its last value
        double precision, intent(in),    dimension(:)   :: baby_point
        !> The live data array to be inserted into
        double precision, intent(inout), dimension(:,:) :: live_data
        !> The index to sort by (in this case it's M%l0)
        integer,intent(in)          :: loglike_pos   

        double precision :: baby_loglike  !loglikelihood of the baby point
        integer          :: nlive         !number of live points
        integer          :: baby_position !where to insert the baby point in the array


        nlive       = size(live_data,2)    ! size of the array

        baby_loglike = baby_point(loglike_pos) ! loglikelihood of the baby point

        ! search for the position with a binary search algorithm
        baby_position =  binary_search(1,nlive,baby_loglike,loglike_pos,live_data)

        ! Delete the lowest likelihood point by shifting the points in the data
        ! array from baby_position-1 and below down by one and add the new point
        ! at baby_position-1
        live_data(:,:baby_position-1) = eoshift(live_data(:,:baby_position-1),dim=2,shift=+1,boundary=baby_point)

    end subroutine insert_point_into_live_data

    !> Insert the new point into the posterior array (pun very much intended)
    !!
    !! Since the data array is already sorted, one can insert a new point using
    !! binary search insertion algorithm
    subroutine insert_into_posterior(point,posterior_data)
        !> The point to be inserted by order of its last value
        double precision, intent(in),    dimension(:)   :: point
        !> The live data array to be inserted into
        double precision, intent(inout), dimension(:,:) :: posterior_data
        !> The index to sort by (in this case it's M%l0)

        integer          :: nposterior     !number of live points
        integer          :: point_position !where to insert the baby point in the array


        nposterior       = size(posterior_data,2)    ! size of the array

        ! search for the position with a binary search algorithm
        ! (note the minus signs so as to interface with the binary search, since
        ! we want the order in terms of highest to lowest)
        point_position =  binary_search(1,nposterior,-point(1),1,-posterior_data)

        ! shift the points up one (since we are sorting in order of highest to lowest)
        posterior_data(:,point_position:) = eoshift(posterior_data(:,point_position:),dim=2,shift=-1,boundary=point)

    end subroutine insert_into_posterior

    !> Insert the new point into the incubating stack
    !!
    !! Since the data array is already sorted, one can insert a new point using
    !! binary search insertion algorithm
    subroutine insert_point_into_incubating_data(incubating_point,incubating_data,loglike_pos)
        !> The point to be inserted by order of its last value
        double precision, intent(in),    dimension(:)   :: incubating_point
        !> The live data array to be inserted into
        double precision, intent(inout), dimension(:,:) :: incubating_data
        !> The index to sort by (in this case its the M%l1
        integer,intent(in)          :: loglike_pos   

        double precision :: incubating_loglike  !loglikelihood of the incubating point
        integer          :: nincubating         !number of live points
        integer          :: incubating_position !where to insert the incubating point in the array


        nincubating = size(incubating_data,2) ! size of the array

        incubating_loglike = incubating_point(loglike_pos) ! loglikelihood of the incubating point

        ! search for the position with a binary search algorithm
        incubating_position = binary_search(1,nincubating,incubating_loglike,loglike_pos,incubating_data)

        ! Shift the stack up 
        incubating_data(:,incubating_position:) = eoshift(incubating_data(:,incubating_position:),dim=2,shift=-1,boundary=incubating_point)

    end subroutine insert_point_into_incubating_data




    !> Binary search algorithm
    !!
    !! Assuming the data are ordered from lowest to highest, this algorithm
    !! takes in loglike, and states at which index this point belongs at.
    !!
    !! e.g. if the array is:
    !!
    !! index | 1 | 2 | 3 | 4 | 5 | 6 |
    !! ------|---|---|---|---|---|---|
    !! value | 3 | 10| 11| 15| 16| 17|
    !!
    !! then  14 would have index 4, 2 would have index 1 and 20 would have index 7.
    recursive function binary_search(imin,imax,value,array_index,data_array) result(point_pos)
        implicit none
        !> Lower bound of search
        integer, intent(in) :: imin
        !> Upper bound of search
        integer, intent(in) :: imax
        !> inquiry value
        double precision, intent(in) :: value

        !> The data array to be consulted
        double precision, intent(in), dimension(:,:) :: data_array

        !> The index of the 2D array which we are sorting by
        integer, intent(in) :: array_index
        integer :: point_pos

        if(value < data_array(array_index,imin)) then
            ! beneath the bounds of the array, has index of the minimum
            point_pos = imin
        else if (value > data_array(array_index,imax)) then
            ! above the bounds of the array, has index one above the maximum
            point_pos = imax+1
        else if (imin+1==imax) then
            ! directly between imin and imax, so set it to be equal to imax
            point_pos = imax
        else
            ! calculate the lower value of the midpoint of the two limits
            ! e.g. imax=6,imin=1 gives 3 
            point_pos = (imin + imax)/2

            if( value < data_array(array_index,point_pos) ) then
                ! point in lower subset
                point_pos = binary_search(imin,point_pos,value,array_index,data_array)
            else 
                ! point in upper subset
                point_pos = binary_search(point_pos,imax,value,array_index,data_array)
            end if
        end if

    end function binary_search















end module nested_sampling_parallel_module
