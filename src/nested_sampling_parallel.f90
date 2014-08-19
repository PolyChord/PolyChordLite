module nested_sampling_parallel_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSamplingP(M,settings)
        use mpi_module
        use model_module,    only: model
        use utils_module,    only: logzero,loginf
        use settings_module, only: program_settings
        use utils_module,    only: logsumexp
        use random_module,   only: random_integers
        use feedback_module

        implicit none
        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings



        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data

        ! The new-born baby point
        double precision,    dimension(M%nTotal)   :: baby_point
        ! The recently dead point
        double precision,    dimension(M%nTotal)   :: late_point

        ! Point to seed a new one from
        double precision,    dimension(M%nTotal)   :: seed_point

        ! temp variable for getting a random integer
        integer, dimension(1) :: point_number 

        double precision :: baby_likelihood
        double precision :: late_likelihood

        double precision, dimension(6) :: evidence_vec

        ! Means to be calculated
        double precision :: mean_likelihood_calls

        logical :: more_samples_needed

        integer :: ndead

        double precision,    dimension(M%nDims,2)   :: min_max_array

        double precision, dimension(M%nTotal,settings%nlive) :: incubating_data

        double precision, allocatable, dimension(:,:) :: live_data_local
        integer :: nprocs
        integer :: myrank
        integer :: nlive_local

        integer :: i_live
        integer :: nslaves

        integer,parameter :: RUNTAG=0
        integer,parameter :: ENDTAG=1

        integer mpi_status(MPI_STATUS_SIZE)

        double precision,allocatable, dimension(:) :: running_likelihoods

        double precision :: loglikelihood_bound







        nprocs = mpi_size()  ! Get the number of MPI procedures
        myrank = mpi_rank()  ! Get the MPI label of the current processor



        if(myrank==0) then 
            call write_opening_statement(M,settings)
            call write_started_generating(settings%feedback)
        end if


        !======= 1) Initialisation =====================================
        ! (i)   On all nodes generate initial live points in parallel by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables for the master node
        ! (iii) Send out the first nprocs-1 messages to the slaves



        !~~~ (i) Generate Live Points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !   Create initial live points on all processors, and then merge them onto
        !   the root with MPI_GATHER

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
        !  (e) running_likelihoods    | Array of likelihood contours that are running at the moment
        !  (f) min_max_arrray         | Array of maximums and minimums for each coordinate - allows rescaling

        if(myrank==0) then 

            ! (a) Compute the average loglikelihood and initialise the evidence vector accordingly
            evidence_vec = logzero
            evidence_vec(6) = logsumexp(live_data(M%l0,:)) - log(settings%nlive+0d0)

            ! (b) initialise the mean number of likelihood calls to 1
            mean_likelihood_calls = 1d0

            ! (c) record the first dead point 
            late_point = live_data(:,1)
            ndead = 1

            ! Get the likelihood contour
            late_likelihood = late_point(M%l0)

            ! If we're saving the dead points, open the relevant file
            if (settings%save_dead) open(unit=222, file='dead_points.dat')

            ! (d) Initialise the incubating stack
            incubating_data = 0d0
            incubating_data(M%l1,:) = loginf

            ! (e) Initialise the running likelihood stack
            allocate(running_likelihoods(nprocs-1))
            running_likelihoods=loginf

            ! (f) calculate the minimums and maximums of the live data
            min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
            min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

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
                    ! get a random number in [1,nlive]
                    point_number = random_integers(1,settings%nlive) 
                    ! get this point from live_data 
                    seed_point = live_data(:,point_number(1))
                end do

                ! Record the likelihood bound which this seed will generate from
                seed_point(M%l1) = loglikelihood_bound

                ! Record that this likelihood is now active
                running_likelihoods(nslaves) = loglikelihood_bound


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








        ! definitely more samples needed than this
        more_samples_needed = .true.


        do while ( more_samples_needed )

            if(myrank == 0) then
                ! The master nodes tasks:
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


                ! (1) Find the lowest likelihood which is neither running nor incubating
                i_live=1
                loglikelihood_bound = late_likelihood
                do while( loglikelihood_bound<late_likelihood .or. &
                        any(loglikelihood_bound==running_likelihoods(:)) .or. &
                        any(loglikelihood_bound==incubating_data(M%l1,:)) )

                    loglikelihood_bound= live_data(M%l0,i_live)
                    i_live=i_live+1

                end do

                !write(*,'("live likelihoods:    ", <settings%nlive>E13.4)') live_data(M%l0,:)
                !write(*,'("incubating_data:     ", <settings%nlive>E13.4)') incubating_data(M%l1,:)
                !write(*,'("running_likelihoods: ", <nprocs>E13.4)')         running_likelihoods
                !write(*,'("loglikelihood_bound: ", E13.4)')                 loglikelihood_bound
                !write(*,*) '------------------------------------------'

                ! Listen for a signal from any waiting slave
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

                ! Select a seed point for the generator
                !  -excluding the points which have likelihoods equal to the
                !   loglikelihood bound
                seed_point(M%l0)=loglikelihood_bound
                do while ( seed_point(M%l0)<=loglikelihood_bound )
                    ! get a random number in [1,nlive]
                    point_number = random_integers(1,settings%nlive) 
                    ! get this point from live_data 
                    seed_point = live_data(:,point_number(1))
                end do

                ! Record the likelihood bound which this seed will generate from
                seed_point(M%l1) = loglikelihood_bound

                ! Record that this likelihood is now running
                running_likelihoods(mpi_status(MPI_SOURCE)) = loglikelihood_bound



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


                ! transfer from incubating stack to live_data if necessary
                do while( incubating_data(M%l1,1) == late_likelihood )

                    ! birth the new point from the incubator
                    baby_point = incubating_data(:,1)             
                    baby_likelihood  = baby_point(M%l0)
                    ! pop the stack
                    incubating_data(:,:settings%nlive-1) = incubating_data(:,2:)

                    ! Insert the new point
                    call insert_point_into_live_data(baby_point,live_data,M%l0)

                    ! Calculate the new evidence (and check to see if we're accurate enough)
                    call settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed,evidence_vec)

                    ! update the mean number of likelihood calls
                    mean_likelihood_calls = mean_likelihood_calls + (baby_point(M%d0) - late_point(M%d0) ) / (settings%nlive + 0d0)


                    ! record the next point thats about to die
                    late_point = live_data(:,1)
                    ndead = ndead + 1

                    ! Get the new likelihood contour
                    late_likelihood = late_point(M%l0)


                    ! Feedback to command line every nlive iterations
                    if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                        write(*,'("ndead     = ", I20                  )') ndead
                        write(*,'("efficiency= ", F20.2                )') mean_likelihood_calls
                        write(*,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                        write(*,*)
                    end if

                end do

                ! Halt if we've reached the desired maximum iterations
                if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

                ! Write the dead points to a file if desired
                if (settings%save_dead) write(222,'(<M%nTotal+1>E17.9)') exp(ndead*log(settings%nlive/(settings%nlive+1d0))), late_point

                ! update the minimum and maximum values of the live points
                min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
                min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)


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

                    ! close the dead points file if we're done
                    if (settings%save_dead) close(222)

                    call write_final_results(M,evidence_vec,ndead,settings%feedback)
                end if



            else
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
        !> The index to sort by (in this case its M%l0
        integer,intent(in)          :: loglike_pos   

        double precision :: baby_loglike  !loglikelihood of the baby point
        integer          :: nlive         !number of live points
        integer          :: baby_position !where to insert the baby point in the array


        nlive       = size(live_data,2)    ! size of the array

        baby_loglike = baby_point(loglike_pos) ! loglikelihood of the baby point

        ! search for the position with a binary search algorithm
        baby_position =  binary_search(1,nlive,baby_loglike,loglike_pos,live_data)

        ! Delete the lowest likelihood point by shifting the points in the data
        ! array from baby_position and below down by one
        live_data(:,1:baby_position-2) = live_data(:,2:baby_position-1)

        ! Add the new point
        live_data(:,baby_position-1) = baby_point(:)

    end subroutine insert_point_into_live_data


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
        incubating_data(:,incubating_position+1:) = incubating_data(:,incubating_position:)

        ! Add the new point
        incubating_data(:,incubating_position) = incubating_point(:)

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
        double precision, intent(inout), dimension(:,:) :: data_array

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
