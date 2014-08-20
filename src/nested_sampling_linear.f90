module nested_sampling_linear_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSamplingL(M,settings)
        use model_module,      only: model
        use utils_module,      only: logzero,loginf,DBL_FMT
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp
        use random_module,     only: random_integers
        use read_write_module, only: write_resume_file
        use feedback_module

        implicit none
        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings



        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data

        double precision,    dimension(M%nDims,2)   :: min_max_array

        logical :: more_samples_needed

        ! The new-born baby point
        double precision,    dimension(M%nTotal)   :: baby_point
        double precision :: baby_likelihood

        ! The recently dead point
        double precision,    dimension(M%nTotal)   :: late_point
        double precision :: late_likelihood

        ! Point to seed a new one from
        double precision,    dimension(M%nTotal)   :: seed_point
        ! temp variable for getting a random integer
        integer, dimension(1) :: point_number 


        ! Evidence info
        double precision, dimension(6) :: evidence_vec


        logical :: resume=.false.
        ! Means to be calculated
        double precision :: mean_likelihood_calls

        integer :: ndead



        call write_opening_statement(M,settings) 
        inquire(file=trim(settings%file_root)//'.resume',exist=resume)
        if(resume) write(*,*) "Resuming from previous run"
        !======= 1) Initialisation =====================================
        ! (i)   generate initial live points by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables

        ! Create initial live points
        if(resume) then
            ! If there is a resume file present, then load the live points from that
            write(*,*) "Reading live data"
            open(1000,file=trim(settings%file_root)//'.resume',action='read')
            read(1000,'(<M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
        else
            call write_started_generating(settings%feedback)

            ! Otherwise generate them anew
            live_data = GenerateLivePoints(M,settings%nlive)
            ! Sort them in order of likelihood, first point lowest, last point highest
            call quick_sort(live_data,M%l0)

            call write_finished_generating(settings%feedback)
        end if





        !~~~ (ii) Initialise all variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! There are several variables used throughout the rest of the
        ! algorithm that need to be initialised here
        !  (a) evidence_vec           | Vector containing the evidence, its error, and any other 
        !                             |  things that need to be accumulated over the run.
        !                             |  we need to initialise its sixth argument.
        !  (b) mean_likelihood_calls  | Mean number of likelihood calls over the past nlive iterations
        !  (c) ndead                  | Number of iterations/number of dead points
        !  (d) min_max_array          | Array of maximums and minimums for each coordinate - allows rescaling

        ! (a)
        if(resume) then
            ! If resuming, get the accumulated stats to calculate the
            ! evidence from the resume file
            read(1000,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
        else
            ! Otherwise compute the average loglikelihood and initialise the evidence vector accordingly
            evidence_vec = logzero
            evidence_vec(6) = logsumexp(live_data(M%l0,:)) - log(settings%nlive+0d0)
        end if

        ! (b) initialise the mean number of likelihood calls to 1
        mean_likelihood_calls = sum(live_data(M%d0,:))/settings%nlive

        ! (c) get number of dead points
        if(resume) then
            ! If resuming, then get the number of dead points from the resume file
            write(*,*) "Reading ndead"
            read(1000,'(I)') ndead
            close(1000)
        else
            ! Otherwise no dead points originally
            ndead = 0
        end if

        ! If we're saving the dead points, open the relevant file
        if (settings%save_dead) then
            if(resume) then
                open(unit=222, file='dead_points.dat',position='append',action='write')
            else
                open(unit=222, file='dead_points.dat',action='write')
            end if
        end if

        ! (d) calculate the minimums and maximums of the live data
        min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
        min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

        call write_started_sampling(settings%feedback)

        ! definitely more samples needed than this
        more_samples_needed = .true.

        do while ( more_samples_needed )


            ! Record the point that's just died
            late_point = live_data(:,1)
            ndead = ndead + 1

            ! Get the likelihood contour
            late_likelihood = late_point(M%l0)

            ! Select a seed point for the generator
            !  -excluding the points which have likelihoods equal to the
            !   loglikelihood bound
            seed_point(M%l0)=late_likelihood
            do while (seed_point(M%l0)<=late_likelihood)
                ! get a random number in [1,nlive]
                point_number = random_integers(1,settings%nlive) 
                ! get this point from live_data 
                seed_point = live_data(:,point_number(1))
            end do

            ! Record the likelihood bound which this seed will generate from
            seed_point(M%l1) = late_likelihood

            ! Generate a new point within the likelihood bound of the late point
            baby_point = settings%sampler(seed_point, min_max_array, M)
            baby_likelihood  = baby_point(M%l0)

            ! update the mean number of likelihood calls
            mean_likelihood_calls = mean_likelihood_calls + (baby_point(M%d0) - late_point(M%d0) ) / (settings%nlive + 0d0)

            ! Insert the new point
            call insert_point_into_live_data(baby_point,live_data,M%l0)

            ! Calculate the new evidence (and check to see if we're accurate enough)
            call settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed,evidence_vec)

            ! Halt if we've reached the desired maximum iterations
            if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

            ! Write the dead points to a file if desired
            if (settings%save_dead) write(222,'(<M%nTotal+1>E17.9)') exp(ndead*log(settings%nlive/(settings%nlive+1d0))), late_point

            ! Feedback to command line every nlive iterations
            if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                write(*,'("ndead     = ", I20                  )') ndead
                write(*,'("efficiency= ", F20.2                )') mean_likelihood_calls
                write(*,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                write(*,*)
                call write_resume_file(settings,M,live_data,evidence_vec,ndead) 
            end if

        end do

        ! update the minimum and maximum values of the live points
        min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
        min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

        ! close the dead points file if we're done
        if (settings%save_dead) close(222)

        call write_final_results(M,evidence_vec,ndead,settings%feedback)

    end subroutine NestedSamplingL




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
        ! array from baby_position and below down by one and add the new point
        live_data(:,:baby_position-1) = eoshift(live_data(:,:baby_position-1),dim=2,shift=+1,boundary=baby_point)

    end subroutine insert_point_into_live_data


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















end module nested_sampling_linear_module
