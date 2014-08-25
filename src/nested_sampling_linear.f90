module nested_sampling_linear_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSamplingL(M,settings)
        use model_module,      only: model
        use utils_module,      only: logzero,loginf,DBL_FMT,read_resume_unit,stdout_unit
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

        double precision,    dimension(M%nDims,2)   :: min_max_array

        double precision, allocatable, dimension(:,:) :: posterior_array
        double precision, dimension(M%nDims+2) :: posterior_point
        integer :: nposterior
        integer :: insertion_index(1)

        logical :: more_samples_needed

        ! The new-born baby point
        double precision,    dimension(M%nTotal)   :: baby_point
        double precision :: baby_likelihood

        ! The recently dead point
        double precision,    dimension(M%nTotal)   :: late_point
        double precision :: late_likelihood
        double precision :: late_logweight
        integer :: late_index

        ! Point to seed a new one from
        double precision,    dimension(M%nTotal)   :: seed_point


        ! Evidence info
        double precision, dimension(6) :: evidence_vec


        logical :: resume=.false.
        ! Means to be calculated
        double precision :: mean_likelihood_calls
        integer :: total_likelihood_calls

        integer :: ndead

        double precision :: lognlive 
        double precision :: lognlivep1 
        double precision :: logminimumweight

        call write_opening_statement(M,settings) 
        inquire(file=trim(settings%file_root)//'.resume',exist=resume)
        if(resume) write(stdout_unit,'("Resuming from previous run")')
        !======= 1) Initialisation =====================================
        ! (i)   generate initial live points by sampling
        !       randomly from the prior (i.e. unit hypercube)
        ! (ii)  Initialise all variables

        ! Create initial live points
        if(resume) then
            ! If there is a resume file present, then load the live points from that
            write(stdout_unit,'("Reading live data")')
            open(read_resume_unit,file=trim(settings%file_root)//'.resume',action='read')
            ! Read the index of the late point
            read(read_resume_unit,'(I)') late_index
            ! Read the live data
            read(read_resume_unit,'(<M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
        else
            call write_started_generating(settings%feedback)

            ! Otherwise generate them anew
            live_data = GenerateLivePoints(M,settings%nlive)

            ! Sort them in order of likelihood, first point lowest, last point highest
            late_index = create_linked_list(M,live_data,settings%nlive)

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
        !  (e) posterior_array        | Array of weighted posterior points

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

        ! (c) get number of dead points
        if(resume) then
            ! If resuming, then get the number of dead points from the resume file
            write(stdout_unit,'("Reading ndead")')
            read(read_resume_unit,'(I)') ndead
        else
            ! Otherwise no dead points originally
            ndead = 0
        end if

        ! (d) calculate the minimums and maximums of the live data
        min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
        min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

        ! Calculate these global variables so we don't need to again
        lognlive         = log(settings%nlive+0d0)
        lognlivep1       = log(settings%nlive+1d0)
        logminimumweight = log(settings%minimum_weight)


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
            read(read_resume_unit,'(<M%nDims+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)
        end if



        ! Close the resume file if we've openend it
        if(resume) close(read_resume_unit)













        call write_started_sampling(settings%feedback)

        ! definitely more samples needed than this
        more_samples_needed = .true.

        do while ( more_samples_needed )


            ! Record the point that's just died
            late_point = live_data(:,late_index)
            ndead = ndead + 1

            ! Get the likelihood contour
            late_likelihood = late_point(M%l0)

            ! Calculate the late logweight
            late_logweight = (ndead-1)*lognlive - ndead*lognlivep1 

            ! Select a seed point for the generator
            !  -excluding the points which have likelihoods equal to the
            !   loglikelihood bound
            seed_point(M%l0)=late_likelihood
            do while (seed_point(M%l0)<=late_likelihood)
                ! get a random number in [1,nlive]
                ! get this point from live_data 
                seed_point = live_data(:,random_integer(settings%nlive))
            end do

            ! Record the likelihood bound which this seed will generate from
            seed_point(M%l1) = late_likelihood

            ! Generate a new point within the likelihood bound of the late point
            baby_point = settings%sampler(seed_point, min_max_array, M)
            baby_likelihood  = baby_point(M%l0)

            ! update the mean number of likelihood calls
            mean_likelihood_calls = mean_likelihood_calls + (baby_point(M%d0) - late_point(M%d0) ) / (settings%nlive + 0d0)

            ! update the total number of likelihood calls
            total_likelihood_calls = total_likelihood_calls + baby_point(M%d0)

            ! Insert the new point

            call insert_into_live(M,settings%nlive,baby_point,live_data,late_index)

            ! Calculate the new evidence (and check to see if we're accurate enough)
            call settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed,evidence_vec)

            ! Halt if we've reached the desired maximum iterations
            if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

            ! Update the set of weighted posteriors
            if( settings%calculate_posterior .and. late_point(M%l0) + late_logweight - evidence_vec(1) > logminimumweight ) then

                ! calculate a new point for insertion
                posterior_point(1) = late_point(M%l0) + late_logweight
                posterior_point(2) = late_point(M%l0)
                posterior_point(3:) = late_point(M%p0:M%p1)

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


            ! Feedback to command line every nlive iterations
            if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                write(stdout_unit,'("ndead     = ", I20                  )') ndead
                write(stdout_unit,'("efficiency= ", F20.2                )') mean_likelihood_calls
                write(stdout_unit,'("log(Z)    = ", F20.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                write(stdout_unit,'("")')
            end if

            ! update the minimum and maximum values of the live points
            min_max_array(:,1) = minval(live_data(M%h0:M%h1,:),2)
            min_max_array(:,2) = maxval(live_data(M%h0:M%h1,:),2)

            ! Update the resume and posterior files every update_resume iterations, or at program termination
            if (mod(ndead,settings%update_resume) .eq. 0 .or.  more_samples_needed==.false.)  then
                if(settings%write_resume) call write_resume_file(settings,M,late_index,live_data,evidence_vec,ndead,nposterior,posterior_array) 
                if(settings%calculate_posterior) call write_posterior_file(settings,M,posterior_array,evidence_vec(1),nposterior)  
            end if

        end do

        call write_final_results(M,evidence_vec,ndead,total_likelihood_calls,settings%feedback)

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


    end function GenerateLivePoints



    function create_linked_list(M,live_data,nlive)  result(lowest_index)
        use model_module,      only: model
        implicit none

        !> The model details 
        type(model),            intent(in) :: M
        !> The number of live points
        integer, intent(in) :: nlive
        !> The live_data array
        double precision, intent(inout), dimension(M%nTotal,nlive) :: live_data



        double precision, dimension(nlive) :: loglikes

        ! The lowest index point
        integer :: lowest_index

        ! Temporary variables
        integer :: prev_index(1)
        integer :: next_index(1)
        double precision :: prev_loglike
        double precision :: next_loglike
        
        ! Loop variable
        integer :: i_live



        loglikes=live_data(M%l0,:)

        ! (1) Initialise the first point

        ! Find the lowest loglikelihood in live_data...
        next_loglike = minval(loglikes) 
        ! ... and its index (note that this is the function output)
        next_index   = minloc(loglikes)
        lowest_index = next_index(1)

        ! The pointer to the previous point for the lowest index points to nothing
        live_data(M%prevlive,lowest_index) = 0d0


        ! (2) Loop over the remaining live points.
        ! At each stage, we assign the previous point the address of the next
        ! point, and the next point the address of the previous point
        do i_live=2,nlive

            ! Catch the case of equal loglikelihoods (only really occurs for
            ! logzero i.e. unphysical points)
            ! In the event of several lowest values, minloc returns the lowest
            ! index. We can therefore check to see if there are multiple points
            ! with the same contour by searching in the index above the last
            ! one.
            if(count(loglikes(next_index(1)+1:)==next_loglike)>0) then

                ! find the index of the next loglikelihood
                prev_index   = next_index
                next_index   = minloc(loglikes(next_index(1)+1:),loglikes(next_index(1)+1:)==next_loglike)+next_index(1)

                ! Give the previous point the address of the next point
                live_data(M%nextlive,prev_index(1)) = next_index(1)
                ! Give the next point the address of the previous point
                live_data(M%prevlive,next_index(1)) = prev_index(1)


            else
                ! Find the next loglikelihood...
                prev_loglike = next_loglike
                next_loglike = minval(loglikes,loglikes>prev_loglike) 

                ! ... and the index of that likelihood
                prev_index   = next_index
                next_index   = minloc(loglikes,loglikes>prev_loglike) 

                ! Give the previous point the address of the next point
                live_data(M%nextlive,prev_index(1)) = next_index(1)
                ! Give the next point the address of the previous point
                live_data(M%prevlive,next_index(1)) = prev_index(1)
            end if
            
        end do

        ! Assign the highest point a pointer to say that it is at the end
        live_data(M%nextlive,next_index) = -1d0

    end function create_linked_list


    subroutine insert_into_live(M,nlive,baby_point,live_data,lowest_index)
        use model_module,      only: model
        implicit none
        !> The model details
        type(model),            intent(in) :: M
        !> The number of live points
        integer,intent(in) :: nlive
        !> The point to be inserted 
        double precision, intent(in),    dimension(M%nTotal)   :: baby_point
        !> The live data array to be inserted into
        double precision, intent(inout), dimension(M%nTotal,nlive) :: live_data
        !> The index to sort by (in this case it's M%l0)
        integer,intent(inout)          :: lowest_index

        ! The index of baby_point (set to be the old lowest_index)
        integer :: baby_index

        ! The index directly above baby_point
        integer :: next_index
        ! The index directly below baby_point
        integer :: prev_index

        ! The new lowest index once the dead point has been deleted
        integer :: new_lowest_index

        ! The position of baby_index in live_data will be the position of the
        ! late point, since we're about to delete it
        baby_index = lowest_index

        ! Initialise the insert index at the start point
        next_index = lowest_index

        ! Find the position to insert the baby_point by searching sequentially through the array
        do while( live_data(M%l0,next_index) < baby_point(M%l0) )
            prev_index = next_index
            next_index = nint(live_data(M%nextlive,next_index)) 
            if(next_index==-1) exit
        end do
        ! We have now returned the indices of the points surrounding baby_point


        if(prev_index/=lowest_index) then
            ! If the baby point is above the lowest index:

            ! Find the new lowest_index
            new_lowest_index=nint(live_data(M%nextlive,lowest_index))

            ! 'Delete' the late point by re-pointing the new lowest_index point to 0
            live_data( M%prevlive , new_lowest_index ) = 0d0

            ! Pass on the new lowest index
            lowest_index=new_lowest_index

        else
            ! Point the previous index at 0
            prev_index = 0d0
            ! tell the algorithm that the new lowest index is the baby index
            lowest_index=baby_index
        end if

        ! Insert baby_point at the same place as the late point
        live_data(:,baby_index) = baby_point


        ! Update the links of baby point
        live_data(M%nextlive,baby_index) = next_index
        live_data(M%prevlive,baby_index) = prev_index

        ! Update the links of the two points surrounding baby_point
        !  - point the point above baby_point to baby_point if baby_index is not
        !    the highest
        if(next_index/=-1) live_data(M%prevlive, next_index ) = baby_index
        !  - point the point below baby_point to baby_point if baby_index is not
        !    the lowest
        if(prev_index/=0) live_data(M%nextlive, prev_index) = baby_index


    end subroutine insert_into_live
        







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















end module nested_sampling_linear_module
