module nested_sampling_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSampling(M,settings)
        use model_module,    only: model
        use utils_module,    only: logzero
        use settings_module, only: program_settings
        use utils_module,    only: logsumexp
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

        double precision :: likelihood_bound
        double precision :: baby_likelihood
        double precision :: late_likelihood

        double precision, dimension(2) :: evidence_vec

        ! Means to be calculated
        double precision :: mean_likelihood_calls

        logical :: more_samples_needed


        integer :: ndead


        !------- 1) Initialisation ---------------------------------
        ! Need to initialise
        !  * live_data
        !  * mean log likelihood to be passed to the evidence calculator
        !  * more_samples_needed
        !  * mean_likelihood_calls
        !  * ndead

        ! Create initial live points
        call write_started_generating(settings%feedback)
        live_data = GenerateLivePoints(M,settings)
        call write_finished_generating(settings%feedback)

        ! Set the number of likelihood calls for each point to 1
        live_data(M%d0,:) = 1

        ! Definitely need more samples than this
        more_samples_needed=.true.

        ! Compute the average loglikelihood and hand it to the evidence calculator via the second argument,
        ! and pass the maximum likelihood value in the first argument
        ! (note the negative number of dead points in the third argument to trigger initialisation). 
        evidence_vec = settings%evidence_calculator(                                       &
                                  live_data(M%l0,settings%nlive),                          &
                                  logsumexp(live_data(M%l0,:)) - log(settings%nlive+0d0) , &
                                  -1,more_samples_needed)

        ! initialise the mean number of likelihood calls to 1
        mean_likelihood_calls = 1d0

        ! Set the initial trial values of the chords as the diagonal of the hypercube
        live_data(M%d0+1,:) = sqrt(M%nDims+0d0)

        ! no dead points originally
        ndead = 0


        if (settings%save_dead) open(unit=222, file='dead_points.dat')

        call write_started_sampling(settings%feedback)

        do while (more_samples_needed)

            ! Get the likelihood contour
            likelihood_bound = live_data(M%l0,1)

            ! Select a seed point

            ! Generate a new point within the likelihood bounds
            baby_point = settings%sampler(live_data, likelihood_bound, M)
            late_point = live_data(:,1)

            ndead = ndead + 1
            baby_likelihood  = baby_point(M%l0)
            late_likelihood = late_point(M%l0)

            ! Calculate the new evidence
            evidence_vec =  settings%evidence_calculator(baby_likelihood,late_likelihood,ndead,more_samples_needed)
            ! update the mean number of likelihood calls
            mean_likelihood_calls = mean_likelihood_calls + (baby_point(M%d0) - late_point(M%d0) ) / (settings%nlive + 0d0)

            ! Insert the new point
            call insert_baby_point(baby_point,live_data)


            if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then

                write(*,'("ndead     = ", I12                  )') ndead
                write(*,'("efficiency= ", F12.2                )') mean_likelihood_calls
                write(*,'("log(Z)    = ", F12.5, " +/- ", F12.5)') evidence_vec(1), exp(0.5*evidence_vec(2)-evidence_vec(1)) 
                write(*,*)
            end if


            if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

            if (settings%save_dead) write(222,'(<M%nTotal+1>E17.9)') exp(ndead*log(settings%nlive/(settings%nlive+1d0))), late_point

        end do

        !-----------------------------------------------
        call write_final_results(M,evidence_vec,ndead,settings%feedback)
        !-----------------------------------------------
        if (settings%save_dead) close(222)

    end subroutine NestedSampling




    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePoints(M,settings) result(live_data)
        use model_module,    only: model, calculate_point
        use random_module,   only: random_reals
        use settings_module, only: program_settings
        use feedback_module, only: write_generating_live_points 
        implicit none

        !> The model details (loglikelihood, priors, ndims etc...)
        type(model), intent(in) :: M

        !> The program settings
        type(program_settings), intent(in) :: settings

        !live_data(:,i) constitutes the information in the ith live point in the unit hypercube: 
        ! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data

        ! Loop variable
        integer i_live


        ! Generate nlive points
        do i_live=1, settings%nlive

            ! Generate a random coordinate in the first nDims rows of live_data
            live_data(:,i_live) = random_reals(M%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, live_data(:,i_live) )

            call write_generating_live_points(settings%feedback,i_live,settings%nlive)


        end do

        ! Sort them in order of likelihood, first argument lowest, last argument highest
        call quick_sort(live_data)



    end function GenerateLivePoints



    !> Sort nlive vectors in data_array(:,nlive) by the last element of each
    !! vector
    recursive subroutine quick_sort(data_array)

        !> an array of nlive  (coord,derived,like) to be sorted by the
        !! likelihood value at the end of each vector
        double precision, dimension(:,:), intent(inout) :: data_array

        integer :: iq

        if( size(data_array, 2 ) > 1) then
            call partition(data_array,iq)
            call quick_sort( data_array(:,:iq-1) )
            call quick_sort( data_array(:,iq:)   )
        endif
        contains

        !> @todo comment this quick_sort partition
        subroutine partition(data_array, marker)

            double precision, intent(inout), dimension(:,:) :: data_array
            integer, intent(out) :: marker
            integer :: i, j
            double precision, dimension(size(data_array,1)) :: temp,x
            integer :: loglike_pos

            loglike_pos = size(data_array,1)

            x(:) = data_array(:,1)
            i= 0
            j= size( data_array, 2 ) + 1

            do
                j = j-1
                do
                    if (data_array(loglike_pos,j) <= x(loglike_pos)) exit
                    j = j-1
                end do

                i = i+1
                do
                    if (data_array(loglike_pos,i) >= x(loglike_pos)) exit
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





    !> Insert the new point into the data array
    !!
    !! Since the data array is already sorted, one can insert a new point using
    !! binary search insertion algorithm, which is contained in here
    subroutine insert_baby_point(baby_point,data_array)
        double precision, intent(in),    dimension(:)   :: baby_point
        double precision, intent(inout), dimension(:,:) :: data_array

        double precision :: loglike 
        integer          :: loglike_pos
        integer          :: nlive

        integer :: baby_position

        loglike_pos = size(data_array,1)
        nlive       = size(data_array,2)

        loglike = baby_point(loglike_pos)

        if( loglike > data_array(loglike_pos,nlive) ) then
            baby_position = nlive
        else
            baby_position =  binary_search(1,nlive)
        end if

        ! Delete the lowest likelihood point by shifting the points in the data
        ! array from baby_position and below down by one
        data_array(:,1:baby_position-1) = data_array(:,2:baby_position)

        ! Add the new point
        data_array(:,baby_position) = baby_point(:)

        contains

        !> Binary search algorithm
        recursive function binary_search(imin,imax) result(baby_point_pos)
            integer, intent(in) :: imin
            integer, intent(in) :: imax

            integer :: baby_point_pos

            baby_point_pos = (imin + imax)/2

            if ( baby_point_pos == imin ) then
                return
            else if( data_array(loglike_pos,baby_point_pos) > loglike ) then
                ! baby_point in lower subset
                baby_point_pos = binary_search(imin,baby_point_pos)
            else 
                ! baby_point in upper subset
                baby_point_pos = binary_search(baby_point_pos,imax)
            end if

        end function binary_search

    end subroutine insert_baby_point















end module nested_sampling_module
