module nested_sampling_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSampling(M,settings)
        use model_module,    only: model
        use settings_module, only: program_settings
        use galileo_module,  only: GalileanSampling
        use feedback_module

        implicit none
        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings


        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data

        ! The new-born point
        double precision,    dimension(M%nTotal)   :: new_point

        double precision :: likelihood_bound,new_likelihood,late_likelihood
        double precision, dimension(2) :: evidence_vec

        double precision :: mean_loglike

        logical :: more_samples_needed


        integer :: ndead


        call write_started_generating(settings%feedback)

        ! Create initial live points
        live_data = GenerateLivePoints(M,settings)

        call write_finished_generating(settings%feedback)

        ! Compute the average loglikelihood and hand it to the evidence calculator
        mean_loglike = sum(exp(live_data(M%l0,:)))/settings%nlive
        evidence_vec = settings%evidence_calculator(mean_loglike,mean_loglike,-1,more_samples_needed)


        ! Count the number of dead points
        ndead = 0

        ! Definately need more samples than this
        more_samples_needed=.true.


            !write(*,'(41F9.6 F16.8)') live_data
            !write(*,*) '----------------------------------------'


        call write_started_sampling(settings%feedback)

        do while (more_samples_needed)

            likelihood_bound = live_data(M%l0,1)

            ! Generate a new point within the likelihood bounds
            new_point = settings%sampler(live_data, likelihood_bound, M)

            ndead = ndead + 1
            new_likelihood  = new_point(M%l0)
            late_likelihood = live_data(M%l0,1)

            ! Calculate the new evidence
            evidence_vec =  settings%evidence_calculator(new_likelihood,late_likelihood,ndead,more_samples_needed)

            if (settings%feedback>=1 .and. mod(ndead,settings%nlive) .eq.0 ) then
                if (settings%feedback>=2) then
                    write(*,'("new_point: (", <M%nDims>F10.5, ") ->", F12.5 )') new_point(M%p0:M%p1), new_point(M%l0)
                end if
                write(*,'("ndead   = ", I12                  )') ndead
                write(*,'("Z       = ", E12.5, " +/- ", E12.5)') evidence_vec(1:2)
                if (evidence_vec(1) > 0 ) then
                    write(*,'("log(Z)  = ", F12.5, " +/- ", F12.5)') log(evidence_vec(1)), evidence_vec(2)/evidence_vec(1) 
                end if
            end if

            ! Insert the new point
            call insert_new_point(new_point,live_data)

            if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.

        end do

        call write_final_results(evidence_vec,ndead,settings%feedback)

    end subroutine NestedSampling




    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    function GenerateLivePoints(M,settings) result(live_data)
        use model_module,    only: model, calculate_point
        use random_module,   only: random_hypercube_point
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
            live_data(:,i_live) = random_hypercube_point(M%nDims)

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
    subroutine insert_new_point(new_point,data_array)
        double precision, intent(in),    dimension(:)   :: new_point
        double precision, intent(inout), dimension(:,:) :: data_array

        double precision :: loglike 
        integer          :: loglike_pos
        integer          :: nlive

        integer :: new_position

        loglike_pos = size(data_array,1)
        nlive       = size(data_array,2)

        loglike = new_point(loglike_pos)

        if( loglike > data_array(loglike_pos,nlive) ) then
            new_position = nlive
        else
            new_position =  binary_search(1,nlive)
        end if

        ! Delete the lowest likelihood point by shifting the points in the data
        ! array from new_position and below down by one
        data_array(:,1:new_position-1) = data_array(:,2:new_position)

        ! Add the new point
        data_array(:,new_position) = new_point(:)

        contains

        !> Binary search algorithm
        recursive function binary_search(imin,imax) result(newpoint_pos)
            integer, intent(in) :: imin
            integer, intent(in) :: imax

            integer :: newpoint_pos

            newpoint_pos = (imin + imax)/2

            if ( newpoint_pos == imin ) then
                return
            else if( data_array(loglike_pos,newpoint_pos) > loglike ) then
                ! new_point in lower subset
                newpoint_pos = binary_search(imin,newpoint_pos)
            else 
                ! new_point in upper subset
                newpoint_pos = binary_search(newpoint_pos,imax)
            end if

        end function binary_search

    end subroutine insert_new_point















end module nested_sampling_module
