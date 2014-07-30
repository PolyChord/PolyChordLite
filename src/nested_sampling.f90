module nested_sampling_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSampling(M,settings)
        use model_module,    only: model
        use settings_module, only: program_settings
        use galileo_module,  only: GalileanSample

        implicit none
        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings


        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-physical coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data

        ! The new-born point
        double precision,    dimension(M%nTotal)   :: new_point

        double precision :: likelihood_bound


        ! Create initial live points
        write(*,*) 'generating live points'
        call GenerateLivePoints(live_data,M)

            !write(*,'(41F9.6 F16.8)') live_data
            !write(*,*) '----------------------------------------'

        write(*,*) 'started sampling'
        do while (.true.)

            likelihood_bound = live_data(M%l0,1)

            ! Generate a new point within the likelihood bounds
            call settings%sampler(new_point, live_data, likelihood_bound, M)

            ! Insert the point
            call insert_new_point(new_point,live_data)

            write(*,*) likelihood_bound

        end do

    end subroutine NestedSampling




    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    subroutine GenerateLivePoints(live_data,M)
        use model_module
        use random_module, only: random_coordinate
        implicit none

        type(model), intent(in) :: M     !> The model details (loglikelihood, priors, ndims etc...)

        !>live_data(:,i) constitutes the information in the ith live point in the unit hypercube: 
        !! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(:,:), intent(out) :: live_data

        ! Loop variable
        integer i_live

        integer nlive

        nlive = size(live_data,2)


        ! Generate nlive points
        do i_live=1, nlive

            ! Generate a random coordinate in the first nDims rows of live_data
            call random_coordinate( live_data(M%h0:M%h1,i_live) )

            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, live_data(:,i_live) )

            ! Calculate the likelihood and store it in the last index
            live_data(M%l0,i_live) = M%loglikelihood( live_data(M%p0:M%p1,i_live))

            ! Calculate the derived parameters
            call calculate_derived_parameters( M, live_data(:,i_live) )

            ! Calculate the derived parameters
            !> @todo Need to calculate the derived parameters

        end do

        ! Sort them in order of likelihood, first argument lowest, last argument highest
        call quick_sort(live_data)



    end subroutine 



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
        else if( loglike < data_array(loglike_pos,1) ) then
            new_position = 1
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
