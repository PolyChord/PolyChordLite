module nested_sampling_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSampling(M,settings)
        use model_module, only: model
        use settings_module,  only: program_settings

        implicit none
        type(model),            intent(in) :: M
        type(program_settings), intent(in) :: settings


        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,settings%nlive) :: live_data

        ! The new-born point
        !double precision, dimension(M%nTotal) :: new_point

        ! Create initial live points
        call GenerateLivePoints(live_data,M,settings%nlive)


        ! Generate a new point within the likelihood bounds
        !call GallileanSample(new_point, live_data, M, likelihood_bound)





    end subroutine NestedSampling















    !> Generate an initial set of live points distributed uniformly in the unit hypercube
    subroutine GenerateLivePoints(live_data,M,nlive)
        use model_module, only: model
        use random_module, only: random_coordinate
        implicit none

        type(model), intent(in) :: M     !> The model details (loglikelihood, priors, ndims etc...)
        integer,     intent(in) :: nlive !> The number of live points to be generated

        !>live_data(:,i) constitutes the information in the ith live point in the unit hypercube: 
        !! ( <-hypercube coordinates->, <-derived parameters->, likelihood)
        double precision, dimension(M%nTotal,nlive), intent(out) :: live_data

        ! Loop variable
        integer i_live


        ! Generate nlive points
        do i_live=1, nlive

            ! Generate a random coordinate in the first nDims rows of live_data
            call random_coordinate( live_data(1:M%nDims,i_live), M%nDims)

            ! Calculate the likelihood and store it in the last index
            live_data(M%nTotal,i_live) = M%loglikelihood( live_data(1:M%nDims,i_live) )

            ! Calculate the derived parameters
            !> @todo Need to calculate the derived parameters

        end do

        ! Sort them in order of likelihood, first argument lowest, last argument
        ! highest
        call quick_sort(live_data)

        !write(*,'(2F9.6, F10.2)') live_data


    end subroutine 



    !> Sort nlive vectors in data_array(:,nlive) by the last element of each
    !! vector
    recursive subroutine quick_sort(data_array)

        double precision, dimension(:,:), intent(inout) :: data_array

        integer :: iq

        if( size(data_array, 2 ) > 1) then
            call quick_sort_partition(data_array,iq)
            call quick_sort( data_array(:,:iq-1) )
            call quick_sort( data_array(:,iq:)   )
        endif

        

    end subroutine quick_sort

    subroutine quick_sort_partition(data_array, marker)

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
                ! exchange A(i) and A(j)
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

    end subroutine quick_sort_partition


















end module nested_sampling_module
