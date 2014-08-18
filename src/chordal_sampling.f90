module chordal_module
    implicit none


    contains

    function ChordalSampling(settings,seed_point, loglikelihood_bound,min_max_array, M,feedback)  result(baby_point)
        use settings_module, only: program_settings
        use random_module, only: random_skewed_direction,random_direction,random_reals,random_integers
        use model_module,  only: model, calculate_point
        use utils_module, only: logzero

        implicit none

        ! ------- Inputs -------
        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M

        !> The seed point
        double precision, intent(in), dimension(M%nTotal)   :: seed_point

        !> The current loglikelihood bound
        double precision, intent(in) :: loglikelihood_bound

        !> The minimum and maximum values from each of the live points
        double precision, intent(in),    dimension(:,:)   :: min_max_array

        !> Optional argument to cause the sampler to print out relevent information
        integer, intent(in), optional :: feedback

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: baby_point


        ! ------- Local Variables -------
        double precision,    dimension(M%nDims)   :: nhat

        
        integer :: i


        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Sampler    : Chordal" )')
                write(*,'( "  num chords = ",I4 )') settings%num_chords
            end if
            return
        end if


        ! Start the baby point at the seed point
        baby_point = seed_point

        ! Set the number of likelihood evaluations to zero
        baby_point(M%d0) = 0

        ! Re-scale the unit hypercube so that min->0, max->1 of min_max_array
        call re_scale(baby_point(M%h0:M%h1),min_max_array)

        do i=1,settings%num_chords

            ! Get a new random direction
            nhat = random_direction(M%nDims) 

            ! Generate a new random point along the chord defined by baby_point and nhat
            baby_point = random_chordal_point( nhat, baby_point, loglikelihood_bound,min_max_array, M)
        end do

        ! de-scale the unit hypercube so that 0->min, 1->max of min_max_array
        call de_scale(baby_point(M%h0:M%h1),min_max_array)

    end function ChordalSampling



    function random_chordal_point(nhat,seed_point,loglikelihood_bound,min_max_array,M) result(baby_point)
        use model_module,  only: model, calculate_point
        use utils_module,  only: logzero, distance
        use random_module, only: random_reals
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: seed_point
        !> The root value to find
        double precision, intent(in) :: loglikelihood_bound
        !> The minimum and maximum values from each of the live points
        double precision, intent(in),    dimension(M%nDims,2)   :: min_max_array

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: baby_point

        ! The upper bound
        double precision,    dimension(M%nTotal)   :: u_bound
        ! The lower bound
        double precision,    dimension(M%nTotal)   :: l_bound

        double precision :: trial_chord_length


        ! record the number of likelihood calls
        u_bound(M%d0) = seed_point(M%d0)
        l_bound(M%d0) = 0

        ! set the likelihoods of start bounds so that the loops below are entered
        u_bound(M%l0) = loglikelihood_bound
        l_bound(M%l0) = loglikelihood_bound

        ! set the trial chord length at half the bound of the old length
        trial_chord_length = seed_point(M%d0+1)/2d0 

        ! Expand the region by doubling until u_bound has a loglikelihood less
        ! than the log likelihood bound
        do while(u_bound(M%l0) >= loglikelihood_bound )
            trial_chord_length = 2d0*trial_chord_length
            u_bound(M%h0:M%h1) = seed_point(M%h0:M%h1) + nhat * trial_chord_length
            call de_scale(u_bound(M%h0:M%h1),min_max_array)
            call calculate_point(M,u_bound)
            call re_scale(u_bound(M%h0:M%h1),min_max_array)
        end do

        ! repeat for the l_bound
        trial_chord_length = trial_chord_length/2
        do while(l_bound(M%l0) >= loglikelihood_bound )
            trial_chord_length = 2*trial_chord_length
            l_bound(M%h0:M%h1) = seed_point(M%h0:M%h1) - nhat * trial_chord_length
            call de_scale(l_bound(M%h0:M%h1),min_max_array)
            call calculate_point(M,l_bound)
            call re_scale(l_bound(M%h0:M%h1),min_max_array)
        end do

        baby_point = find_positive_within(l_bound,u_bound)

        ! Store the new length
        baby_point(M%d0+1) = max(&
            distance(u_bound(M%h0:M%h1),seed_point(M%h0:M%h1)),&
            distance(l_bound(M%h0:M%h1),seed_point(M%h0:M%h1))&
            ) 


        contains

        recursive function find_positive_within(l_bound,u_bound) result(finish_point)
            implicit none
            !> The upper bound
            double precision, intent(inout), dimension(M%nTotal)   :: u_bound
            !> The lower bound
            double precision, intent(inout), dimension(M%nTotal)   :: l_bound

            ! The output finish point
            double precision,    dimension(M%nTotal)   :: finish_point

            double precision,dimension(1) :: random_temp

            ! Draw a random point within l_bound and u_bound
            random_temp =random_reals(1)
            finish_point(M%h0:M%h1) = l_bound(M%h0:M%h1)*(1d0-random_temp(1)) + random_temp(1) * u_bound(M%h0:M%h1)

            ! Pass on the number of likelihood calls that have been made
            finish_point(M%d0) = l_bound(M%d0) + u_bound(M%d0)
            ! zero the likelihood calls for l_bound and u_bound, as these are
            ! now stored in point
            l_bound(M%d0) = 0
            u_bound(M%d0) = 0


            ! calculate the likelihood 
            call de_scale(finish_point(M%h0:M%h1),min_max_array)
            call calculate_point(M,finish_point)
            call re_scale(finish_point(M%h0:M%h1),min_max_array)

            ! If we're not within the likelihood bound then we need to sample further
            if( finish_point(M%l0) < loglikelihood_bound ) then

                if ( dot_product(finish_point(M%h0:M%h1)-seed_point(M%h0:M%h1),nhat) > 0d0 ) then
                    ! If finish_point is on the u_bound side of seed_point, then
                    ! contract u_bound
                    u_bound = finish_point
                else
                    ! If finish_point is on the l_bound side of seed_point, then
                    ! contract l_bound
                    l_bound = finish_point
                end if

                ! Call the function again
                finish_point = find_positive_within(l_bound,u_bound)

            end if
            ! otherwise finish_point is returned

        end function find_positive_within


    end function random_chordal_point


    !> [min,max] -> [0,1]
    subroutine re_scale(hypercube_coord,min_max_array)
        implicit none
        !> The hypercube coordinate to scale down by min_max_array
        double precision, intent(inout),    dimension(:)   :: hypercube_coord
        !> The minimum and maximum values from each of the live points
        double precision, intent(in),    dimension(:,:)   :: min_max_array

        hypercube_coord(:) = (hypercube_coord(:)- min_max_array(:,1)) / (min_max_array(:,2)-min_max_array(:,1))
    end subroutine re_scale

    subroutine de_scale(re_scaled_hypercube_coord,min_max_array)
        implicit none
        !> The hypercube coordinate to scale down by min_max_array
        double precision, intent(inout),    dimension(:)   :: re_scaled_hypercube_coord
        !> The minimum and maximum values from each of the live points
        double precision, intent(in),    dimension(:,:)   :: min_max_array

        re_scaled_hypercube_coord(:) = min_max_array(:,1) + (min_max_array(:,2)-min_max_array(:,1))*re_scaled_hypercube_coord(:)
    end subroutine de_scale

end module chordal_module

