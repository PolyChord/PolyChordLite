module chordal_module
    implicit none


    contains

    function ChordalSampling(settings,seed_point,min_max_array, M,feedback)  result(baby_point)
        use settings_module, only: program_settings
        use random_module, only: random_gaussian
        use model_module,  only: model, calculate_point
        use utils_module, only: logzero,stdout_unit

        implicit none

        ! ------- Inputs -------
        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M

        !> The seed point
        double precision, intent(in), dimension(M%nTotal)   :: seed_point

        !> The minimum and maximum values from each of the live points
        double precision, intent(in),    dimension(:,:)   :: min_max_array

        !> Optional argument to cause the sampler to print out relevent information
        integer, intent(in), optional :: feedback

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(M%nTotal)   :: baby_point


        ! ------- Local Variables -------
        double precision,    dimension(M%nDims,product(settings%num_chords))   :: nhats

        double precision  :: max_chord

        double precision :: step_length

        integer :: i


        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(stdout_unit,'( "Sampler    : Chordal" )')
                do i=1,maxval(M%grade)
                    write(stdout_unit,'( "  num chords(",I4,") = ",I8 )') i, settings%num_chords(i)
                end do
            end if
            return
        end if


        ! Start the baby point at the seed point
        baby_point = seed_point

        ! Set the number of likelihood evaluations to zero
        baby_point(M%nlike) = 0

        ! Re-scale the unit hypercube so that min->0, max->1 of min_max_array
        call re_scale(baby_point(M%h0:M%h1),min_max_array)

        ! Record the step length
        step_length = seed_point(M%last_chord)

        ! Initialise max_chord at 0
        max_chord = 0

        ! Generate the set of nhats to use
        nhats=0
        do i=1,M%nDims
            nhats(i,: product(settings%num_chords) : product(settings%num_chords)/product(settings%num_chords(:M%grade(i))) ) &
                = random_gaussian(product(settings%num_chords(:M%grade(i))))
        end do
        do i=1,product(settings%num_chords)
            nhats(:,i) = nhats(:,i)/sqrt(dot_product(nhats(:,i),nhats(:,i)))
        end do

        do i=1,product(settings%num_chords)
            ! Give the baby point the step length
            baby_point(M%last_chord) = step_length

            ! Generate a new random point along the chord defined by baby_point and nhat
            baby_point = random_chordal_point( nhats(:,i), baby_point,min_max_array, M)

            ! keep track of the largest chord
            max_chord = max(max_chord,baby_point(M%last_chord))
        end do

#ifdef MPI
        ! Make sure to hand back any incubator information which has likely been
        ! overwritten
        baby_point(M%incubator) = seed_point(M%incubator)
#endif

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_point(M%last_chord) = max_chord

        ! de-scale the unit hypercube so that 0->min, 1->max of min_max_array
        call de_scale(baby_point(M%h0:M%h1),min_max_array)

    end function ChordalSampling



    function random_chordal_point(nhat,seed_point,min_max_array,M) result(baby_point)
        use model_module,  only: model, calculate_point
        use utils_module,  only: logzero, distance
        use random_module, only: random_real
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: seed_point
        !> The minimum and maximum values from each of the live points
        double precision, intent(in),    dimension(M%nDims,2)   :: min_max_array

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: baby_point

        ! The upper bound
        double precision,    dimension(M%nTotal)   :: u_bound
        ! The lower bound
        double precision,    dimension(M%nTotal)   :: l_bound

        double precision :: trial_chord_length

        ! estimate at an appropriate chord
        trial_chord_length = seed_point(M%last_chord)

        ! record the number of likelihood calls
        u_bound(M%nlike) = seed_point(M%nlike)
        l_bound(M%nlike) = 0


        ! Select initial start and end points
        l_bound(M%h0:M%h1) = seed_point(M%h0:M%h1) - random_real() * trial_chord_length * nhat 
        u_bound(M%h0:M%h1) = l_bound(M%h0:M%h1) + trial_chord_length * nhat 

        ! Calculate initial likelihoods
        call de_scale(u_bound(M%h0:M%h1),min_max_array)
        call de_scale(l_bound(M%h0:M%h1),min_max_array)
        call calculate_point(M,u_bound)
        call calculate_point(M,l_bound)
        call re_scale(u_bound(M%h0:M%h1),min_max_array)
        call re_scale(l_bound(M%h0:M%h1),min_max_array)

        ! expand u_bound until it's outside the likelihood region
        do while(u_bound(M%l0) >= seed_point(M%l1) )
            u_bound(M%h0:M%h1) = u_bound(M%h0:M%h1) + nhat * trial_chord_length
            call de_scale(u_bound(M%h0:M%h1),min_max_array)
            call calculate_point(M,u_bound)
            call re_scale(u_bound(M%h0:M%h1),min_max_array)
        end do

        ! expand l_bound until it's outside the likelihood region
        do while(l_bound(M%l0) >= seed_point(M%l1) )
            l_bound(M%h0:M%h1) = l_bound(M%h0:M%h1) - nhat * trial_chord_length
            call de_scale(l_bound(M%h0:M%h1),min_max_array)
            call calculate_point(M,l_bound)
            call re_scale(l_bound(M%h0:M%h1),min_max_array)
        end do

        ! Sample within this bound
        baby_point = find_positive_within(l_bound,u_bound)

        ! Pass on the loglikelihood bound
        baby_point(M%l1) = seed_point(M%l1)

        ! Estimate the next appropriate chord
        baby_point(M%last_chord) = distance( u_bound(M%h0:M%h1),l_bound(M%h0:M%h1) )!distance( baby_point(M%h0:M%h1),seed_point(M%h0:M%h1) )

        contains

        recursive function find_positive_within(l_bound,u_bound) result(finish_point)
            implicit none
            !> The upper bound
            double precision, intent(inout), dimension(M%nTotal)   :: u_bound
            !> The lower bound
            double precision, intent(inout), dimension(M%nTotal)   :: l_bound

            ! The output finish point
            double precision,    dimension(M%nTotal)   :: finish_point

            double precision :: random_temp

            ! Draw a random point within l_bound and u_bound
            random_temp =random_real()
            finish_point(M%h0:M%h1) = l_bound(M%h0:M%h1)*(1d0-random_temp) + random_temp * u_bound(M%h0:M%h1)

            ! Pass on the number of likelihood calls that have been made
            finish_point(M%nlike) = l_bound(M%nlike) + u_bound(M%nlike)
            ! zero the likelihood calls for l_bound and u_bound, as these are
            ! now stored in point
            l_bound(M%nlike) = 0
            u_bound(M%nlike) = 0


            ! calculate the likelihood 
            call de_scale(finish_point(M%h0:M%h1),min_max_array)
            call calculate_point(M,finish_point)
            call re_scale(finish_point(M%h0:M%h1),min_max_array)

            ! If we're not within the likelihood bound then we need to sample further
            if( finish_point(M%l0) < seed_point(M%l1) ) then

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

