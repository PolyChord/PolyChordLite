module chordal_module
    implicit none

    contains

    function SliceSampling(loglikelihood,priors,settings,live_data,seed_point)  result(baby_point)
        use priors_module, only: prior
        use settings_module, only: program_settings

        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        ! ------- Inputs -------
        !> The prior information
        type(prior), dimension(:), intent(in) :: priors

        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> The seed point
        double precision, intent(in), dimension(:)   :: seed_point

        !> The directions of the chords
        double precision, intent(in), allocatable, dimension(:,:) :: live_data

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(size(seed_point))   :: baby_point


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat

        double precision  :: max_chord

        double precision :: step_length

        integer :: i_chords


        ! Start the baby point at the seed point
        baby_point = seed_point

        ! Set the number of likelihood evaluations to zero
        baby_point(settings%nlike) = 0

        ! Record the step length
        step_length = seed_point(settings%last_chord)

        ! Initialise max_chord at 0
        max_chord = 0

        do i_chords=1,settings%chain_length
            ! Give the baby point the step length
            baby_point(settings%last_chord) = step_length

            ! Get a new random direction
            call settings%get_nhat(live_data,nhat)

            ! Generate a new random point along the chord defined by baby_point and nhat
            baby_point = slice_sample(loglikelihood,priors, nhat, baby_point, settings)

            ! keep track of the largest chord
            max_chord = max(max_chord,baby_point(settings%last_chord))
        end do

        ! Make sure to hand back any incubator information which has likely been
        ! overwritten (this is only relevent in parallel mode)
        baby_point(settings%daughter) = seed_point(settings%daughter)

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_point(settings%last_chord) = max_chord

    end function SliceSampling

    subroutine HitAndRun(settings,live_data,nhat)
        use settings_module, only: program_settings
        use random_module, only: random_direction
        implicit none

        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> Any data from the live points which is needed
        double precision, intent(in), allocatable, dimension(:,:) :: live_data

        ! ------- Outputs -------
        !> The newly generated point
        double precision, intent(out),   dimension(:)     :: nhat

        ! Get a new isotropic random direction
        nhat = random_direction(settings%nDims)

    end subroutine HitAndRun

    subroutine Adaptive_Parallel(settings,live_data,nhat)
        use settings_module, only: program_settings
        use random_module, only: random_distinct_integers
        use utils_module, only: mod2
        implicit none

        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> Any data from the live points which is needed
        double precision, intent(in), allocatable, dimension(:,:) :: live_data

        ! ------- Outputs -------
        !> The newly generated point
        double precision, intent(out),   dimension(:)     :: nhat
        integer,             dimension(2)                :: nhat_indices(2)
        integer :: nlive

        nlive = count( nint(live_data(settings%nDims+1,:)) == 1 )

        ! Get two distinct indices 
        nhat_indices = random_distinct_integers(nlive,2)

        ! Define the direction as the difference between these
        nhat = live_data(:,nhat_indices(1)) - live_data(:,nhat_indices(2))

        ! Normalise nhat
        nhat = nhat/sqrt(mod2(nhat))


    end subroutine Adaptive_Parallel



    function SliceSampling_Adaptive_Graded(loglikelihood,priors,settings,live_data,seed_point)  result(baby_point)
        use priors_module, only: prior
        use settings_module, only: program_settings
        use random_module, only: random_direction

        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        ! ------- Inputs -------
        !> The prior information
        type(prior), dimension(:), intent(in) :: priors

        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> The seed point
        double precision, intent(in), dimension(:)   :: seed_point

        !> The directions of the chords
        double precision, intent(in), allocatable, dimension(:,:) :: live_data

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(size(seed_point))   :: baby_point


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat

        double precision  :: max_chord

        double precision :: step_length

        integer :: i_chords


        ! Start the baby point at the seed point
        baby_point = seed_point

        ! Set the number of likelihood evaluations to zero
        baby_point(settings%nlike) = 0

        baby_point = run_sub_chain(loglikelihood,priors,settings,live_data,seed_point,1)


    end function SliceSampling_Adaptive_Graded

    recursive function run_sub_chain(loglikelihood,priors,settings,live_data,seed_point,grade)  result(baby_point)
        use priors_module, only: prior
        use settings_module, only: program_settings
        use random_module, only: random_gaussian,random_distinct_integers
        use utils_module, only: mod2

        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        ! ------- Inputs -------
        !> The prior information
        type(prior), dimension(:), intent(in) :: priors

        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> The seed point
        double precision, intent(in), dimension(:)   :: seed_point

        !> The directions of the chords
        double precision, intent(in), allocatable, dimension(:,:) :: live_data

        integer, intent(in) :: grade

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(size(seed_point))   :: baby_point


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat
        integer,             dimension(2)                :: nhat_indices(2)

        double precision  :: max_chord

        double precision :: step_length

        integer :: max_grade
        integer :: num_steps_taken

        integer :: i_chords

        integer :: nlive


        ! Find the maximum grade
        max_grade = maxval(settings%grade)

        ! Start the baby point at the seed point
        baby_point = seed_point

        ! Record the step length
        step_length = seed_point(settings%last_chord)

        ! Initialise max_chord at 0
        max_chord = 0

        ! We've initially taken no steps
        num_steps_taken = 0

        do while(.true.)

            ! If we're not at the maximum grade, then recurse down one more
            if(grade<max_grade) then
                baby_point = run_sub_chain(loglikelihood,priors,settings,live_data,baby_point,grade+1)
            end if
            ! keep track of the largest chord
            max_chord = max(max_chord,baby_point(settings%last_chord))

            ! Exit if we've taken enough steps
            if(num_steps_taken >= settings%chain_lengths(grade)) exit

            ! Give the baby point the step length
            baby_point(settings%last_chord) = step_length

            if(grade==1) then
                nlive = count( nint(live_data(settings%nDims+1,:)) == 1 )
                ! Get two distinct indices 
                nhat_indices = random_distinct_integers(nlive,2)

                ! Define the direction as the difference between these
                nhat = live_data(:,nhat_indices(1)) - live_data(:,nhat_indices(2))

                ! Normalise nhat
                nhat = nhat/sqrt(mod2(nhat))
            else
                ! Get a set of nDims gaussian random variables
                nhat = random_gaussian(settings%nDims)
                ! Zero out the unused dimensions
                where( settings%grade <grade )  nhat=0
                ! Normalise nhat
                nhat = nhat/sqrt(mod2(nhat))
            end if

            ! Generate a new random point along the chord defined by baby_point and nhat
            baby_point = slice_sample(loglikelihood,priors, nhat, baby_point, settings)

            ! Iterate the number of steps taken
            num_steps_taken = num_steps_taken+1

            ! keep track of the largest chord
            max_chord = max(max_chord,baby_point(settings%last_chord))
        end do

        ! Make sure to hand back any incubator information which has likely been
        ! overwritten (this is only relevent in parallel mode)
        baby_point(settings%daughter) = seed_point(settings%daughter)

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_point(settings%last_chord) = max_chord

    end function run_sub_chain




    !> Slice sample along the direction defined by nhat, starting at x0.
    !!
    !! We loosely follow the notation found in Neal's landmark paper:
    !! [Neal 2003](http://projecteuclid.org/download/pdf_1/euclid.aos/1056562461).
    !!
    !! We use the 'stepping out proceedure' to expand the bounds R and L until
    !! they lie outside the iso-likelihood contour, and then contract inwards
    !! using the  'shrinkage' procedure
    !!
    !! Each seed point x0 contains an initial estimate of the width w.
    !!
    function slice_sample(loglikelihood,priors,nhat,x0,S) result(baby_point)
        use settings_module, only: program_settings
        use priors_module, only: prior
        use utils_module,  only: logzero, distance
        use random_module, only: random_real
        use calculate_module, only: calculate_point
        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        !> The prior information
        type(prior), dimension(:), intent(in) :: priors
        !> program settings
        type(program_settings), intent(in) :: S
        !> The direction to search for the root in
        double precision, intent(in),    dimension(S%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(S%nTotal)   :: x0

        ! The output finish point
        double precision,    dimension(S%nTotal)   :: baby_point

        ! The upper bound
        double precision,    dimension(S%nTotal)   :: R
        ! The lower bound
        double precision,    dimension(S%nTotal)   :: L

        double precision :: w

        ! estimate at an appropriate chord
        w = x0(S%last_chord)

        ! record the number of likelihood calls
        R(S%nlike) = x0(S%nlike)
        L(S%nlike) = 0


        ! Select initial start and end points
        L(S%h0:S%h1) = x0(S%h0:S%h1) - random_real() * w * nhat 
        R(S%h0:S%h1) = L(S%h0:S%h1) + w * nhat 

        ! Calculate initial likelihoods
        call calculate_point(loglikelihood,priors,R,S)
        call calculate_point(loglikelihood,priors,L,S)

        ! expand R until it's outside the likelihood region
        do while(R(S%l0) >= x0(S%l1) .and. R(S%l0) > logzero )
            R(S%h0:S%h1) = R(S%h0:S%h1) + nhat * w
            call calculate_point(loglikelihood,priors,R,S)
        end do

        ! expand L until it's outside the likelihood region
        do while(L(S%l0) >= x0(S%l1) .and. L(S%l0) > logzero )
            L(S%h0:S%h1) = L(S%h0:S%h1) - nhat * w
            call calculate_point(loglikelihood,priors,L,S)
        end do

        ! Sample within this bound
        baby_point = find_positive_within(L,R)

        ! Pass on the loglikelihood bound
        baby_point(S%l1) = x0(S%l1)

        ! Estimate the next appropriate chord
        baby_point(S%last_chord) = distance( R(S%h0:S%h1),L(S%h0:S%h1) )/2d0

        contains

        recursive function find_positive_within(L,R) result(x1)
            implicit none
            !> The upper bound
            double precision, intent(inout), dimension(S%nTotal)   :: R
            !> The lower bound
            double precision, intent(inout), dimension(S%nTotal)   :: L

            ! The output finish point
            double precision,    dimension(S%nTotal)   :: x1

            double precision :: random_temp

            ! Draw a random point within L and R
            random_temp =random_real()
            x1(S%h0:S%h1) = L(S%h0:S%h1)*(1d0-random_temp) + random_temp * R(S%h0:S%h1)

            ! Pass on the number of likelihood calls that have been made
            x1(S%nlike) = L(S%nlike) + R(S%nlike)
            ! zero the likelihood calls for L and R, as these are
            ! now stored in point
            L(S%nlike) = 0
            R(S%nlike) = 0


            ! calculate the likelihood 
            call calculate_point(loglikelihood,priors,x1,S)

            ! If we're not within the likelihood bound then we need to sample further
            if( x1(S%l0) < x0(S%l1) .or. x1(S%l0) <= logzero ) then

                if ( dot_product(x1(S%h0:S%h1)-x0(S%h0:S%h1),nhat) > 0d0 ) then
                    ! If x1 is on the R side of x0, then
                    ! contract R
                    R = x1
                else
                    ! If x1 is on the L side of x0, then
                    ! contract L
                    L = x1
                end if

                ! Call the function again
                x1 = find_positive_within(L,R)

            end if
            ! otherwise x1 is returned

        end function find_positive_within


    end function slice_sample


    ! Information extractors



    !> Blanck extractor
    !!
    !! Use this if your sampling method doesn't need any information from the
    !! live points
    subroutine no_processing(settings,live_points,live_data,loglikelihood_bound)
        use settings_module, only: program_settings
        implicit none

        ! ------- Inputs -------
        !> program settings 
        class(program_settings), intent(in) :: settings

        !> The live points
        double precision, intent(in), dimension(:,:) :: live_points

        !> The loglikelihood bound to define the live points
        double precision, intent(in) :: loglikelihood_bound

        ! ------ Result -----------
        !> The processed data
        double precision, intent(out), allocatable, dimension(:,:) :: live_data

        if(.not. allocated(live_data)) allocate(live_data(0,0))


    end subroutine no_processing


    subroutine get_live_coordinates(settings,live_points,live_data,loglikelihood_bound)
        use utils_module, only: flag_waiting
        use settings_module, only: program_settings
        implicit none

        ! ------- Inputs -------
        !> program settings 
        class(program_settings), intent(in) :: settings

        !> The live points
        double precision, intent(in), dimension(:,:) :: live_points

        !> The loglikelihood bound to define the live points
        double precision, intent(in) :: loglikelihood_bound

        ! ------ Result -----------
        !> The processed data
        double precision, intent(out), allocatable, dimension(:,:) :: live_data

        integer :: i_data
        integer :: i_live

        if(.not. allocated(live_data)) allocate(live_data(settings%nDims+1,settings%nlive))

        live_data(:,settings%nDims+1) = 0

        i_data=1
        do i_live=1,size(live_points,2)
            if( live_points(settings%l1,i_live)<=loglikelihood_bound .and. live_points(settings%daughter,i_live)>=flag_waiting )  then
                ! Extract the coordinates
                live_data(:settings%nDims,i_data) = live_points(settings%h0:settings%h1,i_live) 
                ! Recort that this is a point to be drawn from
                live_data(settings%nDims+1,i_data) = 1
                i_data = i_data+1
            end if
            if(i_data>settings%nlive) exit
        end do

    end subroutine get_live_coordinates


end module chordal_module

