module chordal_module
    implicit none

    contains

    function SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)  result(baby_points)
        use priors_module, only: prior
        use settings_module, only: program_settings,phantom_type,live_type
        use random_module, only: random_orthonormal_basis
        use utils_module, only: distance

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
        double precision, intent(in), dimension(settings%nTotal)   :: seed_point

        !> The directions of the chords
        double precision, intent(in), dimension(settings%nDims,settings%nDims) :: cholesky

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(settings%nTotal,settings%nDims*settings%num_repeats)   :: baby_points

        
        double precision, dimension(settings%nTotal)   :: previous_point

        double precision, dimension(settings%nDims,settings%nDims) :: basis


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat

        double precision  :: max_chord

        double precision :: step_length
        double precision :: modulus

        integer :: i_chords
        integer :: i_repeat


        ! Start the baby point at the seed point
        previous_point = seed_point

        ! Set the number of likelihood evaluations to zero
        previous_point(settings%nlike) = 0

        ! Record the step length
        step_length = seed_point(settings%last_chord)

        ! Give the start point the step length
        previous_point(settings%last_chord) = step_length

        ! Initialise max_chord at 0
        max_chord = 0


        do i_repeat=1,settings%num_repeats
            ! Generate a random orthonormal set of vectors
            basis = random_orthonormal_basis(settings%nDims)

            do i_chords=(i_repeat-1)*settings%nDims+1,(i_repeat-1)*settings%nDims+settings%nDims
                ! Get a new random direction
                nhat =basis(:,1+mod(i_chords-1,settings%nDims) ) 
                nhat = matmul(cholesky,nhat)
                modulus = sqrt(dot_product(nhat,nhat))
                nhat = nhat/modulus

                ! Generate a new random point along the chord defined by the previous point and nhat
                baby_points(:,i_chords) = slice_sample(loglikelihood,priors, nhat, previous_point, settings)

                ! Set this one to be a phantom point
                baby_points(settings%point_type,i_chords) = phantom_type

                ! keep track of the largest chord
                max_chord = max(max_chord,baby_points(settings%last_chord,i_chords))

                ! Save this for the next loop
                previous_point = baby_points(:,i_chords)

                ! Give the previous point the step length
                previous_point(settings%last_chord) = step_length

                ! Zero the likelihood calls
                previous_point(settings%nlike) = 0
            end do

        end do

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_points(settings%last_chord,:) = max_chord

        ! Set the last one to be a live type
        baby_points(settings%point_type,size(baby_points,2)) = live_type

    end function SliceSampling


    function GradedSliceSampling(loglikelihood,priors,settings,choleskys,seed_point)  result(baby_points)
        use priors_module, only: prior
        use settings_module, only: program_settings,phantom_type,live_type
        use random_module, only: random_skewed_direction
        use utils_module, only: distance

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
        double precision, intent(in), dimension(settings%nTotal)   :: seed_point

        !> The directions of the chords
        double precision, intent(in), dimension(settings%nDims,settings%nDims,settings%grades%min_grade:settings%grades%max_grade) :: choleskys

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(settings%nTotal,settings%nDims*settings%num_repeats)   :: baby_points

        
        double precision, dimension(settings%nTotal)   :: previous_point


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat

        double precision  :: max_chord

        double precision :: step_length

        integer :: i_chords


        baby_points=0d0
        ! Start the baby point at the seed point
        previous_point = seed_point

        ! Set the number of likelihood evaluations to zero
        previous_point(settings%nlike) = 0

        ! Record the step length
        step_length = seed_point(settings%last_chord)

        ! Give the start point the step length
        previous_point(settings%last_chord) = step_length

        ! Initialise max_chord at 0
        max_chord = 0
        i_chords=0

        call run_sub_chain(settings%grades%min_grade)

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_points(settings%last_chord,:) = max_chord

        ! Set the last one to be a live type
        baby_points(settings%point_type,size(baby_points,2)) = live_type


        contains
        recursive subroutine run_sub_chain(grade)
            implicit none
            integer, intent(in) :: grade

            integer :: i_step

            ! Skip over unused parameters
            if(settings%grades%chain_lengths(grade)==0 .and. grade<settings%grades%max_grade) call run_sub_chain(grade+1)

            do i_step=1,settings%grades%chain_lengths(grade)

                ! Increment i_chords
                i_chords = i_chords+1

                ! Take a step
                ! Get a new random direction
                nhat = random_skewed_direction(settings%nDims,choleskys(:,:,grade))

                ! Generate a new random point along the chord defined by the previous point and nhat
                baby_points(:,i_chords) = slice_sample(loglikelihood,priors, nhat, previous_point, settings)

                ! Set this one to be a phantom point
                baby_points(settings%point_type,i_chords) = phantom_type

                ! Zero the likelihood calls if its not the slowest likelihood
                if(grade/=settings%grades%min_grade) baby_points(settings%nlike,i_chords) = 0

                ! keep track of the largest chord
                max_chord = max(max_chord,baby_points(settings%last_chord,i_chords))

                ! Save this for the next loop
                previous_point = baby_points(:,i_chords)

                ! Give the previous point the step length
                previous_point(settings%last_chord) = step_length

                ! Zero the likelihood calls
                previous_point(settings%nlike) = 0

                ! If we're not at the maximum grade, then recurse down one more
                if(grade<settings%grades%max_grade) then
                    call run_sub_chain(grade+1)
                end if

            end do

        end subroutine run_sub_chain



    end function GradedSliceSampling





    function AdaptiveParallelSliceSampling(loglikelihood,priors,settings,live_points,seed_point)  result(baby_points)
        use priors_module, only: prior
        use settings_module, only: program_settings,phantom_type,live_type
        use random_module, only: random_distinct_integers
        use utils_module, only: distance

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
        double precision, intent(in), dimension(settings%nTotal)   :: seed_point

        !> The directions of the chords
        double precision, intent(in), dimension(:,:) :: live_points

        ! ------- Outputs -------
        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision, dimension(settings%nTotal,settings%nDims*settings%num_repeats)   :: baby_points

        
        double precision, dimension(settings%nTotal)   :: previous_point


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat

        double precision  :: max_chord

        double precision :: step_length

        integer :: i_chords

        integer :: random_pair(2)


        ! Start the baby point at the seed point
        previous_point = seed_point

        ! Set the number of likelihood evaluations to zero
        previous_point(settings%nlike) = 0

        ! Record the step length
        step_length = seed_point(settings%last_chord)

        ! Give the start point the step length
        previous_point(settings%last_chord) = step_length

        ! Initialise max_chord at 0
        max_chord = 0

        do i_chords=1,settings%nDims*settings%num_repeats

            random_pair=random_distinct_integers(size(live_points,2),2)

            ! Get a new random direction
            nhat = live_points(settings%h0:settings%h1,random_pair(1))-live_points(settings%h0:settings%h1,random_pair(2))
            nhat = nhat/sqrt(dot_product(nhat,nhat))

            ! Generate a new random point along the chord defined by the previous point and nhat
            baby_points(:,i_chords) = slice_sample(loglikelihood,priors, nhat, previous_point, settings)

            ! Set this one to be a phantom point
            baby_points(settings%point_type,i_chords) = phantom_type

            ! keep track of the largest chord
            max_chord = max(max_chord,baby_points(settings%last_chord,i_chords))

            ! Save this for the next loop
            previous_point = baby_points(:,i_chords)

            ! Give the previous point the step length
            previous_point(settings%last_chord) = step_length

            ! Zero the likelihood calls
            previous_point(settings%nlike) = 0

        end do

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_points(settings%last_chord,:) = max_chord

        ! Set the last one to be a live type
        baby_points(settings%point_type,size(baby_points,2)) = live_type

    end function AdaptiveParallelSliceSampling





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


end module chordal_module

