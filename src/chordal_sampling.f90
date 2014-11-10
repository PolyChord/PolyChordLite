module chordal_module
    implicit none

    contains

    function SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)  result(baby_points)
        use priors_module, only: prior
        use settings_module, only: program_settings,phantom_type,live_type
        use random_module, only: random_orthonormal_basis,random_real
        use utils_module, only: logzero

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
        double precision,    dimension(settings%nTotal,settings%num_babies)   :: baby_points

        
        double precision, dimension(settings%nTotal)   :: previous_point


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat
        double precision,    dimension(settings%nDims,settings%num_babies)   :: nhats

        double precision  :: max_chord

        double precision :: step_length

        integer :: i_babies

        double precision,dimension(settings%grades%min_grade:settings%grades%max_grade) :: timers
        double precision :: time0,time1
        integer, dimension(settings%num_babies) :: grade_order

        logical do_timing


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

        if(settings%do_timing) then
            do_timing = random_real()<1d0/settings%print_timing
        else
            do_timing=.false.
        end if

        ! Generate a choice of nhats in the orthogonalised space
        if(settings%do_grades) then
            nhats = generate_nhats(settings,grade_order)
            if(do_timing) timers=0
        else
            nhats = generate_nhats(settings)
        end if

        ! Transform to the unit hypercube
        nhats = matmul(cholesky,nhats)



        do i_babies=1,settings%num_babies
            ! Zero the timer
            if(do_timing) call cpu_time(time0)

            ! Get a new random direction
            nhat = nhats(:,i_babies)

            ! Normalise it
            nhat = nhat/sqrt(dot_product(nhat,nhat))

            ! Generate a new random point along the chord defined by the previous point and nhat
            baby_points(:,i_babies) = slice_sample(loglikelihood,priors, nhat, previous_point, settings)

            if(baby_points(settings%l0,i_babies)<=logzero) then
                baby_points(settings%l0,settings%num_babies)=logzero
                return
            end if

            ! Set this one to be a phantom point
            baby_points(settings%point_type,i_babies) = phantom_type

            if(settings%do_grades) then
                if(grade_order(i_babies)/=settings%grades%min_grade) baby_points(settings%nlike,i_babies) = 0
            end if

            ! keep track of the largest chord
            max_chord = max(max_chord,baby_points(settings%last_chord,i_babies))

            ! Save this for the next loop
            previous_point = baby_points(:,i_babies)

            ! Give the previous point the step length
            previous_point(settings%last_chord) = step_length

            ! Zero the likelihood calls
            previous_point(settings%nlike) = 0

            ! Note the time for this grade
            if(do_timing) then
                call cpu_time(time1)
                if(settings%do_grades) timers(grade_order(i_babies)) = timers(grade_order(i_babies)) + time1-time0
            end if
        end do
        if(do_timing) then 
            write(*, '( <settings%grades%max_grade-settings%grades%min_grade+1>(F8.2, "% "), "(Total time:", F8.2, "seconds)" )') timers/sum(timers)*100, sum(timers)
        end if

        ! Hand back the maximum chord this time to be used as the step length
        ! next time this point is drawn
        baby_points(settings%last_chord,:) = max_chord

        ! Set the last one to be a live type
        baby_points(settings%point_type,settings%num_babies) = live_type

    end function SliceSampling

    function generate_nhats(settings,grade_order) result(nhats)
        use settings_module, only: program_settings
        use random_module, only: random_orthonormal_basis,shuffle_deck
        implicit none
        class(program_settings), intent(in) :: settings

        integer, dimension(settings%num_babies), intent(out),optional :: grade_order

        double precision,    dimension(settings%nDims,settings%num_babies)   :: nhats


        integer :: i_grade
        integer :: grade_nDims
        integer :: grade_index
        integer :: num_repeats
        integer :: i_repeat

        integer :: i_babies

        integer, dimension(settings%num_babies) :: deck

        ! Initialise at 0
        nhats=0d0

        i_babies=1
        if(present(grade_order)) then

            ! Generate a sequence of random bases
            do i_grade=settings%grades%min_grade,settings%grades%max_grade

                grade_nDims = settings%grades%grade_nDims(i_grade)
                grade_index = settings%grades%grade_index(i_grade) 
                num_repeats = settings%grades%num_repeats(i_grade) 

                do i_repeat=1, num_repeats
                    nhats(grade_index:grade_index+grade_nDims-1,i_babies:i_babies+grade_nDims-1) = random_orthonormal_basis(grade_nDims)
                    grade_order(i_babies:i_babies+grade_nDims-1)=i_grade
                    i_babies = i_babies + grade_nDims
                end do

            end do
        else
            ! Generate a sequence of random bases
            do i_repeat=1, settings%num_repeats
                nhats(:,i_babies:i_babies+settings%nDims-1) = random_orthonormal_basis(settings%nDims)
                i_babies = i_babies + settings%nDims
            end do
        end if

        ! Create a shuffled deck
        deck = (/ (i_babies,i_babies=1,settings%num_babies) /)
        call shuffle_deck(deck(2:settings%num_babies))

        ! Shuffle the nhats
        nhats = nhats(:,deck)

        ! Shuffle the grade order
        if(present(grade_order)) grade_order = grade_order(deck)


    end function generate_nhats




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

        double precision :: temp_random

        integer :: i_step

        ! estimate at an appropriate chord
        w = x0(S%last_chord)

        ! record the number of likelihood calls
        R(S%nlike) = x0(S%nlike)
        L(S%nlike) = 0


        ! Select initial start and end points
        temp_random = random_real()
        L(S%h0:S%h1) = x0(S%h0:S%h1) -   temp_random   * w * nhat 
        R(S%h0:S%h1) = x0(S%h0:S%h1) + (1-temp_random) * w * nhat 

        ! Calculate initial likelihoods
        call calculate_point(loglikelihood,priors,R,S)
        call calculate_point(loglikelihood,priors,L,S)

        ! expand R until it's outside the likelihood region
        i_step=0
        do while(R(S%l0) >= x0(S%l1) .and. R(S%l0) > logzero )
            i_step=i_step+1
            R(S%h0:S%h1) = x0(S%h0:S%h1) + nhat * w * i_step
            call calculate_point(loglikelihood,priors,R,S)
            if(i_step>100) write(*,*) ' too many R steps '
        end do

        ! expand L until it's outside the likelihood region
        i_step=0
        do while(L(S%l0) >= x0(S%l1) .and. L(S%l0) > logzero )
            i_step=i_step+1
            L(S%h0:S%h1) = x0(S%h0:S%h1) - nhat * w * i_step
            call calculate_point(loglikelihood,priors,L,S)
            if(i_step>100) write(*,*) ' too many L steps '
        end do

        ! Sample within this bound
        i_step=0
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

            double precision :: x0Rd
            double precision :: x0Ld

            i_step=i_step+1
            if (i_step>100) then
                write(*,*) 'Polychord Warning: Non deterministic loglikelihood'
                x1(S%l0) = logzero
                return
            end if
            
            ! Find the distance between x0 and L 
            x0Ld= distance(x0(S%h0:S%h1),L(S%h0:S%h1))
            ! Find the distance between x0 and R 
            x0Rd= distance(x0(S%h0:S%h1),R(S%h0:S%h1))

            ! Draw a random point within L and R
            x1(S%h0:S%h1) = x0(S%h0:S%h1)+ (random_real() * (x0Rd+x0Ld) - x0Ld) * nhat 

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

