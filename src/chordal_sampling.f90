module chordal_module
    implicit none

    contains

    function SliceSampling(loglikelihood,priors,settings,logL,seed_point,cholesky,nlike,num_repeats)  result(baby_points)
        use priors_module, only: prior
        use settings_module, only: program_settings
        use random_module, only: random_orthonormal_basis,random_real

        implicit none
        interface
            function loglikelihood(theta,phi)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                double precision :: loglikelihood
            end function
        end interface

        ! ------- Inputs -------
        !> The prior information
        type(prior), dimension(:), intent(in) :: priors

        !> program settings (mostly useful to pass on the number of live points)
        type(program_settings), intent(in) :: settings

        !> The seed point
        double precision, intent(in), dimension(settings%nTotal)   :: seed_point

        !> The loglikelihood contour to sample on
        double precision, intent(in)   :: logL

        !> The directions of the chords
        double precision, intent(in), dimension(settings%nDims,settings%nDims) :: cholesky

        !> The number of babies within each grade to produce
        integer, dimension(:), intent(in) :: num_repeats

        ! ------- Outputs -------
        !> The number of likelihood evaluations
        integer, dimension(:), intent(out) :: nlike

        !> The newly generated point, plus the loglikelihood bound that
        !! generated it
        double precision,    dimension(settings%nTotal,sum(num_repeats))   :: baby_points


        ! ------- Local Variables -------
        double precision,    dimension(settings%nDims)   :: nhat
        double precision,    dimension(settings%nDims,sum(num_repeats))   :: nhats
        integer,   dimension(sum(num_repeats))   :: speeds ! The speed of each nhat

        double precision, dimension(settings%nTotal)   :: previous_point

        integer :: i_babies

        double precision :: w

        ! Start the baby point at the seed point
        previous_point = seed_point

        ! Initialise the likelihood counter at 0
        nlike = 0

        ! Generate a choice of nhats in the orthogonalised space
        nhats = generate_nhats(settings,num_repeats,speeds)

        ! Transform to the unit hypercube
        nhats = matmul(cholesky,nhats)

        do i_babies=1,size(nhats,2)
            ! Get a new random direction
            nhat = nhats(:,i_babies)

            ! Normalise it
            w = sqrt(dot_product(nhat,nhat))
            nhat = nhat/w
            w = w * 3d0 !* exp( lgamma(0.5d0 * settings%nDims) - lgamma(1.5d0 + 0.5d0 * settings%nDims) ) * settings%nDims

            ! Generate a new random point along the chord defined by the previous point and nhat
            baby_points(:,i_babies) = slice_sample(loglikelihood,priors,logL, nhat, previous_point, w, settings,nlike(speeds(i_babies)))

            ! Save this for the next loop
            previous_point = baby_points(:,i_babies)

        end do

    end function SliceSampling

    function generate_nhats(settings,num_repeats,speeds) result(nhats)
        use settings_module, only: program_settings
        use random_module, only: random_orthonormal_bases,shuffle_deck
        implicit none
        type(program_settings), intent(in) :: settings

        integer, dimension(:), intent(in) :: num_repeats !> The number of babies within each grade to produce
        integer,  intent(out), dimension(sum(num_repeats))   :: speeds !> The speed of each nhat
        double precision,    dimension(settings%nDims,sum(num_repeats))   :: nhats


        integer :: i_grade
        integer :: grade_nDims
        integer :: grade_index

        integer :: i_babies

        integer, dimension(size(nhats,2)) :: deck

        ! Initialise at 0
        nhats=0d0

        i_babies=1

        ! Generate a sequence of random bases
        do i_grade=1,size(num_repeats)

            grade_nDims   = sum(settings%grade_dims(i_grade:))
            grade_index   = sum(settings%grade_dims(:i_grade-1))+1

            nhats(grade_index:,i_babies:i_babies+num_repeats(i_grade)-1) = random_orthonormal_bases(grade_nDims,num_repeats(i_grade))

            speeds(i_babies:i_babies+num_repeats(i_grade)-1) = i_grade

            i_babies = i_babies + num_repeats(i_grade)

        end do

        ! Create a shuffled deck
        ! N.B. we want the first one to be a slow likelihood evaluation, as this
        ! will be the first one we're doing this run
        deck = [ (i_babies,i_babies=1,size(nhats,2)) ]
        call shuffle_deck(deck(2:))

        ! Shuffle the nhats
        nhats = nhats(:,deck)

        ! Shuffle the speeds
        speeds = speeds(deck)


    end function generate_nhats






    !> Slice sample at slice logL, along the direction defined by nhat, starting at x0.
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
    function slice_sample(loglikelihood,priors,logL,nhat,x0,w,S,n) result(baby_point)
        use settings_module, only: program_settings
        use priors_module, only: prior
        use utils_module,  only: logzero, distance
        use random_module, only: random_real
        use calculate_module, only: calculate_point
        implicit none
        interface
            function loglikelihood(theta,phi)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                double precision :: loglikelihood
            end function
        end interface

        !> The prior information
        type(prior), dimension(:), intent(in) :: priors
        !> program settings
        type(program_settings), intent(in) :: S
        !> The Loglikelihood bound
        double precision, intent(in)                           :: logL
        !> The direction to search for the root in
        double precision, intent(in),    dimension(S%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(S%nTotal)   :: x0
        !> The initial width
        double precision, intent(in) :: w
        !> The number of likelihood calls
        integer, intent(inout) :: n

        ! The output finish point
        double precision,    dimension(S%nTotal)   :: baby_point

        ! The upper bound
        double precision,    dimension(S%nTotal)   :: R
        ! The lower bound
        double precision,    dimension(S%nTotal)   :: L


        double precision :: temp_random

        integer :: i_step

        ! Select initial start and end points
        temp_random = random_real()
        L(S%h0:S%h1) = x0(S%h0:S%h1) -   temp_random   * w * nhat 
        R(S%h0:S%h1) = x0(S%h0:S%h1) + (1-temp_random) * w * nhat 

        ! Calculate initial likelihoods
        call calculate_point(loglikelihood,priors,R,S,n)
        call calculate_point(loglikelihood,priors,L,S,n)

        ! expand R until it's outside the likelihood region
        i_step=0
        do while( R(S%l0) >= logL .and. R(S%l0) > logzero )
            i_step=i_step+1
            R(S%h0:S%h1) = x0(S%h0:S%h1) + nhat * w * i_step
            call calculate_point(loglikelihood,priors,R,S,n)
            n=n+1
        end do
        if(i_step>100) write(*,'(" too many R steps (",I10,")")') i_step

        ! expand L until it's outside the likelihood region
        i_step=0
        do while(L(S%l0) >= logL      .and. L(S%l0) > logzero )
            i_step=i_step+1
            L(S%h0:S%h1) = x0(S%h0:S%h1) - nhat * w * i_step
            call calculate_point(loglikelihood,priors,L,S,n)
        end do
        if(i_step>100) write(*,'(" too many L steps (",I10,")")') i_step

        ! Sample within this bound
        i_step=0
        baby_point = find_positive_within(L,R)

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
                write(*,'("Polychord Warning: Non deterministic loglikelihood")')
                x1(S%l0) = logzero
                return
            end if

            ! Find the distance between x0 and L 
            x0Ld= distance(x0(S%h0:S%h1),L(S%h0:S%h1))
            ! Find the distance between x0 and R 
            x0Rd= distance(x0(S%h0:S%h1),R(S%h0:S%h1))

            ! Draw a random point within L and R
            x1(S%h0:S%h1) = x0(S%h0:S%h1)+ (random_real() * (x0Rd+x0Ld) - x0Ld) * nhat 

            ! calculate the likelihood 
            call calculate_point(loglikelihood,priors,x1,S,n)

            ! If we're not within the likelihood bound then we need to sample further
            if( x1(S%l0) < logL .or. x1(S%l0) <= logzero ) then

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

