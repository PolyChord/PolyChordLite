
!> Module containing utilities to generate random numbers

module random_module

#ifdef MPI
    use mpi_module
#endif

    implicit none           

    contains


    ! ===========================================================================================


    !> Initialise the random number generators
    !! 
    !! Takes an optional argument in case you want to set the seed. Otherwise it
    !! just uses the system time. 
    !!
    !! @todo Write seed to log file

    subroutine initialise_random(seed_input)
        implicit none
        !> Optional seed input.
        !! If this isn't included, then the system time is used
        integer, optional, intent(in) :: seed_input

        integer,allocatable,dimension(:) :: seed ! vector to be passed to random_seed
#ifdef MPI
        integer :: mpierror
#endif

        integer :: myrank

        integer :: size_seed
        integer :: dt(8)
        integer :: t
        integer :: i



        ! Get the global ranking
#ifdef MPI
        call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierror)
#else
        myrank = 0
#endif

        call random_seed(size = size_seed)
        allocate(seed(size_seed))

        if (present(seed_input)) then
            ! If the seed argument is present, initialise stream with this
            t = seed_input
        else

            ! Seed from the system time (hopefully one of these will work on your machine)
            call system_clock(t)
            if (t == 0) then
                call date_and_time(values=dt)
                t = &
                    (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 &
                    + dt(7) * 1000 &
                    + dt(8)
            end if

        end if

#ifdef MPI
        ! Broadcast the system time to all nodes
        call MPI_BCAST(t,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)      

        ! Augment the seed on each node by adding 1 to it
        t = t+myrank
#endif

        ! set up the seeds for the better generator
        seed = [ ( floor(basic_random(t)*huge(0)) , i=1,size_seed ) ]

        ! Seed the better generator
        call random_seed(put=seed)

    end subroutine initialise_random

    function basic_random(seed)
        implicit none
        real :: basic_random
        integer :: seed
        integer :: oldseed = 0
        integer, parameter :: c1 = 19423
        integer, parameter :: c2 = 811
        save oldseed

        if (oldseed .eq. 0) oldseed = seed
        oldseed = mod(c1 * oldseed, c2)
        basic_random = 1.0 * oldseed / c2

    end function basic_random

    ! ===========================================================================================


    !>  Random unit hypercube coordinate
    !!
    !! Generate a randomly directed unit vector in the unit hypercube

    function random_reals(nDims) result(reals)
        implicit none

        !> Size of coordinate vector
        integer,intent(in) :: nDims 

        ! The output nDims coordinate
        double precision, dimension(nDims) :: reals

        call random_number(reals)


    end function random_reals

    ! ===========================================================================================

    !>  Random array of logicals, produced with probability p
    function random_logicals(nDims,p)
        implicit none

        !> Size of coordinate vector
        integer,intent(in) :: nDims 

        !> probability of true
        double precision,intent(in) :: p 

        ! The output nDims coordinate
        logical, dimension(nDims) :: random_logicals

        random_logicals = random_reals(nDims)<p


    end function random_logicals

    ! ===========================================================================================


    !>  Single random real number

    function random_real()
        implicit none

        ! The output nDims coordinate
        double precision :: random_real
        double precision :: random_real_vec(1)

        random_real_vec = random_reals(1)
        random_real = random_real_vec(1)

    end function random_real

    ! ===========================================================================================


    !>  Bernoulli trials, probability p, or odds ratio p:q for
    !!
    !! http://en.wikipedia.org/wiki/Bernoulli_trial
    !! 
    !! http://en.wikipedia.org/wiki/Odds#Gambling_odds_versus_probabilities
    !!
    !! Returns true with probability p
    !! Returns true with probability p/(p+q)

    function bernoulli_trials(nDims,p,q)
        implicit none

        ! The output nDims coordinate
        integer, intent(in) :: nDims
        double precision,intent(in) :: p
        double precision,intent(in),optional :: q

        logical, dimension(nDims) :: bernoulli_trials

        if(present(q)) then
            bernoulli_trials = random_reals(nDims) < p/(p+q)
        else
            bernoulli_trials = random_reals(nDims) < p
        end if


    end function bernoulli_trials

    ! ===========================================================================================


    !>  Bernoulli trial (see above)

    function bernoulli_trial(p,q)
        implicit none

        ! The output nDims coordinate
        double precision,intent(in) :: p
        double precision,intent(in),optional :: q

        logical :: bernoulli_trial
        logical :: bernoulli_trial_vec(1)

        if(present(q)) then
            bernoulli_trial_vec = bernoulli_trials(1,p,q)
        else
            bernoulli_trial_vec = bernoulli_trials(1,p)
        end if

        bernoulli_trial = bernoulli_trial_vec(1)

    end function bernoulli_trial

    ! ===========================================================================================


    !>  Random integer between 1 and nmax (inclusive)

    function random_integer(nmax)
        implicit none

        !> Maximum integer to generate
        integer,intent(in) :: nmax

        ! The output nDims coordinate
        integer :: random_integer

        random_integer = ceiling(random_real()*nmax)


    end function random_integer

    ! ===========================================================================================

    !>  Random Gaussian vector
    !!
    !! Generate nDims random gausian numbers with mean 0 and variance 1

    function random_gaussian(nDims)
        use utils_module, only: inv_normal_cdf
        implicit none

        !> Size of vector to be generated
        integer,intent(in) :: nDims

        ! The output unit vector
        double precision, dimension(nDims) :: random_gaussian

        random_gaussian = inv_normal_cdf(random_reals(nDims))

    end function random_gaussian

    ! ===========================================================================================


    !>  Random direction vector
    !!
    !! Generate a randomly directed unit vector in nDims dimensional space
    !!
    !! This is done by generating nDims gaussian random numbers with mean 0 and
    !! variance 1. These are spherically symetrically distributed, so
    !! normalising this give a random direction.

    function random_direction(nDims)
        implicit none

        !> Size of vector to be generated
        integer,intent(in) :: nDims

        ! The output unit vector
        double precision, dimension(nDims) :: random_direction
        double precision :: random_direction2

        random_direction2=0
        do while(random_direction2<=0)

            ! Generate nDims gaussian random numbers
            random_direction = random_gaussian(nDims)
            ! Calculate the modulus squared
            random_direction2 = dot_product(random_direction,random_direction)
        end do

        ! normalise the vector
        random_direction = random_direction / sqrt(random_direction2)

    end function random_direction

    ! ===========================================================================================


    !>  Random subspace vector
    !!
    !! Generate a randomly directed unit vector in the nDims-1 dimensional subspace defined by rhat
    !!
    !! This is done by generating nDims gaussian random numbers with mean 0 and
    !! variance 1. These are spherically symetrically distributed, so
    !! normalising this give a random direction.

    function random_subdirection(nDims,rhat)
        implicit none

        !> Size of vector to be generated
        integer,intent(in) :: nDims
        double precision, dimension(nDims),intent(in) :: rhat

        ! The output unit vector
        double precision, dimension(nDims) :: random_subdirection
        double precision :: random_subdirection2

        random_subdirection2=0
        do while (random_subdirection2<=0)

            ! Generate nDims gaussian random numbers
            random_subdirection = random_gaussian(nDims)

            ! Project into the rhat space
            random_subdirection = random_subdirection - rhat * dot_product(rhat,random_subdirection) / dot_product(rhat,rhat)

            ! Find the modulus squared
            random_subdirection2 = dot_product(random_subdirection,random_subdirection)

        end do
        ! normalise the vector
        random_subdirection = random_subdirection / sqrt(random_subdirection2)

    end function random_subdirection


    ! ===========================================================================================


    !>  Random point in unit sphere
    !!
    !! Generate a uniformly distributed point in the unit n-sphere
    !!
    !! This is done by generating a random direction, and then multiplying it by
    !! a radius distributed as \f$ P(r) ~ r^{d-1}\f$. where \f$d\f$ is the
    !! dimension of the input vector

    function random_point_in_sphere(nDims)
        implicit none

        !> Dimension of sphere
        integer,intent(in) :: nDims 

        ! The output vector
        double precision, dimension(nDims) :: random_point_in_sphere

        ! Temporary variable for storing a random radius
        double precision :: rand_rad


        ! generate a random direction
        random_point_in_sphere = random_direction(nDims)

        ! generate a single random number distributed as ~ r^{d-1}
        rand_rad = random_real()**(1d0/nDims)

        ! Create a point in the sphere
        random_point_in_sphere = random_point_in_sphere * rand_rad


    end function random_point_in_sphere

    !> Construct a randomly oriented ndimensional orthogonormal basis by using 
    !! [Gram-Schmidt Orthogonalisation](http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process). 
    !!
    !! Outputs a matrix where the ith basis vector is stored in basis(:,i)
    function random_orthonormal_basis(nDims) result(basis)
        !> Dimensionality of the basis
        integer, intent(in) ::  nDims

        ! Set of vectors, the ith vector is in basis(:,i)
        double precision, dimension(nDims,nDims) :: basis

        ! Iterators
        integer :: i,j

        do i=1,nDims
            ! Generate a randomly directed vector 
            basis(:,i) = random_direction(nDims)
            ! Othogonalise it with respect to the other vectors using Gram
            ! Schmidt orthogonalisation
            do j= 1,i-1
                basis(:,i) = basis(:,i)- dot_product(basis(:,i),basis(:,j)) * basis(:,j)
            end do
            ! normalise the vector
            basis(:,i) = basis(:,i)/sqrt(dot_product(basis(:,i),basis(:,i)))
        end do

    end function random_orthonormal_basis


    !> Construct a sequence of nhats composed of several random orthonormal bases
    !!
    !! Outputs a matrix where the ith basis vector is stored in basis(:,i)
    function random_orthonormal_bases(nDims,num_nhats) result(nhats)
        !> Dimensionality of the basis
        integer, intent(in) ::  nDims
        integer, intent(in) ::  num_nhats

        ! Set of vectors, the ith vector is in basis(:,i)
        double precision, dimension(nDims,num_nhats) :: nhats

        double precision, dimension(nDims,nDims) :: basis

        integer :: lower_index,upper_index

        ! Fill up the first whole sections with nDims
        lower_index = 1
        upper_index = nDims
        do while(upper_index<num_nhats)
            nhats(:,lower_index:upper_index) = random_orthonormal_basis(nDims)
            lower_index = lower_index+nDims
            upper_index = upper_index+nDims
        end do

        basis = random_orthonormal_basis(nDims)

        nhats(:,lower_index:num_nhats) = basis(:,1:1+num_nhats-lower_index)
        



    end function random_orthonormal_bases

    ! ===========================================================================================
    !> Generate a set of k distinct random integers between [1,m]
    !! There are two ways of doing this:
    !! # Generate successive random numbers, rejecting any which have already
    !!   been found
    !! # Generate a shuffled deck of random numbers, and pick the first k of
    !!   them
    !!
    !! Which of these is the more efficient depends on the size difference of
    !! k and m. In general, the first method will require one to draw an
    !! expected number of random numbers evaluating to:
    !! \f[ m\sum\limits_{i=1}^k \frac{1}{m-(i-1)} \approx m \log\left(\frac{m}{m-k}\right) \f]
    !! whereas the second requires m random numbers to be drawn. The cutoff
    !! moment is therefore when
    !! \f[ k \sim \frac{e-1}{e} n \sim 0.63 n \f]
    !!
    !! We thus use the first method if k is less than 0.63n
    function random_distinct_integers(m,k) result(integers)
        implicit none
        double precision, parameter :: cutoff=0.6321205588285576784d0

        !> The upper bound of integers to be generated
        integer, intent(in) :: m

        !> The number of distinct integers to be generated
        !! If k is not present, then we generate a shuffled deck of integers
        integer,intent(in) :: k

        ! The returning array
        integer, dimension(k) :: integers

        integer :: i
        integer, dimension(m) :: deck 

        if(k<cutoff*m) then
            ! Case 1, generate by rejection

            i=1
            do while (i<=k)
                ! Draw a random integer
                integers(i) = random_integer(m) 

                ! Move to the next index if it's unique
                if(all(integers(i)/=integers(:i-1))) i=i+1
            end do

        else
            ! Case 2, generate by shuffling

            ! Generate a randomly shuffled deck
            deck = [(i,i=1,m)]
            call shuffle_deck(deck)

            ! output the first k of them
            integers = deck

        end if

    end function random_distinct_integers


    !> Shuffle a 'deck' of integers
    !!
    !! This uses a [Fischer-Yates](http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle)
    !! shuffle, which is \f$\mathcal{O}(n)\f$

    subroutine shuffle_deck(deck)
        implicit none

        ! The input deck to be shuffled
        integer, intent(inout), dimension(:) :: deck

        ! Temporary variable for swapping
        integer :: temp

        ! Indices
        integer :: i,j

        ! Size of deck
        integer :: n

        n = size(deck)

        do i=n,1,-1
            ! pick a random integer in [1,i]
            j=random_integer(i)

            ! swap elements i and j
            temp    = deck(i) 
            deck(i) = deck(j)
            deck(j) = temp
        end do

    end subroutine shuffle_deck


    !> This function generates a random integer with variable probability
    !!
    !! The probability is defined by the input array, and need not be
    !! normalised.
    !!
    !! E.g: If one reieves an array (/ 50.0, 120.0, 30.0 /), this function will
    !! return:
    !! |---|-----------------|
    !! | 1 | 25% of the time |
    !! | 2 | 60% of the time |
    !! | 3 | 15% of the time |
    !! |---|-----------------|
    !!
    function random_integer_P(probabilities)
        implicit none
        double precision, dimension(:),intent(in) :: probabilities

        integer :: random_integer_P

        double precision :: norm
        double precision :: cdf
        double precision :: rand


        ! Calculate normalisation constant
        norm = sum(probabilities)

        ! Generate a random real number
        rand = random_real()

        ! Initialise the cumulative probability at 0
        cdf=0d0

        ! iterate through the array, returning when the cdf is bigger than the
        ! random real.
        do random_integer_P=1,size(probabilities)
            cdf = cdf + probabilities(random_integer_P)/norm
            if(rand<cdf) return
        end do


    end function random_integer_P


    !> Generate a random inverse covariance matrix with unit determinant
    !! 
    subroutine random_inverse_covmat(invcovmat,logdetcovmat,sigma,nDims)
        implicit none
        integer,          intent(in)                         :: nDims
        double precision, intent(out),dimension(nDims,nDims) :: invcovmat
        double precision, intent(out)                        :: logdetcovmat
        double precision, intent(in)                         :: sigma

        double precision, dimension(nDims)       :: eigenvalues
        double precision, dimension(nDims,nDims) :: eigenvectors
        integer :: j
        double precision, parameter :: rng=1d-2

        ! Generate a random basis for the eigenvectors
        eigenvectors = random_orthonormal_basis(nDims)
        ! Generate the eigenvalues logarithmically in [rng,1] * sigma
        do j=1,nDims
            eigenvalues(j)  = sigma * rng**((j-1d0)/(nDims-1d0))
        end do
        ! Sort them lowest to highest for consistency
        !call dlasrt('D',nDims,eigenvalues,info)

        ! Create the inverse covariance matrix in the eigenbasis
        invcovmat = 0d0
        do j=1,nDims
            invcovmat(j,j) = 1d0/eigenvalues(j)**2
        end do

        ! Rotate the matrix into the coordinate basis
        invcovmat = matmul(eigenvectors,matmul(invcovmat,transpose(eigenvectors)))

        ! sum up the logs of the eigenvalues to get the log of the determinant
        logdetcovmat = 2 * sum(log(eigenvalues))

    end subroutine random_inverse_covmat

end module random_module
