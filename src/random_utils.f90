! Include the Intel Math Kernel Library - Vector Statistical Library
! This is a library used for generating random numbers optimally on intel cores 
include 'mkl_vsl.f90'  

!> Module containing utilities to generate random numbers
!!
!! For details on the Intel Math Kernel Library - Vector Statistical Library (MKL_VSL)
!! consult the [excellent manual](
!! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/)
!!
!! @todo Save stream state for exact replication if necessary

module random_module

    use mkl_vsl_type
    use mkl_vsl

    use mpi

    implicit none           


    !> The stream state for vector random number generators
    type (vsl_stream_state) :: rng_stream


    contains


    ! ===========================================================================================
    

    !> Initialise the random number generators
    !! 
    !! Takes an optional argument in case you want to set the seed. Otherwise it
    !! just uses the system time. 
    !!
    !! The intel command vslnewstream initialises the module variable
    !! rng_stream, which is then used throughout the code to generate random
    !! numbers from. More can be found on streams [here](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#7_3_2_Creating_and_Initializing_Random_Streams.htm)
    !! We have opted to use a fast implementation of the Mersenne Twister as our
    !! [basic random number generator.](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#7_2_Basic_Generators.htm)
    !!
    !! @todo Error code processing
    !!
    !! @todo Write seed to log file
    !!
    !! @todo parallel stream initialisation : check that they're independent

    subroutine initialise_random(seed_input)
        implicit none
        !> Optional seed input.
        !! If this isn't included, then the system time is used
        integer, optional, intent(in) :: seed_input

        integer :: errcode ! Error code
        integer :: seed    ! seed to be generated from system time
        integer :: mpierror

        ! The choice of random number generator
        ! In this case it is a fast mersenne twister
        integer,parameter       :: basic_rng_type = VSL_BRNG_SFMT19937

        integer :: myrank

        ! Get the global ranking
        call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierror)

        if (present(seed_input)) then
            ! If the seed argument is present, initialise stream with this
            errcode=vslnewstream( rng_stream, basic_rng_type,  seed_input + myrank)
        else
            ! Otherwise initialise it with the system time
            call system_clock(seed)
            ! Broadcast the same seed to everybody from 'root' (0)
            call MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierror)      

            errcode=vslnewstream( rng_stream, basic_rng_type,  seed + myrank )
        end if

    end subroutine initialise_random

    ! ===========================================================================================


    !>  Random unit hypercube coordinate
    !!
    !! Generate a randomly directed unit vector in the unit hypercube
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm) 
    !! intel method to generate uniform random numbers using the 
    !! [vdrnguniform](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-D7AD317E-34EC-4789-8027-01D0E194FAC1.htm)
    !!
    !! Note that [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm)
    !! can come in an accurate form.

    function random_reals(nDims)
        implicit none

        !> Size of coordinate vector
        integer,intent(in) :: nDims 

        ! The output nDims coordinate
        double precision, dimension(nDims) :: random_reals

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        double precision, parameter :: l_bound  = 0.0 ! generate random numbers between 0 and 1
        double precision, parameter :: u_bound  = 1.0 ! 

        integer :: errcode ! Error code

        ! Generate nDims random numbers, stored in the output vector 'random_coordinate'
        ! (v=vector,d=double,rng,uniform)
        errcode=vdrnguniform( method, rng_stream, nDims, random_reals, l_bound, u_bound)


    end function random_reals

    ! ===========================================================================================


    !>  Single random real number
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm) 
    !! intel method to generate uniform random numbers using the 
    !! [vdrnguniform](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-D7AD317E-34EC-4789-8027-01D0E194FAC1.htm)
    !!
    !! Note that [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm)
    !! can come in an accurate form.

    function random_real()
        implicit none

        ! The output nDims coordinate
        double precision :: random_real

        double precision :: random_reals(1)

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        double precision, parameter :: l_bound  = 0.0 ! generate random numbers between 0 and 1
        double precision, parameter :: u_bound  = 1.0 ! 

        integer :: errcode ! Error code

        ! Generate nDims random numbers, stored in the output vector 'random_coordinate'
        ! (v=vector,d=double,rng,uniform)
        errcode=vdrnguniform( method, rng_stream, 1, random_reals, l_bound, u_bound)

        random_real = random_reals(1)


    end function random_real

    ! ===========================================================================================


    !>  Random integer between 1 and nmax (inclusive)
    !!
    !! Generate a random integer  between 1 and nmax inclusive
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm) 
    !! intel method to generate uniform random integers using the 
    !! [virnguniform](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-46FE7EF1-CE00-42CB-891D-17FD96062255.htm)
    !!
    !! Note that [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm)
    !! can come in an accurate form.

    function random_integer(nmax)
        implicit none

        !> Maximum integer to generate
        integer,intent(in) :: nmax

        ! The output nDims coordinate
        integer :: random_integer
        integer :: random_integers(1)

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        integer, parameter :: l_bound  = 1 ! generate random numbers between 0 and 1

        integer :: errcode ! Error code

        ! Generate nDims random numbers, stored in the output vector 'random_coordinate'
        ! (v=vector,i=integer,rng,uniform)
        ! This generates a vector of random integers between [1,nmax+1) = [1,nmax]
        errcode=virnguniform( method, rng_stream, 1, random_integers, l_bound, nmax+1)

        random_integer = random_integers(1)


    end function random_integer

    ! ===========================================================================================

    !>  Random Gaussian vector
    !!
    !! Generate nDims random gausian numbers with mean 0 and variance 1
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_GAUSSIAN_ICDF](https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_4_Gaussian_VSL_RNG_METHOD_GAUSSIAN_ICDF.htm) 
    !! intel method to generate Gaussian random numbers using the
    !! [vdrnggaussian](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-63196F25-5013-4038-8BCD-2613C4EF3DE4.htm) function.

    function random_gaussian(nDims)
        implicit none

        !> Size of vector to be generated
        integer,intent(in) :: nDims

        ! The output unit vector
        double precision, dimension(nDims) :: random_gaussian

        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_GAUSSIAN_ICDF

        double precision, parameter :: sigma = 1.0 ! Standard deviation of gaussian random numbers
        double precision, parameter :: mean  = 0.0 ! Mean of gaussian random numbers

        integer :: errcode ! Error code


        ! Generate nDims gaussian random numbers, stored in a vector
        ! random_gaussian
        ! with mean=0 variance=1
        ! (v=vector,d=double,rng,gaussian)
        errcode=vdrnggaussian( method, rng_stream, nDims, random_gaussian, mean, sigma)

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
        do while(random_direction2==0)

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
        do while (random_subdirection2==0)

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


    !>  Random skewed direction vector
    !!
    !! Generate a randomly directed normalised vector in nDims dimensional space
    !! according to the covariance matrix recieved. The vector is normalised
    !! with respect to the [Mahalanobis_distance](http://en.wikipedia.org/wiki/Mahalanobis_distance)
    !! defined by the variance-covariance matrix.
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_GAUSSIANMV_ICDF](https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_7_GaussianMV_VSL_RNG_METHOD_GAUSSIANMV_ICDF.htm#9_3_7_GaussianMV_VSL_RNG) 
    !! intel method to generate multivariate gaussian random numbers using the 
    !! [vdrnggaussianmv](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-1595CFFA-4878-4BCD-9C1D-2034731C1F4F.htm#GUID-1595CFFA-4878-4BCD-9C1D-2034731C1F4F) function.

    function random_skewed_direction(nDims,cholesky) result(nhat)
        use utils_module, only: distance
        implicit none

        !> Size of vector to be generated
        integer,intent(in) :: nDims

        !> Cholesky decomposition of the covariance matrix
        double precision, dimension(nDims,nDims) :: cholesky

        ! The output normalised vector
        double precision, dimension(nDims) :: nhat

        ! Temporary length 1 vector of nDims vectors for passing to the routine
        double precision, dimension(nDims,1) :: nhat_temp

        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_GAUSSIANMV_ICDF
        ! The type of matrix storage (this indicates 'full' storage, i.e. the
        ! easy to understand kind)
        integer,parameter       :: mstorage=VSL_MATRIX_STORAGE_FULL

        ! Mean of the vector produced
        double precision, dimension(nDims) :: mean

        ! Error code for output of intel routine
        integer :: errcode

        ! The modulus squared of the vector
        double precision ::nhat2

        mean =0
        nhat2=0d0
        do while(nhat2==0d0)
            ! Compute a random gaussian vector subject to the matrix
            errcode=vdrnggaussianmv( method, rng_stream, 1,nhat_temp, nDims, mstorage, mean, cholesky)

            nhat2 = dot_product(nhat_temp(:,1),nhat_temp(:,1))
        end do

        nhat = nhat_temp(:,1)/sqrt(nhat2)
        

    end function random_skewed_direction


    ! ===========================================================================================

    ! ===========================================================================================



    !> Random point in the [beta distribution](http://en.wikipedia.org/wiki/Beta_distribution)
    !!
    !! Generate an nDims-vector of points with each coordinate \f$x_i\f$ distributed as: 
    !! \f$ P(x_i;p,q) ~ x_i^{p-1}(1-x_i)^{q-1}\f$
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_BETA_CJA](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_17_Beta_VSL_RNG_METHOD_BETA_CJA.htm#9_3_16_Beta_VSL_RNG) 
    !! intel method to generate points. Note that this can come in an accurate form                               

    function random_beta_vec(p,q,nDims)
        implicit none

        !> beta shape parameter 1
        double precision, intent(in) :: p
        !> beta shape parameter 2
        double precision, intent(in) :: q
        !> Size of coordinate vector
        integer :: nDims

        ! The output vector
        double precision, dimension(nDims) :: random_beta_vec


        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_BETA_CJA

        double precision, parameter :: offset = 0.0 ! Mean of the beta distribution
        double precision, parameter :: scale_factor  = 1.0 ! width of the beta distribution

        integer :: errcode ! Error code



        ! call the intel procedure to generate the point
        errcode=vdrngbeta( method, rng_stream, nDims, random_beta_vec, p, q, offset, scale_factor )


    end function random_beta_vec

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

        ! generate a single random number distributed as beta(r,1) ~ r^{d-1}
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


    !> De-initialise the random number generators

    subroutine deinitialise_random

        implicit none
        integer :: errcode ! Error code

        errcode=vsldeletestream( rng_stream )

    end subroutine deinitialise_random


end module random_module
