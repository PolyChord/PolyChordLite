include 'mkl_vsl.f90'  ! Include the Intel Math Kernel Library - Vector Statistical Library
                       ! This is a library used for generating random numbers optimally on intel cores

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
    !! @todo parallel stream initialisation

    subroutine initialise_random(seed_input)
        implicit none
        !> Optional seed input.
        !! If this isn't included, then the system time is used
        integer, optional, intent(in) :: seed_input

        integer :: errcode ! Error code
        integer :: seed    ! seed to be generated from system time

        ! The choice of random number generator
        ! In this case it is a fast mersenne twister
        integer,parameter       :: basic_rng_type = VSL_BRNG_SFMT19937


        if (present(seed_input)) then
            ! If the seed argument is present, initialise stream with this
            errcode=vslnewstream( rng_stream, basic_rng_type,  seed_input )
        else
            ! Otherwise initialise it with the system time
            call system_clock(seed)
            errcode=vslnewstream( rng_stream, basic_rng_type,  seed )
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
    !! intel method to generate uniform random numbers. Note that this can come in an accurate form

    function random_hypercube_point(nDims)
        implicit none

        !> Size of coordinate vector
        integer :: nDims 

        ! The output nDims coordinate
        double precision, dimension(nDims) :: random_hypercube_point

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        double precision, parameter :: u_bound  = 0.0 ! generate random numbers between 0 and 1
        double precision, parameter :: l_bound  = 1.0 ! 

        integer :: errcode ! Error code

        ! Generate nDims random numbers, stored in the output vector 'random_coordinate'
        ! (v=vector,d=double,rng,uniform)
        errcode=vdrnguniform( method, rng_stream, nDims, random_hypercube_point, u_bound, l_bound)


    end function random_hypercube_point

    ! ===========================================================================================


    !>  Random real number
    !!
    !! Generate a single random real number in the range [0,1]
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm) 
    !! intel method to generate uniform random numbers. Note that this can come in an accurate form

    function random_real()
        implicit none

        !> Size of coordinate vector
        integer :: nDims 

        ! output random real number
        double precision :: random_real

        ! temporary vector to pass to vdrnguniform
        double precision, dimension(1) :: random_real_vec

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        double precision, parameter :: u_bound  = 0.0 ! generate random numbers between 0 and 1
        double precision, parameter :: l_bound  = 1.0 ! 

        integer :: errcode ! Error code

        ! Generate nDims random numbers, stored in the output vector 'random_coordinate'
        ! (v=vector,d=double,rng,uniform)
        errcode=vdrnguniform( method, rng_stream, 1, random_real_vec, u_bound, l_bound)

        random_real = random_real_vec(1)


    end function random_real
    ! ===========================================================================================


    !>  Random integers between 1 and nmax (inclusive)
    !!
    !! Generate nDims integers between 1 and nmax inclusive
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_UNIFORM_STD](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_1_Uniform_VSL_RNG_METHOD_UNIFORM_STD.htm) 
    !! intel method to generate uniform random numbers. Note that this can come in an accurate form

    function random_integer(nDims,nmax)
        implicit none

        !> Size of coordinate vector
        integer :: nDims 

        !> Maximum integer to generate
        integer :: nmax

        ! The output nDims coordinate
        integer, dimension(nDims) :: random_integer

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        integer, parameter :: u_bound  = 1 ! generate random numbers between 0 and 1

        integer :: errcode ! Error code

        ! Generate nDims random numbers, stored in the output vector 'random_coordinate'
        ! (v=vector,i=integer,rng,uniform)
        ! This generates a vector of random integers between [1,nmax+1) = [1,nmax]
        errcode=virnguniform( method, rng_stream, nDims, random_integer, u_bound, nmax+1)


    end function random_integer

    ! ===========================================================================================


    !>  Random direction vector
    !!
    !! Generate a randomly directed unit vector in nDims dimensional space
    !!
    !! This is done by generating nDims gaussian random numbers with mean 0 and
    !! variance 1. These are spherically symetrically distributed, so
    !! normalising this give a random direction.
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_GAUSSIAN_ICDF](https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_4_Gaussian_VSL_RNG_METHOD_GAUSSIAN_ICDF.htm) 
    !! intel method to generate Gaussian random numbers. 

    function random_direction(nDims)
        implicit none

        !> Size of vector to be generated
        integer :: nDims

        ! The output unit vector
        double precision, dimension(nDims) :: random_direction

        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_GAUSSIAN_ICDF

        double precision, parameter :: sigma = 1.0 ! Standard deviation of gaussian random numbers
        double precision, parameter :: mean  = 0.0 ! Mean of gaussian random numbers

        integer :: errcode ! Error code


        ! Generate nDims gaussian random numbers, stored in a vector
        ! random_direction
        ! with mean=0 variance=1
        ! (v=vector,d=double,rng,gaussian)
        errcode=vdrnggaussian( method, rng_stream, nDims, random_direction, mean, sigma)

        ! normalise the vector
        random_direction = random_direction / sqrt( dot_product(random_direction,random_direction) )

    end function random_direction

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
        integer :: nDims 

        ! The output vector
        double precision, dimension(nDims) :: random_point_in_sphere

        ! Temporary variable for storing a random radius
        double precision :: rand_rad(1)


        ! generate a random direction
        random_point_in_sphere = random_direction(nDims)

        ! generate a single random number distributed as beta(r,1) ~ r^{d-1}
        rand_rad = random_beta_vec(nDims+0d0,1d0,1)

        ! Create a point in the sphere
        random_point_in_sphere = random_point_in_sphere * rand_rad(1)


    end function random_point_in_sphere


    ! ===========================================================================================

    !> De-initialise the random number generators

    subroutine deinitialise_random

        implicit none
        integer :: errcode ! Error code

        errcode=vsldeletestream( rng_stream )

    end subroutine deinitialise_random


end module random_module
