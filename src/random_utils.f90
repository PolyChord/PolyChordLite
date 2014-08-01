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
    !! intel method to generate Gaussian random numbers. Note that this can come
    !! in an accurate form

    subroutine random_coordinate(coordinates)
        implicit none

        !> The output nDims coordinate
        double precision, intent(out), dimension(:) :: coordinates

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        double precision, parameter :: u_bound  = 0.0 ! generate random numbers between 0 and 1
        double precision, parameter :: l_bound  = 1.0 ! 

        integer :: errcode ! Error code

        integer :: nDims ! Size of coordinate vector


        ! Get the dimensionality of the vector to be generated
        nDims = size(coordinates)

        ! Generate nDims random numbers, stored in a vector 'coordinates'
        ! (v=vector,d=double,rng,uniform)
        errcode=vdrnguniform( method, rng_stream, nDims, coordinates, u_bound, l_bound)


    end subroutine random_coordinate

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

    subroutine random_direction(nhat)
        implicit none

        !> The output unit vector
        double precision, intent(out), dimension(:) :: nhat

        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_GAUSSIAN_ICDF

        double precision, parameter :: sigma = 1.0 ! Standard deviation of gaussian random numbers
        double precision, parameter :: mean  = 0.0 ! Mean of gaussian random numbers

        integer :: errcode ! Error code

        integer :: nDims ! Size of coordinate vector


        ! Get the dimensionality of the vector to be generated
        nDims = size(nhat)

        ! Generate nDims gaussian random numbers, stored in a vector nhat
        ! with mean=0 variance=1
        ! (v=vector,d=double,rng,gaussian)
        errcode=vdrnggaussian( method, rng_stream, nDims, nhat, mean, sigma)

        ! normalise the vector
        nhat = nhat / sqrt( dot_product(nhat,nhat) )

        ! nhat is now a random direction

    end subroutine random_direction

    ! ===========================================================================================


    !>  Random point in unit sphere
    !!
    !! Generate a uniformly distributed point in the unit n-sphere
    !!
    !! This is done by generating a random direction, and then multiplying it by
    !! a radius distributed as \f$ P(r) ~ r^{d-1}\f$. where \f$d$\f is the
    !! dimension of the input vector

    subroutine random_point_in_sphere(point)
        implicit none

        !> The output unit vector
        double precision, intent(out), dimension(:) :: point

        integer :: nDims ! Size of coordinate vector

        double precision :: rand_rad(1)


        ! Get the dimensionality of the vector to be generated
        nDims = size(point)

        ! generate a random direction
        call random_direction(point)

        ! generate a random radius
        call random_beta_vec(rand_rad, nDims+0d0,1d0)

        ! Create a point in the sphere
        point = point * rand_rad(1)


    end subroutine random_point_in_sphere

    ! ===========================================================================================


    !> Random point in beta distribution
    !!
    !! Generate a vector of points distributed as \f$ P(r;p,q) ~ r^{p-1}(1-r)^{q-1}\f$
    !!
    !! We use the 
    !! [VSL_RNG_METHOD_BETA_CJA](
    !! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_17_Beta_VSL_RNG_METHOD_BETA_CJA.htm#9_3_16_Beta_VSL_RNG) 
    !! intel method to generate points. Note that this can come in an accurate form                               

    subroutine random_beta_vec(point,p,q)
        implicit none

        !> The output unit vector
        double precision, intent(out), dimension(:) :: point

        !> The input beta parameters
        double precision, intent(in) :: p
        double precision, intent(in) :: q

        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_BETA_CJA

        double precision, parameter :: offset = 0.0 ! Mean of the beta distribution
        double precision, parameter :: scale_factor  = 1.0 ! width of the beta distribution

        integer :: errcode ! Error code

        integer :: nDims ! Size of coordinate vector



        ! Get the dimensionality of the vector to be generated
        nDims = size(point)

        ! call the intel procedure to generate the point
        errcode=vdrngbeta( method, rng_stream, nDims, point, p, q, offset, scale_factor )


    end subroutine random_beta_vec

    ! ===========================================================================================

    !> De-initialise the random number generators

    subroutine deinitialise_random

        implicit none
        integer :: errcode ! Error code

        errcode=vsldeletestream( rng_stream )

    end subroutine deinitialise_random


end module random_module
