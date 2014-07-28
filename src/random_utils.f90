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

    subroutine random_coordinate(coordinates,nDims)
        implicit none

        !> The output nDims coordinate
        double precision, intent(out), dimension(nDims) :: coordinates
        !> Dimensions of the unit vector
        integer,          intent(in)                    :: nDims

        ! Method to generate random numbers 
        ! (This can be upgraded to VSL_RNG_METHOD_UNIFORM_STD_ACCURATE)
        integer,parameter       :: method=VSL_RNG_METHOD_UNIFORM_STD

        double precision, parameter :: u_bound  = 0.0 ! generate random numbers between 0 and 1
        double precision, parameter :: l_bound  = 1.0 ! 

        integer :: errcode ! Error code


        ! Generate nDims random numbers, stored in a vector 'coordinates'
        ! (v=vector,d=double,rng,uniform)
        errcode=vdrnguniform( method, rng_stream, ndims, coordinates, u_bound, l_bound)


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
    !! [VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2](https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/vslnotes/hh_goto.htm#9_3_3_Gaussian_VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2.htm) 
    !! intel method to generate Gaussian random numbers. 

    subroutine random_direction(nhat,nDims)
        implicit none

        !> The output unit vector
        double precision, intent(out), dimension(nDims) :: nhat
        !> Dimensions of the unit vector
        integer,          intent(in)                    :: nDims

        ! Method to generate random numbers 
        integer,parameter       :: method=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2

        double precision, parameter :: sigma = 1.0 ! Standard deviation of gaussian random numbers
        double precision, parameter :: mean  = 0.0 ! Mean of gaussian random numbers

        integer :: errcode ! Error code


        ! Generate nDims gaussian random numbers, stored in a vector nhat
        ! with mean=0 variance=1
        ! (v=vector,d=double,rng,gaussian)
        errcode=vdrnggaussian( method, rng_stream, nDims, nhat, mean, sigma)

        ! normalise the vector
        nhat = nhat / sqrt( dot_product(nhat,nhat) )

        ! nhat is now a random direction

    end subroutine random_direction

    ! ===========================================================================================
    

    !> De-initialise the random number generators

    subroutine deinitialise_random

        implicit none
        integer :: errcode ! Error code

        errcode=vsldeletestream( rng_stream )

    end subroutine deinitialise_random


end module random_module
