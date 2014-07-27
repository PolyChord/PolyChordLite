include 'mkl_vsl.f90'

module random_module

    use mkl_vsl_type
    use mkl_vsl

    implicit none

    type (vsl_stream_state) :: rng_stream
    integer,parameter       :: basic_rng_type = VSL_BRNG_MT19937
    integer                 :: seed


    contains



    subroutine initialise_random(seed_input)
        implicit none
        integer, optional, intent(in) :: seed_input

        integer :: errcode


        if (present(seed_input)) then
            ! If the seed argument is present, initialise stream with this
            errcode=vslnewstream( rng_stream, basic_rng_type,  seed_input )
        else
            ! Otherwise initialise it with the system time
            call system_clock(seed)
            errcode=vslnewstream( rng_stream, brng,  seed )
        end if


    end subroutine initialise_random


    !>  Random direction vector generator
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
        errcode=vdrnggaussian( method, rng_stream, nDims, nhat, mean, sigma)

        ! normalise the vector
        nhat = nhat / sqrt( dot_product(nhat,nhat) )

        ! nhat is now a random direction



    end subroutine random_direction





    subroutine deinitialise_random

        implicit none
        integer :: errcode

        errcode=vsldeletestream( rng_stream )

    end subroutine deinitialise_random


end module random_module
