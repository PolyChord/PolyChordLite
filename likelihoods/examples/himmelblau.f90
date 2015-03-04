module loglikelihood_module

    contains

    !> Upside down [Himmelblau function](http://en.wikipedia.org/wiki/Himmelblaus_function).
    !!
    !! \f[ -\log L(x, y) = (x^2+y-11)^2 + (x+y^2-7)^2 \f]
    !!
    !! Himmelblau's function is a multi-modal function, used to test the performance of 
    !! optimization  algorithms. 
    !! It has one local maximum at <math>x = -0.270845 \,</math> and <math>y =
    !! -0.923039 \,</math> where <math>f(x,y) = 181.617 \,</math>, and four
    !! identical local minima:
    !! * \f$ f(3.0, 2.0) = 0.0,            \f$
    !! * \f$ f(-2.805118, 3.131312) = 0.0, \f$
    !! * \f$ f(-3.779310, -3.283186) = 0.0,\f$
    !! * \f$ f(3.584428, -1.848126) = 0.0. \f$
    !!
    !! The locations of all the Maxima and minima can be found
    !! analytically. However, because they are roots of cubic polynomials, when
    !! written in terms of radicals, the expressions are somewhat
    !! complicated.
    !!
    !! The function is named after David Mautner Himmelblau (1924-2011), who
    !! introduced it. 
    !! 
    !! Note that this is a 2D function, and should hence only be used for settings%nDims=2
    !!
    function loglikelihood(theta,phi)
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision                            :: loglikelihood ! loglikelihood value to output

        ! Normalisation
        loglikelihood = -log(0.4071069421432255d0) 

        ! Evaluate
        loglikelihood =  loglikelihood  -  (theta(1)**2 + theta(2)- 11d0)**2    -   (theta(1)+theta(2)**2-7d0)**2 

    end function loglikelihood

    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

    end subroutine setup_loglikelihood

end module loglikelihood_module
