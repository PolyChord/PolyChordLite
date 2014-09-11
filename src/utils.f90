module utils_module

    !> The effective value of \f$ log(0) \f$
    double precision, parameter :: logzero = -huge(0d0)
    !> The effective value of \f$ log(\inf) \f$
    double precision, parameter :: loginf = +huge(0d0) 

    !> The maximum character length
    integer, parameter :: STR_LENGTH = 100

    !> \f$ 2\pi \f$ in double precision
    double precision, parameter :: TwoPi = 8d0*atan(1d0)

    !> The default double format
    !!
    !! should have write statements along the lines of 
    !! write(*,'(E<DBL_FMT(1)>.<DBL_FMT(2)>)')
    integer, parameter, dimension(2) :: DBL_FMT=(/17,8/)

    !> unit for stdout
    integer, parameter :: stdout_unit = 6
    !> unit for reading from resume file
    integer, parameter :: read_resume_unit = 10
    !> unit for writing to resum file
    integer, parameter :: write_resume_unit = 11
    !> unit for writing to txt file
    integer, parameter :: write_txt_unit = 12
    !> unit for writing dead file
    integer, parameter :: write_dead_unit = 13
    !> unit for reading covariance matrices
    integer, parameter :: read_covmat_unit = 14
    !> unit for reading covariance matrices
    integer, parameter :: write_phys_unit = 15
    !> unit for writing evidence distribution
    integer, parameter :: write_ev_unit = 16

    ! All series used to approximate F are computed with relative
    ! tolerance:
    double precision eps
    parameter( eps = 1d-15 )
    ! which means that we neglect all terms smaller than eps times the
    ! current sum

    contains



    !> Euclidean distance of two coordinates
    !!
    !! returns \f$\sqrt{\sum_i (a_i-b_i)^2 } \f$
    function distance(a,b)
        implicit none
        !> First vector
        double precision, dimension(:) :: a
        !> Second vector
        double precision, dimension(:) :: b

        double precision :: distance

        distance = sqrt( dot_product(a-b,a-b) )

    end function distance


    !> Double comparison
    function dbleq(a,b)
        implicit none
        double precision :: a,b
        logical :: dbleq
        double precision, parameter :: eps = 1d-7

        dbleq =  abs(a-b) < eps * max(abs(a),abs(b)) 

    end function dbleq







    !> How to actually calculate sums from logs.
    !!
    !! i.e. if one has a set of logarithms \f$\{\log(L_i)\}\f$, how should one
    !! calculate \f$ \log(\sum_i L_i)\f$ without underflow?
    !!
    !! One does it with the 'log-sum-exp' trick, by subtracting off the maximum
    !! value so that at least the maxmimum value doesn't underflow, and then add
    !! it back on at the end:
    function logsumexp(vector)
        implicit none
        !> vector of log(w
        double precision, dimension(:),intent(in) :: vector

        double precision :: logsumexp
        double precision :: maximumlog

        maximumlog = maxval(vector)

        logsumexp =  maximumlog + log(sum(exp(vector - maximumlog)))

    end function logsumexp


    function logaddexp(loga,logb)
        implicit none
        double precision :: loga
        double precision :: logb
        double precision :: logaddexp

        if (loga>logb) then
            logaddexp = loga + log(exp(logb-loga) + 1)
        else
            logaddexp = logb + log(exp(loga-logb) + 1)
        end if

    end function logaddexp

    function logsubexp(loga,logb)
        implicit none
        double precision :: loga
        double precision :: logb
        double precision :: logsubexp

        if(loga>logb) then
            logsubexp = loga + log(1-exp(logb-loga))
        else 
            logsubexp = logzero
        end if

    end function logsubexp


    !> Hypergeometric function 1F1
    !!
    !!
    function Hypergeometric1F1(a,b,z)
        implicit none
        double precision, intent(in)  :: a,b,z
        double precision              :: Hypergeometric1F1
        integer n
        double precision change

        ! This computes the hypergeometric 1F1 function using a
        ! truncated power series, stopping when the relative change
        ! due to higher terms is less than epsilon
        !
        ! (note that this is very similar to 2F1, but without the c)

        ! http://en.wikipedia.org/wiki/Hypergeometric_function#The_hypergeometric_series

        Hypergeometric1F1=0d0
        n=0

        do while( abovetol(change,Hypergeometric1F1) )
           change = Pochhammer(a,n) * Pochhammer(b,n) * z**n / gamma(1d0+n)

           Hypergeometric1F1 = Hypergeometric1F1 + change
           n=n+1
        enddo

    end function Hypergeometric1F1




    !> Hypergeometric function 
    function Hypergeometric2F1(a,b,c,z)
        implicit none
        double precision, intent(in)  :: a,b,c,z
        double precision              :: Hypergeometric2F1
        integer n
        double precision change

        ! This computes the hypergeometric 2F1 function using a
        ! truncated power series, stopping when the relative change
        ! due to higher terms is less than epsilon

        ! http://en.wikipedia.org/wiki/Hypergeometric_function#The_hypergeometric_series

        Hypergeometric2F1=0d0
        n=0

        do while( abovetol(change,Hypergeometric2F1) )
           change = Pochhammer(a,n) * Pochhammer(b,n) / Pochhammer(c,n) * z**n / gamma(1d0+n)

           Hypergeometric2F1 = Hypergeometric2F1 + change
           n=n+1
        enddo

    end function Hypergeometric2F1




    recursive function Pochhammer (x,n) result (xn)
        ! This function computes the rising factorial x^(n):
        ! for a non-negative integer n and real x
        !
        ! http://en.wikipedia.org/wiki/Pochhammer_symbol

        implicit none
        double precision, intent(in)  :: x
        integer,          intent(in)  :: n
        double precision              :: xn

        if (n<=0) then
            xn = 1
        else
            xn = Pochhammer(x,n-1)*(x+n-1)
        endif

    end function Pochhammer

    function abovetol (change,current_sum)
        ! This function outputs true if change is outside of the
        ! tolerance epsilon from current_sum
        implicit none
        double precision, intent(in) :: change,current_sum
        logical abovetol

        if (current_sum == 0d0) then
            ! this check is useful for entering loops
            abovetol = .true.
        else
            abovetol = abs( change/current_sum ) > eps
        endif

    end function abovetol



end module utils_module
