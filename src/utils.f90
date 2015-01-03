module utils_module

    implicit none

    !> The effective value of \f$ log(0) \f$
    double precision, parameter :: logzero = -sqrt(huge(0d0))
    !> The effective value of \f$ log(\inf) \f$
    double precision, parameter :: loginf = +sqrt(huge(0d0))

    !> The maximum character length
    integer, parameter :: STR_LENGTH = 100

    !> \f$ 2\pi \f$ in double precision
    double precision, parameter :: TwoPi = 8d0*atan(1d0)

    !> The default formats
    integer, parameter :: fmt_len = 200
    character(7) :: DB_FMT='E17.8E3'
    character(4) :: FLT_FMT='F8.2'
    character(3) :: INT_FMT='I20'

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
    !> unit for writing live points
    integer, parameter :: write_phys_unit = 15
    !> unit for writing clusters of live points
    integer, parameter :: write_phys_cluster_unit = 16
    !> unit for writing evidence distribution
    integer, parameter :: write_ev_unit = 17
    !> unit for writing evidence distribution
    integer, parameter :: write_stats_unit = 18
    !> unit for reading unnormalised posterior files
    integer, parameter :: read_untxt_unit = 19
    !> unit for writing unnormalised posterior files
    integer, parameter :: write_untxt_unit = 20


    !> Log[1/2 Erfc[j/Sqrt[2]]]
    double precision, parameter,dimension(20) :: logsigma = (/-1.84102, -3.78318, -6.60773, -10.3601, -15.065, -20.7368, -27.3843, -35.0134, -43.6281, -53.2313, -63.8249, -75.4107, -87.9897, -101.563, -116.131, -131.695, -148.256, -165.812, -184.366, -203.917 /)

    ! All series used to approximate F are computed with relative
    ! tolerance:
    double precision eps
    parameter( eps = 1d-15 )
    ! which means that we neglect all terms smaller than eps times the
    ! current sum


    integer,parameter :: flag_blank     = -2
    integer,parameter :: flag_gestating = -1
    integer,parameter :: flag_waiting   = 0


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

        distance = sqrt(distance2(a,b))

    end function distance

    !> Euclidean distance squared of two coordinates
    !!
    !! returns \f$\sum_i (a_i-b_i)^2 \f$
    function distance2(a,b)
        implicit none
        !> First vector
        double precision, dimension(:) :: a
        !> Second vector
        double precision, dimension(:) :: b

        double precision :: distance2

        distance2 = dot_product(a-b,a-b) 

    end function distance2


    !> Mutual proximity 
    !!
    function MP(a,b,live_data)
        implicit none
        double precision, dimension(:)   :: a
        double precision, dimension(:)   :: b
        double precision, dimension(:,:) :: live_data

        double precision :: MP

        double precision :: dab2

        integer i

        MP=0d0

        dab2 = distance2(a,b)

        do i=1,size(live_data,2)
            if(distance2(live_data(:,i),a) > dab2 .and. distance2(live_data(:,i),b) > dab2 ) MP = MP+1d0
        end do

        MP = MP/(size(live_data,2)+0d0)

    end function MP

    function MP2(seed,baby,live_data)
        implicit none
        double precision, dimension(:)   :: seed
        double precision, dimension(:)   :: baby
        double precision, dimension(:,:) :: live_data

        double precision :: MP2

        double precision :: dab2

        integer i

        MP2=0d0

        dab2 = distance2(seed,baby)

        do i=1,size(live_data,2)
            if(distance2(live_data(:,i),seed) > dab2 ) MP2 = MP2+1d0
        end do

        MP2 = MP2/(size(live_data,2)+0d0)

    end function MP2

    !> Modulus squared of a vector
    !!
    !! returns \f$\sum_i (a_i)^2 \f$
    function mod2(a)
        implicit none
        !> First vector
        double precision, dimension(:) :: a

        double precision :: mod2

        mod2 = dot_product(a,a)

    end function mod2


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

    !> This function increases loga by logb (and by log c if present)
    subroutine logincexp(loga,logb,logc)
        implicit none
        double precision,intent(inout)       :: loga
        double precision,intent(in)          :: logb
        double precision,intent(in),optional :: logc

        if (loga>logb) then
            loga = loga + log(exp(logb-loga) + 1)
        else
            loga = logb + log(exp(loga-logb) + 1)
        end if

        if(present(logc)) then
            if (loga>logc) then
                loga = loga + log(exp(logc-loga) + 1)
            else
                loga = logc + log(exp(loga-logc) + 1)
            end if
        end if

    end subroutine logincexp

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

        if (current_sum <= 0d0) then
            ! this check is useful for entering loops
            abovetol = .true.
        else
            abovetol = abs( change/current_sum ) > eps
        endif

    end function abovetol

    function calc_cholesky(a) result(L)
        implicit none
        double precision, intent(in),dimension(:,:) :: a
        double precision, dimension(size(a,1),size(a,2)) :: L
        integer :: i,j

        ! Set it all to zero to begin with
        L = 0

        ! Zero out the upper half
        do i=1,size(a,1)

            L(i,i)= sqrt( a(i,i) - sum(L(i,:i-1)**2) )

            do j=i+1,size(a,1)
                L(j,i) = (a(i,j) - sum(L(i,:i-1)*L(j,:i-1)))/L(i,i)
            end do

        end do

    end function calc_cholesky

    function calc_covmat(lives,phantoms) result(covmat)
        implicit none
        double precision, intent(in), dimension(:,:) :: lives
        double precision, intent(in), dimension(:,:) :: phantoms

        double precision, dimension(size(lives,1),size(lives,1)) :: covmat

        double precision, dimension(size(lives,1)) :: mean

        integer :: nDims,nlive,nphantom

        nDims = size(lives,1)
        nlive = size(lives,2)
        nphantom = size(phantoms,2)

        ! Compute the mean 
        mean = ( sum(lives,dim=2) + sum(phantoms,dim=2) ) 
        mean = mean / (nlive + nphantom)

        ! Compute the covariance matrix
        covmat = matmul(lives    - spread(mean,dim=2,ncopies=nlive) ,    transpose(lives    - spread(mean,dim=2,ncopies=nlive)    ) ) &
            +    matmul(phantoms - spread(mean,dim=2,ncopies=nphantom) , transpose(phantoms - spread(mean,dim=2,ncopies=nphantom) ) )   
        covmat = covmat / (nlive + nphantom -1)

    end function calc_covmat

    !> This function computes the similarity matrix of an array of data.
    !!
    !! Assume that the data_array can be considered an indexed array of vectors
    !! V = ( v_i : i=1,n )
    !!
    !! The similarity matrix can be expressed very neatly as
    !! d_ij = (v_i-v_j) . (v_i-v_j)
    !!      = v_i.v_i + v_j.v_j - 2 v_i.v_j
    !!
    !! The final term can be written as a data_array^T data_array, and the first
    !! two are easy to write. We can therefore calculate this in two lines with
    !! instrisic functions
    function calc_similarity_matrix(data_array) result(similarity_matrix)

        double precision, intent(in), dimension(:,:) :: data_array

        double precision, dimension(size(data_array,2),size(data_array,2)) :: similarity_matrix

        integer :: i



        similarity_matrix = spread( [( dot_product(data_array(:,i),data_array(:,i)),i=1,size(data_array,2) )], dim=2,ncopies=size(data_array,2) )

        similarity_matrix = similarity_matrix + transpose(similarity_matrix) - 2d0 * matmul( transpose(data_array),data_array )

    end function calc_similarity_matrix




    function normal_cdf(x)
        implicit none
        double precision, intent(in),dimension(:) :: x
        double precision,dimension(size(x)) :: normal_cdf
        double precision, parameter :: sqrt2 = sqrt(2d0)

        normal_cdf = 0.5d0 * ( 1d0 + erf(x/sqrt2))

    end function normal_cdf

    function inv_normal_cdf(x)
        implicit none
        double precision, intent(in),dimension(:) :: x
        double precision,dimension(size(x)) :: inv_normal_cdf
        integer :: i

        inv_normal_cdf = [( r8_normal_01_cdf_inverse(x(i)), i=1,size(x) )]

    end function inv_normal_cdf


    ! This was downloaded from:
    ! http://people.sc.fsu.edu/~jburkardt/f_src/asa241/asa241.f90

    function r8_normal_01_cdf_inverse ( p )

        !*****************************************************************************80
        !
        !! R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
        !
        !  Discussion:
        !
        !    The result is accurate to about 1 part in 10**16.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    27 December 2004
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Michael Wichura.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Michael Wichura,
        !    The Percentage Points of the Normal Distribution,
        !    Algorithm AS 241,
        !    Applied Statistics,
        !    Volume 37, Number 3, pages 477-484, 1988.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) P, the value of the cumulative probability 
        !    densitity function.  0 < P < 1.  If P is outside this range,
        !    an "infinite" value will be returned.
        !
        !    Output, real ( kind = 8 ) D_NORMAL_01_CDF_INVERSE, the normal deviate 
        !    value with the property that the probability of a standard normal 
        !    deviate being less than or equal to the value is P.
        !
        implicit none

        real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
            3.3871328727963666080D+00, &
            1.3314166789178437745D+02, &
            1.9715909503065514427D+03, &
            1.3731693765509461125D+04, &
            4.5921953931549871457D+04, &
            6.7265770927008700853D+04, &
            3.3430575583588128105D+04, &
            2.5090809287301226727D+03 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
            1.0D+00, &
            4.2313330701600911252D+01, &
            6.8718700749205790830D+02, &
            5.3941960214247511077D+03, &
            2.1213794301586595867D+04, &
            3.9307895800092710610D+04, &
            2.8729085735721942674D+04, &
            5.2264952788528545610D+03 /)
        real   ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
            1.42343711074968357734D+00, &
            4.63033784615654529590D+00, &
            5.76949722146069140550D+00, &
            3.64784832476320460504D+00, &
            1.27045825245236838258D+00, &
            2.41780725177450611770D-01, &
            2.27238449892691845833D-02, &
            7.74545014278341407640D-04 /)
        real ( kind = 8 ), parameter :: const1 = 0.180625D+00
        real ( kind = 8 ), parameter :: const2 = 1.6D+00
        real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
            1.0D+00, &
            2.05319162663775882187D+00, &
            1.67638483018380384940D+00, &
            6.89767334985100004550D-01, &
            1.48103976427480074590D-01, &
            1.51986665636164571966D-02, &
            5.47593808499534494600D-04, &
            1.05075007164441684324D-09 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
            6.65790464350110377720D+00, &
            5.46378491116411436990D+00, &
            1.78482653991729133580D+00, &
            2.96560571828504891230D-01, &
            2.65321895265761230930D-02, &
            1.24266094738807843860D-03, &
            2.71155556874348757815D-05, &
            2.01033439929228813265D-07 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
            1.0D+00, &
            5.99832206555887937690D-01, &
            1.36929880922735805310D-01, &
            1.48753612908506148525D-02, &
            7.86869131145613259100D-04, &
            1.84631831751005468180D-05, &
            1.42151175831644588870D-07, &
            2.04426310338993978564D-15 /)
        real ( kind = 8 ) p
        real ( kind = 8 ) q
        real ( kind = 8 ) r
        real ( kind = 8 ) r8_normal_01_cdf_inverse 
        real ( kind = 8 ), parameter :: split1 = 0.425D+00
        real ( kind = 8 ), parameter :: split2 = 5.0D+00

        if ( p <= 0.0D+00 ) then
            r8_normal_01_cdf_inverse = - huge ( p )
            return
        end if

        if ( 1.0D+00 <= p ) then
            r8_normal_01_cdf_inverse = huge ( p )
            return
        end if

        q = p - 0.5D+00

        if ( abs ( q ) <= split1 ) then

            r = const1 - q * q
            r8_normal_01_cdf_inverse = q * r8poly_value ( 8, a, r ) &
                / r8poly_value ( 8, b, r )

        else

            if ( q < 0.0D+00 ) then
                r = p
            else
                r = 1.0D+00 - p
            end if

            if ( r <= 0.0D+00 ) then
                r8_normal_01_cdf_inverse = - 1.0D+00
                stop
            end if

            r = sqrt ( -log ( r ) )

            if ( r <= split2 ) then

                r = r - const2
                r8_normal_01_cdf_inverse = r8poly_value ( 8, c, r ) &
                    / r8poly_value ( 8, d, r )

            else

                r = r - split2
                r8_normal_01_cdf_inverse = r8poly_value ( 8, e, r ) &
                    / r8poly_value ( 8, f, r )

            end if

            if ( q < 0.0D+00 ) then
                r8_normal_01_cdf_inverse = - r8_normal_01_cdf_inverse
            end if

        end if

        return
    end function r8_normal_01_cdf_inverse


    function r8poly_value ( n, a, x )

        !*****************************************************************************80
        !
        !! R8POLY_VALUE evaluates an R8POLY
        !
        !  Discussion:
        !
        !    For sanity's sake, the value of N indicates the NUMBER of 
        !    coefficients, or more precisely, the ORDER of the polynomial,
        !    rather than the DEGREE of the polynomial.  The two quantities
        !    differ by 1, but cause a great deal of confusion.
        !
        !    Given N and A, the form of the polynomial is:
        !
        !      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    13 August 2004
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) N, the order of the polynomial.
        !
        !    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
        !    A(1) is the constant term.
        !
        !    Input, real ( kind = 8 ) X, the point at which the polynomial is 
        !    to be evaluated.
        !
        !    Output, real ( kind = 8 ) R8POLY_VALUE, the value of the polynomial at X.
        !
        implicit none

        integer ( kind = 4 ) n

        real ( kind = 8 ) a(n)
        integer ( kind = 4 ) i
        real ( kind = 8 ) r8poly_value
        real ( kind = 8 ) x

        r8poly_value = 0.0D+00
        do i = n, 1, -1
            r8poly_value = r8poly_value * x + a(i)
        end do

        return
    end function r8poly_value



end module utils_module
