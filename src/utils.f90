
module utils_module

    implicit none
include 'mkl_vml.f90'  

    !> The effective value of \f$ log(0) \f$
    double precision, parameter :: logzero = -sqrt(huge(0d0))
    !> The effective value of \f$ log(\inf) \f$
    double precision, parameter :: loginf = +sqrt(huge(0d0))

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
    !> unit for writing live points
    integer, parameter :: write_phys_unit = 15
    !> unit for writing clusters of live points
    integer, parameter :: write_phys_cluster_unit = 16
    !> unit for writing evidence distribution
    integer, parameter :: write_ev_unit = 17

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

    !> Wrapper for the intel vml lgamma function
    !!
    !! This function returns log(gamma(x))
    function lgamma(x)
        implicit none
        double precision :: x

        double precision :: lgamma

        double precision, dimension(1) :: input
        double precision, dimension(1) :: output

        input(1) = x
        call vdlgamma(1,input,output)
        lgamma = output(1)

    end function lgamma






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

        if (current_sum == 0d0) then
            ! this check is useful for entering loops
            abovetol = .true.
        else
            abovetol = abs( change/current_sum ) > eps
        endif

    end function abovetol

    function calc_cholesky(matrix) result(cholesky)
        implicit none
        double precision, intent(in),dimension(:,:) :: matrix
        double precision, dimension(size(matrix,1),size(matrix,2)) :: cholesky
        integer :: errcode
        integer :: i,j

        cholesky = matrix
        call dpotrf ('L',size(matrix,1),cholesky,size(matrix,1),errcode)

        ! Zero out the upper half
        do i=1,size(matrix,1)
            do j=i+1,size(matrix,2)
                cholesky(i,j)=0d0
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



end module utils_module
