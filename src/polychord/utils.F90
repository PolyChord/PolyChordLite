
module utils_module
    use, intrinsic :: ieee_arithmetic
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, parameter :: dp = kind(1.d0)

    !> The maximum character length
    integer, parameter :: STR_LENGTH = 300

    !> \f$ 2\pi \f$ in real(dp)
    real(dp), parameter :: TwoPi = 8d0*atan(1d0)
    !> \f$ \log(2\pi) \f$ in real(dp)
    real(dp), parameter :: logTwoPi = log(8d0*atan(1d0))

    !> The default writing formats
    integer, parameter :: fmt_len = 200
    character(8) :: DB_FMT='E24.15E3'
    character(4) :: FLT_FMT='F8.2'
    character(3) :: INT_FMT='I12'

    !> Feedback levels
    integer, parameter :: title_fb = 0
    integer, parameter :: normal_fb = 1
    integer, parameter :: fancy_fb = 2
    integer, parameter :: verbose_fb = 3

    !> unit for stderr
    integer, parameter :: stderr_unit = 0
    !> unit for stdin
    integer, parameter :: stdin_unit = 5
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
    !> unit for writing equal weights files
    integer, parameter :: write_equals_unit = 19
    !> unit for writing weights files
    integer, parameter :: write_posterior_unit = 19

    !> unit for deleting generic files
    integer, parameter :: delete_unit = 20


    !> unit for params file
    integer, parameter :: params_unit = 21
    !> unit for paramnames file
    integer, parameter :: paramnames_unit = 22

    !> unit for prior file
    integer, parameter :: write_prior_unit = 23

    !> unit for writing dead birth file
    integer, parameter :: write_dead_birth_unit = 24

    !> unit for writing live birth file
    integer, parameter :: write_live_birth_unit = 25

    !> unit for writing max file
    integer, parameter :: write_max_unit = 26

    !> unit for properties file
    integer, parameter :: properties_unit = 27

    ! All series used to approximate F are computed with relative
    ! tolerance:
    real(dp) eps
    parameter( eps = 1d-15 )
    ! which means that we neglect all terms smaller than eps times the
    ! current sum

    !> Minimum variance for regularization of covariance matrices
    !! Used when clusters have insufficient points or near-zero variance
    real(dp), parameter :: min_variance = 1.d-20

    !> Relative threshold for detecting zero eigenvalues in regularization
    real(dp), parameter :: eigenvalue_rel_threshold = 1.d-10
    !> Absolute floor threshold for eigenvalues (prevents issues when max eigenvalue is tiny)
    real(dp), parameter :: eigenvalue_abs_threshold = 1.d-20

#ifdef HAVE_LAPACK
    !> LAPACK interface for symmetric eigenvalue decomposition
    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            integer, intent(out) :: info
            real(dp), intent(inout) :: a(lda,*)
            real(dp), intent(out) :: w(*), work(*)
        end subroutine dsyev
    end interface
#endif

    integer,parameter :: flag_blank     = -2
    integer,parameter :: flag_gestating = -1
    integer,parameter :: flag_waiting   = 0


    contains
    
    ! Write format for n integers
    function integer_format(n)
        implicit none
        integer, intent(in) :: n
        character(len=fmt_len) :: integer_format

        write(integer_format,'("(",I0,A,")")') n,INT_FMT   ! define the integer format

    end function integer_format

    ! Write format for n doubles
    function double_format(n)
        implicit none
        integer, intent(in) :: n
        character(len=fmt_len) :: double_format

        write(double_format,'("(",I0,A,")")') n,DB_FMT   ! define the integer format

    end function double_format

    !> Swaps two integers via a temporary variable
    subroutine swap_integers(a,b)
        implicit none
        integer,intent(inout) :: a,b
        integer :: temp
        temp=a
        a=b
        b=temp
    end subroutine swap_integers

    !> location of minimum in an array of doubles
    !!
    !! This is just a wrapper around minloc, but gets around the annoying ideosyncracy
    !! of fortran that minloc must return a length 1 array
    function minpos(a)
        implicit none
        real(dp), intent(in), dimension(:) :: a
        integer :: minpos
        integer :: minpos_vec(1)

        minpos_vec = minloc(a)
        minpos = minpos_vec(1)
    end function minpos



    function cyc(iterator,cycle_size)
        implicit none
        integer, intent(in) :: iterator
        integer, intent(in) :: cycle_size
        logical :: cyc

        if(cycle_size<=0) then
            cyc = .false.
        else
            cyc = mod(iterator,cycle_size)==0
        end if

    end function cyc



    !> Euclidean distance of two coordinates
    !!
    !! returns \f$\sqrt{\sum_i (a_i-b_i)^2 } \f$
    function distance(a,b,w)
        implicit none
        !> First vector
        real(dp), dimension(:) :: a
        !> Second vector
        real(dp), dimension(:) :: b
        !> Wraparound parameters
        logical, dimension(:) :: w

        real(dp) :: distance

        distance = sqrt(distance2(a,b,w))

    end function distance

    !> Euclidean distance squared of two coordinates
    !!
    !! returns \f$\sum_i (a_i-b_i)^2 \f$
    function distance2(a,b,w)
        implicit none
        !> First vector
        real(dp), dimension(:) :: a
        !> Second vector
        real(dp), dimension(:) :: b
        !> Wraparound parameters
        logical, dimension(:) :: w

        real(dp), dimension(size(a)) :: dx

        real(dp) :: distance2

        dx = a-b
        where(w) dx = dx - nint(dx)
        distance2 = dot_product(dx, dx) 
    end function distance2

    function loggamma(n)
        use iso_c_binding
        implicit none
        real(dp)            :: loggamma
        real(dp),intent(in) :: n
        real(c_double) :: c_n
        interface 
            function lgamma (y) bind(c)
                use iso_c_binding
                implicit none
                real(c_double)        :: lgamma
                real(c_double), value :: y
            end function
        end interface


        c_n = n

        loggamma = lgamma(c_n)


    end function loggamma


    !> Mutual proximity 
    !!
    function MP(a,b,live_data,w)
        implicit none
        real(dp), dimension(:)   :: a
        real(dp), dimension(:)   :: b
        real(dp), dimension(:,:) :: live_data
        !> Wraparound parameters
        logical, dimension(:) :: w

        real(dp) :: MP

        real(dp) :: dab2

        integer i

        MP=0d0

        dab2 = distance2(a,b,w)

        do i=1,size(live_data,2)
            if(distance2(live_data(:,i),a,w) > dab2 .and. distance2(live_data(:,i),b,w) > dab2 ) MP = MP+1d0
        end do

        MP = MP/(size(live_data,2)+0d0)

    end function MP

    function MP2(seed,baby,live_data,w)
        implicit none
        real(dp), dimension(:)   :: seed
        real(dp), dimension(:)   :: baby
        real(dp), dimension(:,:) :: live_data
        !> Wraparound parameters
        logical, dimension(:) :: w

        real(dp) :: MP2

        real(dp) :: dab2

        integer i

        MP2=0d0

        dab2 = distance2(seed,baby,w)

        do i=1,size(live_data,2)
            if(distance2(live_data(:,i),seed,w) > dab2 ) MP2 = MP2+1d0
        end do

        MP2 = MP2/(size(live_data,2)+0d0)

    end function MP2

    !> Modulus squared of a vector
    !!
    !! returns \f$\sum_i (a_i)^2 \f$
    function mod2(a)
        implicit none
        !> First vector
        real(dp), dimension(:) :: a

        real(dp) :: mod2

        mod2 = dot_product(a,a)

    end function mod2


    !> Double comparison
    function dbleq(a,b)
        implicit none
        real(dp) :: a,b
        logical :: dbleq
        real(dp), parameter :: eps = 1d-7

        dbleq =  abs(a-b) < eps * max(abs(a),abs(b)) 

    end function dbleq



    !> Identity matrix ( nDims x nDims )
    function identity_matrix(nDims)
        implicit none
        !> dimensionality of the identity matrix
        integer,intent(in) :: nDims
        !> The identity matrix to be returned
        real(dp), dimension(nDims,nDims) :: identity_matrix

        integer :: i_dims ! iterator over dimensions

        identity_matrix=0d0
        do i_dims=1,nDims
            identity_matrix(i_dims,i_dims) = 1d0
        end do


    end function identity_matrix

    !> Trace of a matrix
    function trace(a)
        implicit none
        !> The identity matrix to be returned
        real(dp), dimension(:,:),intent(in) :: a

        real(dp) :: trace

        integer :: i ! iterator over dimensions

        trace= sum ( [( a(i,i), i=1,min(size(a,1),size(a,2)) )] ) 

    end function trace


    function delete_file(file_name,feedback) result(deleted)
        implicit none
        character(STR_LENGTH),intent(in) :: file_name
        logical, optional, intent(in) :: feedback

        logical :: deleted ! whether or not there was a file to be deleted

        ! Check that file exists:
        inquire( file=trim(file_name), exist=deleted)

        if(deleted) then
            if(present(feedback)) then
                if(feedback) write(stdout_unit,'("Deleting file: ", A)') trim(file_name)
            end if
            ! open the file
            open(delete_unit,file=trim(file_name)) 
            ! Delete it if it exists
            close(delete_unit,status='delete')
        end if


    end function delete_file






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
        real(dp), dimension(:),intent(in) :: vector

        real(dp) :: logsumexp
        real(dp) :: maximumlog

        maximumlog = maxval(vector)
        logsumexp  =  maximumlog + log(sum(exp(vector - maximumlog)))

    end function logsumexp


    function logaddexp(loga,logb)
        implicit none
        real(dp) :: loga
        real(dp) :: logb
        real(dp) :: logaddexp

        if (loga>logb) then
            logaddexp = loga + log(exp(logb-loga) + 1)
        else
            logaddexp = logb + log(exp(loga-logb) + 1)
        end if

    end function logaddexp

    function logaddexp_(loga,logb)
        implicit none
        real(dp), dimension(:) :: loga
        real(dp), dimension(size(loga)) :: logb
        real(dp), dimension(size(loga)) :: logaddexp_

        where (loga>logb)
            logaddexp_ = loga + log(exp(logb-loga) + 1)
        else where
            logaddexp_ = logb + log(exp(loga-logb) + 1)
        end where

    end function logaddexp_

    function logsubexp(loga,logb)
        implicit none
        real(dp) :: loga
        real(dp) :: logb
        real(dp) :: logsubexp

        if(loga>logb) then
            logsubexp = loga + log(1-exp(logb-loga))
        else 
            logsubexp = -huge(1d0)
        end if

    end function logsubexp

    !> This function increases loga by logb (and by log c if present)
    subroutine logincexp(loga,logb,logc)
        implicit none
        real(dp),intent(inout)       :: loga
        real(dp),intent(in)          :: logb
        real(dp),intent(in),optional :: logc

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


    !> Return the sorted indices of an array of doubles
    !!
    function sort_doubles(a) result(k)
        implicit none
        real(dp),intent(in), dimension(:) :: a
        real(dp), dimension(size(a)) :: b
        integer, dimension(size(a)) :: i
        integer, dimension(size(a)) :: k

        integer :: j

        b = a
        i = [ (j,j=1,size(a)) ]

        call quicksort(b,i)

        k=i
        !do j=1,size(a)
            !k(i(j)) = j
        !end do

    end function sort_doubles

    recursive subroutine quicksort(A,Ai)
        real(dp), intent(inout), dimension(:) :: A
        integer, intent(inout), dimension(:) :: Ai
        integer :: iq

        if(size(A) > 1) then
            call Partition(A,Ai,iq)
            call quicksort(A(:iq-1),Ai(:iq-1))
            call quicksort(A(iq:),Ai(iq:))
        endif
    end subroutine quicksort

    subroutine Partition(A,Ai,marker)
        real(dp), intent(inout), dimension(:) :: A
        integer, intent(inout), dimension(:) :: Ai
        integer, intent(out) :: marker
        integer :: i, j
        real(dp) :: temp
        integer :: tempi
        real(dp) :: x      ! pivot point
        x = A(1)
        i= 0
        j= size(A) + 1

        do
            j = j-1
            do
                if (A(j) <= x) exit
                j = j-1
            end do
            i = i+1
            do
                if (A(i) >= x) exit
                i = i+1
            end do
            if (i < j) then
                ! exchange A(i) and A(j)
                temp = A(i)
                A(i) = A(j)
                A(j) = temp
                tempi = Ai(i)
                Ai(i) = Ai(j)
                Ai(j) = tempi

            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do

    end subroutine Partition





    !> Hypergeometric function 1F1
    !!
    !!
    function Hypergeometric1F1(a,b,z)
        implicit none
        real(dp), intent(in)  :: a,b,z
        real(dp)              :: Hypergeometric1F1
        integer n
        real(dp) change

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
        real(dp), intent(in)  :: a,b,c,z
        real(dp)              :: Hypergeometric2F1
        integer n
        real(dp) change

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
    real(dp), intent(in)  :: x
    integer,          intent(in)  :: n
    real(dp)              :: xn

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
        real(dp), intent(in) :: change,current_sum
        logical abovetol

        if (current_sum <= 0d0) then
            ! this check is useful for entering loops
            abovetol = .true.
        else
            abovetol = abs( change/current_sum ) > eps
        endif

    end function abovetol

    !> Compute Cholesky decomposition of a covariance matrix
    !!
    !! If n_points is provided and n_points < nDims + 1, the covariance matrix
    !! is rank-deficient. In this case, we first apply eigendecomposition-based
    !! regularization to make it positive-definite before computing Cholesky.
    !!
    !! @param[in] a         Input covariance matrix
    !! @param[in] n_points  Optional: number of points used to compute covariance
    !! @return    L         Lower-triangular Cholesky factor
    function calc_cholesky(a, n_points) result(L)
        implicit none
        real(dp), intent(in),dimension(:,:) :: a
        integer, intent(in), optional :: n_points
        real(dp), dimension(size(a,1),size(a,2)) :: a_reg
        real(dp), dimension(size(a,1),size(a,2)) :: L
        integer :: i,j, nDims
        logical :: needs_regularization

        nDims = size(a, 1)

        ! Determine if eigendecomposition-based regularization is needed
        ! This is required when the covariance matrix is rank-deficient (n < D+1)
        needs_regularization = .false.
        if (present(n_points)) then
            if (n_points < nDims + 1) then
                needs_regularization = .true.
            end if
        end if

        if (needs_regularization) then
            ! Use eigendecomposition-based regularization for rank-deficient matrices
            ! This replaces zero eigenvalues with 1.0 (unit variance in null-space)
            a_reg = regularize_covmat(a)
        else
            ! === FIX 3: Pre-regularize the matrix to avoid numerical issues ===
            ! Add min_variance to diagonal before attempting decomposition
            ! This prevents failure for near-singular matrices from tight clusters
            a_reg = a + identity_matrix(nDims) * min_variance
        end if

        ! Set it all to zero to begin with
        L = 0

        ! Zero out the upper half
        do i=1,size(a,1)

            L(i,i)= a_reg(i,i) - sum(L(i,:i-1)**2)
            if (L(i,i).le.0d0) then
                ! If the cholesky decomposition does not exist, then set it to
                ! be a re-scaled identity matrix (original PolyChord behavior)
                L = identity_matrix(nDims) * sqrt(trace(a_reg))
                return
            else
                L(i,i)=sqrt(L(i,i))
            end if

            do j=i+1,size(a_reg,1)
                L(j,i) = (a_reg(i,j) - sum(L(i,:i-1)*L(j,:i-1)))/L(i,i)
            end do

        end do

    end function calc_cholesky

    !> Compute eigenvalues and eigenvectors of a symmetric matrix using Jacobi method
    !!
    !! This is a pure-Fortran fallback for systems without LAPACK.
    !! Uses the cyclic Jacobi algorithm which iteratively applies Givens rotations
    !! to diagonalize the matrix.
    !!
    !! @param[inout] a  On entry: symmetric matrix. On exit: destroyed.
    !! @param[out]   w  Eigenvalues in ascending order
    !! @param[out]   v  Eigenvectors as columns (v(:,i) is eigenvector for w(i))
    !! @param[out]   info  0 on success, >0 if failed to converge
    subroutine jacobi_eigen(a, w, v, info)
        implicit none
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: w
        real(dp), intent(out), dimension(:,:) :: v
        integer, intent(out) :: info

        integer :: n, i, j, k, p, q, sweep
        real(dp) :: threshold, off_diag, c, s, t, tau, h, g, temp
        real(dp) :: tol, sum_off
        integer, parameter :: max_sweeps = 50

        n = size(a, 1)
        info = 0
        tol = 1.0d-15

        ! Initialize eigenvector matrix to identity
        v = 0.0d0
        do i = 1, n
            v(i,i) = 1.0d0
        end do

        ! Main Jacobi iteration
        do sweep = 1, max_sweeps
            ! Compute sum of off-diagonal elements
            sum_off = 0.0d0
            do i = 1, n-1
                do j = i+1, n
                    sum_off = sum_off + abs(a(i,j))
                end do
            end do

            ! Check for convergence
            if (sum_off < tol * n * n) exit

            ! Threshold for this sweep
            if (sweep < 4) then
                threshold = 0.2d0 * sum_off / (n * n)
            else
                threshold = 0.0d0
            end if

            ! Sweep through all off-diagonal elements
            do p = 1, n-1
                do q = p+1, n
                    off_diag = abs(a(p,q))

                    ! Skip if element is small enough
                    if (sweep > 4 .and. off_diag < tol * abs(a(p,p)) .and. &
                        off_diag < tol * abs(a(q,q))) then
                        a(p,q) = 0.0d0
                        a(q,p) = 0.0d0
                        cycle
                    end if

                    if (off_diag > threshold) then
                        h = a(q,q) - a(p,p)

                        if (abs(h) < tol * off_diag) then
                            t = 1.0d0
                            if (a(p,q) < 0.0d0) t = -1.0d0
                        else
                            tau = h / (2.0d0 * a(p,q))
                            if (tau >= 0.0d0) then
                                t = 1.0d0 / (tau + sqrt(1.0d0 + tau*tau))
                            else
                                t = -1.0d0 / (-tau + sqrt(1.0d0 + tau*tau))
                            end if
                        end if

                        c = 1.0d0 / sqrt(1.0d0 + t*t)
                        s = t * c
                        tau = s / (1.0d0 + c)
                        h = t * a(p,q)

                        ! Update diagonal elements
                        a(p,p) = a(p,p) - h
                        a(q,q) = a(q,q) + h
                        a(p,q) = 0.0d0
                        a(q,p) = 0.0d0

                        ! Update rest of row/column p and q
                        do k = 1, p-1
                            g = a(k,p)
                            h = a(k,q)
                            a(k,p) = g - s * (h + g * tau)
                            a(k,q) = h + s * (g - h * tau)
                        end do
                        do k = p+1, q-1
                            g = a(p,k)
                            h = a(k,q)
                            a(p,k) = g - s * (h + g * tau)
                            a(k,q) = h + s * (g - h * tau)
                        end do
                        do k = q+1, n
                            g = a(p,k)
                            h = a(q,k)
                            a(p,k) = g - s * (h + g * tau)
                            a(q,k) = h + s * (g - h * tau)
                        end do

                        ! Update eigenvector matrix
                        do k = 1, n
                            g = v(k,p)
                            h = v(k,q)
                            v(k,p) = g - s * (h + g * tau)
                            v(k,q) = h + s * (g - h * tau)
                        end do
                    end if
                end do
            end do
        end do

        ! Check convergence
        if (sweep > max_sweeps) then
            info = 1
            write(*,'(A)') 'PolyChord WARNING (jacobi_eigen): Failed to converge'
        end if

        ! Extract eigenvalues from diagonal
        do i = 1, n
            w(i) = a(i,i)
        end do

        ! Sort eigenvalues in ascending order (simple bubble sort - n is small)
        do i = 1, n-1
            do j = i+1, n
                if (w(j) < w(i)) then
                    ! Swap eigenvalues
                    temp = w(i)
                    w(i) = w(j)
                    w(j) = temp
                    ! Swap eigenvector columns
                    do k = 1, n
                        temp = v(k,i)
                        v(k,i) = v(k,j)
                        v(k,j) = temp
                    end do
                end if
            end do
        end do

    end subroutine jacobi_eigen


    !> Compute eigenvalues and eigenvectors of a symmetric matrix
    !!
    !! Uses LAPACK DSYEV if available, otherwise falls back to Jacobi method.
    !!
    !! @param[in]   a_in  Symmetric matrix (not modified)
    !! @param[out]  w     Eigenvalues in ascending order
    !! @param[out]  v     Eigenvectors as columns
    !! @param[out]  info  0 on success, >0 on failure
    subroutine sym_eigen(a_in, w, v, info)
        implicit none
        real(dp), intent(in), dimension(:,:) :: a_in
        real(dp), intent(out), dimension(:) :: w
        real(dp), intent(out), dimension(:,:) :: v
        integer, intent(out) :: info

        integer :: n, lwork
        real(dp), allocatable :: work(:)
#ifdef HAVE_LAPACK
        real(dp), allocatable :: a_copy(:,:)
#endif

        n = size(a_in, 1)

#ifdef HAVE_LAPACK
        ! Use LAPACK DSYEV
        allocate(a_copy(n,n))
        a_copy = a_in

        ! Query optimal workspace size
        allocate(work(1))
        lwork = -1
        call dsyev('V', 'U', n, a_copy, n, w, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Compute eigenvalues and eigenvectors
        call dsyev('V', 'U', n, a_copy, n, w, work, lwork, info)

        ! DSYEV stores eigenvectors in a_copy
        v = a_copy

        deallocate(work, a_copy)
#else
        ! Use pure-Fortran Jacobi fallback
        ! jacobi_eigen destroys the input matrix, so we work on a copy
        ! and store eigenvectors separately
        allocate(work(n*n))
        work = reshape(a_in, [n*n])
        v = reshape(work, [n,n])  ! v is used as working copy for matrix
        deallocate(work)

        ! jacobi_eigen expects separate output for eigenvectors
        block
            real(dp), dimension(n,n) :: a_work, v_work
            a_work = a_in
            call jacobi_eigen(a_work, w, v_work, info)
            v = v_work
        end block
#endif

    end subroutine sym_eigen


    !> Regularize a covariance matrix to ensure positive-definiteness
    !!
    !! For rank-deficient covariance matrices (from clusters with n < D+1 points),
    !! this function replaces zero eigenvalues with 1.0 (unit variance in null-space).
    !!
    !! Algorithm:
    !! 1. Eigendecomposition: Σ = Q Λ Q^T
    !! 2. Replace eigenvalues below threshold with 1.0
    !! 3. Reconstruct: Σ_reg = Q Λ_reg Q^T
    !!
    !! @param[in]  covmat  Input covariance matrix (may be singular)
    !! @return     covmat_reg  Regularized positive-definite covariance matrix
    function regularize_covmat(covmat) result(covmat_reg)
        implicit none
        real(dp), intent(in), dimension(:,:) :: covmat
        real(dp), dimension(size(covmat,1), size(covmat,2)) :: covmat_reg

        integer :: n, i, j, info
        real(dp), allocatable :: w(:), v(:,:), w_reg(:)
        real(dp) :: max_eigenvalue, threshold

        n = size(covmat, 1)
        allocate(w(n), v(n,n), w_reg(n))

        ! Step 1: Eigendecomposition
        call sym_eigen(covmat, w, v, info)

        if (info /= 0) then
            write(*,'(A,I3)') 'PolyChord WARNING (regularize_covmat): Eigendecomposition failed with info=', info
            write(*,'(A)') '                   Returning identity matrix'
            covmat_reg = identity_matrix(n)
            deallocate(w, v, w_reg)
            return
        end if

        ! Check for NaN in eigenvalues
        if (any(ieee_is_nan(w))) then
            write(*,'(A)') 'PolyChord WARNING (regularize_covmat): NaN in eigenvalues, returning identity matrix'
            covmat_reg = identity_matrix(n)
            deallocate(w, v, w_reg)
            return
        end if

        ! Step 2: Handle negative eigenvalues (numerical error) and compute threshold
        w_reg = max(0.0d0, w)
        max_eigenvalue = maxval(w_reg)

        ! Hybrid threshold: relative + absolute floor
        threshold = max(max_eigenvalue * eigenvalue_rel_threshold, eigenvalue_abs_threshold)

        ! Step 3: Replace small/zero eigenvalues with 1.0 (unit variance in null-space)
        ! Also handle rank-0 case (all eigenvalues near zero)
        if (max_eigenvalue < eigenvalue_abs_threshold) then
            ! All eigenvalues effectively zero - return identity
            write(*,'(A)') 'PolyChord INFO (regularize_covmat): All eigenvalues near zero, using identity matrix'
            covmat_reg = identity_matrix(n)
            deallocate(w, v, w_reg)
            return
        end if

        do i = 1, n
            if (w_reg(i) < threshold) then
                w_reg(i) = 1.0d0
            end if
        end do

        ! Step 4: Reconstruct covariance matrix: Σ_reg = V * diag(w_reg) * V^T
        ! Use efficient matmul instead of nested loops for BLAS optimization
        block
            real(dp), dimension(n,n) :: temp_mat
            ! Calculate V * Lambda (scale columns of V by eigenvalues)
            do j = 1, n
                temp_mat(:,j) = v(:,j) * w_reg(j)
            end do
            ! Calculate (V * Lambda) * V^T
            covmat_reg = matmul(temp_mat, transpose(v))
        end block

        ! Sanity check
        if (any(ieee_is_nan(covmat_reg))) then
            write(*,'(A)') 'PolyChord WARNING (regularize_covmat): NaN in reconstructed matrix, returning identity'
            covmat_reg = identity_matrix(n)
        end if

        deallocate(w, v, w_reg)

    end function regularize_covmat

    function calc_covmat(x,wraparound) result(covmat)
        implicit none
        real(dp), intent(in), dimension(:,:) :: x
        logical, intent(in), dimension(size(x,1)) :: wraparound

        real(dp), dimension(size(x,1),size(x,1)) :: covmat

        real(dp), dimension(size(x,1),size(x,2)) :: dx
        real(dp), dimension(size(x,1)) :: mu, circle_mu

        integer :: nDims,n,i
        real(dp) :: sum_s, sum_c  ! For wraparound checks

        nDims = size(x,1)
        n = size(x,2)

        ! Handle insufficient points (n < 2)
        if (n < 2) then
            ! Return unit identity matrix to allow local exploration
            ! This is mathematically sound: represents zero correlation and unit variance
            covmat = identity_matrix(nDims)
            return
        endif

        ! === FIX 2: Wraparound atan2(0,0) protection ===
        ! Compute the circle mean with protection against atan2(0,0)
        circle_mu = 0d0
        if (any(wraparound)) then
            do i = 1, nDims
                if (wraparound(i)) then
                    sum_s = sum(sin(x(i,:)*TwoPi))
                    sum_c = sum(cos(x(i,:)*TwoPi))
                    ! Check if magnitude is effectively zero (all points identical or antipodal)
                    if (sqrt(sum_s**2 + sum_c**2) < 1.d-15) then
                        ! Mean is undefined, use first point's value
                        circle_mu(i) = x(i,1)
                    else
                        ! Normal case: compute angle
                        circle_mu(i) = atan2(sum_s, sum_c) / TwoPi
                    endif
                endif
            end do
        endif

        ! Compute the mean
        dx = x - spread(circle_mu,dim=2,ncopies=n)
        where(spread(wraparound,dim=2,ncopies=n)) dx = dx - nint(dx)
        mu = modulo(sum(dx,dim=2)/n + circle_mu, 1d0)

        ! Compute the covariance matrix
        dx = x - spread(mu,dim=2,ncopies=n)
        where(spread(wraparound,dim=2,ncopies=n)) dx = dx - nint(dx)

        covmat = matmul(dx,transpose(dx))/(n-1)

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

        real(dp), intent(in), dimension(:,:) :: data_array

        real(dp), dimension(size(data_array,2),size(data_array,2)) :: similarity_matrix

        integer :: i



        similarity_matrix = spread( [( dot_product(data_array(:,i),data_array(:,i)),i=1,size(data_array,2) )], dim=2,ncopies=size(data_array,2) )

        similarity_matrix = similarity_matrix + transpose(similarity_matrix) - 2d0 * matmul( transpose(data_array),data_array )

    end function calc_similarity_matrix




    !> This function relabels an array with more sensible indices
    !!
    !! We do this by constructing a mapping of each original label to a new label
    !!
    !! mapping :  new labels ---> old labels
    function relabel(array,num_labels) result(array_relabel)
        implicit none
        integer,intent(in),dimension(:)  :: array
        integer,intent(out)              :: num_labels

        integer,dimension(size(array)) :: array_relabel

        integer,dimension(size(array)) :: mapping

        integer :: npoints
        integer :: i_point
        integer :: i_label

        ! Find the number of points
        npoints = size(array)

        ! We will re-label the array type in array(1) with the integer 1
        mapping(1) = array(1)
        num_labels = 1

        do i_point=1,npoints
            ! If the array type for i_point is not already included in the
            ! array, then add it
            if( all(array(i_point)/=mapping(1:num_labels)) ) then
                num_labels=num_labels+1
                mapping(num_labels) = array(i_point)
            end if
        end do

        ! mapping now contains the random integers that are found in array

        ! We now relabel according to the inverse mapping
        do i_label=1,num_labels
            where(array==mapping(i_label)) array_relabel=i_label
        end do

    end function



    !> The volume of a unit n-sphere
    function Vn(nDims)
        implicit none
        integer,intent(in) :: nDims
        real(dp) :: Vn
        real(dp), parameter :: sqrtpi = sqrt(4d0*atan(1d0))
        Vn = sqrtpi**nDims /gamma(1d0+nDims/2d0)
    end function Vn






    function normal_cdf(x)
        implicit none
        real(dp), intent(in),dimension(:) :: x
        real(dp),dimension(size(x)) :: normal_cdf
        real(dp), parameter :: sqrt2 = sqrt(2d0)

        normal_cdf = 0.5d0 * ( 1d0 + erf(x/sqrt2))

    end function normal_cdf

    function inv_normal_cdf(x)
        implicit none
        real(dp), intent(in),dimension(:) :: x
        real(dp),dimension(size(x)) :: inv_normal_cdf
        integer :: i

        inv_normal_cdf = [( r8_normal_01_cdf_inverse(x(i)), i=1,size(x) )]

    end function inv_normal_cdf

    function convert_c_string(c_string) result(string)
        use iso_c_binding
        implicit none
        character(len=1,kind=c_char), intent(in), dimension(:) :: c_string
        character(len=size(c_string)) :: string

        integer :: i

        string = ' '
        do i=1,len(string)
            if( c_string(i) == c_null_char ) exit
            string(i:i) = c_string(i)
        end do

    end function convert_c_string

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

    !> Compute the loglikelihood of a multivariate gaussian
    function log_gauss(theta,mean,invcovmat,logdetcovmat)
        implicit none
        !> The input vector
        real(dp), intent(in), dimension(:) :: theta
        !> The mean
        real(dp), intent(in), dimension(:) :: mean
        !> The precomputed inverse covariance matrix
        real(dp), intent(in), dimension(:,:) :: invcovmat
        !> The precomputed logarithm of the determinant
        real(dp), intent(in) :: logdetcovmat


        ! The output
        real(dp) :: log_gauss

        ! Gaussian normalisation
        log_gauss = - ( size(theta) * logTwoPi + logdetcovmat )/2d0 

        log_gauss = log_gauss - dot_product(theta-mean,matmul(invcovmat,theta-mean))/2d0

    end function log_gauss

    !> This gets the wallclock timer from the mpi library
    function time() 
        implicit none
#ifdef MPI
        include 'mpif.h'
#endif
        real(dp) :: time

#ifdef MPI
        time = MPI_Wtime()
#else 
        call cpu_time(time)
#endif
      
    end function time


end module utils_module
