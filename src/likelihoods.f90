module example_likelihoods
    use model_module , only: model,logzero

    contains

    !> Basic Gaussian likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, and all sigmas at 0.01

    function gaussian_loglikelihood(M,theta,feedback)
        implicit none
        class(model),     intent(in)               :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: gaussian_loglikelihood

        double precision, dimension(M%nDims) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(M%nDims) :: mu    ! Mean

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision

        
        ! Initialise the mean and standard deviation
        mu    = 5d-1   ! mean in the center
        sigma = 1d-2  ! all sigma set relatively small

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : Gaussian" )')
            end if
            if(feedback>=2) then
                write(*,'( "     mean: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') mu
                write(*,'( "     sigma: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') sigma
            end if
            gaussian_loglikelihood = logzero
            return
        end if


        ! Gaussian normalisation
        gaussian_loglikelihood = - M%nDims/2d0 * log( TwoPi ) - sum( log( sigma ) ) 

        ! theta dependence
        gaussian_loglikelihood = gaussian_loglikelihood - sum( ( ( theta - mu ) / sigma ) ** 2d0 ) / 2d0


    end function gaussian_loglikelihood

    !> Correlated gaussian loglikelihood
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function gaussian_loglikelihood_corr(M,theta,feedback)
        use random_module, only: random_real
        implicit none
        class(model),     intent(in)               :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: gaussian_loglikelihood_corr

        double precision, allocatable, dimension(:,:), save :: invcovmat ! covariance matrix
        double precision, allocatable, dimension(:,:), save :: covmat ! Standard deviation (uncorrelated) 
        double precision, allocatable, dimension(:),   save :: mu    ! Mean
        double precision, allocatable, dimension(:),   save :: eigenvalues
        double precision, save :: logdetcovmat

        logical,save :: initialised=.false.

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision

        integer :: info
        integer :: i

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : Correlated Gaussian" )')
            end if
            if(feedback>=2) then
                write(*,'( "     mean: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') mu
                write(*,'( "     sigma: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') covmat
            end if
            gaussian_loglikelihood_corr = logzero
            return
        end if



        if(.not. initialised) then
            allocate(invcovmat(M%nDims,M%nDims),covmat(M%nDims,M%nDims),mu(M%nDims),eigenvalues(max(3*M%nDims-1,1)))
            ! Generate a random mean vector
            mu = random_real(M%nDims)

            ! generate a 'random' positive semi-definite matrix
            !    start by generating a random vector of length ndims*ndims and
            !    reshaping it into a square matrix
            covmat = reshape(random_real(M%nDims*M%nDims),(/M%nDims,M%nDims/))
            ! symmetrise this to make it positive semi-definate
            covmat = matmul(transpose(covmat),covmat)

            covmat = 0d0
            do i=1,M%nDims
                covmat(i,i) = 0.01
            end do

            write(*,'("/------ covmat -------------\")')
            write(*,'("|", <M%nDims>F9.5,"|")') covmat
            write(*,'("\---------------------------/")')

            ! Invert the symmetric matrix
            invcovmat = covmat      !initialise it with the covmat
            !    Perform scalar decomposition
            !    (d=double,po=symmetric positive definite,trf=triangular matrix factorization)
            !   'U' = compute upper half of matrix
            call potrf(invcovmat, 'U', info)
            write(*,'("/------ invcovmat ----------\")')
            write(*,'("|", <M%nDims>F9.5,"|")') invcovmat
            write(*,'("\---------------------------/")')
            !    Invert matrix
            !    (d=double,po=symmetric positive definite,trf=inverse matrix using the factorization)
            !   'U' = compute upper half of matrix
            call potri(invcovmat,'U', info)
            write(*,'("/------ invcovmat ----------\")')
            write(*,'("|", <M%nDims>F9.5,"|")') invcovmat
            write(*,'("\---------------------------/")')

            ! Compute the eigenvalues with ssyev
            !   (d=double, sy=symmetric, ev=eigenvalues)
            !   'N' = don't compute eigenvectors
            !   'U' = compute upper half of matrix
            call syev(covmat,eigenvalues,'N','U',info)
            write(*,*) '------ eigenvalues --------'
            write(*,'("|", <M%nDims>F9.5,"|")') eigenvalues
            write(*,*) '---------------------------'
            ! sum up the logs of the eigenvalues to get the determinant
            logdetcovmat = sum(log(eigenvalues))






            initialised=.true.
        end if
        
        ! Compute log likelihood
        gaussian_loglikelihood_corr = log_gauss(theta,mu,invcovmat,logdetcovmat)


    end function gaussian_loglikelihood_corr



    !> Compute the loglikelihood of a multivariate gaussian
    function log_gauss(theta,mean,invcovmat,logdetcovmat)
        implicit none
        !> The input vector
        double precision, intent(in), dimension(:) :: theta
        !> The mean
        double precision, intent(in), dimension(:) :: mean
        !> The precomputed inverse covariance matrix
        double precision, intent(in), dimension(:,:) :: invcovmat
        !> The precomputed logarithm of the determinant
        double precision, intent(in) :: logdetcovmat

        ! temp variable
        double precision, dimension(size(theta)) :: Etheta

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision

        ! The output
        double precision :: log_gauss

        ! Gaussian normalisation
        log_gauss = - size(theta)/2d0 * ( log( TwoPi ) + logdetcovmat )

        ! note the use of blas routine dsymv to avoid having to symmetrise the matrix
        !   (d=double,sy=symmetric,mv=matrix-vector product)
        call dsymv(invcovmat,theta-mean,Etheta,'U')
        log_gauss = log_gauss - dot_product(theta-mean,Etheta)/2d0

    end function log_gauss





    !> Pyramidal likelihood centered on 0.5.
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function pyramidal_loglikelihood(M,theta,feedback)
        implicit none
        class(model),intent(in)                    :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: pyramidal_loglikelihood

        double precision, dimension(M%nDims) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(M%nDims) :: center    ! Mean
        
        center    = 5d-1   ! mean in the center
        sigma = 1d-2 

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : Pyramidal" )')
            end if
            if(feedback>=2) then
                write(*,'( "     center: ")')
                write(*,'( " [", <M%nDims>F15.9 ,"]")') center
            end if
            return
        end if



        ! normalisation
        pyramidal_loglikelihood =   -(M%nDims)*log(2d0) - log(gamma(1d0+M%nDims/2d0)) -sum(log(sigma))

        ! theta dependence
        pyramidal_loglikelihood = pyramidal_loglikelihood - maxval(abs(theta-center)/sigma)**2

    end function pyramidal_loglikelihood



end module example_likelihoods
