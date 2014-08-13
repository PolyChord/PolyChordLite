module example_likelihoods
    use model_module,    only: model
    use utils_module,    only: logzero

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
        use random_module, only: random_reals, random_skewed_direction
        implicit none
        class(model),     intent(in)               :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: gaussian_loglikelihood_corr

        double precision, allocatable, dimension(:,:), save :: invcovmat ! covariance matrix
        double precision, allocatable, dimension(:),   save :: mu    ! Mean
        double precision, save :: logdetcovmat

        logical,save :: initialised=.false.

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision
        double precision, parameter :: sigma = 0.01 ! width of peak

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : Correlated Gaussian" )')
                write(*,'( "  sigma     = ",E15.7 )') sigma
            end if
            gaussian_loglikelihood_corr = logzero

            if(.not. initialised) then
                allocate(invcovmat(M%nDims,M%nDims), &
                         mu(M%nDims))

                ! Generate a random mean vector, localised around the center 
                mu = 0.5d0
                ! Generate a random covariance vector and its inverse and logdet
                call generate_covariance(invcovmat,logdetcovmat,sigma,M%nDims)

                initialised=.true.
            end if

            return
        end if


        
        ! Compute log likelihood
        gaussian_loglikelihood_corr = log_gauss(theta,mu,invcovmat,logdetcovmat)


    end function gaussian_loglikelihood_corr

    !> Cluster of correlated gaussians
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function gaussian_loglikelihood_cluster(M,theta,feedback)
        use random_module, only: random_reals
        use utils_module,  only: logsumexp
        implicit none
        class(model),     intent(in)               :: M
        double precision, intent(in), dimension(:) :: theta
        integer,optional, intent(in)               :: feedback


        double precision :: gaussian_loglikelihood_cluster

        double precision, allocatable, dimension(:,:,:), save :: invcovmat ! list of covariance matrices
        double precision, allocatable, dimension(:,:),   save :: mu    ! list of means
        double precision, allocatable, dimension(:),     save :: logdetcovmat !list of log(determinants)

        double precision, allocatable, dimension(:),     save :: log_likelihoods

        logical,save :: initialised=.false.

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision
        double precision, parameter :: sigma = 0.01 ! width of peak
        integer, parameter :: num_peaks = 10
        integer :: i !iterator


        ! Feedback if requested
        if(present(feedback)) then

            if(feedback>=0) then
                write(*,'( "Likelihood : Clustered Gaussian" )')
                write(*,'( "  num_peaks = ",I4 )') num_peaks
                write(*,'( "  sigma     = ",E15.7 )') sigma
            end if

            if(.not. initialised) then
                allocate(invcovmat(M%nDims,M%nDims,num_peaks),&
                         mu(M%nDims,num_peaks),               &
                         logdetcovmat(num_peaks),             &
                         log_likelihoods(num_peaks)           &
                     )

                ! Generate num_peaks random mean vectors, localised around the center 
                do i=1,num_peaks
                    mu(:,i) = 0.5d0 + 10*sigma*(2d0*random_reals(M%nDims) -1d0)
                end do

                ! Set them all on the center to create a 'boo-bah' structure
                mu(:,:) = 0.5d0



                ! Generate a i random covariance matricesi, their inverses and logdets
                do i=1,num_peaks
                    call generate_covariance(invcovmat(:,:,i),logdetcovmat(i),sigma,M%nDims)
                end do
                write(*,'(<M%nDims>E12.5)') invcovmat(:,:,1)
                initialised=.true.
            end if

            gaussian_loglikelihood_cluster = logzero
            return
        end if


        ! Compute log likelihood
        gaussian_loglikelihood_cluster = logzero
        do i =1,num_peaks
            log_likelihoods(i) = log_gauss(theta(:),mu(:,i),invcovmat(:,:,i),logdetcovmat(i))
        end do
        gaussian_loglikelihood_cluster = logsumexp(log_likelihoods)
        gaussian_loglikelihood_cluster = gaussian_loglikelihood_cluster - log(num_peaks+0d0)


    end function gaussian_loglikelihood_cluster





    subroutine generate_covariance(invcovmat,logdetcovmat,sigma,nDims)
        use random_module, only: random_reals, random_orthonormal_basis
        implicit none
        double precision, intent(out),dimension(nDims,nDims) :: invcovmat
        double precision, intent(out)                        :: logdetcovmat
        double precision, intent(in)                         :: sigma
        integer,          intent(in)                         :: nDims

        double precision, dimension(nDims)       :: eigenvalues
        double precision, dimension(nDims,nDims) :: eigenvectors
        integer :: j
        double precision, parameter :: rng=5e-1

        ! Generate a random basis for the eigenvectors
        eigenvectors = random_orthonormal_basis(nDims)
        ! Generate the eigenvalues logarithmically in [rng,1] * sigma
        eigenvalues  = sigma *( rng**random_reals(nDims))

        ! Create the inverse covariance matrix in the eigenbasis
        invcovmat = 0d0
        do j=1,nDims
            invcovmat(j,j) = 1d0/eigenvalues(j)**2
        end do

        ! Rotate the matrix into the coordinate basis
        invcovmat = matmul(eigenvectors,matmul(invcovmat,transpose(eigenvectors)))

        ! sum up the logs of the eigenvalues to get the log of the determinant
        logdetcovmat = 2 * sum(log(eigenvalues))

    end subroutine generate_covariance


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

        double precision, parameter :: TwoPi = 8d0*atan(1d0) ! 2\pi in double precision

        ! The output
        double precision :: log_gauss

        ! Gaussian normalisation
        log_gauss = - ( size(theta) * log( TwoPi ) + logdetcovmat )/2d0 

        log_gauss = log_gauss - dot_product(theta-mean,matmul(invcovmat,theta-mean))/2d0

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

    !> Upside down [Rosenbrock function](http://en.wikipedia.org/wiki/Rosenbrock_function).
    !!
    !! \f[ -\log\mathcal{L}(\theta) = \sum_{i=1}^{N-1}   (a-\theta_i)^2+ b (\theta_{i+1} -\theta_i^2 )^2 \f]
    !!
    !! This is the industry standard 'Banana'. It useful for testing whether the
    !! algorithm is capable of navigating a curving degenerate space and able to
    !! find the global maximum. As is conventional, we choose \f$a=1\f$ and
    !! \f$b=100\f$
    !!
    !! It is only valid in nDims>2.
    !!
    !! To be precise, we choose our log likelihood as the negative of the 'true'
    !! Rosenbrock function, as our algorithm is a maximiser.
    !!
    !! In addition, we add an offset to normalise the loglikelihood so that in
    !! 2D it should integrate to 1. (There is no analytic formula for ND).
    !! 
    !! The global maximum is atop a long, narrow, parabolic shaped ridge.
    !!
    !! dimension | extrema
    !! ----------|------
    !! \f$ 2 \f$ | One maximum at \f$(1,1)\f$
    !! \f$ 3 \f$ | One maximum at \f$(1,1,1)\f$ 
    !! \f$4-7\f$ | One global maximum at \f$(1,\ldots,1)\f$, one local near \f$(-1,1,\ldots,1)\f$
    !!
    !! \f$(1,\ldots,1)\f$ is always the
    !! global maximum, but it is unclear analytically how many (if any) local maxima there
    !! are elsewhere in dimensions higher than 7.
    !!
    function rosenbrock_loglikelihood(M,theta,feedback) result(loglikelihood)
        implicit none
        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        class(model),intent(in)                    :: M
        !> The input parameters
        double precision, intent(in), dimension(:) :: theta
        !> Optional argument. If provided, then the function simply prints out a few details
        integer,optional, intent(in)               :: feedback

        double precision, parameter :: a = 1d0
        double precision, parameter :: b = 1d2
        double precision, parameter :: pi = 4d0*atan(1d0) ! \pi in double precision

        ! The return value
        double precision :: loglikelihood

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Likelihood : RosenBrock" )')
            end if
            return
        end if

        ! Normalisation for 2D
        loglikelihood = -log(pi/sqrt(b))
        
        ! Sum expressed with fortran intrinsics
        loglikelihood =  loglikelihood  - sum( (a-theta(1:M%nDims-1))**2 + b*(theta(2:M%nDims) - theta(1:M%nDims-1)**2)**2 )

        
    end function rosenbrock_loglikelihood




end module example_likelihoods
