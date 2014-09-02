module example_likelihoods
    use model_module,    only: model
    use utils_module,    only: logzero,TwoPi,stdout_unit
#ifdef MPI
    use mpi_module
#endif

    contains

    !> Basic Gaussian likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, and all sigmas at 0.01

    function gaussian_loglikelihood(theta,phi,context) result(loglikelihood)
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        double precision :: loglikelihood

        double precision, dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(size(theta)) :: mu    ! Mean

        
        ! Initialise the mean and standard deviation
        mu    = 5d-1   ! mean in the center
        sigma = 1d-2  ! all sigma set relatively small

        ! Gaussian normalisation
        loglikelihood = - sum( log( sigma ) + log( TwoPi )/2d0 ) 

        ! theta dependence
        loglikelihood = loglikelihood - sum( ( ( theta - mu ) / sigma ) ** 2d0 ) / 2d0

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if


    end function gaussian_loglikelihood



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
    function rosenbrock_loglikelihood(theta,phi,context) result(loglikelihood)
        use utils_module, only: logzero
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        ! The return value
        double precision :: loglikelihood

        ! Parameters of the rosenbrock function
        double precision, parameter :: a = 1d0
        double precision, parameter :: b = 1d2
        double precision, parameter :: pi = 4d0*atan(1d0) ! \pi in double precision

        ! number of dimensions
        integer :: nDims

        ! Get the number of dimensions
        nDims = size(theta)

        ! Normalisation for 2D
        loglikelihood = -log(pi/sqrt(b))

        ! Sum expressed with fortran intrinsics
        loglikelihood =  loglikelihood  - sum( (a-theta(1:nDims-1))**2 + b*(theta(2:nDims) - theta(1:nDims-1)**2)**2 )

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

    end function rosenbrock_loglikelihood



    !> Upside down [Himmelblau function](http://en.wikipedia.org/wiki/Himmelblau's_function).
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
    !! Note that this is a 2D function, and should hence only be used for M%nDims=2
    !!
    function himmelblau_loglikelihood(theta,phi,context) result(loglikelihood)
        use utils_module, only: logzero
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        ! The return value
        double precision :: loglikelihood

        ! Normalisation
        loglikelihood = -log(0.4071069421432255d0) 

        ! Evaluate
        loglikelihood =  loglikelihood  -  (theta(1)**2 + theta(2)- 11d0)**2    -   (theta(1)+theta(2)**2-7d0)**2 

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

    end function himmelblau_loglikelihood


    !> Upside down [Rastrigin function](http://en.wikipedia.org/wiki/Rastrigin_function).
    !!
    !! \f[ -\log L(\theta) = f(\mathbf{x}) = A n + \sum_{i=1}^n \left[\theta_i^2 - A\cos(2 \pi \theta_i)\right] \f]
    !!
    !! where  \f$ A=10\f$  and \f$ \theta_i\in[-5.12,5.12] \f$ . It has a global minimum 
    !! at \f$ \theta = \mathbf{0}\f$  where \f$ f(\theta)=0\f$ .
    !!
    !! In mathematical optimization, the Rastrigin function is a
    !! non-convex function used as a performance test problem for optimization
    !! algorithms. It is a typical example of non-linear multimodal function. It
    !! was first proposed by Rastrigin as a 2-dimensional function and has been generalized
    !! by Mühlenbein et al. Finding the minimum of this function is a
    !! fairly difficult problem due to its large search space and its large number
    !! of local minima.
    !!
    function rastrigin_loglikelihood(theta,phi,context) result(loglikelihood)
        use utils_module, only: logzero,TwoPi
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        double precision, parameter :: A=10d0

        ! The return value
        double precision :: loglikelihood

        loglikelihood =  loglikelihood  -  sum( log(4991.217507308888d0) + theta**2 - A*cos(TwoPi*theta) )

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

    end function rastrigin_loglikelihood




    !> Eggbox likelihood
    !!
    !! \f[ -\log L(x, y) = (2+ \prod_i^n \cos(\theta_i/2) )^5 \f]
    !!
    function eggbox_loglikelihood(theta,phi,context) result(loglikelihood)
        use utils_module, only: logzero
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        ! The return value
        double precision :: loglikelihood

        ! No normalisation implemented yet
        loglikelihood = 0

        ! calculate
        loglikelihood =  loglikelihood  -  (2 + product(cos(theta / 2d0) ) )**5

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

    end function eggbox_loglikelihood















    !> Correlated gaussian loglikelihood
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function gaussian_loglikelihood_corr(theta,phi,context)
        use random_module, only: random_reals
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context


        double precision :: gaussian_loglikelihood_corr

        integer :: nDims

        double precision, allocatable, dimension(:,:), save :: invcovmat ! covariance matrix
        double precision, allocatable, dimension(:),   save :: mu    ! Mean
        double precision, save :: logdetcovmat

        logical,save :: initialised=.false.

        double precision, parameter :: sigma = 0.01 ! width of peak


        nDims = size(theta)

        if(.not. initialised) then
            allocate(invcovmat(nDims,nDims), &
                mu(nDims))

            ! Create a mean vector at the center of the space
            mu = 0.5d0
#ifdef MPI
            ! Generate a random covariance matrix, its inverse and logdet on the root node
            if(mpi_rank()==0) call generate_covariance(invcovmat,logdetcovmat,sigma,nDims)
            ! Broadcast the covariance matrix and normalisation data to the
            ! rest of the nodes
            ! Covariance matrix:
            call MPI_BCAST(           &
                invcovmat,            & ! inverse covariance matrix data to be broadcast
                nDims*nDims,      & ! size of the data
                MPI_DOUBLE_PRECISION, & ! type of data
                0,                    & ! root node id
                MPI_COMM_WORLD,       & ! communication info
                mpierror)               ! error (from mpiutils)
            ! normalisation logdet covmat
            call MPI_BCAST(           &
                logdetcovmat,         & ! log(determinant of inverse covariance matrix) data to be broadcast
                1,                    & ! size of the data
                MPI_DOUBLE_PRECISION, & ! type of data
                0,                    & ! root node id
                MPI_COMM_WORLD,       & ! communication info
                mpierror)               ! error (from mpiutils)
#else
            ! Generate a random covariance matrix, its inverse and logdet
            call generate_covariance(invcovmat,logdetcovmat,sigma,nDims)
#endif

            initialised=.true.
        end if


        
        ! Compute log likelihood
        gaussian_loglikelihood_corr = log_gauss(theta,mu,invcovmat,logdetcovmat)


        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

    end function gaussian_loglikelihood_corr

    !> Cluster of correlated gaussians
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function gaussian_loglikelihood_cluster(theta,phi,context)
        use random_module, only: random_reals
        use utils_module,  only: logsumexp
#ifdef MPI
        use mpi_module
#endif
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context


        double precision :: gaussian_loglikelihood_cluster

        integer :: nDims

        double precision, allocatable, dimension(:,:,:), save :: invcovmat ! list of covariance matrices
        double precision, allocatable, dimension(:,:),   save :: mu    ! list of means
        double precision, allocatable, dimension(:),     save :: logdetcovmat !list of log(determinants)

        double precision, allocatable, dimension(:),     save :: log_likelihoods

        logical,save :: initialised=.false.

        double precision, parameter :: sigma = 0.01 ! width of peak
        integer, parameter :: num_peaks = 10
        integer :: i !iterator

        nDims=size(theta)

        if(.not. initialised) then
            allocate(invcovmat(nDims,nDims,num_peaks),&
                mu(nDims,num_peaks),               &
                logdetcovmat(num_peaks),             &
                log_likelihoods(num_peaks)           &
                )

#ifdef MPI

            if(mpi_rank()==0) then
                ! Generate num_peaks random mean vectors, localised around the center on the root node
                do i=1,num_peaks
                    mu(:,i) = 0.5d0 + 10*sigma*(2d0*random_reals(nDims) -1d0)
                end do
                ! Generate num_peaks random covariance matrices, their inverses and logdets on the root node
                do i=1,num_peaks
                    call generate_covariance(invcovmat(:,:,i),logdetcovmat(i),sigma,nDims)
                end do
            end if
            ! Broadcast the means, covariances and normalisation data to the
            ! rest of the nodes
            ! Means:
            call MPI_BCAST(                &
                mu,                        & ! inverse covariance matrix data to be broadcast
                nDims*num_peaks,         & ! size of the data
                MPI_DOUBLE_PRECISION,      & ! type of data
                0,                         & ! root node id
                MPI_COMM_WORLD,            & ! communication info
                mpierror)                    ! error (from mpiutils)

            ! Inverse Covariance matrices
            call MPI_BCAST(                &
                invcovmat,                 & ! inverse covariance matrix data to be broadcast
                nDims*nDims*num_peaks, & ! size of the data
                MPI_DOUBLE_PRECISION,      & ! type of data
                0,                         & ! root node id
                MPI_COMM_WORLD,            & ! communication info
                mpierror)                    ! error (from mpiutils)

            ! normalisation logdetcovmat
            call MPI_BCAST(                &
                logdetcovmat,              & ! log(determinant of inverse covariance matrix) data to be broadcast
                num_peaks,                 & ! size of the data
                MPI_DOUBLE_PRECISION,      & ! type of data
                0,                         & ! root node id
                MPI_COMM_WORLD,            & ! communication info
                mpierror)                    ! error (from mpiutils)
#else
            ! Generate num_peaks random mean vectors, localised around the center 
            do i=1,num_peaks
                mu(:,i) = 0.5d0 + 10*sigma*(2d0*random_reals(nDims) -1d0)
            end do

            ! Generate num_peaks random covariance matrices, their inverses and logdets
            do i=1,num_peaks
                call generate_covariance(invcovmat(:,:,i),logdetcovmat(i),sigma,nDims)
            end do
#endif

            initialised=.true.
        end if


        ! Compute log likelihood
        gaussian_loglikelihood_cluster = logzero
        do i =1,num_peaks
            log_likelihoods(i) = log_gauss(theta(:),mu(:,i),invcovmat(:,:,i),logdetcovmat(i))
        end do
        gaussian_loglikelihood_cluster = logsumexp(log_likelihoods)
        gaussian_loglikelihood_cluster = gaussian_loglikelihood_cluster - log(num_peaks+0d0)


        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

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
        double precision, parameter :: rng=1e-1

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


        ! The output
        double precision :: log_gauss

        ! Gaussian normalisation
        log_gauss = - ( size(theta) * log( TwoPi ) + logdetcovmat )/2d0 

        log_gauss = log_gauss - dot_product(theta-mean,matmul(invcovmat,theta-mean))/2d0

    end function log_gauss


end module example_likelihoods
