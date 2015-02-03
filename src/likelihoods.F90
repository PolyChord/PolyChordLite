!> This file contains various likelihoods ready to be used in main.f90

module example_likelihoods
    use utils_module,    only: logzero,TwoPi,stdout_unit,Hypergeometric1F1,Hypergeometric2F1,Pochhammer 
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

    !> Twin gaussian peaks likelihood with mean mu(:) and an uncorrelated covariance sigma(:).
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.
    !!
    !! The mean is set at 0.5 by default, apart from the x directions, where
    !! they are separated by 6 sigma widths and all sigmas at 0.01

    function twin_gaussian_loglikelihood(theta,phi,context) result(loglikelihood)
        use utils_module, only: logaddexp
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context

        double precision :: loglikelihood 
        double precision :: loglikelihood1
        double precision :: loglikelihood2

        double precision, dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(size(theta)) :: mu1   ! Mean
        double precision, dimension(size(theta)) :: mu2   ! Mean


        ! Initialise the mean and standard deviation
        sigma = 1d-2  ! all sigma set relatively small
        mu1    = 5d-1 ! mean in the center
        mu1(1) = 5d-1 + 10*sigma(1)
        mu2    = 5d-1 ! mean in the center
        mu2(1) = 5d-1 - 10*sigma(1)

        ! Gaussian normalisation
        loglikelihood1 = - sum( log( sigma ) + log( TwoPi )/2d0 ) 
        loglikelihood2 = loglikelihood1

        ! theta dependence
        loglikelihood1 = loglikelihood1 - sum( ( ( theta - mu1 ) / sigma ) ** 2d0 ) / 2d0
        loglikelihood2 = loglikelihood2 - sum( ( ( theta - mu2 ) / sigma ) ** 2d0 ) / 2d0

        loglikelihood = logaddexp(loglikelihood1,loglikelihood2) - log(2d0)

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if


    end function twin_gaussian_loglikelihood



    !> Upside down [Rosenbrock function](http://en.wikipedia.org/wiki/Rosenbrock_function).
    !!
    !! \f[ -\log\mathcal{L}(\theta) = \sum_{i=1}^{N-1}   (a-\theta_i)^2+ b (\theta_{i+1} -\theta_i^2 )^2 \f]
    !!
    !! This is the industry standard 'Banana'. It is useful for testing whether the
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
    function himmelblau_loglikelihood(theta,phi,context) result(loglikelihood)
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
        use utils_module, only: TwoPi
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

        loglikelihood =  - sum( log(4991.21750d0) + theta**2 - A*cos(TwoPi*theta) )

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

    !> Model of the planck loglikelihood
    !!
    !! This reads in a covariance matrix from data/planck_covmat.dat and
    !! produces a correlated gaussian from it. This effectively models many of
    !! the difficulties associated with the planck loglikelihood
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function planck_loglikelihood(theta,phi,context)
        use random_module, only: random_reals
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                 :: context


        double precision :: planck_loglikelihood

        integer :: nDims

        double precision, allocatable, dimension(:,:), save :: invcovmat ! covariance matrix
        double precision, allocatable, dimension(:),   save :: mu    ! Mean
        double precision, save :: logdetcovmat

        logical,save :: initialised=.false.


        nDims = size(theta)

        if(.not. initialised) then
            allocate(invcovmat(nDims,nDims), &
                mu(nDims))

            ! create a rough mean vector for planck
            mu( 1) =  0.2207254E-01  !omegabh2    
            mu( 2) =  0.1196086E+00  !omegach2    
            mu( 3) =  0.1041323E+01  !theta       
            mu( 4) =  0.9667156E-01  !tau         
            mu( 5) =  0.3102562E+01  !logA         
            mu( 6) =  0.9616422E+00  !ns           
            mu( 7) =  0.1685974E+03  !aps100      
            mu( 8) =  0.5352338E+02  !aps143      
            mu( 9) =  0.1070724E+03  !aps217      
            mu(10) =  0.7815372E+01  !acib143     
            mu(11) =  0.2861451E+02  !acib217     
            mu(12) =  0.5130602E+01  !asz143      
            mu(13) =  0.8821973E+00  !psr         
            mu(14) =  0.4167502E+00  !cibr        
            mu(15) =  0.5319894E+00  !ncib        
            mu(16) =  0.1000581E+01  !cal0        
            mu(17) =  0.9964250E+00  !cal2        
            mu(18) =  0.4724667E+00  !xi          
            mu(19) =  0.4426342E+01  !aksz        
            mu(20) =  0.5281187E+00  !bm_1_1      
  
  
            ! Generate a planck covariance matrix from file and compute inverses
            ! and logdets 
            call read_covariance(invcovmat,logdetcovmat,nDims)

            initialised=.true.
        end if



        ! Compute log likelihood
        planck_loglikelihood = log_gauss(theta,mu,invcovmat,logdetcovmat)


        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if

    end function planck_loglikelihood




    !> Gaussian Needle loglikelihood











    !> Random Correlated gaussian loglikelihood
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

#ifdef MPI
        integer :: mpierror
#endif


        nDims = size(theta)

        if(.not. initialised) then
            allocate(invcovmat(nDims,nDims), &
                mu(nDims))

            ! Create a mean vector at the center of the space
            mu = 0.5d0
            ! Generate a random covariance matrix, its inverse and logdet on the root node
            call generate_covariance(invcovmat,logdetcovmat,sigma,nDims)
#ifdef MPI
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
        integer, parameter :: num_peaks = 3
        integer :: i !iterator
#ifdef MPI
        integer :: mpierror
#endif

        nDims=size(theta)

        if(.not. initialised) then
            allocate(invcovmat(nDims,nDims,num_peaks),&
                mu(nDims,num_peaks),               &
                logdetcovmat(num_peaks),             &
                log_likelihoods(num_peaks)           &
                )

            ! Generate num_peaks random mean vectors, localised around the center on the root node
            do i=1,num_peaks
                mu(:,i) = 0.5d0 + 10*sigma*(2d0*random_reals(nDims) -1d0)
            end do
            ! Generate num_peaks random covariance matrices, their inverses and logdets on the root node
            do i=1,num_peaks
                call generate_covariance(invcovmat(:,:,i),logdetcovmat(i),sigma,nDims)
            end do

#ifdef MPI
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
        integer,          intent(in)                         :: nDims
        double precision, intent(out),dimension(nDims,nDims) :: invcovmat
        double precision, intent(out)                        :: logdetcovmat
        double precision, intent(in)                         :: sigma

        double precision, dimension(nDims)       :: eigenvalues
        double precision, dimension(nDims,nDims) :: eigenvectors
        integer :: j
        double precision, parameter :: rng=1e-2

        ! Generate a random basis for the eigenvectors
        eigenvectors = random_orthonormal_basis(nDims)
        ! Generate the eigenvalues logarithmically in [rng,1] * sigma
        do j=1,nDims
            eigenvalues(j)  = sigma * rng**((j-1d0)/(nDims-1d0))
        end do
        ! Sort them lowest to highest for consistency
        !call dlasrt('D',nDims,eigenvalues,info)

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


    subroutine read_covariance(invcovmat,logdetcovmat,nDims)
        use random_module, only: random_reals, random_orthonormal_basis
        use utils_module, only: read_covmat_unit
        implicit none
        integer,          intent(in)                         :: nDims
        double precision, intent(out),dimension(nDims,nDims) :: invcovmat
        double precision, intent(out)                        :: logdetcovmat

        !double precision, dimension(nDims)       :: eigenvalues
        !double precision, dimension(nDims,nDims) :: eigenvectors
        integer :: j,k
        !double precision, parameter :: rng=1e-2
        double precision, dimension(nDims,nDims) :: covmat

        ! Read in the covariance matrix
        open(read_covmat_unit,file="data/planck_covmat.dat",action='read')
        read(read_covmat_unit,*) covmat
        close(read_covmat_unit)

        ! calculate the cholesky decomposition
        invcovmat=covmat
        !call dpotrf('U',nDims,invcovmat,nDims,info)

        ! Calculate the determinant
        logdetcovmat=0
        do j=1,nDims
            logdetcovmat = logdetcovmat+log(invcovmat(j,j))*2
        end do

        ! calculate the inverse function
        !call dpotri('U',nDims,invcovmat,nDims,info)
        ! symmetrise
        do j=1,nDims
            do k=j,nDims
                invcovmat(j,k) = invcovmat(j,k)
                invcovmat(k,j) = invcovmat(j,k)
            end do
        end do
    end subroutine read_covariance



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

    !> Pyramidal likelihood centered on 0.5.
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function pyramidal_loglikelihood(theta,phi,context) result(loglikelihood)
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)  :: theta
        !> Output derived parameters
        double precision, intent(out), dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                :: context
                            
        double precision :: loglikelihood

        double precision, dimension(size(theta)) :: sigma ! Standard deviation (uncorrelated) 
        double precision, dimension(size(theta)) :: center    ! Mean
        
        center= 5d-1   ! mean in the center
        sigma = 1d-3 

        ! normalisation
        loglikelihood =   - log(gamma(1d0+size(theta)/2d0)) -sum(log(2d0*sigma))

        ! theta dependence
        loglikelihood = loglikelihood - maxval(abs(theta-center)/sigma)**2

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if


    end function pyramidal_loglikelihood


    !> Gaussian shell defined as,
    !!
    !!  \f$ latex here \f$
    !!
    !! where the parameters used are,
    !!      - w = width of the guassian shell 
    !!      - R = Radius of the gaussian shell (ie: distance from the center)
    !!      - in progress...
    function gaussian_shell(theta,phi,context) result(loglikelihood)
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)  :: theta
        !> Output derived parameters
        double precision, intent(out), dimension(:) :: phi
        !> Pointer to any additional information
        integer,          intent(in)                :: context

        !The output of the function 
        double precision :: loglikelihood

        !Additional parameters that are needed
        double precision, parameter :: pi       = 4d0*atan(1d0) ! \pi in double precision
        double precision, parameter :: width    = 0.1           ! width of shell peak
        double precision, parameter :: radius   = 1d0           ! radius to the shell
        double precision            :: distance                 ! storing distance to center of theta
        double precision, dimension(size(theta))    :: center   ! center of the gaussian shell
        integer                     :: dim                      ! dimensions of theta/the space

        !For reducing the computational load by storing the normalisation:
        logical, save               :: stored   = .false.
        double precision, save      :: normalisation
        
        !Getting the dimension and storing in dim
        dim      = size(theta)
        !Define where the center of the gaussian shell is
        center = 0d0
        !Calculate the distance between the parameter space points theta and center
        distance = sqrt( sum( (center-theta)**2 ) )

        !Computing the normalisation upon first call to this function:
        if(.NOT. stored) then
            !For the first call to this function, we calculate the normalisation
            !and change the flag to say this is completed:
            !   choose from the options below 
            !   (Mathematic normalisation is most exact)

            !No normalisation:
            !normalisation = 0

            !Mathematica normalisation:
            !normalisation = - (width**dim)*( ((2*pi)**(dim/2d0)) * Hypergeometric1F1( (1d0-dim)/2d0,1d0/2d0,-(radius*radius)/(2d0*width*width) ) + ( (2**((1d0+dim)/2d0)) * (pi**(dim/2d0)) * radius * GAMMA( (1d0+dim)/2d0 ) * Hypergeometric1F1( 1d0-(dim/2d0),3d0/2d0,-(radius*radius)/(2d0*width*width) ))/(width * GAMMA( dim/2d0 )) )

            !Temporary normalisation (when you know the numerical value):
            !normalisation = -0.45

            !Will's approximate normalisation based on assumptions,
            !   1. width*dimension << radius 
            !Basically, over the region of interest the integral in polar coords (a gaussian and the nD volume element) 
            !simplify nicely and allow the limits of the integral to go to -inf. 
            normalisation = log_gamma(1d0+(dim/2d0)) - (1d0/2d0)*log(2d0) - log(dim+0d0) - ((1d0+dim)/2d0)*log(pi) - (dim-1d0)*log(radius) - log(width)

            !Change the flag:
            stored = .true.
        end if

        !Use the stored value as the normalisation:
        loglikelihood = normalisation

        !Calculating the loglikelihood dependant on theta here:
        loglikelihood = loglikelihood - ( (distance-radius)**2d0 /(2d0*width*width) )

        ! Use up these parameters to stop irritating warnings
        if(size(phi)>0) then
            phi= context
            phi=0d0
        end if


    end function gaussian_shell

end module example_likelihoods
