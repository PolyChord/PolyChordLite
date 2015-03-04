module loglikelihood_module


    double precision, parameter :: sigma = 0.1  ! width of peak

    double precision, allocatable, dimension(:,:)  :: invcovmat    ! inverse covariance matrix
    double precision, allocatable, dimension(:)    :: mu           ! Mean
    double precision                               :: logdetcovmat ! log(det(covariance matrix))

    contains
    !> Random Correlated gaussian loglikelihood
    !! 
    !! It is normalised so that it should output an evidence of 1.0 for
    !! effectively infinite priors.

    function loglikelihood(theta,phi)
        use utils_module, only: log_gauss
        implicit none
        !> Input parameters
        double precision, intent(in), dimension(:)   :: theta
        !> Output derived parameters
        double precision, intent(out),  dimension(:) :: phi

        double precision :: loglikelihood

        ! Compute log likelihood
        loglikelihood = log_gauss(theta,mu,invcovmat,logdetcovmat)

    end function loglikelihood



    subroutine setup_loglikelihood(settings,mpi_communicator)
#ifdef MPI
        use mpi
#endif
        use settings_module,   only: program_settings
        use random_module,     only: random_inverse_covmat
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

        integer :: nDims
        integer :: mpierror

        ! Get the dimensionality from settings
        nDims = settings%nDims

        ! Allocate the mean vector and inverse covariance matrix
        allocate(mu(nDims),invcovmat(nDims,nDims))

        ! Create a mean vector at the center of the space
        mu = 0.5d0

        ! Generate a random covariance matrix, its inverse and logdet on the root node
        call random_inverse_covmat(invcovmat,logdetcovmat,sigma,nDims)

#ifdef MPI
        ! Broadcast the covariance matrix and normalisation data to the
        ! rest of the nodes
        ! Covariance matrix:
        call MPI_BCAST(           &
            invcovmat,            & ! inverse covariance matrix data to be broadcast
            nDims*nDims,          & ! size of the data
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

    end subroutine setup_loglikelihood

end module loglikelihood_module
