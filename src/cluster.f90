module cluster_module
    use model_module, only: model
    implicit none

    !> Type to contain all of the cluster information, including a pointer to
    !! the clustering algorithm
    type :: cluster_info

        !> The maximum number of clusters
        integer :: ncluster_max = 1

        !> The number of detected clusters
        integer :: ncluster = 1

        !> The number of live points
        double precision :: nlive

        !> The dimensionality
        integer :: nDims

        !> The set of means of the clusters
        double precision, allocatable, dimension(:,:) :: mean

        !> The covariance matrices
        double precision, allocatable, dimension(:,:,:) :: covmat
        double precision, allocatable, dimension(:,:,:) :: invcovmat

        !> The cholesky representation of the inverse covariance matrices
        double precision, allocatable, dimension(:,:,:) :: cholesky


        !> Pointer to the cluster detection procedure.
        !!
        !! e.g: DENCLUE, DBSCAN, X-MEANS 
        !! 
        procedure(clus), pass(cluster),       pointer :: detect_clusters

    end type cluster_info

    interface
        !> Interface to the sampling procedure
        !!
        !! The sampling procedure takes the details of the model (M) and the current
        !! set of live points (live_data) in order to generate a new_point
        !! uniformly sampled from within the loglikelihood contour specifed by
        !! loglikelihood_bound
        subroutine clus(cluster,M,baby_point,late_point) 

            import :: model
            import :: cluster_info
            implicit none

            ! ------- Inputs -------
            !> cluster details to be modified
            class(cluster_info), intent(inout) :: cluster

            !> Model details
            type(model), intent(in) :: M
            !> The new-born baby point
            double precision,    dimension(M%nTotal)   :: baby_point
            !> The recently dead point
            double precision,    dimension(M%nTotal)   :: late_point
        end subroutine clus
    end interface


    contains

    !> Update the stored value of the covariance matrix
    !!
    !! Note that we make extensive use of the intel function 
    !! [dsyr](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/hh_goto.htm#GUID-38EB4174-468F-4045-B2E0-A38C4FB861D0.htm)(n,alpha,x,1,A,n),
    !! which computes:
    !! \f[ A + alpha * x x^T \f]
    !! where A is an n by n symmetric matrix, x is a n-vector alpha is a scalar
    !!
    !! We also use [dpotrf](https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/GUID-15A9E907-2600-430B-BF11-DFA0A375711B.htm#GUID-15A9E907-2600-430B-BF11-DFA0A375711B)
    !! which computes the upper triangular cholesky factor. This is required to
    !! pass to the random direction generator
    subroutine update_global_covariance(cluster,M,baby_point,late_point) 
        use utils_module, only: logzero
        use model_module, only: model
        implicit none

        ! ------- Inputs -------
        !> cluster details to be modified
        class(cluster_info), intent(inout) :: cluster

        !> Model details
        type(model), intent(in) :: M
        !> The new-born baby point
        double precision,    dimension(M%nTotal)   :: baby_point
        !> The recently dead point
        double precision,    dimension(M%nTotal)   :: late_point

        integer :: info

        integer :: nDims

        nDims = cluster%nDims


        ! Initially
        ! covmat = <x_i x_j> - <x_i><x_j>


        ! Remove the old mean from the covariance matrix to start with
        ! covmat  -> covmat + <x_i><x_j>
        ! covmat = <x_i x_j>
        call dsyr('U',nDims,+1d0,cluster%mean(:,1),1,cluster%covmat(:,:,1),nDims) 


        ! Add the new point to the mean
        ! mean = <x_i>
        cluster%mean(:,1) = cluster%mean(:,1) + baby_point(M%h0:M%h1) / cluster%nlive
        ! Remove the old point from mean
        cluster%mean(:,1) = cluster%mean(:,1) - late_point(M%h0:M%h1) / cluster%nlive

        ! Add the new point to the sum of the squares (now stored in covmat)
        ! covmat = <x_i x_j>
        call dsyr('U',nDims,1d0/cluster%nlive,baby_point(M%h0:M%h1),1,cluster%covmat(:,:,1),nDims)
        ! Remove the old point
        call dsyr('U',nDims,-1d0/cluster%nlive,late_point(M%h0:M%h1),1,cluster%covmat(:,:,1),nDims)

        ! add the new mean back in
        ! covmat -> covmat - <x_i><x_j>
        ! covmat = <x_i x_j> - <x_i><x_j>
        call dsyr('U',nDims,-1d0,cluster%mean(:,1),1,cluster%covmat(:,:,1),nDims) 


        ! calculate the cholesky factorisation of the covariance matrix
        cluster%cholesky = cluster%covmat
        call dpotrf('U',nDims,cluster%cholesky(:,:,1),nDims,info)

        ! calculate the inverse covariance matrix
        cluster%invcovmat = cluster%cholesky
        call dpotri('U',nDims,cluster%invcovmat(:,:,1),nDims,info)


    end subroutine update_global_covariance


    !> Do all of the necessary initialisation and allocation
    subroutine initialise_cluster(cluster)
        implicit none
        type(cluster_info)     :: cluster


        
        ! allocate the arrays
        allocate(&
            cluster%cholesky(cluster%nDims,cluster%nDims,cluster%ncluster_max), &
            cluster%covmat(cluster%nDims,cluster%nDims,cluster%ncluster_max), &
            cluster%invcovmat(cluster%nDims,cluster%nDims,cluster%ncluster_max), &
            cluster%mean(cluster%nDims,cluster%ncluster_max)&
            )
        cluster%cholesky = 0
        cluster%covmat = 0
        cluster%mean = 0

    end subroutine initialise_cluster


    function logsumexp(a,b)
        implicit none
        double precision :: a
        double precision :: b
        double precision :: logsumexp

        logsumexp= max(a,b) + log(exp(a-max(a,b))+exp(b-max(a,b)))

    end function logsumexp



end module cluster_module
