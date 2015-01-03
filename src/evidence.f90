!> This module contains tools to calculate evidence estimators
module evidence_module
    implicit none

    type :: run_time_info

        ! Number of active clusters 
        ! These are clusters that are currently evolving
        integer :: ncluster_A
        ! Number of passive clusters
        ! These are clusters that have 'died'
        integer :: ncluster_P
        ! Total number of clusters, including pieces of clusters before
        ! splitting
        integer :: ncluster_T

        ! Global evidence, and evidence squared
        ! log(Z)
        double precision :: logevidence
        ! log(Z^2)
        double precision :: logevidence2

        ! local cluster variables
        ! n_i - live points in each cluster
        integer, allocatable, dimension(:)            :: n
        ! log(L_i) - log likelihood contour of each cluster
        double precision, allocatable, dimension(:)   :: logL
        ! log(X_i) - log volume for the cluster (since splitting)
        double precision, allocatable, dimension(:)   :: logX
        ! log(Z_i) - log evidence for the cluster (since splitting)
        double precision, allocatable, dimension(:)   :: logZ
        ! log(Z_i^2) - log evidence for the cluster (since splitting)
        double precision, allocatable, dimension(:)   :: logZ2

        ! local cluster correlations
        ! log(X_iX_j)
        double precision, allocatable, dimension(:,:) :: logXX
        ! log(Z_iX_j)
        double precision, allocatable, dimension(:,:) :: logZX
        
    end type run_time_info


    contains

    subroutine allocate_run_time_info(settings,info)
        use utils_module,    only: logzero
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),intent(out) :: info

        allocate(                                               &
            info%n(settings%ncluster),                          &
            info%logL(settings%ncluster),                       &
            info%logX(settings%ncluster),                       &
            info%logZ(settings%ncluster*2),                     &
            info%logZ2(settings%ncluster*2),                    &
            info%logXX(settings%ncluster,settings%ncluster),    &
            info%logZX(settings%ncluster*2,settings%ncluster)   &
            )

        ! Initially there is one active cluster, and no inactive ones
        info%ncluster_A = 1
        info%ncluster_P = 0
        info%ncluster_T = info%ncluster_A

        ! All live points in the first cluster to start with
        info%n=0
        info%n(1) = settings%nlive

        ! All evidences set to logzero
        info%logevidence=logzero
        info%logevidence2=logzero
        info%logZ=logzero
        info%logZ2=logzero

        ! All volumes set to logzero, apart from first cluster where they're set to log(1)=0
        info%logX       = logzero
        info%logXX      = logzero
        info%logX(1)    = 0
        info%logXX(1,1) = 0

        ! Cross correlation set to logzero
        info%logZX = logzero

        ! Initial error bounds set to logzero
        info%logL = logzero

    end subroutine allocate_run_time_info


    subroutine update_evidence(r,i,newloglike)
        use utils_module, only: logsumexp,logincexp
        implicit none

        ! The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: r

        !> The cluster index to update
        integer,intent(in) :: i
        !> The loglikelihood to update
        double precision,intent(in) :: newloglike

        ! Iterator
        integer :: j

        ! Temporary variables for notational ease
        double precision :: log2
        double precision :: logni
        double precision :: logni1
        double precision :: logni2
        double precision :: logLi
        double precision :: logXi
        double precision :: logXi2
        double precision :: logZiXi
        double precision :: logZXi

        log2  = log(2d0)
        logni = log( r%n(i) +0d0 )
        logni1= log( r%n(i) +1d0 )
        logni2= log( r%n(i) +2d0 )

        logLi = r%logL(i) 
        logXi = r%logX(i) 
        logXi2= r%logXX(i,i) 

        logZiXi = r%logZX(i,i)
        logZXi =  logsumexp( r%logZX(:r%ncluster_T,i) )

        ! Local evidence
        call logincexp( r%logZ(i), logXi+logLi-logni1  )
        ! Global evidence
        call logincexp( r%logevidence   , logXi+logLi-logni1  )

        ! Local evidence error
        call logincexp( r%logZ2(i),      log2 + logZiXi + logLi - logni1, log2 + logXi2  + 2*logLi - logni1 - logni2)
        ! Global evidence error
        call logincexp( r%logevidence2 , log2 + logZXi  + logLi - logni1, log2 + logXi2  + 2*logLi - logni1 - logni2)

        ! Update cross correlations logZX for all clusters
        r%logZX(i,i) = r%logZX(i,i) + logni - logni1
        call logincexp( r%logZX(i,i) , logXi2 + logLi + logni - logni1 - logni2 )

        do j=1,r%ncluster_A
            if(j/=i) call logincexp( r%logZX(i,j) , r%logXX(i,j) + logLi - logni1 )
        end do

        do j=1,r%ncluster_T
            if(j/=i) r%logZX(j,i) = r%logZX(j,i) + logni - logni1
        end do


        ! Update local volumes and cross correlations
        r%logX(i)  = r%logX(i) + logni - logni1

        do j=1,r%ncluster_A
            if(j/=i) then
                r%logXX(i,j) = r%logXX(i,j) + logni - logni1
                r%logXX(j,i) = r%logXX(i,j)
            else
                r%logXX(i,i) = r%logXX(i,i) + logni - logni2
            end if
        end do

        ! Update the loglikelihood
        r%logL(i) = newloglike

        ! Decrease the number of live points in cluster i
        r%n(i) = r%n(i) - 1


    end subroutine update_evidence

    subroutine delete_evidence(r,i,deleted) 
        implicit none
        ! The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: r

        !> The cluster index to update
        integer,intent(in) :: i

        !> Whether or not we're bifurcating or deleting. .true. => deleting
        logical, intent(in) :: deleted

        integer :: top_bound

        top_bound = r%ncluster_A

        if(.not. deleted) top_bound = top_bound + r%ncluster_P

        ! move everything to the end by a cyclic shift
        r%n(    i:top_bound)  = cshift(r%n(    i:top_bound),shift=1)
        r%logL( i:top_bound)  = cshift(r%logL( i:top_bound),shift=1)
        r%logX( i:top_bound)  = cshift(r%logX( i:top_bound),shift=1)
        r%logZ( i:top_bound)  = cshift(r%logZ( i:top_bound),shift=1)
        r%logZ2(i:top_bound)  = cshift(r%logZ2(i:top_bound),shift=1)

        r%logXX(i:top_bound,:) = cshift(r%logXX(i:top_bound,:),shift=1,dim=1)
        r%logXX(:,i:top_bound) = cshift(r%logXX(:,i:top_bound),shift=1,dim=2)

        r%logZX(i:top_bound,:) = cshift(r%logZX(i:top_bound,:),shift=1,dim=1)
        r%logZX(:,i:top_bound) = cshift(r%logZX(:,i:top_bound),shift=1,dim=2)

        ! Now decrease the number of active points and increase the number
        ! of passive points
        r%ncluster_A= r%ncluster_A - 1
        if(deleted) r%ncluster_P = r%ncluster_P + 1


    end subroutine delete_evidence


    !> This function takes the evidence info in cluster i and splits it into
    !! several clusters according to the cluster numbers ni.
    !!
    !! It places these new clusters at the end of the active array, but it does
    !! not delete the old cluster. This will be done by the delete_clusters
    !! routine below.
    !!
    subroutine bifurcate_evidence(r,i,ni) 
        use utils_module, only: logzero
        implicit none

        ! The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: r

        !> The cluster index to be split
        integer,intent(in)              :: i
        !> The numbers in each new clusters
        integer,intent(in),dimension(:) :: ni



        ! Iterator
        integer :: j,k
        ! Addresses of new clusters
        integer, dimension(size(ni)) :: ii

        ! Temporary variables for notational ease
        double precision                     :: logn
        double precision                     :: logn1
        double precision,dimension(size(ni)) :: logni
        double precision,dimension(size(ni)) :: logni1
        double precision                     :: logX
        double precision                     :: logX2

        ! shift the existing passive clusters to make room for the size(ni) 
        ! new active clusters.
        r%logZ (r%ncluster_A+1+size(ni):r%ncluster_T+size(ni))   = r%logZ( r%ncluster_A+1:r%ncluster_T) 
        r%logZ2(r%ncluster_A+1+size(ni):r%ncluster_T+size(ni))   = r%logZ2(r%ncluster_A+1:r%ncluster_T) 
        r%logZX(r%ncluster_A+1+size(ni):r%ncluster_T+size(ni),:) = r%logZX(r%ncluster_A+1:r%ncluster_T,:)

        ! For notational convenience, create an array of integers detailing the
        ! addresses of the new clusters
        ii = [ ( j, j=r%ncluster_A+1,r%ncluster_A+size(ni) ) ]

        ! Define some useful shorthands
        logn  = log( sum(ni) + 0d0 ) 
        logn1 = log( sum(ni) + 1d0 ) 
        logni = log( ni + 0d0 )
        logni1= log( ni + 1d0 )
        logX  = r%logX(i) 
        logX2 = r%logXX(i,i) 


        ! Note that we now have new clusters
        r%ncluster_A = r%ncluster_A + size(ni)
        r%ncluster_T = r%ncluster_T + size(ni)

        ! Initialise the default variables
        r%logL(ii) = r%logL(i) 

        r%logZ(ii) = logzero
        r%logZ2(ii) = logzero

        r%n(ii) = ni


        ! Now initialise the volumes 
        ! NOTE: ii here is an array of indices
        r%logX(ii) = logni - logn + logX

        ! initialise the cross correlations 
        do j=1,size(ii)
            do k=1,size(ii)

                if(k==j) then
                    r%logXX(ii(j),ii(j)) = logni(j)+logni1(j)-logn-logn1 + logX2
                else 
                    r%logXX(ii(j),ii(k)) = logni(j)+logni(k) -logn-logn1 + logX2
                end if

            end do
        end do

        ! initialise the cross correlations that are uncorrelated
        do k=1,r%ncluster_A

            if(all(k/=ii)) then

                r%logXX(ii,k) = r%logX(ii) + r%logX(k)
                r%logXX(k,ii) = r%logXX(ii,k)

            end if

        end do

        r%logZX(ii,:) = logzero

        do k=1,r%ncluster_T
            r%logZX(k,ii) = logni - logn + r%logZX(k,i)
        end do

    end subroutine bifurcate_evidence














end module evidence_module
