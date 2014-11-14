!> This module contains tools to calculate evidence estimators
module evidence_module
    implicit none

    type :: run_time_info

        ! Number of active clusters
        integer :: ncluster_A
        ! Number of passive clusters
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
            info%logZ(settings%nclustertot),                    &
            info%logZ2(settings%nclustertot),                   &
            info%logXX(settings%ncluster,settings%ncluster),    &
            info%logZX(settings%nclustertot,settings%ncluster)  &
            )

        ! Initially there is one active cluster, and no inactive ones
        info%ncluster_A = 1
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


    subroutine write_cluster_info(info)
        implicit none
        type(run_time_info), intent(in) :: info


        write(*,'("ncluster_A:", I4)') info%ncluster_A
        write(*,'("ncluster_T:", I4)') info%ncluster_T
        !write(*,'("logevidence:", E17.6)') info%logevidence
        !write(*,'("logevidence2:", E17.6)') info%logevidence2
        write(*,'("n:", <size(info%n)>I4)') info%n
        !write(*,'("logL:", <size(info%logL)>E17.6)') info%logL
        !write(*,'("logX:", <size(info%logX)>E17.6)') info%logX
        !write(*,'("logZ:", <size(info%logZ)>E17.6)') info%logZ
        !write(*,'("logZ2:", <size(info%logZ2)>E17.6)') info%logZ2
        !write(*,'("logXX:", <size(info%logXX,1)>E17.6)') info%logXX
        !write(*,'("logZX:", <size(info%logZX,1)>E17.6)') info%logZX

    end subroutine write_cluster_info





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
                r%logXX(j,i) = r%logXX(j,i) + logni - logni1
            else
                r%logXX(i,i) = r%logXX(i,i) + logni - logni2
            end if
        end do


        ! Decrease the number of live points in cluster i
        r%n(i) = r%n(i) - 1

        ! Update the loglikelihood
        r%logL(i) = newloglike


    end subroutine update_evidence


    subroutine bifurcate_evidence(r,i,ni) 
        use utils_module, only: logzero
        implicit none

        ! The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: r

        !> The new cluster indices
        !! NOTE: this is an array
        integer,intent(in),dimension(:) :: i
        !> The numbers in each new clusters
        integer,intent(in),dimension(:) :: ni

        ! Iterator
        integer :: j,k

        integer old_ncluster_A

        ! Temporary variables for notational ease
        double precision                     :: logn
        double precision                     :: logn1
        double precision,dimension(size(ni)) :: logni
        double precision,dimension(size(ni)) :: logni1
        double precision                     :: logX

        logn  = log( sum(ni) +0d0 ) 
        logn1 = log( sum(ni) +1d0 ) 
        logni = log( ni +0d0 )
        logni1= log( ni +1d0 )
        logX  = r%logX(i(1)) 



        ! Note that we now have new clusters
        old_ncluster_A = r%ncluster_A
        r%ncluster_A = r%ncluster_A + size(ni) -1
        r%ncluster_T = r%ncluster_T + size(ni)


        ! Make room for them in the r variable by pushing the passive cluster
        ! r up a few places (mostly evidence r)
        r%logZX(r%ncluster_A+2:,:) = r%logZX(old_ncluster_A+1:,:)
        r%logZ( r%ncluster_A+2:)   = r%logZ( old_ncluster_A+1:)
        r%logZ2(r%ncluster_A+2:)   = r%logZ2(old_ncluster_A+1:)

        ! Copy up the details in cluster i(1) to r%ncluster_A+1
        r%logZX(r%ncluster_A+1:,:) = r%logZX(i(1):,:)
        r%logZ( r%ncluster_A+1:)   = r%logZ( i(1):)
        r%logZ2(r%ncluster_A+1:)   = r%logZ2(i(1):)

        ! Now update the volumes 
        ! NOTE: i here is an array of indices
        r%logX(i) = logni - logn + logX

        ! Update the cross correlations 
        do j=1,size(i)
            do k=1,size(i)

                if(k==j) then
                    r%logXX(i(j),i(j)) = logni(j)+logni1(j)-logn-logn1 + logX
                else
                    r%logXX(i(j),i(k)) = logni(j)+logni(k)-logn-logn1  + logX
                end if

            end do
        end do

        ! Update the cross correlations that are uncorrelated
        do j=1,size(i)
            do k=1,r%ncluster_A

                if(all(k/=i(:))) then
                    r%logXX(i(j),k) = r%logX(i(j)) * r%logX(k)
                    r%logXX(k,i(j)) = r%logXX(i(j),k)
                end if

            end do
        end do

        ! Initialise the remaining variables
        r%logL(i) = r%logL(i(1)) 

        r%logZ(i) = logzero
        r%logZ2(i) = logzero
        r%logZX(:,i) = logzero

        r%n(i) = ni



    end subroutine bifurcate_evidence














end module evidence_module
