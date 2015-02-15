module run_time_module
    implicit none

    !> The run time information.
    !!
    !! This is what needs to be saved in order to resume a run.
    !! Bundling these all into the same type enables easy passing of data 
    !! from one fuction to another
    type run_time_info

        !> Number of dead points
        integer :: ndead

        !> Total number of likelihood calls
        integer :: nlike

        !> The number currently evolving clusters
        integer :: ncluster
        !> The number of live points in each cluster
        integer, allocatable, dimension(:) :: nlive
        !> The number of phantom points in each cluster
        integer, allocatable, dimension(:) :: nphantom
        !> The number of posterior points in each cluster
        integer, allocatable, dimension(:) :: nposterior

        !> Live points
        double precision, allocatable, dimension(:,:,:) :: live
        !> Phantom points
        double precision, allocatable, dimension(:,:,:) :: phantom
        !> Posterior points
        double precision, allocatable, dimension(:,:,:) :: posterior

        !> Covariance Matrices
        double precision, allocatable, dimension(:,:,:) :: covmat
        !> Cholesky decompositions
        double precision, allocatable, dimension(:,:,:) :: cholesky

        !> Global evidence estimate
        double precision :: logZ
        !> Global evidence^2 estimate
        double precision :: logZ2
        !> Local volume estimate
        double precision, allocatable, dimension(:)   :: logXp
        !> global evidence volume cross correlation
        double precision, allocatable, dimension(:)   :: logZXp
        !> Local evidence estimate
        double precision, allocatable, dimension(:)   :: logZp
        !> Local evidence^2 estimate 
        double precision, allocatable, dimension(:)   :: logZp2
        !> local evidence volume cross correlation
        double precision, allocatable, dimension(:)   :: logZpXp
        !> local volume cross correlation
        double precision, allocatable, dimension(:,:) :: logXpXq

        !> Minimum loglikelihood
        double precision :: logL
        !> Minimum loglikelihoods
        double precision, allocatable, dimension(:) :: logLp
        !> Cluster containing the minimum loglikelihood point
        integer :: p
        !> The minimum loglikelihood point within this cluster
        integer :: i

    end type run_time_info

    contains

    !> This is a self explanatory subroutine.
    !!
    !! It allocates the arrays for a single cluster 
    subroutine initialise_run_time_info(settings,RTI)
        use utils_module,    only: logzero,identity_matrix
        use settings_module, only: program_settings

        implicit none
        !> Program settings
        type(program_settings), intent(in) :: settings
        !> Run time information
        type(run_time_info),intent(out) :: RTI

        ! Allocate all of the arrays with one cluster
        RTI%ncluster = 1
        allocate(                                            &
            RTI%live(settings%nTotal,settings%nlive,1),      &
            RTI%phantom(settings%nTotal,settings%nlive,1),   &
            RTI%posterior(settings%nTotal,settings%nlive,1), &
            RTI%logZp(1),                                    &
            RTI%logXp(1),                                    &
            RTI%logZXp(1),                                   &
            RTI%logZp2(1),                                   &
            RTI%logZpXp(1),                                  &
            RTI%logXpXq(1,1),                                &
            RTI%logLp(1),                                    &
            RTI%nlive(1),                                    &
            RTI%nphantom(1),                                 &
            RTI%nposterior(1),                               &
            RTI%cholesky(settings%nDims,settings%nDims,1),   &
            RTI%covmat(settings%nDims,settings%nDims,1)      &
            )

        ! All evidences set to logzero
        RTI%logZ=logzero
        RTI%logZ2=logzero
        RTI%logZp=logzero
        RTI%logZXp=logzero
        RTI%logZp2=logzero
        RTI%logZpXp=logzero

        ! All volumes set to 1
        RTI%logXp=1d0
        RTI%logXpXq=1d0

        !Initially no live points at all
        RTI%nlive=0
        RTI%nphantom=0
        RTI%nposterior=0

        !No likelihood calls
        RTI%nlike=0

        !No dead points
        RTI%ndead=0

        !Cholesky and covmat set to identity
        RTI%cholesky(:,:,1) = identity_matrix(settings%nDims)
        RTI%covmat(:,:,1)   = identity_matrix(settings%nDims)

        ! Loglikelihoods at zero
        RTI%logL  = logzero
        RTI%logLp = logzero


    end subroutine initialise_run_time_info

    subroutine update_evidence(RTI)
        use utils_module, only: logsumexp,logincexp
        implicit none

        ! The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: RTI

        !> The cluster index to update
        integer :: p
        !> The loglikelihood to update
        double precision :: logL

        ! Iterator
        integer :: q

        ! Temporary variables for notational ease
        double precision,parameter :: log2 = log(2d0)
        double precision :: lognp
        double precision :: lognp1
        double precision :: lognp2


        logL  = RTI%logL
        p     = RTI%p

        lognp = log( RTI%nlive(p) +0d0 )
        lognp1= log( RTI%nlive(p) +1d0 )
        lognp2= log( RTI%nlive(p) +2d0 )


        ! Global evidence
        call logincexp( RTI%logZ, RTI%logXp(p)+logL-lognp1  )
        ! Local evidence
        call logincexp( RTI%logZp(p) , RTI%logXp(p)+logL-lognp1  )
        ! Local volume
        RTI%logXp(p)  = RTI%logXp(p) + lognp - lognp1


        ! Global evidence error
        call logincexp( RTI%logZ2 ,                                 &
            log2 + RTI%logZXp(p)  + logL - lognp1,              &
            log2 + RTI%logXpXq(p,p)  + 2*logL - lognp1 - lognp2 &
            )

        ! global evidence volume cross correlation p=p
        RTI%logZXp(p) = RTI%logZXp(p) + lognp - lognp1
        call logincexp( RTI%logZXp(p), &
            RTI%logXpXq(p,p)+ logL + lognp - lognp1 - lognp2 &
            )

        ! global evidence volume cross correlation p/=q
        do q=1,RTI%ncluster
            if(p/=q) call logincexp( RTI%logZXp(q) , RTI%logXpXq(p,q)+ logL - lognp1 )
        end do


        ! Local evidence error
        call logincexp( RTI%logZp2(p),                             &
            log2 + RTI%logZpXp(p)  + logL - lognp1,            &
            log2 + RTI%logXpXq(p,p)  + 2*logL - lognp1 - lognp2 &
            )


        ! Local evidence volume cross correlation
        RTI%logZpXp(p) = RTI%logZpXp(p) + lognp - lognp1
        call logincexp( RTI%logZpXp(p) , RTI%logXpXq(p,p)+ logL + lognp - lognp1 - lognp2 )


        ! Local volume cross correlation (p=p)
        RTI%logXpXq(p,p) = RTI%logXpXq(p,p) + lognp - lognp2

        ! Local volume cross correlation (p=q)
        do q=1,RTI%ncluster
            if(p/=q) then
                RTI%logXpXq(p,q) = RTI%logXpXq(p,q) + lognp - lognp1
                RTI%logXpXq(q,p) = RTI%logXpXq(q,p) + lognp - lognp1
            end if
        end do

    end subroutine update_evidence

    subroutine calculate_covmats(settings,RTI)
        use settings_module, only: program_settings
        use utils_module, only: calc_cholesky
        implicit none

        type(program_settings), intent(in) :: settings  !> Program settings
        type(run_time_info),intent(inout) :: RTI        !> Run time information

        integer :: i_cluster ! cluster iterator
        double precision, dimension(settings%nDims) :: mean ! The mean of a given cluster

        ! For each cluster:
        do i_cluster = 1,RTI%ncluster
            ! Calculate the mean
            mean = ( sum(RTI%live(settings%h0:settings%h1,1:RTI%nlive(i_cluster),i_cluster),dim=2) &
                + sum(RTI%phantom(settings%h0:settings%h1,1:RTI%nphantom(i_cluster),i_cluster),dim=2) ) &
                / (RTI%nlive(i_cluster) + RTI%nphantom(i_cluster) )

            ! Calculate the covariance by using a matrix multiplication
            RTI%covmat(:,:,i_cluster) = & 
                matmul(&
                RTI%live(settings%h0:settings%h1,1:RTI%nlive(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nlive(i_cluster)) , &
                transpose( RTI%live(settings%h0:settings%h1,1:RTI%nlive(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nlive(i_cluster)) ) &
                )&
                +&
                matmul(&
                RTI%phantom(settings%h0:settings%h1,1:RTI%nphantom(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nphantom(i_cluster)) , &
                transpose( RTI%phantom(settings%h0:settings%h1,1:RTI%nphantom(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nphantom(i_cluster)) ) &
                )

            ! Calculate the cholesky decomposition
            RTI%cholesky(:,:,i_cluster) = calc_cholesky(RTI%covmat(:,:,i_cluster))
        end do


    end subroutine calculate_covmats

    !> Calculate unbiased evidence estimates and errors. 
    !!
    !! The evidences generated by nested sampling are distributed according to a log-normal distribution:
    !! http://en.wikipedia.org/wiki/Log-normal_distribution
    !!
    !! What we accumulate in the routine update_evidence is log(<Z>), and log(<Z^2>).
    !! What we want is <log(Z)>,and 
    subroutine calculate_logZ_estimate(RTI,logZ,sigmalogZ,logZp,sigmalogZp)
        use utils_module, only: logzero
        implicit none

        type(run_time_info),intent(in)                                  :: RTI        !> Run time information
        double precision, intent(out)                                   :: logZ       !>
        double precision, intent(out)                                   :: sigmalogZ  !>
        double precision, intent(out), dimension(RTI%ncluster),optional :: logZp      !>
        double precision, intent(out), dimension(RTI%ncluster),optional :: sigmalogZp !>

        logZ       = max(logzero,2*RTI%logZ - 0.5*RTI%logZ2)
        sigmalogZ  = sqrt(abs(RTI%logZ2 - 2*RTI%logZ))

        if(present(logZp).and.present(sigmalogZp))then
            logZp      = max(logzero,2*RTI%logZp - 0.5*RTI%logZp2)
            sigmalogZp = sqrt(abs(RTI%logZp2 - 2*RTI%logZp))
        end if


    end subroutine calculate_logZ_estimate

    function live_logZ(settings,RTI)
        use utils_module, only: logzero,logsumexp,logincexp
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information

        double precision ::live_logZ ! Amount of evidence remaining in the live points

        integer :: i_cluster ! cluster iterator

        ! Initialise it with no log evidence
        live_logZ = logzero

        ! Sum up over all the clusters mean(likelihood) * volume
        do i_cluster = 1,RTI%ncluster
            call logincexp(live_logZ, &
                logsumexp(RTI%live(settings%l0,:RTI%nlive(i_cluster),i_cluster)) &
                - log(RTI%nlive(i_cluster)+0d0) &
                + RTI%logXp(i_cluster) &
                )
        end do

    end function live_logZ

    function replace_point(settings,RTI,baby_points,cluster_id) result(replaced)
        use utils_module, only: logzero,logsumexp,logincexp
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information
        integer,intent(in) :: cluster_id               !> Cluster identifier
        !> New-born baby points, created by slice sampling routine
        double precision,intent(in),dimension(settings%nTotal,settings%num_babies) :: baby_points

        logical :: replaced ! Have we successfully replaced a point?

        ! live point, last of the baby points
        double precision,dimension(settings%nTotal) :: point

        double precision :: logL ! loglikelihood bound

        integer :: i_point ! point iterator


        ! The loglikelihood contour is defined by the cluster it belongs to
        logL = RTI%logLp(cluster_id)

        ! Assign the points to cluster_id, if they are:
        ! 1) Within the isolikelihood contour of the cluster.
        ! 2) Within the voronoi cell of the cluster.

        do i_point=1,settings%num_babies-1
            ! Assign a temporary variable
            point = baby_points(:,i_point)

            if( point(settings%l0) > logL ) then ! (1)
                if( identify_cluster(settings,RTI,point) == cluster_id) then !(2)
                    call add_point(point,RTI%phantom,RTI%nphantom,cluster_id)  ! Add the new phantom point
                end if
            end if

        end do

        ! Now assign the live point
        point = baby_points(:,i_point)
        if( point(settings%l0) > logL ) then ! (1)
            if( identify_cluster(settings,RTI,point) == cluster_id) then !(2)

                replaced = .true.                                   ! Mark this as a replaced live point
                call kill_outermost_live_point(settings,RTI)        ! Kill the outermost live point
                call add_point(point,RTI%live,RTI%nlive,cluster_id) ! Add the new live point

            else
                replaced = .false.                                  ! We haven't killed of any points
            end if
        end if

    end function replace_point

    subroutine add_point(point,array,narray,cluster_id)
        implicit none
        !> Point to be added to end of array
        double precision, dimension(:), intent(in) :: point      
        !> Array to be added to
        double precision, dimension(:,:,:), allocatable, intent(inout) :: array
        !> number of points in array (second index)
        integer,dimension(:), allocatable, intent(inout) :: narray
        !> cluster identity (third index)
        integer, intent(in) :: cluster_id


        ! Increase the number of points in the cluster
        narray(cluster_id) = narray(cluster_id) + 1

        ! If this takes us over the size of the array, then double it
        if(narray(cluster_id) > size(array,2) ) call reallocate_3(array, u2=size(array,2)*2 )

        ! Add the point to the end position
        array(:,narray(cluster_id),cluster_id) = point

    end subroutine add_point

    subroutine delete_point(i_point,array,narray,cluster_id)
        implicit none
        !> Position of point to be deleted from array
        integer, intent(in) :: i_point
        !> Array to be delete from
        double precision, dimension(:,:,:), allocatable, intent(inout) :: array
        !> number of points in array (second index)
        integer,dimension(:), allocatable, intent(inout) :: narray
        !> cluster identity (third index)
        integer, intent(in) :: cluster_id

        ! delete the point by overwriting it with the point at the end 
        array(:,i_point,cluster_id) = array(:,narray(cluster_id),cluster_id)

        ! reduce the number of points in the cluster
        narray(cluster_id) = narray(cluster_id) - 1

    end subroutine delete_point


    subroutine kill_outermost_live_point(settings,RTI) 
        use settings_module,   only: program_settings
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: RTI

        call update_evidence(RTI)                        ! Update the evidence value
        call delete_point(RTI%i,RTI%live,RTI%nlive,RTI%p)! Delete the live point from the array
        call find_global_min(settings,RTI)               ! Find the new global minimum

    end subroutine kill_outermost_live_point

    subroutine find_global_min(settings,RTI) 
        use settings_module,   only: program_settings
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: RTI


    end subroutine find_global_min




    function identify_cluster(settings,RTI,point) result(cluster)
        use settings_module,   only: program_settings
        use utils_module,      only: loginf,distance2
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(in) :: RTI

        double precision, dimension(settings%nTotal),intent(in)   :: point

        integer :: cluster

        integer :: i_cluster
        integer :: i_live

        double precision :: temp_distance2
        double precision :: closest_distance2

        if( RTI%ncluster == 1) then
            cluster=1
            return
        end if

        closest_distance2=loginf

        ! Find the cluster this point is nearest to
        do i_cluster=1,RTI%ncluster
            do i_live=1,RTI%nlive(i_cluster)
                temp_distance2 = distance2(point(settings%h0:settings%h1),RTI%live(settings%h0:settings%h1,i_live,i_cluster) )
                if(temp_distance2 < closest_distance2) then
                    cluster = i_cluster
                    closest_distance2 = temp_distance2
                end if
            end do
        end do

    end function identify_cluster




end module
