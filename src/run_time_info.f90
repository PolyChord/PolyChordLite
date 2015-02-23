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

        !> Live points
        double precision, allocatable, dimension(:,:,:) :: live
        !> The number of live points in each cluster
        integer, allocatable, dimension(:) :: nlive
        !> Phantom points
        double precision, allocatable, dimension(:,:,:) :: phantom
        !> The number of phantom points in each cluster
        integer, allocatable, dimension(:) :: nphantom
        !> Posterior points
        double precision, allocatable, dimension(:,:,:) :: posterior
        !> The number of posterior points in each cluster
        integer, allocatable, dimension(:) :: nposterior


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

        !> Minimum loglikelihoods
        double precision, allocatable, dimension(:) :: logLp
        !> The minimum loglikelihood point within each cluster
        integer,allocatable, dimension(:)           :: i

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
        allocate(                                               &
            RTI%live(settings%nTotal,settings%nlive,1),         &
            RTI%phantom(settings%nTotal,settings%nlive,1),      &
            RTI%posterior(settings%nposterior,settings%nlive,1),&
            RTI%logZp(1),                                       &
            RTI%logXp(1),                                       &
            RTI%logZXp(1),                                      &
            RTI%logZp2(1),                                      &
            RTI%logZpXp(1),                                     &
            RTI%logXpXq(1,1),                                   &
            RTI%logLp(1),                                       &
            RTI%i(1),                                           &
            RTI%nlive(1),                                       &
            RTI%nphantom(1),                                    &
            RTI%nposterior(1),                                  &
            RTI%cholesky(settings%nDims,settings%nDims,1),      &
            RTI%covmat(settings%nDims,settings%nDims,1)         &
            )

        ! All evidences set to logzero
        RTI%logZ=logzero
        RTI%logZ2=logzero
        RTI%logZp=logzero
        RTI%logZXp=logzero
        RTI%logZp2=logzero
        RTI%logZpXp=logzero

        ! All volumes set to 1
        RTI%logXp=0d0
        RTI%logXpXq=0d0

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
        RTI%logLp = logzero
        ! First position default lowest
        RTI%i     = 0


    end subroutine initialise_run_time_info

    function update_evidence(RTI,p) result(logweight)
        use utils_module, only: logsumexp,logincexp
        implicit none

        !> The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: RTI

        !> The cluster index to update
        integer :: p

        ! The loglikelihood to update
        double precision :: logL

        ! The logweight of the deleted point
        double precision :: logweight

        ! Iterator
        integer :: q

        ! Temporary variables for notational ease
        double precision,parameter :: log2 = log(2d0)
        double precision :: lognp
        double precision :: lognp1
        double precision :: lognp2

        logL  = RTI%logLp(p)

        lognp = log( RTI%nlive(p) +0d0 )
        lognp1= log( RTI%nlive(p) +1d0 )
        lognp2= log( RTI%nlive(p) +2d0 )

        ! Output the logweight
        logweight =  RTI%logXp(p) - lognp1

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

        ! Update the number of dead points
        RTI%ndead = RTI%ndead+1

    end function update_evidence

    !> This function takes the evidence info in cluster i and splits it into
    !! several clusters according to the cluster numbers ni.
    !!
    !! It places these new clusters at the end of the active array, and deletes the old cluster
    !!
    subroutine add_cluster(settings,RTI,p,cluster_list,num_new_clusters) 
        use settings_module, only: program_settings
        use utils_module, only: logzero,logsumexp,logaddexp
        use array_module, only: reallocate_3_d,reallocate_2_d,reallocate_1_d,reallocate_1_i,add_point
        implicit none

        type(program_settings), intent(in) :: settings  !> Program settings
        !> The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: RTI

        !> The cluster index to be split
        integer,intent(in)              :: p
        !> The numbers in each new clusters
        integer,intent(in),dimension(:) :: cluster_list
        !> The number of new clusters
        integer,intent(in) :: num_new_clusters


        !Iterators
        integer :: i_live
        integer :: i_phantom,j_phantom
        integer :: i_cluster

        ! Constructor
        integer :: i


        integer, dimension(RTI%ncluster-1)   :: old_save
        integer, dimension(RTI%ncluster-1)   :: old_target
        integer, dimension(num_new_clusters) :: new_target
        integer                              :: num_old_clusters

        double precision, dimension(size(RTI%live,1),size(cluster_list))     :: old_live
        double precision, dimension(size(RTI%phantom,1),size(RTI%phantom,2),size(RTI%phantom,3)) :: old_phantom
        integer, dimension(size(RTI%nphantom)) :: old_nphantom

        double precision, dimension(num_new_clusters) :: logni
        double precision, dimension(num_new_clusters) :: logni1
        double precision :: logn
        double precision :: logn1
        double precision :: logXp
        double precision, dimension(RTI%ncluster-1) :: logXpXq
        double precision :: logXp2
        double precision :: logZXp

        ! 1) Save the old points as necessary
        old_live  = RTI%live(:,:RTI%nlive(p),p)  ! Save the old live points
        old_phantom = RTI%phantom                ! Save the old phantom points
        old_nphantom= RTI%nphantom               ! Save the old numbers of phantom points

        old_save   = [(i,i=1,p-1),(i,i=p+1,RTI%ncluster)]    ! The indices of the old clusters to save
        num_old_clusters = RTI%ncluster-1                    ! The number of old clusters
        RTI%ncluster = RTI%ncluster+num_new_clusters-1       ! The total number of new clusters

        old_target = [(i,i=1,num_old_clusters)]              ! Where the old clusters will be
        new_target = [(i,i=num_old_clusters+1,RTI%ncluster)] ! Where we're going to insert the new clusters

        ! Define some useful variables
        logXp  = RTI%logXp(p)
        logXp2 = RTI%logXpXq(p,p)
        logZXp = RTI%logZXp(p)
        logXpXq= [ RTI%logXpXq(p,:p-1) , RTI%logXpXq(p,p+1:) ]



        ! 2) Reallocate the arrays

        ! Reallocate the live,phantom and posterior points
        call reallocate_3_d(RTI%live,      new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nlive,     new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)
        call reallocate_3_d(RTI%phantom,   new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nphantom,  new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)
        call reallocate_3_d(RTI%posterior, new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nposterior,new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)

        ! Reallocate the cholesky matrices
        call reallocate_3_d(RTI%cholesky, new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_3_d(RTI%covmat,   new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)

        ! Reallocate the evidence arrays 
        call reallocate_1_d(RTI%logXp,   RTI%ncluster,old_save,old_target)
        call reallocate_1_d(RTI%logZXp,  RTI%ncluster,old_save,old_target)
        call reallocate_1_d(RTI%logZp,   RTI%ncluster,old_save,old_target)
        call reallocate_1_d(RTI%logZp2,  RTI%ncluster,old_save,old_target)
        call reallocate_1_d(RTI%logZpXp, RTI%ncluster,old_save,old_target)
        call reallocate_2_d(RTI%logXpXq, RTI%ncluster,RTI%ncluster,old_save,old_save,old_target,old_target)

        call reallocate_1_d(RTI%logLp,   RTI%ncluster,old_save,old_target) 
        call reallocate_1_i(RTI%i,       RTI%ncluster,old_save,old_target)



        ! 3) Assign the new live points to their new clusters

        RTI%nlive(new_target) = 0 ! Zero the number of live points in the new clusters
        do i_live=1,size(cluster_list)
            ! Insert the new points in the correct positions
            call add_point(old_live(:,i_live),RTI%live,RTI%nlive,new_target(cluster_list(i_live)))
        end do

        ! Find the new minimum loglikelihoods
        call find_min_loglikelihoods(settings,RTI) 

        ! 4) Reassign all the phantom points 
        RTI%nphantom = 0
        RTI%nposterior = 0
        do i_cluster=1,size(old_nphantom)
            do i_phantom=1,old_nphantom(i_cluster)
                ! Reallocate all of the phantom points
                j_phantom = identify_cluster(settings,RTI,old_phantom(:,i_phantom,i_cluster))
                if(old_phantom(settings%l0,i_phantom,i_cluster) > RTI%logLp(j_phantom) ) &
                    call add_point(old_phantom(:,i_phantom,i_cluster),RTI%phantom,RTI%nphantom,j_phantom)
            end do
        end do

        ! 5) Initialise the new evidences and volumes
        ! Re-define where we're going to move the old clusters
        ! so that we don't include the deleted cluster

        ! Find the number of live+phantom points in each cluster
        logni  = log(RTI%nlive(new_target) + RTI%nphantom(new_target) + 0d0)
        logni1 = log(RTI%nlive(new_target) + RTI%nphantom(new_target) + 1d0)
        logn   = logsumexp(logni)
        logn1  = logaddexp(logn,0d0)

        ! Initialise the new volumes
        RTI%logXp(new_target) = logXp + logni - logn
        ! Initialise the new global evidence -local volume cross correlation
        RTI%logZXp(new_target) = logZXp + logni - logn 
        ! Initialise local evidences at 0
        RTI%logZp(new_target) = logzero
        RTI%logZp2(new_target) = logzero
        RTI%logZpXp(new_target) = logzero


        ! Initialise the volume cross correlations
        if(num_old_clusters>0) then
            RTI%logXpXq(new_target,old_target) = spread(logXpXq,1,num_new_clusters) + transpose(spread(logni,1,num_old_clusters)) - logn
            RTI%logXpXq(old_target,new_target) = transpose(RTI%logXpXq(new_target,old_target))
        end if

        ! When they're both new clusters
        RTI%logXpXq(new_target,new_target) &
            = logXp2 + spread(logni,1,num_new_clusters) + transpose(spread(logni,1,num_new_clusters)) -logn-logn1

        ! the squared one is slightly different
        do i_cluster=1,num_new_clusters

            RTI%logXpXq(new_target(i_cluster),new_target(i_cluster)) &
                = logXp2 + logni(i_cluster)+ logni1(i_cluster)-logn-logn1

        end do

    end subroutine add_cluster

    subroutine delete_cluster(settings,RTI) 
        use settings_module, only: program_settings
        use array_module, only: reallocate_3_d,reallocate_2_d,reallocate_1_d,reallocate_1_i
        implicit none
        type(program_settings), intent(in) :: settings  !> Program settings
        !> The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: RTI

        !The cluster index to be deleted
        integer            :: p(1)

        ! new indices of clusters
        integer,dimension(RTI%ncluster-1) :: indices

        ! Contstructor iterator
        integer :: i


        if(any(RTI%nlive==0)) then

            p=minloc(RTI%nlive,RTI%nlive==0)

            ! Get the positions of the clusters to be saved
            indices = [(i,i=1,p(1)-1),(i,i=p(1)+1,RTI%ncluster)]

            ! Reduce the number of clusters
            RTI%ncluster=RTI%ncluster-1

            ! Reallocate the live,phantom and posterior points
            call reallocate_3_d(RTI%live,      new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nlive,     new_size1=RTI%ncluster, save_indices1=indices)
            call reallocate_3_d(RTI%phantom,   new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nphantom,  new_size1=RTI%ncluster, save_indices1=indices)
            call reallocate_3_d(RTI%posterior, new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nposterior,new_size1=RTI%ncluster, save_indices1=indices)

            ! Reallocate the cholesky matrices
            call reallocate_3_d(RTI%cholesky, new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_3_d(RTI%covmat,   new_size3=RTI%ncluster, save_indices3=indices)

            ! Reallocate the evidence arrays 
            call reallocate_1_d(RTI%logXp,   new_size1=RTI%ncluster,save_indices1=indices)
            call reallocate_1_d(RTI%logZXp,  new_size1=RTI%ncluster,save_indices1=indices)
            call reallocate_1_d(RTI%logZp,   new_size1=RTI%ncluster,save_indices1=indices)
            call reallocate_1_d(RTI%logZp2,  new_size1=RTI%ncluster,save_indices1=indices)
            call reallocate_1_d(RTI%logZpXp, new_size1=RTI%ncluster,save_indices1=indices)
            call reallocate_2_d(RTI%logXpXq, new_size1=RTI%ncluster,new_size2=RTI%ncluster,save_indices1=indices,save_indices2=indices)

            call reallocate_1_d(RTI%logLp,   new_size1=RTI%ncluster,save_indices1=indices) 
            call reallocate_1_i(RTI%i,       new_size1=RTI%ncluster,save_indices1=indices)
        end if


    end subroutine delete_cluster


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
            mean = ( sum(RTI%live(settings%h0:settings%h1,:RTI%nlive(i_cluster),i_cluster),dim=2) &
                + sum(RTI%phantom(settings%h0:settings%h1,:RTI%nphantom(i_cluster),i_cluster),dim=2) ) &
                / (RTI%nlive(i_cluster) + RTI%nphantom(i_cluster) )

            ! Calculate the covariance by using a matrix multiplication
            RTI%covmat(:,:,i_cluster) =( & 
                matmul(&
                RTI%live(settings%h0:settings%h1,:RTI%nlive(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nlive(i_cluster)) , &
                transpose( RTI%live(settings%h0:settings%h1,:RTI%nlive(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nlive(i_cluster)) ) &
                )&
                +&
                matmul(&
                RTI%phantom(settings%h0:settings%h1,:RTI%nphantom(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nphantom(i_cluster)) , &
                transpose( RTI%phantom(settings%h0:settings%h1,:RTI%nphantom(i_cluster),i_cluster) &
                - spread(mean,dim=2,ncopies=RTI%nphantom(i_cluster)) ) &
                ) &
                )/ (RTI%nlive(i_cluster) + RTI%nphantom(i_cluster) ) 

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
            if(RTI%nlive(i_cluster)>0) then
                call logincexp(live_logZ, &
                        logsumexp(RTI%live(settings%l0,:RTI%nlive(i_cluster),i_cluster)) &
                        - log(RTI%nlive(i_cluster)+0d0) &
                        + RTI%logXp(i_cluster) &
                        )
            end if
        end do

    end function live_logZ






    function replace_point(settings,RTI,baby_points,cluster_add) result(replaced)
        use utils_module, only: logsumexp,logincexp,minpos
        use settings_module, only: program_settings
        use calculate_module, only: calculate_posterior_point
        use random_module, only: bernoulli_trial
        use array_module, only: add_point,delete_point

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information
        integer,intent(in) :: cluster_add              !> Cluster to add to
        !> New-born baby points, created by slice sampling routine
        double precision,intent(in),dimension(settings%nTotal,settings%num_babies) :: baby_points

        logical :: replaced ! Have we successfully replaced a point?

        ! live point, last of the baby points
        double precision,dimension(settings%nTotal) :: point

        double precision :: logL ! loglikelihood bound

        integer :: i_baby ! point iterator

        integer                                     :: cluster_del     ! cluster to delete from
        double precision,dimension(settings%nTotal) :: deleted_point   ! point we have just deleted
        double precision                            :: logweight       ! The log weighting of this point
        
        integer :: i_phantom ! phantom iterator


        ! The loglikelihood contour is defined by the cluster it belongs to
        logL = RTI%logLp(cluster_add)

        ! Assign the phantom points to cluster_add, if they are:
        ! (1) Within the isolikelihood contour of the cluster.
        ! (2) Within the voronoi cell of the cluster.

        do i_baby=1,settings%num_babies-1
            ! Assign a temporary variable
            point = baby_points(:,i_baby)

            if( point(settings%l0) > logL ) then ! (1)
                if( identify_cluster(settings,RTI,point) == cluster_add) then !(2)
                    call add_point(point,RTI%phantom,RTI%nphantom,cluster_add)  ! Add the new phantom point
                end if
            end if

        end do

        ! Now assign the live point
        point = baby_points(:,i_baby)

        if( point(settings%l0) > logL ) then ! (1)
            if( identify_cluster(settings,RTI,point) == cluster_add) then !(2)

                replaced = .true.  ! Mark this as a replaced live point

                cluster_del   = minpos(RTI%logLp)                                                ! find the cluster we're deleting from
                logweight     = update_evidence(RTI,cluster_del)                                 ! Update the evidence value
                deleted_point = delete_point(RTI%i(cluster_del),RTI%live,RTI%nlive,cluster_del)  ! Delete the live point from the array
                call add_point(point,RTI%live,RTI%nlive,cluster_add)                             ! Add the new live point
                call find_min_loglikelihoods(settings,RTI)                                       ! Find the new minimum likelihoods


                ! Calculate the posterior point and add it to the array
                if(settings%calculate_posterior .and.  bernoulli_trial(settings%thin_posterior)) &
                    call add_point(&
                    calculate_posterior_point(settings,deleted_point,logweight,RTI%logZ,logsumexp(RTI%logXp)),&
                    RTI%posterior,RTI%nposterior,cluster_del )


                ! Now we delete the phantoms
                i_phantom = 1
                do while(i_phantom<=RTI%nphantom(cluster_del))

                    ! Delete points lower than the new loglikelihood bound
                    if ( RTI%phantom(settings%l0,i_phantom,cluster_del) < RTI%logLp(cluster_del) ) then

                        ! Delete this point
                        deleted_point = delete_point(i_phantom,RTI%phantom,RTI%nphantom,cluster_del)

                        ! Calculate the posterior point and add it to the array
                        if(settings%calculate_posterior .and. bernoulli_trial(settings%thin_posterior)) &
                            call add_point(&
                            calculate_posterior_point(settings,deleted_point,logweight,RTI%logZp(cluster_del),RTI%logXp(cluster_del)),&
                            RTI%posterior,RTI%nposterior,cluster_del )

                    else
                        i_phantom = i_phantom+1
                    end if

                end do
            else
                replaced = .false.                                  ! We haven't killed of any points
            end if
        else
            replaced = .false.
        end if

    end function replace_point




    subroutine find_min_loglikelihoods(settings,RTI)
        use utils_module, only: minpos,loginf
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information
        
        integer :: i_cluster     ! cluster iterator

        ! Iterate through each cluster
        do i_cluster=1,RTI%ncluster

            ! Find the position of the lowest point in this cluster
            RTI%i(i_cluster)     = minpos(RTI%live(settings%l0,:RTI%nlive(i_cluster),i_cluster))

            if(RTI%i(i_cluster) == 0) then
                ! If the cluster is empty, we need to signal to delete all points
                RTI%logLp(i_cluster) = loginf
            else
                ! Find the likelihood of the lowest point in this cluster
                RTI%logLp(i_cluster) = RTI%live(settings%l0,RTI%i(i_cluster),i_cluster) 
            end if

        end do

    end subroutine find_min_loglikelihoods



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
