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

        !> The number currently evolving clusters
        integer :: ncluster
        !> The number of dead clusters
        integer :: ncluster_dead

        !> Total number of likelihood calls
        integer,allocatable,dimension(:) :: nlike

        !> The number of repeats within each parameter speed to do
        integer,allocatable, dimension(:)           :: num_repeats

        !> The number of live points in each cluster
        integer, allocatable, dimension(:) :: nlive 
        !> The number of phantom points in each cluster
        integer, allocatable, dimension(:) :: nphantom
        !> The number of weighted posterior points in each cluster
        integer, allocatable, dimension(:) :: nposterior
        !> The number of equally weighted posterior points in each cluster
        integer, allocatable, dimension(:) :: nequals


        integer, allocatable, dimension(:) :: nposterior_dead
        integer, allocatable, dimension(:) :: nposterior_global
        integer, allocatable, dimension(:) :: nequals_dead
        integer, allocatable, dimension(:) :: nequals_global

        !> Live points
        double precision, allocatable, dimension(:,:,:) :: live
        !> Phantom points
        double precision, allocatable, dimension(:,:,:) :: phantom
        !> Posterior stack
        double precision, allocatable, dimension(:,:,:) :: posterior_stack
        !> The number of posterior points in each cluster in the posterior stack
        integer, allocatable, dimension(:) :: nposterior_stack


        !> weighted posterior points
        double precision, allocatable, dimension(:,:,:) :: posterior
        double precision, allocatable, dimension(:,:,:) :: posterior_dead
        double precision, allocatable, dimension(:,:,:) :: posterior_global

        !> Equally weighted posterior points
        double precision, allocatable, dimension(:,:,:) :: equals
        double precision, allocatable, dimension(:,:,:) :: equals_dead
        double precision, allocatable, dimension(:,:,:) :: equals_global


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
        double precision, allocatable, dimension(:)   :: logZp_dead
        !> Local evidence^2 estimate 
        double precision, allocatable, dimension(:)   :: logZp2
        double precision, allocatable, dimension(:)   :: logZp2_dead
        !> local evidence volume cross correlation
        double precision, allocatable, dimension(:)   :: logZpXp
        !> local volume cross correlation
        double precision, allocatable, dimension(:,:) :: logXpXq

        !> Minimum loglikelihoods
        double precision, allocatable, dimension(:) :: logLp
        !> The minimum loglikelihood point within each cluster
        integer,allocatable, dimension(:)           :: i

        !> Maximum weight
        double precision, allocatable, dimension(:) :: maxlogweight
        double precision, allocatable, dimension(:) :: maxlogweight_dead
        double precision                            :: maxlogweight_global

        !> what to thin the posterior by
        double precision :: thin_posterior

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
        RTI%ncluster_dead = 0
        allocate(                                                       &
            RTI%live(settings%nTotal,settings%nlive,1),                 &
            RTI%phantom(settings%nTotal,settings%nlive,1),              &
            RTI%posterior_stack(settings%nposterior,settings%nlive,1),  &
            RTI%posterior(settings%nposterior,settings%nlive,1),        &
            RTI%posterior_dead(settings%nposterior,settings%nlive,0),   &
            RTI%posterior_global(settings%nposterior,settings%nlive,1), &
            RTI%equals(settings%np,settings%nlive,1),                   &
            RTI%equals_dead(settings%np,settings%nlive,0),              &
            RTI%equals_global(settings%np,settings%nlive,1),            &
            RTI%logZp(1),                                               &
            RTI%logZp_dead(0),                                          &
            RTI%logXp(1),                                               &
            RTI%logZXp(1),                                              &
            RTI%logZp2(1),                                              &
            RTI%logZp2_dead(0),                                         &
            RTI%logZpXp(1),                                             &
            RTI%logXpXq(1,1),                                           &
            RTI%logLp(1),                                               &
            RTI%i(1),                                                   &
            RTI%nlive(1),                                               &
            RTI%nphantom(1),                                            &
            RTI%nposterior_stack(1),                                    &
            RTI%nposterior(1),                                          &
            RTI%nposterior_dead(0),                                     &
            RTI%nposterior_global(1),                                   &
            RTI%nequals(1),                                             &
            RTI%nequals_dead(0),                                        &
            RTI%nequals_global(1),                                      &
            RTI%maxlogweight(1),                                        &
            RTI%maxlogweight_dead(0),                                   &
            RTI%cholesky(settings%nDims,settings%nDims,1),              &
            RTI%covmat(settings%nDims,settings%nDims,1),                &
            RTI%num_repeats(size(settings%grade_dims)),                 &
            RTI%nlike(size(settings%grade_dims))                        &
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
        RTI%nposterior_stack=0
        RTI%nposterior=0
        RTI%nequals=0
        RTI%nposterior_global=0
        RTI%nequals_global=0

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

        ! maximum logweight
        RTI%maxlogweight = logzero
        RTI%maxlogweight_dead = logzero
        RTI%maxlogweight_global = logzero

        RTI%thin_posterior = 0d0


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

        lognp = log( RTI%nlive(p) + 0d0 )
        lognp1= log( RTI%nlive(p) + 1d0 )
        lognp2= log( RTI%nlive(p) + 2d0 )

        ! Output the logweight
        logweight =  RTI%logXp(p) - lognp1

        ! Global evidence
        call logincexp( RTI%logZ,      RTI%logXp(p)+logL-lognp1  )
        ! Local evidence
        call logincexp( RTI%logZp(p) , RTI%logXp(p)+logL-lognp1  )
        ! Local volume
        RTI%logXp(p)  = RTI%logXp(p) + lognp - lognp1


        ! Global evidence error
        call logincexp( RTI%logZ2 ,                             &
            log2 + RTI%logZXp(p)     +   logL - lognp1,         &
            log2 + RTI%logXpXq(p,p)  + 2*logL - lognp1 - lognp2 &
            )

        ! global evidence volume cross correlation p=p
        RTI%logZXp(p) = RTI%logZXp(p) + lognp - lognp1
        call logincexp( RTI%logZXp(p),                       &
            RTI%logXpXq(p,p)+ logL + lognp - lognp1 - lognp2 &
            )

        ! global evidence volume cross correlation p/=q
        do q=1,RTI%ncluster
            if(p/=q) call logincexp( RTI%logZXp(q) , RTI%logXpXq(p,q)+ logL - lognp1 )
        end do


        ! Local evidence error
        call logincexp( RTI%logZp2(p),                          &
            log2 + RTI%logZpXp(p)    +   logL - lognp1,         &
            log2 + RTI%logXpXq(p,p)  + 2*logL - lognp1 - lognp2 &
            )


        ! Local evidence volume cross correlation
        RTI%logZpXp(p) = RTI%logZpXp(p) + lognp - lognp1
        call logincexp( RTI%logZpXp(p) ,                     &
            RTI%logXpXq(p,p)+ logL + lognp - lognp1 - lognp2 &
            )


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

    !> This function takes the evidence info in cluster p and splits it into
    !! several clusters according to the cluster numbers ni.
    !!
    !! It places these new clusters at the end of the active array, and deletes the old cluster
    !!
    subroutine add_cluster(settings,RTI,p,cluster_list,num_new_clusters) 
        use settings_module, only: program_settings
        use utils_module, only: logsumexp,logaddexp
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

        double precision, dimension(size(RTI%live,1),RTI%nlive(p))     :: old_live
        integer :: old_nlive
        double precision, dimension(size(RTI%phantom,1),size(RTI%phantom,2),size(RTI%phantom,3)) :: old_phantom
        integer, dimension(size(RTI%nphantom)) :: old_nphantom

        double precision, dimension(size(RTI%posterior,1),RTI%nposterior(p)) :: old_posterior
        integer :: old_nposterior

        double precision, dimension(size(RTI%equals,1),RTI%nequals(p)) :: old_equals
        integer :: old_nequals

        double precision :: old_maxlogweight

        double precision, dimension(num_new_clusters) :: logni
        double precision, dimension(num_new_clusters) :: logni1
        double precision :: logn
        double precision :: logn1
        double precision :: logXp
        double precision :: logZp
        double precision :: logZp2
        double precision, dimension(RTI%ncluster-1) :: logXpXq
        double precision :: logXp2
        double precision :: logZXp
        double precision :: logZpXp

        ! 1) Save the old points as necessary
        old_live  = RTI%live(:,:RTI%nlive(p),p)                 ! Save the old live points
        old_nlive = RTI%nlive(p)                                ! Save the old number of live points 
        old_posterior = RTI%posterior(:,:RTI%nposterior(p),p)   ! Save the old posterior points
        old_nposterior = RTI%nposterior(p)                      ! Save the old number of posterior points  
        old_equals = RTI%equals(:,:RTI%nequals(p),p)            ! Save the old equals points
        old_nequals = RTI%nequals(p)                            ! Save the old number of equals points  
        old_phantom = RTI%phantom                               ! Save the old phantom points
        old_nphantom= RTI%nphantom                              ! Save the old numbers of phantom points
        old_maxlogweight = RTI%maxlogweight(p)                  ! save the old max log weight for the relevant cluster

        old_save   = [(i,i=1,p-1),(i,i=p+1,RTI%ncluster)]    ! The indices of the old clusters to save
        num_old_clusters = RTI%ncluster-1                    ! The number of old clusters
        RTI%ncluster = RTI%ncluster+num_new_clusters-1       ! The total number of new clusters

        old_target = [(i,i=1,num_old_clusters)]              ! Where the old clusters will be
        new_target = [(i,i=num_old_clusters+1,RTI%ncluster)] ! Where we're going to insert the new clusters

        ! Define some useful variables
        logXp  = RTI%logXp(p)
        logXp2 = RTI%logXpXq(p,p)
        logZp  = RTI%logZp(p)
        logZp2 = RTI%logZp2(p)
        logZXp = RTI%logZXp(p)
        logZpXp= RTI%logZpXp(p)
        logXpXq= [ RTI%logXpXq(p,:p-1) , RTI%logXpXq(p,p+1:) ]



        ! 2) Reallocate the arrays

        ! Reallocate the live,phantom and posterior points
        call reallocate_3_d(RTI%live,            new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nlive,           new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)
        call reallocate_3_d(RTI%phantom,         new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nphantom,        new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)
        call reallocate_3_d(RTI%posterior_stack, new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nposterior_stack,new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)
        call reallocate_3_d(RTI%posterior,       new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nposterior,      new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)
        call reallocate_3_d(RTI%equals,          new_size3=RTI%ncluster, save_indices3=old_save,target_indices3=old_target)
        call reallocate_1_i(RTI%nequals,         new_size1=RTI%ncluster, save_indices1=old_save,target_indices1=old_target)

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

        call reallocate_1_d(RTI%maxlogweight,   RTI%ncluster,old_save,old_target) 


        ! 3) Assign the new live points to their new clusters

        RTI%nlive(new_target) = 0 ! Zero the number of live points in the new clusters
        do i_live=1,old_nlive
            ! Insert the new points in the correct positions
            call add_point(old_live(:,i_live),RTI%live,RTI%nlive,new_target(cluster_list(i_live)))
        end do

        ! Find the new minimum loglikelihoods
        call find_min_loglikelihoods(settings,RTI) 

        ! 4) Reassign the posterior points
        RTI%nposterior_stack(new_target) = 0

        RTI%nposterior(new_target) = old_nposterior
        RTI%posterior(:,:old_nposterior,new_target) = spread(old_posterior,dim=3,ncopies=num_new_clusters)

        RTI%nequals(new_target) = old_nequals
        RTI%equals(:,:old_nequals,new_target) = spread(old_equals,dim=3,ncopies=num_new_clusters)

        RTI%maxlogweight(new_target) = old_maxlogweight


        ! 4) Reassign all the phantom points 
        RTI%nphantom = 0
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
        ! Initialise local evidences 
        RTI%logZp(new_target) = logZp + logni - logn 
        RTI%logZp2(new_target) = logZp2 + logni + logni1 - logn - logn1
        RTI%logZpXp(new_target) = logZpXp + logni + logni1 - logn - logn1 


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




        ! reduce the logweighting of posterior points by a factor proportional to the evidence
        do i_cluster=1,num_new_clusters
            RTI%equals(settings%pos_l,:RTI%nequals(new_target(i_cluster)),new_target(i_cluster)) = RTI%equals(settings%pos_l,:RTI%nequals(new_target(i_cluster)),new_target(i_cluster)) + RTI%logZp(new_target(i_cluster)) -logZp
            RTI%posterior(settings%pos_l,:RTI%nposterior(new_target(i_cluster)),new_target(i_cluster)) = RTI%posterior(settings%pos_l,:RTI%nposterior(new_target(i_cluster)),new_target(i_cluster)) + RTI%logZp(new_target(i_cluster)) -logZp
        end do

    end subroutine add_cluster

    function delete_cluster(settings,RTI) 
        use settings_module, only: program_settings
        use array_module, only: reallocate_3_d,reallocate_2_d,reallocate_1_d,reallocate_1_i
        implicit none
        !> The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: RTI
        type(program_settings), intent(in) :: settings
        logical :: delete_cluster

        !The cluster index to be deleted
        integer            :: p(1)

        integer :: size_n

        ! new indices of clusters
        integer,dimension(RTI%ncluster-1) :: indices

        ! Contstructor iterator
        integer :: i

        delete_cluster=.false.

        if(any(RTI%nlive==0)) then

            delete_cluster=.true.

            ! Update the posterior arrays
            call update_posteriors(settings,RTI) 

            p=minloc(RTI%nlive,RTI%nlive==0)

            ! Get the positions of the clusters to be saved
            indices = [(i,i=1,p(1)-1),(i,i=p(1)+1,RTI%ncluster)]

            ! Reduce the number of clusters
            RTI%ncluster=RTI%ncluster-1
            ! Increase the number of dead clusters
            RTI%ncluster_dead = RTI%ncluster_dead + 1

            ! Reallocate the dead arrays
            size_n = max(size(RTI%posterior_dead,2),RTI%nposterior(p(1)))
            call reallocate_3_d(RTI%posterior_dead, new_size2=size_n,new_size3=RTI%ncluster_dead)
            call reallocate_1_i(RTI%nposterior_dead,new_size1=RTI%ncluster_dead)
            size_n = max(size(RTI%equals_dead,2),RTI%nequals(p(1)))
            call reallocate_3_d(RTI%equals_dead,    new_size2=size_n,new_size3=RTI%ncluster_dead)
            call reallocate_1_i(RTI%nequals_dead,   new_size1=RTI%ncluster_dead)
            call reallocate_1_d(RTI%logZp_dead,     new_size1=RTI%ncluster_dead)
            call reallocate_1_d(RTI%logZp2_dead,    new_size1=RTI%ncluster_dead)
            call reallocate_1_d(RTI%maxlogweight_dead,new_size1=RTI%ncluster_dead)

            ! Place the newly dead cluster into this
            RTI%nposterior_dead(RTI%ncluster_dead) = RTI%nposterior(p(1))
            RTI%posterior_dead(:,:RTI%nposterior_dead(RTI%ncluster_dead),RTI%ncluster_dead) = RTI%posterior(:,:RTI%nposterior(p(1)),p(1))
            RTI%nequals_dead(RTI%ncluster_dead) = RTI%nequals(p(1))
            RTI%equals_dead(:,:RTI%nequals_dead(RTI%ncluster_dead),RTI%ncluster_dead) = RTI%equals(:,:RTI%nequals(p(1)),p(1))
            RTI%logZp_dead(RTI%ncluster_dead)  = RTI%logZp(p(1))
            RTI%logZp2_dead(RTI%ncluster_dead) = RTI%logZp2(p(1))
            RTI%maxlogweight_dead(RTI%ncluster_dead) = RTI%maxlogweight(p(1))


            ! Reallocate the live,phantom and posterior points
            call reallocate_3_d(RTI%live,            new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nlive,           new_size1=RTI%ncluster, save_indices1=indices)
            call reallocate_3_d(RTI%phantom,         new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nphantom,        new_size1=RTI%ncluster, save_indices1=indices)
            call reallocate_3_d(RTI%posterior_stack, new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nposterior_stack,new_size1=RTI%ncluster, save_indices1=indices)
            call reallocate_3_d(RTI%posterior,       new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nposterior,      new_size1=RTI%ncluster, save_indices1=indices)
            call reallocate_3_d(RTI%equals,          new_size3=RTI%ncluster, save_indices3=indices)
            call reallocate_1_i(RTI%nequals,         new_size1=RTI%ncluster, save_indices1=indices)

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

            call reallocate_1_d(RTI%maxlogweight,   new_size1=RTI%ncluster,save_indices1=indices) 
        end if


    end function delete_cluster


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
    subroutine calculate_logZ_estimate(RTI,logZ,varlogZ,logZp,varlogZp,logZp_dead,varlogZp_dead)
        use utils_module, only: logzero
        implicit none

        type(run_time_info),intent(in)                                  :: RTI        !> Run time information
        double precision, intent(out)                                   :: logZ       !>
        double precision, intent(out)                                   :: varlogZ  !>
        double precision, intent(out), dimension(RTI%ncluster),optional :: logZp      !>
        double precision, intent(out), dimension(RTI%ncluster),optional :: varlogZp !>
        double precision, intent(out), dimension(RTI%ncluster_dead),optional :: logZp_dead      !>
        double precision, intent(out), dimension(RTI%ncluster_dead),optional :: varlogZp_dead !>

        logZ       = max(logzero,2*RTI%logZ - 0.5*RTI%logZ2)
        varlogZ    = RTI%logZ2 - 2*RTI%logZ

        if(present(logZp).and.present(varlogZp))then
            logZp      = max(logzero,2*RTI%logZp - 0.5*RTI%logZp2)
            varlogZp = RTI%logZp2 - 2*RTI%logZp
        end if

        if(present(logZp_dead).and.present(varlogZp_dead))then
            logZp_dead      = max(logzero,2*RTI%logZp_dead - 0.5*RTI%logZp2_dead)
            varlogZp_dead   = RTI%logZp2_dead - 2*RTI%logZp_dead
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
        double precision,intent(in),dimension(:,:),allocatable :: baby_points

        logical :: replaced ! Have we successfully replaced a point?

        ! live point, last of the baby points
        double precision,dimension(settings%nTotal) :: point

        double precision :: logL ! loglikelihood bound

        integer :: i_baby ! point iterator

        integer                                     :: cluster_del     ! cluster to delete from
        double precision,dimension(settings%nTotal) :: deleted_point   ! point we have just deleted
        double precision,dimension(settings%nposterior) :: posterior_point   ! temporary posterior point
        double precision                            :: logweight       ! The log weighting of this point
        

        ! The loglikelihood contour is defined by the cluster it belongs to
        logL = RTI%logLp(cluster_add)

        ! Assign the phantom points to cluster_add, if they are:
        ! (1) Within the isolikelihood contour of the cluster.
        ! (2) Within the voronoi cell of the cluster.

        do i_baby=1,size(baby_points,2)-1
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


                ! Calculate the posterior point and add it to the posterior stack
                posterior_point = calculate_posterior_point(settings,deleted_point,logweight,RTI%logZ,logsumexp(RTI%logXp))
                call add_point(posterior_point,RTI%posterior_stack,RTI%nposterior_stack,cluster_del )
                RTI%maxlogweight(cluster_del) = max(RTI%maxlogweight(cluster_del),posterior_point(settings%pos_w)+posterior_point(settings%pos_l))
                RTI%maxlogweight_global=max(RTI%maxlogweight_global,RTI%maxlogweight(cluster_del))

            else
                replaced = .false.                                  ! We haven't killed of any points
            end if
        else
            replaced = .false.
        end if

    end function replace_point

    subroutine clean_phantoms(settings,RTI)
        use utils_module, only: logsumexp,logincexp,minpos
        use settings_module, only: program_settings
        use calculate_module, only: calculate_posterior_point
        use random_module, only: bernoulli_trial
        use array_module, only: add_point,delete_point

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information

        double precision,dimension(settings%nTotal) :: deleted_point   ! point we have just deleted
        double precision,dimension(settings%nposterior) :: posterior_point   ! temporary posterior point
        
        integer :: i_cluster  ! cluster iterator
        integer :: i_phantom  ! phantom iterator
        integer :: i_stack(1) ! phantom iterator
        integer,dimension(size(RTI%nposterior_stack)) :: nposterior_stack ! original size of array (just live posteriors)

        nposterior_stack = RTI%nposterior_stack

        ! Now we delete the phantoms
        do i_cluster=1,RTI%ncluster
            i_phantom=1
            do while(i_phantom<=RTI%nphantom(i_cluster))

                ! Find the location of the live point which this point belongs to
                i_stack = minloc( RTI%posterior_stack(settings%pos_l,:nposterior_stack(i_cluster),i_cluster), &
                        RTI%posterior_stack(settings%pos_l,:nposterior_stack(i_cluster),i_cluster) &
                        > RTI%phantom(settings%l0,i_phantom,i_cluster)  )

                if(i_stack(1)<=0) then
                    i_phantom=i_phantom+1
                else

                    ! Delete this point
                    deleted_point = delete_point(i_phantom,RTI%phantom,RTI%nphantom,i_cluster)

                    ! Calculate the posterior point and add it to the array
                    if(settings%equals .or. settings%posteriors ) then
                        if(bernoulli_trial(RTI%thin_posterior)) then
                            posterior_point = calculate_posterior_point(settings,deleted_point,&
                                    RTI%posterior_stack(settings%pos_w,i_stack(1),i_cluster), &!logweight
                                    RTI%posterior_stack(settings%pos_Z,i_stack(1),i_cluster), &!RTI%logZ
                                    RTI%posterior_stack(settings%pos_X,i_stack(1),i_cluster)  &!logsumexp(RTI%logXp))
                                    )

                            call add_point(posterior_point,RTI%posterior_stack,RTI%nposterior_stack,i_cluster )
                            RTI%maxlogweight(i_cluster) = max(RTI%maxlogweight(i_cluster),posterior_point(settings%pos_w)+posterior_point(settings%pos_l))
                            RTI%maxlogweight_global=max(RTI%maxlogweight_global,RTI%maxlogweight(i_cluster))
                        end if
                    end if
                end if

            end do
        end do

    end subroutine clean_phantoms





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





    subroutine update_posteriors(settings,RTI)
        use settings_module, only: program_settings
        use random_module, only: bernoulli_trial
        use array_module, only: delete_point,add_point
        implicit none
        !> Program settings
        type(program_settings), intent(in) :: settings
        !> Run time information
        type(run_time_info),intent(inout) :: RTI


        double precision, dimension(settings%np) :: posterior_point

        integer :: i_post
        integer :: i_cluster


        ! Add in the phantoms from the stack (used to do this at every iteration, but this was the computational bottleneck)
        call clean_phantoms(settings,RTI) 

        if(settings%equals) then
            ! Clean the global equally weighted posteriors
            i_post=1
            do while(i_post<=RTI%nequals_global(1)) 
                ! We don't need to bother with points that have a weight equal to
                ! the max weight, these are automatically accepted
                if(RTI%equals_global(settings%p_w,i_post,1)<RTI%maxlogweight_global) then

                    ! accept with probability p = weight/maxweight
                    if(bernoulli_trial( exp(RTI%equals_global(settings%p_w,i_post,1)-RTI%maxlogweight(1)) )) then
                        ! If accepted, then their new probability is maxlogweight (for
                        ! the next round of stripping)
                        RTI%equals_global(settings%p_w,i_post,1) = RTI%maxlogweight_global
                        ! move on to the next point
                        i_post=i_post+1
                    else
                        ! delete this point if not accepted
                        posterior_point = delete_point(i_post,RTI%equals_global,RTI%nequals_global,1)
                    end if
                else
                    i_post=i_post+1 ! move on to the next point if automatically accepted
                end if
            end do

            ! Clean the local equally weighted posteriors
            if(settings%cluster_posteriors) then
                ! strip out old posterior points from equals
                do i_cluster=1,RTI%ncluster
                    i_post=1
                    do while(i_post<=RTI%nequals(i_cluster)) 
                        ! We don't need to bother with points that have a weight equal to
                        ! the max weight, these are automatically accepted
                        if(RTI%equals(settings%p_w,i_post,i_cluster)<RTI%maxlogweight(i_cluster)) then

                            ! accept with probability p = weight/maxweight
                            if(bernoulli_trial( exp(RTI%equals(settings%p_w,i_post,i_cluster)-RTI%maxlogweight(i_cluster)) )) then
                                ! If accepted, then their new probability is maxlogweight (for
                                ! the next round of stripping)
                                RTI%equals(settings%p_w,i_post,i_cluster) = RTI%maxlogweight(i_cluster) 
                                ! move on to the next point
                                i_post=i_post+1
                            else
                                ! delete this point if not accepted
                                posterior_point = delete_point(i_post,RTI%equals,RTI%nequals,i_cluster)
                            end if
                        else
                            i_post=i_post+1 ! move on to the next point if automatically accepted
                        end if
                    end do
                end do
            end if
        end if

        ! Add new posterior points from the stack
        do i_cluster=1,RTI%ncluster
            do i_post=1,RTI%nposterior_stack(i_cluster)

                if(settings%equals) then

                    ! Add the global equally weighted posteriors
                    if(bernoulli_trial( exp( RTI%posterior_stack(settings%pos_w,i_post,i_cluster) + RTI%posterior_stack(settings%pos_l,i_post,i_cluster) - RTI%maxlogweight_global) )) then
                        posterior_point(settings%p_w) = RTI%maxlogweight_global
                        posterior_point(settings%p_2l) = -2*RTI%posterior_stack(settings%pos_l,i_post,i_cluster)
                        posterior_point(settings%p_p0:settings%p_d1) = RTI%posterior_stack(settings%pos_p0:settings%pos_d1,i_post,i_cluster)
                        call add_point(posterior_point,RTI%equals_global,RTI%nequals_global,1)
                    end if

                    ! Add the clustered equally weighted posteriors
                    if(settings%cluster_posteriors) then
                        ! Test for acceptance of equally weighted posteriors
                        if(bernoulli_trial( exp( RTI%posterior_stack(settings%pos_w,i_post,i_cluster) + RTI%posterior_stack(settings%pos_l,i_post,i_cluster) - RTI%maxlogweight(i_cluster)) )) then
                            posterior_point(settings%p_w) = RTI%maxlogweight(i_cluster)
                            posterior_point(settings%p_2l) = -2*RTI%posterior_stack(settings%pos_l,i_post,i_cluster)
                            posterior_point(settings%p_p0:settings%p_d1) = RTI%posterior_stack(settings%pos_p0:settings%pos_d1,i_post,i_cluster)
                            call add_point(posterior_point,RTI%equals,RTI%nequals,i_cluster)
                        end if
                    end if

                end if

                if(settings%posteriors) then
                    ! Add point to global weighted posterior array
                    call add_point(RTI%posterior_stack(:,i_post,i_cluster),RTI%posterior_global,RTI%nposterior_global,1)
                    ! Add point to cluster weighted posterior array
                    if(settings%cluster_posteriors) call add_point(RTI%posterior_stack(:,i_post,i_cluster),RTI%posterior,RTI%nposterior,i_cluster)
                end if
            end do
        end do

        RTI%nposterior_stack=0 ! reset the stack to zero

    end subroutine update_posteriors



end module
