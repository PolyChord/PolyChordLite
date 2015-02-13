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
        double precision, allocatable :: logL
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

        ! Decrease the number of live points in cluster p
        RTI%nlive(p) = RTI%nlive(p) - 1


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


    subroutine find_min_loglike(settings,RTI)
        use utils_module, only: loginf
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information

        ! The minimum loglikelihood to be returned
        double precision :: logL

        integer :: i_live(1)

        integer :: i_cluster !> cluster iterator

        RTI%logL = loginf

        ! Find the separate minima of all the clusters
        do i_cluster=1,RTI%ncluster

            ! Find the position of the lowest point in this cluster
            i_live = minloc(RTI%live(settings%l0,:RTI%nlive(i_cluster),i_cluster))
            ! Find the likelihood of the lowest point in this cluster
            logL = RTI%live(settings%l0,i_live(1),i_cluster) 

            ! If this is the lowest likelihood we've found 
            if(logL<RTI%logL) then
                RTI%logL = logL      !> Record that we've found a new low
                RTI%p    = i_cluster !> Record the cluster identity
                RTI%i    = i_live(1) !> Record the live point identity
            end if

        end do


    end subroutine find_min_loglike

    function live_logZ(settings,RTI)
        use utils_module, only: logzero,logsumexp,logincexp
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings !> Program settings
        type(run_time_info),intent(inout)  :: RTI      !> Run time information

        double precision ::live_logZ !> Amount of evidence remaining in the live points

        integer :: i_cluster

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

end module
