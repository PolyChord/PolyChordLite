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
            integer :: nclusters
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
            double precision, allocatable, dimension(:,:,:) :: covmats
            !> Cholesky decompositions
            double precision, allocatable, dimension(:,:,:) :: choleskys

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

        end type run_time_info

    contains

        !> This is a self explanatory subroutine.
        !!
        !! It allocates the arrays for a single cluster 
        subroutine allocate_run_time_info(settings,RTI)
            use utils_module,    only: logzero
            use settings_module, only: program_settings

            implicit none
            type(program_settings), intent(in) :: settings
            type(run_time_info),intent(out) :: RTI

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
                RTI%nlive(1),                                    &
                RTI%nphantom(1),                                 &
                RTI%nposterior(1)                                &
                )

            ! All evidences set to logzero
            RTI%logZ=logzero
            RTI%logZ2=logzero
            RTI%logXp=1
            RTI%logXpXq=1
            RTI%logZp=logzero
            RTI%logZXp=logzero
            RTI%logZp2=logzero
            RTI%logZpXp=logzero


        end subroutine allocate_run_time_info

    subroutine update_evidence(RTI,p,logL)
        use utils_module, only: logsumexp,logincexp
        implicit none

        ! The variable containing all of the runtime information
        type(run_time_info), intent(inout) :: RTI

        !> The cluster index to update
        integer,intent(in) :: p
        !> The loglikelihood to update
        double precision,intent(in) :: logL

        ! Iterator
        integer :: q

        ! Temporary variables for notational ease
        double precision,parameter :: log2 = log(2d0)
        double precision :: lognp
        double precision :: lognp1
        double precision :: lognp2


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
        do q=1,RTI%nclusters
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
        do q=1,RTI%nclusters
            if(p/=q) then
                RTI%logXpXq(p,q) = RTI%logXpXq(p,q) + lognp - lognp1
                RTI%logXpXq(q,p) = RTI%logXpXq(q,p) + lognp - lognp1
            end if
        end do

        ! Decrease the number of live points in cluster p
        RTI%nlive(p) = RTI%nlive(p) - 1


    end subroutine update_evidence


end module
