!> This module contains tools to calculate evidence estimators
module evidence_module
    implicit none

    contains


    !> Moments-based Evidence calculator based on Keeton's paper ([arXiv:1102.0996](http://arxiv.org/abs/1102.0996))
    !!
    !! The evidence calculator recieves the relevent information, namely 
    !! * the loglikelihood of the newly created point ( new_loglikelihood )
    !! * the loglikelihood of the dying point  ( old_loglikelihood )
    !! * number of iterations/dead points ( ndead )
    !!
    !! It ouputs 
    !! * a length 6 vector ( evidence_vec ) with the [logevidence, logevidence error] in the value of the first two positions
    !!   (the remainder of the vector contains data which needs to be accumulated in order to make the calculation)
    !! * whether more samples are needed in the logical variable more_samples_needed
    !!
    !! The Keeton evidence calculates using the algorithm outlined in ([arXiv:1102.0996](http://arxiv.org/abs/1102.0996)),
    !! using a moment-based error analysis.
    !!
    !! Let \f$M\f$ be the number of live points and \f$N\f$ be the number of dead points.
    !! Effectively we wish to evaluate the evidence composed of two contributions:
    !! \f[ Z = Z_\mathrm{dead} + Z_\mathrm{live}\f]
    !!
    !! where:
    !!
    !! \f[ Z_\mathrm{dead} = \frac{1}{M} \sum_{i=1}^{N} L_i\left(\frac{M}{M+1}\right)^{i} \qquad
    !!     Z_\mathrm{live} = \frac{1}{M} \left(\frac{M}{M+1}\right)^{N}\sum_{i=1}^{M} L_i^\mathrm{live} = X_N\langle L^\mathrm{live}\rangle \f]
    !!
    !! In addition, we wish to calculate the error, which can be expressed as:
    !! \f[ \sigma_Z^2 =  \sigma_\mathrm{dead}^2 +\sigma_\mathrm{cross}^2 + \sigma_\mathrm{live}^2 \f]
    !! where
    !! \f[ \sigma_\mathrm{dead}^2 = \frac{2}{M(M+1)} \sum_{k=1}^{N} L_k \left(\frac{M}{M+1}\right)^{k} \sum_{i=1}^{k} L_i \left(\frac{M+1}{M+2}\right)^{i} - \frac{1}{M^2} \left[\sum_{i=1}^{N} L_i \left(\frac{M}{M+1}\right)^{i}\right]^2 \f]
    !! \f[ \sigma_\mathrm{cross}^2 = 2 \langle L^\mathrm{live}\rangle \left(\frac{M}{M+1}\right)^N \sum_{i=1}^{N} L_i \left[ \frac{1}{M+1}\left(\frac{M+1}{M+2}\right)^i - \frac{1}{M}\left(\frac{M}{M+1}\right)^i \right]  \f]
    !! \f[ \sigma_\mathrm{live}^2 = \langle L^\mathrm{live}\rangle^2 \left(\frac{M}{M+1}\right)^N \left[ \left(\frac{M+1}{M+2}\right)^N - \left(\frac{M}{M+1}\right)^N \right]  \f]
    !!
    !! These rather impenetrable expressions may be simplified by defining:
    !! \f[ X_k = \left(\frac{M}{M+1}\right)^k \qquad
    !!  Y_k = \left(\frac{M+1}{M+2}\right)^k \qquad
    !!  Z_\mathrm{dead}^\prime  = \frac{1}{M+1} \sum_{i=1}^{N} L_i Y_i \qquad
    !!  \hat{Z}_\mathrm{dead} = \sum_{k=1}^{N} \frac{1}{M}L_k X_k \sum_{i=1}^{k} \frac{1}{M+1}L_i Y_i \f]
    !! so that they become
    !! \f[ \sigma_\mathrm{dead}^2 = 2\hat{Z}_\mathrm{dead} - Z_\mathrm{dead}^2 \f]
    !! \f[ \sigma_\mathrm{cross}^2 = 2 \langle L^\mathrm{live}\rangle X_N (Z_\mathrm{dead}^\prime-Z_\mathrm{dead}) \f]
    !! \f[ \sigma_\mathrm{live}^2 = \langle L^\mathrm{live}\rangle^2 X_N \left[ Y_N - X_N \right]  \f]
    !! 
    !! We can update these accordingly at each iteration:
    !! \f[ Z_\mathrm{live} \rightarrow Z_\mathrm{live} + X_M\frac{L_\mathrm{new}}{M} - X_M\frac{L_\mathrm{old}}{M} \f]
    !! \f[ Z_\mathrm{dead} \rightarrow Z_\mathrm{dead} + \frac{1}{M}L_\mathrm{old}X_N \f]
    !! \f[ Z_\mathrm{dead}^\prime \rightarrow Z_\mathrm{dead}^\prime + \frac{1}{M+1}L_\mathrm{old}Y_N \f]
    !! \f[ \hat{Z}_\mathrm{dead} \rightarrow \hat{Z}_\mathrm{dead} + \frac{1}{M}L_\mathrm{old}X_N Z_\mathrm{dead}^\prime\f]
    !! and ensure that we save both the dead and live evidences for future use. The dead evidence is stored in 
    !! evidence_vec(3), and the live evidence is stored in evidence_vec(6)
    !!
    !! Also note that in order to prevent floating point errors, all of these are done with the logarithmic values
    !!
    subroutine KeetonEvidence(settings,new_loglikelihood,old_loglikelihood,ndead,more_samples_needed, evidence_vec)
        use settings_module
        use utils_module, only: logzero,logaddexp,logsubexp

        implicit none

        ! ------- Inputs ------- 
        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> loglikelihood of the newest created point
        double precision,       intent(in) :: new_loglikelihood

        !> loglikelihood of the most recently dead point
        double precision,       intent(in) :: old_loglikelihood

        !> number of dead points/ number of iterations
        integer,                intent(in) :: ndead

        !> vector containing [logevidence, logevidence error,logZ_dead,logZhat_dead,logZp_dead,logmean_loglike_live]
        double precision, intent(inout), dimension(6) :: evidence_vec

        ! ------- Outputs ------- 
        !> Whether we have obtained enough samples for an accurate evidence
        logical,intent(out) :: more_samples_needed



        ! ------- Local Variables -------

        ! The accumulated evidence associated with the dead points
        ! Z_dead = 1/nlive  sum_k L_k  (nlive/(nlive+1))^k
        double precision :: logZ_dead = logzero

        ! The accumulated evidence squared Z^2
        ! Z2_dead = 1/nlive  sum_k L_k  (nlive/(nlive+1))^k * Z2_deadp_k
        double precision :: logZhat_dead = logzero

        ! part of the accumulated evidence squared
        ! Z2_deadp = 1/(nlive+1)  sum_k L_k  ((nlive+1)/(nlive+2))^k
        double precision :: logZp_dead = logzero

        ! the mean log likelihood of the live points
        ! mean_log_like_live = 1/nlive * sum(L_live)
        double precision :: logmean_loglike_live=0

        ! Two measures of volume at the kth iteration 
        ! logX_k = (  nlive  /(nlive+1))^k
        ! logY_k = ((nlive+1)/(nlive+2))^k
        double precision :: logX_k,logY_k

        ! The temporary mean and variance to eventually be output in evidence_vec 
        double precision :: logmean, logvariance

        ! number of live points in double precision (taken from settings%nlive)
        double precision :: lognlive 
        double precision :: lognlivep1 
        double precision :: lognlivep2 

        double precision,parameter :: log2 = log(2d0)


        ! Get the number of live points (and convert to double precision implicitly)
        lognlive   = log(settings%nlive+0d0)
        lognlivep1 = log(settings%nlive+1d0)
        lognlivep2 = log(settings%nlive+2d0)


        ! Calculate the volumes
        logX_k = ndead * (lognlive  -lognlivep1)
        logY_k = ndead * (lognlivep1-lognlivep2) 

        ! Grab the accumulated values
        logZ_dead           = evidence_vec(3)
        logmean_loglike_live= evidence_vec(4)
        logZhat_dead          = evidence_vec(5)
        logZp_dead         = evidence_vec(6)

        ! Add the evidence contribution of the kth dead point to the accumulated dead evidence
        logZ_dead = logaddexp(logZ_dead , old_loglikelihood+ logX_k -lognlive )

        ! Accumulate part of the dead evidence^2
        logZp_dead = logaddexp(logZp_dead, old_loglikelihood + logY_k -lognlivep1  )

        ! Accumulate the rest of the dead evidence^2
        logZhat_dead = logaddexp(logZhat_dead, old_loglikelihood + logX_k + logZp_dead- lognlive )


        ! Add the contribution of the new_loglikelihood from the mean ...
        logmean_loglike_live = logaddexp(logmean_loglike_live, new_loglikelihood -lognlive )

        ! ... and take away the contribution of the 0d0 from the mean
        logmean_loglike_live = logsubexp(logmean_loglike_live, old_loglikelihood -lognlive )



        ! Calculate the mean evidence as a sum of contributions from the dead
        ! and live evidence
        logmean = logaddexp(logZ_dead, logmean_loglike_live + logX_k) 

        ! Calculate the variance in the evidence
        logvariance = logsubexp(log2 +logZhat_dead,2*logZ_dead)                                         ! dead contribution

        logvariance = logaddexp(logvariance,2*logmean_loglike_live + logX_k+ logsubexp(logY_k,logX_k))    ! live contribution

        logvariance = logaddexp(logvariance, log2 + logmean_loglike_live + logX_k + logZp_dead ) ! cross correlation
        logvariance = logsubexp(logvariance, log2 + logmean_loglike_live + logX_k + logZ_dead       )

        ! Pass the mean and variance to evidence_vec in order to be outputted
        evidence_vec(1) = logmean 
        evidence_vec(2) = logvariance 
        evidence_vec(3) = logZ_dead            
        evidence_vec(4) = logmean_loglike_live 
        evidence_vec(5) = logZhat_dead           
        evidence_vec(6) = logZp_dead          


        ! Test to see whether we have enough samples. If the contribution from
        ! the live evidence is less than a fraction of the dead evidence, then
        ! we have reached convergence of the algorithm
        if (logmean_loglike_live  + logX_k < log(settings%precision_criterion) + logZ_dead) then
            more_samples_needed = .false.
        else
            more_samples_needed = .true.
        end if



    end subroutine KeetonEvidence


    subroutine infer_evidence(settings,loglikelihoods)
        use settings_module, only: program_settings
        use utils_module, only: write_ev_unit,DBL_FMT
#ifdef MPI
        use mpi_module
#endif
        implicit none
        !> Program settings
        type(program_settings), intent(in) :: settings
        !> All of the likelihoods
        double precision, intent(in), dimension(:) :: loglikelihoods

        double precision, dimension(settings%evidence_samples) :: log_evidences

        double precision, dimension(2,settings%evidence_samples-2) :: distribution

        integer :: info,i_err
        integer :: i

#ifdef MPI
        double precision, allocatable, dimension(:) :: log_evidences_local
        integer :: evidence_samples_local
        integer :: nprocs

        nprocs = mpi_size()  ! Get the number of MPI procedures
        write(*,*) 'i am here'

        evidence_samples_local = ceiling(settings%evidence_samples/(nprocs+0d0))

        ! Allocate the arrays for storing log evidences on each node
        allocate(log_evidences_local(evidence_samples_local))
        
        ! Calculate the log evidences on each node
        do i=1,evidence_samples_local
            log_evidences_local(i) = sample_logevidence(loglikelihoods,settings%nlive)
        end do

        ! Gather them onto the root node
        call MPI_GATHER(                 &  
            log_evidences_local,         & ! sending array
            evidence_samples_local,      & ! number of elements to be sent
            MPI_DOUBLE_PRECISION,        & ! type of element to be sent
            log_evidences,               & ! recieving array
            evidence_samples_local,      & ! number of elements to be recieved from each node
            MPI_DOUBLE_PRECISION,        & ! type of element recieved
            0,                           & ! root node address
            MPI_COMM_WORLD,              & ! communication info
            mpierror)                      ! error (from module mpi_module)

#else
        ! Compute evidence_samples volume samples of the log evidence
        do i=1,settings%evidence_samples
            log_evidences(i) = sample_logevidence(loglikelihoods,settings%nlive)
        end do
#endif

        ! Sort them in increasing order 
        call dlasrt('I',settings%evidence_samples,log_evidences,info)

        ! Create a distribution
        distribution(1,:) = log_evidences(2:)
        distribution(2,:) = 2d0/(log_evidences(3:) - log_evidences(:settings%evidence_samples-2))

        ! write to file
        open(write_ev_unit,file=trim(settings%file_root) // '_evidence.txt' , action='write', iostat=i_err) 
        write(write_ev_unit,'(2E<DBL_FMT(1)>.<DBL_FMT(2)>)') distribution
        close(write_ev_unit)


    end subroutine infer_evidence

    function sample_logevidence(loglikelihoods,nlive) result(logevidence)
        use utils_module, only: logsubexp,logsumexp
        use random_module, only: random_reals
        implicit none
        !> All of the likelihoods
        double precision, intent(in), dimension(:) :: loglikelihoods
        !> the number of live points
        integer, intent(in) :: nlive

        ! The result
        double precision :: logevidence

        double precision, dimension(0:size(loglikelihoods)) :: logvolumes
        double precision, dimension(size(loglikelihoods)) :: logweights

        integer :: nlike
        integer :: i

        ! Get the number of likelihoods
        nlike = size(loglikelihoods)

        ! Generate nlike random numbers distributed as P(t) = nlive t^(nlive-1)
        ! take the log to preserve accuracy
        logvolumes(0)=0
        logvolumes(1:) = log(random_reals(nlike))/nlive

        ! shrink each volume by a factor of the previous volume
        do i=1,nlike
            logvolumes(i) = logvolumes(i)+logvolumes(i-1)
        end do

        ! calculate the weights
        do i=1,nlike
            logweights(i) = logsubexp( logvolumes(i-1) , logvolumes(i) )
        end do

        ! Compute the evidence
        logevidence = logsumexp(loglikelihoods+logweights)
        
    end function sample_logevidence


end module evidence_module
