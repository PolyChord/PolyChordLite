!> This module contains tools to calculate evidence estimators
module evidence_module
    implicit none

    contains


    !> Evidence calculator based on Keeton's paper ([arXiv:1102.0996](http://arxiv.org/abs/1102.0996))
    !!
    !! The evidence calculator recieves the relevent information, namely 
    !! * the loglikelihood of the newly created point ( new_loglikelihood )
    !! * the loglikelihood of the dying point  ( old_loglikelihood )
    !! * number of iterations/dead points ( ndead )
    !!
    !! It ouputs 
    !! * a length 2 vector ( evidence_vec ) with the [evidence, evidence error] in the value of the function
    !! * whether more samples are needed in the logical variable more_samples_needed
    !!
    !! The Keeton evidence calculates using the algorithm outlined in ([arXiv:1102.0996](http://arxiv.org/abs/1102.0996)),
    !! using a moment-based error analysis.
    !!
    !! The algorithm terminates when the contribution from the live evidence is
    !! less than a fraction of the contribution from the dead evidence
    function KeetonEvidence(settings,new_loglikelihood,old_loglikelihood,ndead,more_samples_needed) result (evidence_vec)
        use settings_module
        use utils_module, only: logzero

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

        ! ------- Outputs ------- 
        !> Whether we have obtained enough samples for an accurate evidence
        logical,intent(out) :: more_samples_needed

        ! vector containing [evidence, evidence error]
        double precision, dimension(2) :: evidence_vec


        ! ------- Local Variables -------

        ! The accumulated evidence associated with the dead points
        ! Z_dead = 1/nlive  sum_k L_k  (nlive/(nlive+1))^k
        double precision, save :: logZ_dead = logzero

        ! The accumulated evidence squared Z^2
        ! Z2_dead = 1/nlive  sum_k L_k  (nlive/(nlive+1))^k * Z2_dead_part_k
        double precision, save :: logZ2_dead = logzero

        ! part of the accumulated evidence squared
        ! Z2_dead_part = 1/(nlive+1)  sum_k L_k  ((nlive+1)/(nlive+2))^k
        double precision, save :: logZ2_dead_part = logzero

        ! the mean log likelihood of the live points
        ! mean_log_like_live = 1/nlive * sum(L_live)
        double precision, save :: logmean_loglike_live=0

        ! Two measures of volume at the kth iteration 
        ! logX_k = (  nlive  /(nlive+1))^k
        ! logY_k = ((nlive+1)/(nlive+2))^k
        double precision :: logX_k,logY_k

        ! The temporary mean and variance to eventually be output in evidence_vec 
        double precision :: logmean, logvariance

        ! number of live points in double precision (taken from settings%nlive)
        double precision :: nlive 


        ! Throughout the calculation, we measure all likelihoods relative to
        double precision :: max_loglikelihood


        ! If the function is called with a non-positive ndead, then simply store
        ! the mean_log_likelihood that has been calculated and passed via
        ! new_loglikelihood
        if (ndead <= 0) then
            ! Get the current maximum loglikelihood value
            max_loglikelihood = new_loglikelihood
            ! save the average value
            logmean_loglike_live = old_loglikelihood
            more_samples_needed = .true.
            return
        end if

        max_loglikelihood = max(new_loglikelihood,max_loglikelihood)


        ! Get the number of live points (and convert to double precision implicitly)
        nlive = settings%nlive

        logZ_dead =             logZ_dead            - max_loglikelihood
        logZ2_dead =            logZ2_dead           - 2*max_loglikelihood
        logZ2_dead_part =       logZ2_dead_part      - max_loglikelihood
        logmean_loglike_live =  logmean_loglike_live - max_loglikelihood

        ! Calculate the volumes
        logX_k = ndead * log(   nlive  /(nlive+1) )
        logY_k = ndead * log( (nlive+1)/(nlive+2) )


        ! Add the evidence contribution of the kth dead point to the accumulated dead evidence
        logZ_dead = log ( exp(logZ_dead)  + exp( old_loglikelihood-max_loglikelihood + logX_k ) / nlive )

        ! Accumulate part of the dead evidence^2
        logZ2_dead_part = log( exp(logZ2_dead_part) + exp( old_loglikelihood-max_loglikelihood + logY_k ) /(nlive+1)  )

        ! Accumulate the rest of the dead evidence^2
        logZ2_dead = log( exp(logZ2_dead) + exp( old_loglikelihood-max_loglikelihood + logX_k + logZ2_dead_part ) / nlive )


        ! Take away the contribution of the 0d0 from the mean ...
        logmean_loglike_live = log ( exp(logmean_loglike_live) - exp(old_loglikelihood-max_loglikelihood) / nlive )
        ! .. and add the contribution of the new_loglikelihood from the mean
        logmean_loglike_live = log ( exp(logmean_loglike_live) + exp(new_loglikelihood-max_loglikelihood) / nlive )


        ! Calculate the mean evidence as a sum of contributions from the dead
        ! and live evidence
        logmean = log( exp(logZ_dead) + exp(logmean_loglike_live + logX_k) )

        ! Calculate the variance in the evidence
        logvariance = 2*exp(logZ2_dead) - exp(2*logZ_dead)                                      ! dead contribution
        logvariance = logvariance + exp(2*logmean_loglike_live + logX_k+logY_k) - exp(2*logX_k) ! live contribution
        logvariance = log( logvariance + 2*exp( logmean_loglike_live + 2*logX_k) * ( exp(logZ2_dead_part)-exp(logZ_dead)) )  ! cross correlation

        ! Pass the mean and variance to evidence_vec in order to be outputted
        evidence_vec(1) = logmean + old_loglikelihood
        evidence_vec(2) = logvariance + 2*old_loglikelihood

        logZ_dead =             logZ_dead            + max_loglikelihood
        logZ2_dead =            logZ2_dead           + 2*max_loglikelihood
        logZ2_dead_part =       logZ2_dead_part      + max_loglikelihood
        logmean_loglike_live =  logmean_loglike_live + max_loglikelihood


        ! Test to see whether we have enough samples. If the contribution from
        ! the live evidence is less than a fraction of the dead evidence, then
        ! we have reached convergence of the algorithm
        if (logmean_loglike_live  + logX_k < log(settings%precision_criterion) + logZ_dead) then
            more_samples_needed = .false.
        else
            more_samples_needed = .true.
        end if



    end function KeetonEvidence


end module evidence_module
