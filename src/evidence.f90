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
        double precision, save :: Z_dead=0

        ! The accumulated evidence squared Z^2
        ! Z2_dead = 1/nlive  sum_k L_k  (nlive/(nlive+1))^k * Z2_dead_part_k
        double precision, save :: Z2_dead=0

        ! part of the accumulated evidence squared
        ! Z2_dead_part = 1/(nlive+1)  sum_k L_k  ((nlive+1)/(nlive+2))^k
        double precision, save :: Z2_dead_part=0

        ! the mean log likelihood of the live points
        ! mean_log_like_live = 1/nlive * sum(L_live)
        double precision, save :: mean_log_like_live=0

        ! Two measures of volume at the kth iteration 
        ! logX_k = (  nlive  /(nlive+1))^k
        ! logY_k = ((nlive+1)/(nlive+2))^k
        double precision :: logX_k,logY_k

        ! The temporary mean and variance to eventually be output in evidence_vec 
        double precision :: mean, variance

        ! number of live points in double precision (taken from settings%nlive)
        double precision :: nlive 




        ! If the function is called with a non-positive ndead, then simply store
        ! the mean_log_likelihood that has been calculated and passed via
        ! new_loglikelihood
        if (ndead <= 0) then
            ! save the average value
            mean_log_like_live = new_loglikelihood
            more_samples_needed = .true.
            return
        end if

        ! Get the number of live points (and convert to double precision implicitly)
        nlive = settings%nlive

        ! Calculate the volumes
        logX_k = ndead * log(   nlive  /(nlive+1) )
        logY_k = ndead * log( (nlive+1)/(nlive+2) )


        ! Add the evidence contribution of the kth dead point to the accumulated dead evidence
        Z_dead = Z_dead             + exp( old_loglikelihood + logX_k ) / nlive

        ! Accumulate part of the dead evidence^2
        Z2_dead_part = Z2_dead_part + exp( old_loglikelihood + logY_k ) /(nlive+1) 

        ! Accumulate the rest of the dead evidence^2
        Z2_dead = Z2_dead           + exp( old_loglikelihood + logX_k ) / nlive * Z2_dead_part

        ! Take away the contribution of the old_loglikelihood from the mean ...
        mean_log_like_live = mean_log_like_live + - exp(old_loglikelihood) / nlive
        ! .. and add the contribution of the new_loglikelihood from the mean
        mean_log_like_live = mean_log_like_live +   exp(new_loglikelihood) / nlive


        ! Calculate the mean evidence as a sum of contributions from the dead
        ! and live evidence
        mean =        Z_dead 
        mean = mean + mean_log_like_live * exp(logX_k)

        ! Calculate the variance in the evidence
        variance =            2*Z2_dead - Z_dead**2                                          ! dead contribution
        variance = variance + mean_log_like_live**2 * ( exp(logX_k+logY_k) - exp(2*logX_k) ) ! live contribution
        variance = variance + 2*mean_log_like_live * exp(2*logX_k) * ( Z2_dead_part-Z_dead)  ! cross correlation

        ! Pass the mean and variance to evidence_vec in order to be outputted
        evidence_vec(1) = mean
        evidence_vec(2) = sqrt(abs(variance))


        ! Test to see whether we have enough samples. If the contribution from
        ! the live evidence is less than a fraction of the dead evidence, then
        ! we have reached convergence of the algorithm
        if (mean_log_like_live  * exp(logX_k) < settings%precision_criterion * Z_dead) then
            more_samples_needed = .false.
        else
            more_samples_needed = .true.
        end if



    end function KeetonEvidence


end module evidence_module
