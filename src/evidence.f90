!> This module contains tools to calculate evidence estimators
module evidence_module
    implicit none

    contains


    function KeetonEvidence(settings,new_loglikelihood,old_loglikelihood,ndead,evidence_vec) result (more_samples_needed)
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
        !> vector containing [evidence, evidence error]
        double precision, intent(out), dimension(2) :: evidence_vec

        !> Whether we have obtained enough samples for an accurate evidence
        logical :: more_samples_needed

        ! ------- Local Variables -------

        ! The accumulated evidence associated with the dead points
        double precision, save :: Z_dead=0

        ! The accumulated evidence squared Z^2
        double precision, save :: Z2_dead=0

        ! part of the accumulated evidence squared
        double precision, save :: Z2_dead_part=0

        ! the mean log likelihood of the live points
        double precision, save :: mean_log_like_live=0

        double precision :: nlive ! number of live points in double precision (taken from settings%nlive)

        double precision :: logX_k,logY_k

        double precision :: mean, variance


        if (ndead <= 0) then
            ! save the average value
            mean_log_like_live = new_loglikelihood
            more_samples_needed = .true.
            return
        end if
        



        nlive = settings%nlive ! get the number of live points (and convert to double precision implicitly)



        logX_k = ndead * log(   nlive  /(nlive+1) )
        logY_k = ndead * log( (nlive+1)/(nlive+2) )


        Z_dead = Z_dead             + exp( old_loglikelihood + logX_k ) / nlive

        Z2_dead_part = Z2_dead_part + exp( old_loglikelihood + logY_k ) /(nlive+1) 

        Z2_dead = Z2_dead           + exp( old_loglikelihood + logX_k ) / nlive * Z2_dead_part

        mean_log_like_live = mean_log_like_live + ( exp(new_loglikelihood) - exp(old_loglikelihood) )/ nlive


        mean =        Z_dead 
        mean = mean + mean_log_like_live * exp(logX_k)

        variance =            2*Z2_dead - Z_dead**2 
        variance = variance + mean_log_like_live**2 * ( exp(logX_k+logY_k) - exp(2*logX_k) )
        variance = variance + 2*mean_log_like_live * exp(2*logX_k) * ( Z2_dead_part-Z_dead)  



        evidence_vec(1) = mean
        evidence_vec(2) = sqrt(abs(variance))

        


        if (mean_log_like_live  * exp(logX_k) < settings%precision_criterion * Z_dead) then
            more_samples_needed = .false.
        else
            more_samples_needed = .true.
        end if



    end function KeetonEvidence


end module evidence_module
