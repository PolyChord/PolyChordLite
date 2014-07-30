!> This module contains tools to calculate evidence estimators
module evidence_module
    implicit none

    contains


    subroutine KeetonEvidence(settings,new_loglikelihood,old_loglikelihood,ndead,evidence_vec)
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

        ! ------- Local Variables -------
        ! Temporary names for now, will think of better ones shortly
        double precision, save :: evA=0
        double precision, save :: evB=0
        double precision, save :: evC=0
        double precision, save :: evD=0

        integer :: nlive ! number of live points (taken from settings%nlive)





        nlive = settings%nlive ! get the number of live points

        evA = evA + exp( old_loglikelihood + ndead * log(dble(nlive)  /dble(nlive+1)) ) / dble(nlive)

        evB = evB + exp( old_loglikelihood + ndead * log(dble(nlive+1)/dble(nlive+2)) ) / dble(nlive+1) 

        evC = evC + exp( old_loglikelihood + ndead * log(dble(nlive)  /dble(nlive+1)) ) / dble(nlive) * evB

        !evD = evD + ( exp(new_likelihood) - exp(late_likelihood) )/dble(nlive)
        evD = evD + ( exp(new_loglikelihood) - exp(old_loglikelihood) )/dble(nlive)


        evidence_vec(1) = evA + evD * exp( ndead * log(dble(nlive)/dble(nlive+1)) )
        evidence_vec(2) = evD 



    end subroutine KeetonEvidence


end module evidence_module
