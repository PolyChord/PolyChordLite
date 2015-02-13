module calculate_module
    implicit none
    contains

    subroutine calculate_point(loglikelihood,priors,point,settings,nlike)
        use priors_module, only: prior, hypercube_to_physical
        use settings_module, only: program_settings
        use utils_module, only: logzero
        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        type(prior), dimension(:), intent(in) :: priors
        type(program_settings), intent(in) :: settings
        double precision, intent(inout) , dimension(:) :: point
        integer, intent(inout) :: nlike

        if ( any(point(settings%h0:settings%h1)<0d0) .or. any(point(settings%h0:settings%h1)>1d0) )  then
            point(settings%p0:settings%p1) = 0
            point(settings%l0) = logzero
        else
            ! Transform the the hypercube coordinates to the physical coordinates
            point(settings%p0:settings%p1) = hypercube_to_physical( point(settings%h0:settings%h1),priors )

            ! Calculate the likelihood and store it in the last index
            point(settings%l0) = loglikelihood( point(settings%p0:settings%p1), point(settings%d0:settings%d1),settings%context)

            ! accumulate the number of likelihood calls that we've made
            nlike = nlike+1
        end if

    end subroutine calculate_point


end module calculate_module
