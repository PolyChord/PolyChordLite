module calculate_module
    implicit none
    contains

    subroutine calculate_point(loglikelihood,priors,live_data,settings)
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
        double precision, intent(inout) , dimension(:) :: live_data

        if ( any(live_data(settings%h0:settings%h1)<0d0) .or. any(live_data(settings%h0:settings%h1)>1d0) )  then
            live_data(settings%p0:settings%p1) = 0
            live_data(settings%l0) = logzero
        else
            ! Transform the the hypercube coordinates to the physical coordinates
            live_data(settings%p0:settings%p1) = hypercube_to_physical( live_data(settings%h0:settings%h1),priors )

            ! Calculate the likelihood and store it in the last index
            live_data(settings%l0) = loglikelihood( live_data(settings%p0:settings%p1), live_data(settings%d0:settings%d1),settings%context)

            ! accumulate the number of likelihood calls that we've made
            live_data(settings%nlike) = live_data(settings%nlike)+1
        end if

    end subroutine calculate_point

    function gradloglike(loglikelihood,settings,theta,loglike,delta)
        use settings_module, only: program_settings
        implicit none
        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface
        type(program_settings), intent(in) :: settings
        double precision, intent(in), dimension(settings%nDims) :: theta
        double precision, intent(in)                     :: loglike
        double precision, intent(in)                     :: delta

        double precision, dimension(settings%nDims)             :: gradloglike


        double precision, dimension(settings%nDerived)               :: derived
        integer :: context

        double precision, dimension(settings%nDims) :: delta_vec

        integer i

        do i=1,settings%nDims
            delta_vec = 0
            delta_vec(i) = delta
            gradloglike(i) = ( loglikelihood(theta+delta_vec,derived,context) - loglike )/ delta
        end do

    end function gradloglike


end module calculate_module
