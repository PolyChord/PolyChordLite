module calculate_module
    implicit none
    contains

    subroutine calculate_point(loglikelihood,priors,point,settings)
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

        if ( any(point(settings%h0:settings%h1)<0d0) .or. any(point(settings%h0:settings%h1)>1d0) )  then
            point(settings%p0:settings%p1) = 0
            point(settings%l0) = logzero
        else
            ! Transform the the hypercube coordinates to the physical coordinates
            point(settings%p0:settings%p1) = hypercube_to_physical( point(settings%h0:settings%h1),priors )

            ! Calculate the likelihood and store it in the last index
            point(settings%l0) = loglikelihood( point(settings%p0:settings%p1), point(settings%d0:settings%d1),settings%context)

            ! accumulate the number of likelihood calls that we've made
            point(settings%nlike) = point(settings%nlike)+1
        end if

    end subroutine calculate_point

    function calculate_gradloglike(loglikelihood,priors,point,settings,delta) result(gradloglike)
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
        double precision, intent(in) , dimension(:) :: point
        double precision, intent(in)                     :: delta
        double precision, dimension(settings%nDims)             :: gradloglike

        double precision , dimension(settings%nTotal) :: center_point
        double precision , dimension(settings%nTotal) :: outer_point

        integer :: i

        ! Calculate the base point
        center_point=point
        call calculate_point(loglikelihood,priors,center_point,settings)
        if (center_point(settings%l0) <= logzero ) then
            gradloglike=0
            return
        end if

        do i=1,settings%nDims
            ! intialise the outer point at center_point
            outer_point=point
            ! shift it in the ith coordinate by delta
            outer_point(settings%h0+i) = outer_point(settings%h0+i) + delta
            ! calculate the loglikelihood
            call calculate_point(loglikelihood,priors,outer_point,settings)
            if(outer_point(settings%l0)<=logzero) then
                gradloglike=0
                gradloglike(i)=-1d0
                return
            else
                ! grad loglike is the difference between outer and center divided by delta
                gradloglike(i) = (outer_point(settings%l0)-center_point(settings%l0))/delta
            end if
        end do

    end function calculate_gradloglike



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
