module nested_sampling_module
    implicit none

    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    subroutine NestedSampling(loglikelihood,model,settings)
        use model_module, only: model_details
        use settings_module,  only: program_settings

        implicit none
        type(model_details),    intent(in) :: model
        type(program_settings), intent(in) :: settings


        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, likelihood)
        double precision, dimension(model%nDims+1,settings%nlive) :: live_data

        interface
            !> Log(Likelihood) 
            function loglikelihood(theta)
                double precision, intent(in), dimension(:) :: theta
            end function loglikelihood
        end interface


        
        ! Generate initial live points
        !call GenerateLivePoints(live_data)

        ! Calculate likelihoods



        ! Generate a new point within the likelihood bounds
        !new_point = GallileanSample(live_points,likelihood_bound)

    end subroutine NestedSampling





    !> Generate an initial set of live points distributed uniformly in the unit
    !! hypercube
    subroutine GenerateLivePoints(live_data)
        implicit none

        !> This is a very important array. live_data(:,i) constitutes the
        !! information in the ith live point in the unit hypercube:
        !! ( <-hypercube coordinates->, likelihood)
        double precision, intent(out), dimension(:) :: live_data

        ! Generate nlive random numbers
        live_data=1

    end subroutine 


end module nested_sampling_module
