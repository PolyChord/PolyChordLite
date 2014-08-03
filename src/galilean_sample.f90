module galileo_module
    implicit none

    contains

    function GalileanSampling(live_data, loglikelihood_bound, M,feedback)  result(new_point)
        use random_module, only: random_direction,random_hypercube_point
        use model_module,  only: model, calculate_point

        implicit none

        ! ------- Inputs -------
        !> The current set of live points. 2D array:
        !!
        !! First index ranges over ( hypercube coords, physical coords, derived params, loglikelihood),
        !!
        !! Second index ranges over all the live points
        double precision, intent(in), dimension(:,:) :: live_data

        !> The current loglikelihood bound
        double precision, intent(in) :: loglikelihood_bound

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M

        !> Optional argument to cause the sampler to print out relevent information
        integer, intent(in), optional :: feedback

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: new_point

        new_point = live_data(:,1)

        do while(new_point(M%l0)<=loglikelihood_bound)

            ! Generate a random coordinate
            new_point = random_hypercube_point(M%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, new_point )

        end do


    end function GalileanSampling

end module galileo_module

