module galileo_module
    implicit none

    contains

    subroutine GalileanSample(new_point, live_data, likelihood_bound, M) 
        use random_module, only: random_direction,random_hypercube_point
        use model_module,  only: model, calculate_point
        implicit none
        double precision, intent(out),    dimension(:)   :: new_point
        double precision, intent(in), dimension(:,:) :: live_data
        double precision, intent(in) :: likelihood_bound
        type(model),            intent(in) :: M

        new_point = live_data(:,1)

        ! Calculate the likelihood and store it in the last index
        do while(new_point(M%l0)<=likelihood_bound)

            ! Generate a random coordinate
            new_point = random_hypercube_point(M%nDims)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, new_point )

        end do


    end subroutine GalileanSample

end module galileo_module

