module galileo_module
    implicit none

    contains

    subroutine GalileanSample(new_point, live_data, likelihood_bound, M) 
        use random_module, only: random_coordinate
        use model_module,  only: model, hypercube_to_physical
        implicit none
        double precision, intent(out),    dimension(:)   :: new_point
        double precision, intent(in), dimension(:,:) :: live_data
        double precision, intent(in) :: likelihood_bound
        type(model),            intent(in) :: M


        new_point = live_data(:,1)

        ! Calculate the likelihood and store it in the last index
        do while(new_point(M%l0)<=likelihood_bound)

            ! Generate a random coordinate
            call random_coordinate(new_point(:M%nDims))

            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, new_point )


            new_point(M%l0) = M%loglikelihood( new_point(M%p0:M%p1) )
        end do


    end subroutine GalileanSample

end module galileo_module

