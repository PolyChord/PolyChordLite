module test_sampler_module
    implicit none

    contains

    subroutine SphericalCenterSampling(new_point, live_data, likelihood_bound, M) 
        use random_module, only: random_direction,random_coordinate
        use model_module,  only: model, hypercube_to_physical, calculate_derived_parameters
        implicit none
        double precision, intent(out),    dimension(:)   :: new_point
        double precision, intent(in), dimension(:,:) :: live_data
        double precision, intent(in) :: likelihood_bound
        type(model),            intent(in) :: M

        double precision, dimension(1) ::  rad


        new_point = live_data(:,1) -0.5

        call random_coordinate(rad)

        rad = rad * sqrt(dot_product(new_point(M%h0:M%h1),new_point(M%h0:M%h1)))

        new_point = live_data(:,1)

        ! Calculate the likelihood and store it in the last index
        do while(new_point(M%l0)<=likelihood_bound)


            ! Generate a random coordinate
            call random_direction(new_point(:M%nDims))

            new_point(:M%nDims) = 0.5 + rad(1) * new_point(:M%nDims)

            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, new_point )

            new_point(M%l0) = M%loglikelihood( new_point(M%p0:M%p1) )

            ! Calculate the derived parameters
            call calculate_derived_parameters( M, new_point )
        end do


    end subroutine SphericalCenterSampling

end module test_sampler_module

