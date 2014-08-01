module test_sampler_module
    implicit none

    contains

    subroutine BruteForceSampling(new_point, live_data, loglikelihood_bound, M) 

        use random_module, only: random_coordinate
        use model_module,  only: model, hypercube_to_physical, calculate_derived_parameters

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

        ! ------- Outputs -------
        !> The newly generated point
        double precision, intent(out),    dimension(:)   :: new_point



        ! set new_point to be the lowest loglikelihood point for now 
        ! (so the loop below is entered)
        new_point = live_data(:,1)

        do while(new_point(M%l0)<=loglikelihood_bound)

            ! Generate a random point
            call random_coordinate(new_point(:M%nDims))              ! Generate a random coordinate

            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, new_point )

            ! Calculate the loglikelihood and store it in the last index
            new_point(M%l0) = M%loglikelihood( new_point(M%p0:M%p1) )

            ! Calculate the derived parameters
            call calculate_derived_parameters( M, new_point )
        end do


    end subroutine BruteForceSampling

    subroutine SphericalCenterSampling(new_point, live_data, loglikelihood_bound, M) 

        use random_module, only: random_point_in_sphere
        use model_module,  only: model, hypercube_to_physical, calculate_derived_parameters

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

        ! ------- Outputs -------
        !> The newly generated point
        double precision, intent(out),    dimension(:)   :: new_point

        ! ------- Local Variables -------
        double precision, dimension(1) ::  rand_rad

        double precision, parameter :: Vfrac = 3

        double precision :: max_radius
        double precision, dimension(M%nDims) :: center


        ! Set the center of the sampler
        center = 0.5

        ! Calculate the radial distance of the lowest loglikelihood point from the center
        
        max_radius = Vfrac**(1/dble(M%nDims) )*  sqrt( dot_product(  &
            live_data(M%h0:M%h1,1)-center,                           &
            live_data(M%h0:M%h1,1)-center   )) 

        ! reset new_point to be the lowest loglikelihood point for now 
        ! (so the loop below is entered)
        new_point = live_data(:,1)

        do while(new_point(M%l0)<=loglikelihood_bound)

            ! Generate a random point within the unit sphere
            call random_point_in_sphere(new_point(M%h0:M%h1))

            ! scale it to be within max_radius and centered on 0.5
            new_point(M%h0:M%h1) = 0.5 +  max_radius * new_point(M%h0:M%h1)

            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, new_point )

            ! Calculate the loglikelihood and store it in the last index
            new_point(M%l0) = M%loglikelihood( new_point(M%p0:M%p1) )

            ! Calculate the derived parameters
            call calculate_derived_parameters( M, new_point )

        end do


    end subroutine SphericalCenterSampling

end module test_sampler_module

