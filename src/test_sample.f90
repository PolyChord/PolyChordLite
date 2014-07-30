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

        use random_module, only: random_direction,random_coordinate
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
        double precision, dimension(1) ::  rad




        ! Calculate the radial distance of the lowest loglikelihood point from the center
        new_point = live_data(:,1) -0.5  ! Displacement of lowest loglikelihood point from center

        ! create a random radial distance  between 0 and the
        ! displacement of the lowest loglikelihood point from center
        call random_coordinate(rad)      ! generate a random number between 0 and 1
        rad = rad * sqrt(dot_product(new_point(M%h0:M%h1),new_point(M%h0:M%h1))) 


        ! reset new_point to be the lowest loglikelihood point for now 
        ! (so the loop below is entered)
        new_point = live_data(:,1)

        do while(new_point(M%l0)<=loglikelihood_bound)

            ! Generate a random point
            call random_direction(new_point(:M%nDims))               ! Generate a random direction
            new_point(:M%nDims) = 0.5 + rad(1) * new_point(:M%nDims) ! Generate a point a distance rad away 0.5



            ! Transform the the hypercube coordinates to the physical coordinates
            call hypercube_to_physical( M, new_point )

            ! Calculate the loglikelihood and store it in the last index
            new_point(M%l0) = M%loglikelihood( new_point(M%p0:M%p1) )

            ! Calculate the derived parameters
            call calculate_derived_parameters( M, new_point )
        end do


    end subroutine SphericalCenterSampling

end module test_sampler_module

