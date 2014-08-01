module test_sampler_module
    implicit none

    contains

    function BruteForceSampling(live_data, loglikelihood_bound, M) result(new_point)

        use random_module, only: random_hypercube_point
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

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: new_point



        ! set new_point to be the lowest loglikelihood point for now 
        ! (so the loop below is entered)
        new_point = live_data(:,1)

        do while(new_point(M%l0)<=loglikelihood_bound)

            ! Generate a random point
            new_point = random_hypercube_point(M%nDims)              ! Generate a random coordinate

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, new_point )

        end do


    end function BruteForceSampling

    function SphericalCenterSampling(live_data, loglikelihood_bound, M) result(new_point)

        use random_module, only: random_point_in_sphere
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

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: new_point

        ! ------- Local Variables -------
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
            new_point = random_point_in_sphere(M%nDims)

            ! scale it to be within max_radius and centered on 0.5
            new_point(M%h0:M%h1) = center +  max_radius * new_point(M%h0:M%h1)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, new_point )

        end do


    end function SphericalCenterSampling



    function CubicCenterSampling(live_data, loglikelihood_bound, M) result(new_point)

        use random_module, only: random_hypercube_point
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

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: new_point

        ! ------- Local Variables -------
        double precision, parameter :: Vfrac = 3

        double precision :: max_radius
        double precision, dimension(M%nDims) :: center

        integer num_attempts
        integer,save :: max_num_attempts


        num_attempts=0

        ! Set the center of the sampler
        center = 0.5

        ! Calculate the radial distance of the lowest loglikelihood point from the center
        
        max_radius = sqrt( dot_product(  &
            live_data(M%h0:M%h1,1)-center,                           &
            live_data(M%h0:M%h1,1)-center   )) 

        ! set new_point to be the lowest loglikelihood point for now 
        ! (so the loop below is entered)
        new_point = live_data(:,1)

        do while(new_point(M%l0)<=loglikelihood_bound)

            ! Generate a random point within the unit sphere
            new_point = random_hypercube_point(M%nDims)

            ! Re scale it about center with a side length max_radius*2
            new_point(M%h0:M%h1) = center + 2*min(max_radius,1d0) * (new_point(M%h0:M%h1) -center)

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, new_point )

            num_attempts=num_attempts+1
        end do

        if (num_attempts> max_num_attempts) max_num_attempts=num_attempts
        write(*,'("max_num_attempts: ", I9)') max_num_attempts


    end function CubicCenterSampling

end module test_sampler_module

