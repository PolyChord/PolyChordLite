module chordal_module
    implicit none


    contains

    function ChordalSampling(live_data, loglikelihood_bound, M,feedback)  result(new_point)
        use random_module, only: random_direction,random_hypercube_point,random_integer
        use model_module,  only: model, calculate_point, logzero

        implicit none

        ! ------- Inputs -------
        !> The current set of live points. 2D array:
        !!
        !! First index ranges over ( hypercube coords, physical coords, derived params, loglikelihood),
        !!
        !! Second index ranges over all the live points
        double precision, intent(inout), dimension(:,:) :: live_data

        !> The current loglikelihood bound
        double precision, intent(in) :: loglikelihood_bound

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M

        !> Optional argument to cause the sampler to print out relevent information
        integer, intent(in), optional :: feedback

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: new_point


        double precision,    dimension(M%nTotal)   :: random_point

        ! ------- Local Variables -------
        integer, dimension(1) :: point_number 
        integer :: nlive

        double precision,    dimension(M%nDims)   :: nhat

        double precision :: initial_step, acceleration

        initial_step = 0.0001
        acceleration = 2

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Sampler    : Chordal" )')
            end if
            return
        end if


        ! Get the number of live points
        nlive = size(live_data,2)


        random_point(M%d0) = loglikelihood_bound
        ! pick a random point
        !do while(random_point(M%d0) .ge. loglikelihood_bound)
            point_number = random_integer(1,nlive-1)        ! get a random number in [1,nlive-1]
            random_point = live_data(:,1+point_number(1))   ! get this point from live_data 
                                                            ! (excluding the possibility of drawing the late point)
        !end do

        ! get a random direction
        nhat = random_direction(M%nDims)

        ! generate a new random point seeded by this point
        new_point = random_chordal_point( nhat, random_point, initial_step, acceleration, loglikelihood_bound, M)

        ! record the loglikelihood of the created point
        live_data(M%d0,1+point_number(1)) = new_point(M%l0)

        ! mark the new point as ok
        new_point(M%d0) = logzero



    end function ChordalSampling



    function random_chordal_point(nhat,random_point,initial_step,acceleration,loglikelihood_bound,M) result(new_point)
        use model_module,  only: model, calculate_point, logzero
        use random_module, only: random_real
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: random_point
        !> The initial step
        double precision, intent(in) :: initial_step
        !> The acceleration parameter
        double precision, intent(in) :: acceleration
        !> The root value to find
        double precision, intent(in) :: loglikelihood_bound

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: new_point

        double precision,    dimension(M%nTotal)   :: edge_point_1
        double precision,    dimension(M%nTotal)   :: edge_point_2

        double precision random_temp


        ! find an edge in the nhat direction
        edge_point_1 = find_edge( nhat,random_point,initial_step,acceleration,loglikelihood_bound,M)
        ! find an edge in the opposite direction
        edge_point_2 = find_edge(-nhat,random_point,initial_step,acceleration,loglikelihood_bound,M)

        ! Set the loglikelihood of the new point to zero for now
        new_point(M%l0) = logzero

        ! Select a point randomly along the chord (edge_point_1, edge_point_2)
        do while(new_point(M%l0)< loglikelihood_bound)
            random_temp = random_real()
            new_point(M%h0:M%h1) =  random_temp * edge_point_1(M%h0:M%h1)  + (1d0-random_temp)*edge_point_2(M%h0:M%h1)
            call calculate_point(M,new_point)
        end do



    end function random_chordal_point



    function find_edge(nhat,start_point,initial_step,acceleration,loglikelihood_bound,M) result(finish_point)
        use model_module,  only: model, calculate_point
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: start_point
        !> The initial step
        double precision, intent(in) :: initial_step
        !> The acceleration parameter
        double precision, intent(in) :: acceleration
        !> The root value to find
        double precision, intent(in) :: loglikelihood_bound

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: finish_point

        ! The output finish point
        double precision,    dimension(M%nTotal)    :: new_point

        ! saving for root searching
        double precision,    dimension(M%nTotal)    :: old_point

        double precision :: step_length

        step_length=initial_step

        new_point = start_point

        ! Search using an expanding binary search
        do while( new_point(M%l0) > loglikelihood_bound )

            ! Save the old position for later
            old_point = new_point

            ! Take a new step
            new_point = new_point + nhat * step_length

            ! Compute physical coordinates, likelihoods and derived parameters
            call calculate_point( M, new_point )

            ! Increase the step length
            step_length = step_length * acceleration

        end do

        finish_point = new_point

    end function find_edge



end module chordal_module

