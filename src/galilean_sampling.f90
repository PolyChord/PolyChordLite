module galileo_module
    implicit none

    contains

    function GalileanSampling(live_data, loglikelihood_bound, M,feedback)  result(new_point)
        use random_module, only: random_direction,random_hypercube_point,random_integer
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

        ! ------- Local Variables -------
        integer, dimension(1) :: point_number 
        integer :: nlive

        double precision,    dimension(M%nDims)   :: nhat

        double precision :: initial_step, acceleration


        initial_step = 0.01
        acceleration = 2

        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Sampler    : Galilean" )')
            end if
            return
        end if


        ! Get the number of live points
        nlive = size(live_data,2)

        ! pick a random point
        point_number = random_integer(1,nlive)        ! get a random number in [1,nlive]
        new_point = live_data(:,point_number(1))      ! get this point from live_data

        ! get a random direction
        nhat = random_direction(M%nDims)

        write(*,*) new_point(M%h0:M%h1)
        ! find a root in that direction
        new_point = find_edge(nhat,new_point,initial_step,acceleration,loglikelihood_bound,M)

        do while(.true.)
            write(*,*) new_point(M%h0:M%h1)

            ! get a random direction pointing inwards
            nhat = random_inwards_direction(new_point,initial_step,loglikelihood_bound,M)

            new_point = find_edge(nhat,new_point,initial_step,acceleration,loglikelihood_bound,M) 

        end do




    end function GalileanSampling


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

        ! We've found a root, so refine it down to the precision of initial_step
        finish_point =  refine_edge(new_point,old_point,initial_step*1d-2,M,loglikelihood_bound)


    end function find_edge


    function refine_edge(new_point,old_point,prec,M,loglikelihood_bound) result(s)
        use model_module,  only: model, calculate_point
        implicit none
        !> details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> lower bracket
        double precision, intent(in), dimension(M%nTotal)    :: new_point
        !> upper bracket
        double precision, intent(in), dimension(M%nTotal)    :: old_point
        !> initial step
        double precision, intent(in) :: prec
        !> root value to find
        double precision, intent(in) :: loglikelihood_bound

        ! final root
        double precision, dimension(M%nTotal) :: s




        double precision, dimension(M%nTotal)    :: a,b

        a = old_point
        b = new_point

        do while ( dot_product(a(M%h0:M%h1)-b(M%h0:M%h1),a(M%h0:M%h1)-b(M%h0:M%h1)) > prec**2 )

            ! find the midpoint
            s(M%h0:M%h1) =  ( a(M%h0:M%h1) + b(M%h0:M%h1) ) /2d0

            ! calculate the midpoint
            call calculate_point(M,s)

            ! adjust the bounds
            if (s(M%l0) > loglikelihood_bound ) then
                a = s
            else
                b = s
            end if
        end do

        if(a(M%l0) > b(M%l0) ) then
            s=a
        else
            s=b
        end if




    end function refine_edge

    function random_inwards_direction(new_point,delta,loglikelihood_bound,M) result(nhat)
        use random_module, only: random_direction
        use model_module,  only: model, calculate_point,logzero
        implicit none
        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: new_point
        !> The sphere size to search in
        double precision, intent(in) :: delta
        !> The log likelihood bound
        double precision, intent(in) :: loglikelihood_bound

        ! The output random direction pointing inwards
        double precision,    dimension(M%nDims)   :: nhat

        ! temp position
        double precision,    dimension(M%nTotal)   :: temp_point

        ! Initialise at logzero
        temp_point(M%l0) = logzero

        do while( temp_point(M%l0) <loglikelihood_bound)
            nhat=random_direction(M%nDims)
            temp_point = new_point + nhat * delta
            call calculate_point(M, temp_point)
        end do



    end function random_inwards_direction 



end module galileo_module

