module chaotic_chordal_module
    implicit none


    contains

    function ChaoticChordalSampling(live_data, loglikelihood_bound, M,feedback)  result(new_point)
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

        double precision :: max_distance
        integer :: num_chords

        integer :: i

        max_distance = sqrt(M%nDims+0d0)
        num_chords = 4

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
        point_number = random_integer(1,nlive-1)        ! get a random number in [1,nlive-1]
        new_point = live_data(:,1+point_number(1))      ! get this point from live_data 
                                                        ! (excluding the possibility of drawing the late point)

        do i=1,num_chords

            ! get a random direction nhat
            nhat = random_direction(M%nDims)

            ! generate a new random point along the chord defined by new_point and nhat
            new_point = random_chordal_point( nhat, new_point, loglikelihood_bound, max_distance, M)

        end do


    end function ChaoticChordalSampling



    function random_chordal_point(nhat,random_point,loglikelihood_bound,max_distance,M) result(new_point)
        use model_module,  only: model, calculate_point, logzero
        use random_module, only: random_real
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: random_point
        !> The root value to find
        double precision, intent(in) :: loglikelihood_bound
        !> The maximum distance to explore
        double precision, intent(in) :: max_distance

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: new_point

        ! The upper bound
        double precision,    dimension(M%nTotal)   :: u_bound
        ! The lower bound
        double precision,    dimension(M%nTotal)   :: l_bound


        ! Get the start bounds
        u_bound(M%h0:M%h1) = restrict_to_hypercube( random_point(M%h0:M%h1) + nhat * max_distance, nhat, M%nDims)
        l_bound(M%h0:M%h1) = restrict_to_hypercube( random_point(M%h0:M%h1) - nhat * max_distance, nhat, M%nDims)
        call calculate_point(M,u_bound)
        call calculate_point(M,l_bound)

        new_point = find_positive_within(l_bound,u_bound)

        contains

        recursive function find_positive_within(l_bound,u_bound) result(new_point)
            implicit none
            !> The upper bound
            double precision, intent(inout), dimension(M%nTotal)   :: u_bound
            !> The lower bound
            double precision, intent(inout), dimension(M%nTotal)   :: l_bound

            ! The output finish point
            double precision,    dimension(M%nTotal)   :: new_point

            ! Draw a random point within l_bound and u_bound
            new_point(M%h0:M%h1) = l_bound(M%h0:M%h1) + random_real() * (u_bound(M%h0:M%h1) - l_bound(M%h0:M%h1))

            ! calculate the likelihood 
            call calculate_point(M,new_point)

            ! If we're not within the likelihood bound then we need to sample further
            if( new_point(M%l0) < loglikelihood_bound ) then

                if ( dot_product(new_point(M%h0:M%h1)-random_point(M%h0:M%h1),nhat) > 0d0 ) then
                    ! If new_point is on the u_bound side of random_point, then
                    ! contract u_bound
                    u_bound = new_point
                else
                    ! If new_point is on the l_bound side of random_point, then
                    ! contract l_bound
                    l_bound = new_point
                end if

                ! Call the function again
                new_point = find_positive_within(l_bound,u_bound)

            end if
            ! otherwise new_point is returned

        end function find_positive_within


    end function random_chordal_point


    !> Takes a point a, and if it is outside the hyper cube moves it back along
    !! nhat until it is at the edge of the cube
    function restrict_to_hypercube(a,nhat,nDims) result(abar)
        implicit none
        !> Dimensionality of the cube
        integer, intent(in) :: nDims
        !> The direction of the line
        double precision, intent(in),    dimension(nDims)   :: nhat
        !> A point on the line
        double precision, intent(in),    dimension(nDims)   :: a

        !> The restricted point
        double precision,    dimension(nDims)   :: abar

        if( any(a>1 .or. a<0) )  then
            abar = intersect_hypercube_with_line(nhat,a,nDims)
        else
            abar = a
        end if

    end function

    !> Calculates the 'nearest' intersection of a line with the unit hypercube
    !! (where nearest is understood in the sense of the uniform norm - see below)
    !!
    !! The line is parameterised by \f$\lambda\f$ as:
    !! \f[ \mathbf{r}(\lambda) = \mathbf{a} + \lambda \mathbf{\hat{n}}\f]
    !!
    !! We use a clever trick by exploiting the fact that in the [Uniform Norm](http://en.wikipedia.org/wiki/Uniform_norm)
    !! we may define a metric:
    !! \f[ d_\infty(x,y) = \|x-y\|_\infty = \max(\{|x_i-y_i| : i=1\ldots d\}) \f]
    !! such that "spheres" \f$d_\infty(x,0) = const \f$ are in fact hypercubes.
    !!
    !! The unit hypercube is thus: 
    !! \f[ d_\infty(x,\mathbf{0.5}) = 0.5 \f]
    !! where \f$\mathbf{0.5} = (0.5,\ldots,0.5)\f$ is the center of the unit hypercube
    !!
    !! So to find the intersection points of the line with the unit hypercube,
    !! we just normalise any point on the line subject to the measure \f$d_\infty(\mathbf{r},\mathbf{0.5})/0.5\f$
    !!
    !! Note that this will only return the nearest intersection point in terms
    !! of this norm
    !! 
    function intersect_hypercube_with_line(a,nhat,nDims) result(edge_point)
        implicit none
        !> Dimensionality of the cube
        integer, intent(in) :: nDims
        !> A point on the line
        double precision, intent(in),    dimension(nDims)   :: a
        !> The direction of the line
        double precision, intent(in),    dimension(nDims)   :: nhat

        ! The output hypercube points
        integer,    dimension(nDims)   :: edge_point

        ! here we do a lot of things at once:
        ! 1) rescale the origin a is measured from the center of the cube (a->a-0.5)
        ! 2) calculate the d_infty length of this vector
        ! 3) use this length to normalise a-0.5 (note the factor of 2 since ours
        !    is a 'sphere of radius 0.5 in this metric)
        ! 4) move the origin back to 0 by adding 0.5
        ! 5) round to the nearest integer to avoid any rounding errors that may
        !    accumulate and return a vector of integers
        edge_point = nint( 2d0 * (a-0.5d0) / maxval(abs(a-0.5d0))  + 0.5d0 ) 

    end function intersect_hypercube_with_line

end module chaotic_chordal_module

