module chordal_module
    implicit none


    contains

    function ChordalSampling(settings,live_data, loglikelihood_bound, M,feedback)  result(new_point)
        use settings_module, only: program_settings
        use random_module, only: random_direction,random_reals,random_integers
        use model_module,  only: model, calculate_point, logzero

        implicit none

        ! ------- Inputs -------
        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

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

        double precision,    dimension(M%nDims)   :: nhat


        integer :: i


        ! Feedback if requested
        if(present(feedback)) then
            if(feedback>=0) then
                write(*,'( "Sampler    : Chordal" )')
                write(*,'( "  num chords = ",I4 )') settings%num_chords
            end if
            return
        end if



        random_point(M%d0) = loglikelihood_bound
        ! pick a random point
        point_number = random_integers(1,settings%nlive-1) ! get a random number in [1,nlive-1]
        new_point = live_data(:,1+point_number(1))        ! get this point from live_data 
                                                          ! (excluding the possibility of drawing the late point)

        ! Set the number of likelihood evaluations to zero
        new_point(M%d0) = 0

        do i=1,settings%num_chords

            ! get a random direction nhat
            nhat = random_direction(M%nDims)

            ! generate a new random point along the chord defined by new_point and nhat
            new_point = random_chordal_point( nhat, new_point, loglikelihood_bound, M)
        end do


    end function ChordalSampling



    function random_chordal_point(nhat,random_point,loglikelihood_bound,M) result(new_point)
        use model_module,  only: model, calculate_point, logzero
        use random_module, only: random_reals
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: random_point
        !> The root value to find
        double precision, intent(in) :: loglikelihood_bound

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: new_point

        ! The upper bound
        double precision,    dimension(M%nTotal)   :: u_bound
        ! The lower bound
        double precision,    dimension(M%nTotal)   :: l_bound

        double precision :: trial_chord_length


        ! record the number of likelihood calls
        u_bound(M%d0) = random_point(M%d0)
        l_bound(M%d0) = 0

        ! set the likelihoods of start bounds so that the loop below is entered
        u_bound(M%l0) = loglikelihood_bound
        l_bound(M%l0) = loglikelihood_bound

        ! set the trial chord length at half the bound of the old length
        trial_chord_length = random_point(M%d0+1)/2

        do while(u_bound(M%l0) >= loglikelihood_bound )
            trial_chord_length = 2*trial_chord_length
            u_bound(M%h0:M%h1) = random_point(M%h0:M%h1) + nhat * trial_chord_length
            call calculate_point(M,u_bound)
        end do

        trial_chord_length = trial_chord_length/2
        do while(l_bound(M%l0) >= loglikelihood_bound )
            trial_chord_length = 2*trial_chord_length
            l_bound(M%h0:M%h1) = random_point(M%h0:M%h1) - nhat * trial_chord_length
            call calculate_point(M,l_bound)
        end do

        new_point = find_positive_within(l_bound,u_bound)

        ! Store the new length
        new_point(M%d0+1) = max(&
            dot_product(u_bound(M%h0:M%h1)-random_point(M%h0:M%h1),u_bound(M%h0:M%h1)-random_point(M%h0:M%h1)),&
            dot_product(l_bound(M%h0:M%h1)-random_point(M%h0:M%h1),l_bound(M%h0:M%h1)-random_point(M%h0:M%h1)))
        new_point(M%d0+1) = sqrt(new_point(M%d0+1))
        

        contains

        recursive function find_positive_within(l_bound,u_bound) result(point)
            implicit none
            !> The upper bound
            double precision, intent(inout), dimension(M%nTotal)   :: u_bound
            !> The lower bound
            double precision, intent(inout), dimension(M%nTotal)   :: l_bound

            ! The output finish point
            double precision,    dimension(M%nTotal)   :: point

            double precision,dimension(1) :: random_temp

            ! Draw a random point within l_bound and u_bound
            random_temp =random_reals(1)
            point(M%h0:M%h1) = l_bound(M%h0:M%h1)*(1d0-random_temp(1)) + random_temp(1) * u_bound(M%h0:M%h1)

            ! Pass on the number of likelihood calls that have been made
            point(M%d0) = l_bound(M%d0) + u_bound(M%d0)
            ! zero the likelihood calls for l_bound and u_bound, as these are
            ! now stored in point
            l_bound(M%d0) = 0
            u_bound(M%d0) = 0

            ! calculate the likelihood 
            call calculate_point(M,point)

            ! If we're not within the likelihood bound then we need to sample further
            if( point(M%l0) < loglikelihood_bound ) then

                if ( dot_product(point(M%h0:M%h1)-random_point(M%h0:M%h1),nhat) > 0d0 ) then
                    ! If new_point is on the u_bound side of random_point, then
                    ! contract u_bound
                    u_bound = point
                else
                    ! If new_point is on the l_bound side of random_point, then
                    ! contract l_bound
                    l_bound = point
                end if

                ! Call the function again
                point = find_positive_within(l_bound,u_bound)

            end if
            ! otherwise new_point is returned

        end function find_positive_within


    end function random_chordal_point


    !> Takes a point a, and if it is outside the hypercube moves it back along
    !! nhat until it is at the edge of the cube
    !! This function assumes that nhat is pointing away out of the cube
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

        double precision, dimension(nDims) :: lambda

        lambda = nhat

        if( any(a>1 .or. a<0) )  then
            !abar = intersect_hypercube_with_line(nhat,a,nDims)
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

end module chordal_module

