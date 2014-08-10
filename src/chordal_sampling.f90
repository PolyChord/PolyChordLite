module chordal_module
    implicit none


    contains

    function ChordalSampling(settings,seed_point, loglikelihood_bound, cluster, M,feedback)  result(baby_point)
        use settings_module, only: program_settings
        use random_module, only: random_skewed_direction,random_direction,random_reals,random_integers
        use model_module,  only: model, calculate_point
        use utils_module, only: logzero
        use cluster_module, only: cluster_info

        implicit none

        ! ------- Inputs -------
        !> program settings (mostly useful to pass on the number of live points)
        class(program_settings), intent(in) :: settings

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M

        !> The seed point
        double precision, intent(in), dimension(M%nTotal)   :: seed_point

        !> The current loglikelihood bound
        double precision, intent(in) :: loglikelihood_bound

        !> The current clustering information
        type(cluster_info), intent(in) :: cluster

        !> Optional argument to cause the sampler to print out relevent information
        integer, intent(in), optional :: feedback

        ! ------- Outputs -------
        !> The newly generated point
        double precision,    dimension(M%nTotal)   :: baby_point


        ! ------- Local Variables -------
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


        ! Start the baby point at the seed point
        baby_point = seed_point

        ! Set the number of likelihood evaluations to zero
        baby_point(M%d0) = 0

        do i=1,settings%num_chords

            ! get a random direction nhat
            if( settings%do_clustering) then
                nhat = random_skewed_direction(M%nDims,cluster%cholesky(:,:,1),cluster%invcovmat(:,:,1))
            else 
                nhat = random_direction(M%nDims) 
            end if

            ! generate a new random point along the chord defined by baby_point and nhat
            baby_point = random_chordal_point( nhat, baby_point, loglikelihood_bound, M)
        end do


    end function ChordalSampling



    function random_chordal_point(nhat,seed_point,loglikelihood_bound,M) result(baby_point)
        use model_module,  only: model, calculate_point
        use utils_module,  only: logzero
        use random_module, only: random_reals
        implicit none

        !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
        type(model),            intent(in) :: M
        !> The direction to search for the root in
        double precision, intent(in),    dimension(M%nDims)   :: nhat
        !> The start point
        double precision, intent(in),    dimension(M%nTotal)   :: seed_point
        !> The root value to find
        double precision, intent(in) :: loglikelihood_bound

        ! The output finish point
        double precision,    dimension(M%nTotal)   :: baby_point

        ! The upper bound
        double precision,    dimension(M%nTotal)   :: u_bound
        ! The lower bound
        double precision,    dimension(M%nTotal)   :: l_bound

        double precision :: trial_chord_length


        ! record the number of likelihood calls
        u_bound(M%d0) = seed_point(M%d0)
        l_bound(M%d0) = 0

        ! set the likelihoods of start bounds so that the loop below is entered
        u_bound(M%l0) = loglikelihood_bound
        l_bound(M%l0) = loglikelihood_bound

        ! set the trial chord length at half the bound of the old length
        trial_chord_length = seed_point(M%d0+1)/2

        do while(u_bound(M%l0) >= loglikelihood_bound )
            trial_chord_length = 2*trial_chord_length
            u_bound(M%h0:M%h1) = seed_point(M%h0:M%h1) + nhat * trial_chord_length
            call calculate_point(M,u_bound)
        end do

        trial_chord_length = trial_chord_length/2
        do while(l_bound(M%l0) >= loglikelihood_bound )
            trial_chord_length = 2*trial_chord_length
            l_bound(M%h0:M%h1) = seed_point(M%h0:M%h1) - nhat * trial_chord_length
            call calculate_point(M,l_bound)
        end do

        baby_point = find_positive_within(l_bound,u_bound)

        ! Store the new length
        baby_point(M%d0+1) = max(&
            dot_product(u_bound(M%h0:M%h1)-seed_point(M%h0:M%h1),u_bound(M%h0:M%h1)-seed_point(M%h0:M%h1)),&
            dot_product(l_bound(M%h0:M%h1)-seed_point(M%h0:M%h1),l_bound(M%h0:M%h1)-seed_point(M%h0:M%h1)))
        baby_point(M%d0+1) = sqrt(baby_point(M%d0+1))
        

        contains

        recursive function find_positive_within(l_bound,u_bound) result(finish_point)
            implicit none
            !> The upper bound
            double precision, intent(inout), dimension(M%nTotal)   :: u_bound
            !> The lower bound
            double precision, intent(inout), dimension(M%nTotal)   :: l_bound

            ! The output finish point
            double precision,    dimension(M%nTotal)   :: finish_point

            double precision,dimension(1) :: random_temp

            ! Draw a random point within l_bound and u_bound
            random_temp =random_reals(1)
            finish_point(M%h0:M%h1) = l_bound(M%h0:M%h1)*(1d0-random_temp(1)) + random_temp(1) * u_bound(M%h0:M%h1)

            ! Pass on the number of likelihood calls that have been made
            finish_point(M%d0) = l_bound(M%d0) + u_bound(M%d0)
            ! zero the likelihood calls for l_bound and u_bound, as these are
            ! now stored in point
            l_bound(M%d0) = 0
            u_bound(M%d0) = 0

            ! calculate the likelihood 
            call calculate_point(M,finish_point)

            ! If we're not within the likelihood bound then we need to sample further
            if( finish_point(M%l0) < loglikelihood_bound ) then

                if ( dot_product(finish_point(M%h0:M%h1)-seed_point(M%h0:M%h1),nhat) > 0d0 ) then
                    ! If finish_point is on the u_bound side of seed_point, then
                    ! contract u_bound
                    u_bound = finish_point
                else
                    ! If finish_point is on the l_bound side of seed_point, then
                    ! contract l_bound
                    l_bound = finish_point
                end if

                ! Call the function again
                finish_point = find_positive_within(l_bound,u_bound)

            end if
            ! otherwise finish_point is returned

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

