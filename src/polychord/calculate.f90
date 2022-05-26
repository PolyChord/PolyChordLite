module calculate_module
    use utils_module, only: dp
    implicit none
    contains

    subroutine calculate_point(loglikelihood,prior,point,settings,nlike)
        use settings_module, only: program_settings
        implicit none
        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in),   dimension(:) :: theta
                real(dp), intent(out),  dimension(:) :: phi
                real(dp) :: loglikelihood
            end function
        end interface
        interface
            function prior(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function
        end interface

        type(program_settings), intent(in) :: settings
        real(dp), intent(inout) , dimension(:) :: point
        integer, intent(inout) :: nlike

        real(dp),dimension(settings%nDims)    :: cube   ! Hypercube coordinate
        real(dp),dimension(settings%nDims)    :: theta  ! Physical parameters
        real(dp),dimension(settings%nDerived) :: phi    ! derived parameters
        real(dp)                              :: logL

        cube = point(settings%h0:settings%h1)

        if ( any(cube<0) .or. any(cube>1) )  then
            theta = 0
            logL  = settings%logzero
        else
            theta = prior(cube)
            logL  = loglikelihood(theta,phi)
        end if

        if(logL>settings%logzero) nlike = nlike+1

        point(settings%p0:settings%p1) = theta
        point(settings%d0:settings%d1) = phi
        point(settings%l0) = logL

    end subroutine calculate_point

    !> Calculate a posterior point from a live/phantom point
    function calculate_posterior_point(settings,point,logweight,evidence,volume) result(posterior_point)
        use settings_module,   only: program_settings
        use utils_module,      only: logincexp
        implicit none

        type(program_settings), intent(in) :: settings
        real(dp), dimension(settings%nTotal),intent(in) :: point
        real(dp),intent(in) :: logweight
        real(dp),intent(in) :: evidence
        real(dp),intent(in) :: volume
        real(dp), dimension(settings%nposterior) :: posterior_point


        ! Volume
        posterior_point(settings%pos_X)  = volume
        ! Likelihood
        posterior_point(settings%pos_l)  = point(settings%l0)
        ! Un-normalised weighting 
        posterior_point(settings%pos_w)  = logweight
        ! un-normalise cumulative weighting
        posterior_point(settings%pos_Z)  = evidence
        ! Physical parameters
        posterior_point(settings%pos_p0:settings%pos_p1) = point(settings%p0:settings%p1)
        ! Derived parameters
        posterior_point(settings%pos_d0:settings%pos_d1) = point(settings%d0:settings%d1)

    end function calculate_posterior_point


    !> This function computes the distance^2 matrix of an array of data.
    !!
    !! Assume that the points can be considered an indexed array of vectors
    !! V = ( v_i : i=1,n )
    !!
    !! The distance^2 matrix can be expressed very neatly as
    !! d_ij = (v_i-v_j) . (v_i-v_j)
    !!      = v_i.v_i + v_j.v_j - 2 v_i.v_j
    !!
    !! The final term can be written as a points^T points, and the first
    !! two are easy to write. We can therefore calculate this in two lines with
    !! instrisic functions
    function calculate_distance2_matrix(points) result(distance2_matrix)

        real(dp), intent(in), dimension(:,:) :: points

        real(dp), dimension(size(points,1),size(points,1)) ::distance2_matrix 

        integer :: i


        distance2_matrix = spread( &
            [ ( dot_product(points(i,:),points(i,:)), i=1,size(points,1) ) ], &
            dim=2,ncopies=size(points,1) )

        distance2_matrix = distance2_matrix + transpose(distance2_matrix) - 2d0 * matmul(points,transpose(points))

    end function calculate_distance2_matrix








end module calculate_module
