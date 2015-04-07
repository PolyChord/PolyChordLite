module calculate_module
    implicit none
    contains

    subroutine calculate_point(loglikelihood,priors,point,settings,nlike)
        use priors_module, only: prior, hypercube_to_physical
        use settings_module, only: program_settings
        use utils_module, only: logzero
        implicit none
        interface
            function loglikelihood(theta,phi)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                double precision :: loglikelihood
            end function
        end interface

        type(prior), dimension(:), intent(in) :: priors
        type(program_settings), intent(in) :: settings
        double precision, intent(inout) , dimension(:) :: point
        integer, intent(inout) :: nlike

        if ( any(point(settings%h0:settings%h1)<0d0) .or. any(point(settings%h0:settings%h1)>1d0) )  then
            point(settings%p0:settings%p1) = 0
            point(settings%l0) = logzero
        else
            ! Transform the the hypercube coordinates to the physical coordinates
            point(settings%p0:settings%p1) = hypercube_to_physical( point(settings%h0:settings%h1),priors )

            ! Calculate the likelihood and store it in the last index
            point(settings%l0) = loglikelihood( point(settings%p0:settings%p1), point(settings%d0:settings%d1))

            ! accumulate the number of likelihood calls that we've made
            if(point(settings%l0)>logzero) nlike = nlike+1
        end if

    end subroutine calculate_point

    !> Calculate a posterior point from a live/phantom point
    function calculate_posterior_point(settings,point,logweight,evidence,volume) result(posterior_point)
        use settings_module,   only: program_settings
        use utils_module,      only: logincexp
        implicit none

        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nTotal),intent(in) :: point
        double precision,intent(in) :: logweight
        double precision,intent(in) :: evidence
        double precision,intent(in) :: volume
        double precision, dimension(settings%nposterior) :: posterior_point


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


    !> This function computes the similarity matrix of an array of data.
    !!
    !! Assume that the data_array can be considered an indexed array of vectors
    !! V = ( v_i : i=1,n )
    !!
    !! The similarity matrix can be expressed very neatly as
    !! d_ij = (v_i-v_j) . (v_i-v_j)
    !!      = v_i.v_i + v_j.v_j - 2 v_i.v_j
    !!
    !! The final term can be written as a data_array^T data_array, and the first
    !! two are easy to write. We can therefore calculate this in two lines with
    !! instrisic functions
    function calculate_similarity_matrix(data_array) result(similarity_matrix)

        double precision, intent(in), dimension(:,:) :: data_array

        double precision, dimension(size(data_array,2),size(data_array,2)) :: similarity_matrix

        integer :: i


        similarity_matrix = spread( &
            [ ( dot_product(data_array(:,i),data_array(:,i)), i=1,size(data_array,2) ) ], &
            dim=2,ncopies=size(data_array,2) )

        similarity_matrix = similarity_matrix + transpose(similarity_matrix) - 2d0 * matmul( transpose(data_array),data_array )

    end function calculate_similarity_matrix








end module calculate_module
