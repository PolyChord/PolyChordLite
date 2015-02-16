module array_module
    contains



    !> Reallocate a 3D array of doubles
    subroutine reallocate_3(array,u1,u2,u3)
        implicit none
        double precision, dimension(:,:,:),allocatable, intent(inout) :: array
        integer, intent(in),optional :: u1,u2,u3

        integer :: size1,size2,size3

        double precision, dimension(size(array,1),size(array,2),size(array,3)) :: temp_array

        ! Default reallocation size is the size of the original array
        size1 = size(array,1)
        size2 = size(array,2)
        size3 = size(array,3)

        ! If the argument is present, then reallocate to the new size
        if( present(u1) ) size1 = u1
        if( present(u2) ) size2 = u2
        if( present(u3) ) size3 = u3

        temp_array = array                  ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(size1,size2,size3))  ! Re-allocate with new size
        array(1:size(temp_array,1),1:size(temp_array,2),1:size(temp_array,3)) = temp_array

    end subroutine




    subroutine add_point(point,array,narray,cluster_id)
        implicit none
        double precision, dimension(:), intent(in) :: point                     !> Point to be added to end of array
        double precision, dimension(:,:,:), allocatable, intent(inout) :: array !> Array to be added to
        integer,dimension(:), allocatable, intent(inout) :: narray              !> number of points in array (second index)
        integer, intent(in) :: cluster_id                                       !> cluster identity (third index)

        narray(cluster_id) = narray(cluster_id) + 1         ! Increase the number of points in the cluster

        ! If this takes us over the size of the array, then double it
        if(narray(cluster_id) > size(array,2) ) call reallocate_3(array, u2=size(array,2)*2 )

        array(:,narray(cluster_id),cluster_id) = point       ! Add the point to the end position

    end subroutine add_point



    function delete_point(i_point,array,narray,cluster_id) result(point)
        implicit none
        integer, intent(in) :: i_point                                          !> Position of point to be deleted from array
        double precision, dimension(:,:,:), allocatable, intent(inout) :: array !> Array to be delete from
        integer,dimension(:), allocatable, intent(inout) :: narray              !> number of points in array (second index)
        integer, intent(in) :: cluster_id                                       !> cluster identity (third index)
        double precision, dimension(size(array,1)) :: point                     ! The point we have just deleted

        point = array(:,i_point,cluster_id)                                  ! Output the point to be deleted
        array(:,i_point,cluster_id) = array(:,narray(cluster_id),cluster_id) ! delete the point by overwriting it with the point at the end 
        narray(cluster_id) = narray(cluster_id) - 1                          ! reduce the number of points in the cluster

    end function delete_point





end module array_module
