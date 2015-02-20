module array_module
    contains


    !> Reallocate a 1D array of doubles
    subroutine reallocate_1_d(array,u1,p1)
        implicit none
        !> Array to be reallocated
        double precision, dimension(:),allocatable, intent(inout) :: array
        !> New upper bound
        integer, intent(in),optional :: u1
        !> new array configuration after reallocation
        integer, dimension(size(array,1)),intent(in),optional :: p1

        integer :: size1

        integer, dimension(size(array,1)) :: indices1

        integer :: i

        double precision, dimension(size(array,1)) :: temp_array

        ! Default reallocation size is the size of the original array
        size1 = size(array,1)

        ! If the argument is present, then reallocate to the new size
        if( present(u1) ) size1 = u1

        temp_array = array                  ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(size1))  ! Re-allocate with new size

        indices1 = [ (i,i=1,size(temp_array,1)) ]

        ! If the argument is present, define the new positions via the array
        if(present(p1)) indices1=p1

        ! Reassign the old values
        array(indices1) = temp_array

    end subroutine

    !> Reallocate a 2D array of doubles
    subroutine reallocate_2_d(array,u1,u2,p1,p2)
        implicit none
        !> Array to be reallocated
        double precision, dimension(:,:),allocatable, intent(inout) :: array
        !> New upper bound
        integer, intent(in),optional :: u1,u2
        !> new array configuration after reallocation
        integer, dimension(size(array,1)),intent(in),optional :: p1
        integer, dimension(size(array,2)),intent(in),optional :: p2

        integer :: size1,size2

        integer, dimension(size(array,1)) :: indices1
        integer, dimension(size(array,2)) :: indices2

        integer :: i

        double precision, dimension(size(array,1),size(array,2)) :: temp_array

        ! Default reallocation size is the size of the original array
        size1 = size(array,1)
        size2 = size(array,2)

        ! If the argument is present, then reallocate to the new size
        if( present(u1) ) size1 = u1
        if( present(u2) ) size2 = u2

        temp_array = array                  ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(size1,size2))  ! Re-allocate with new size

        indices1 = [ (i,i=1,size(temp_array,1)) ]
        indices2 = [ (i,i=1,size(temp_array,2)) ]

        ! If the argument is present, define the new positions via the array
        if(present(p1)) indices1=p1
        if(present(p2)) indices2=p2

        ! Reassign the old values
        array(indices1,indices2) = temp_array

    end subroutine

    !> Reallocate a 3D array of doubles
    subroutine reallocate_3_d(array,u1,u2,u3,p1,p2,p3)
        implicit none
        !> Array to be reallocated
        double precision, dimension(:,:,:),allocatable, intent(inout) :: array
        !> New upper bound
        integer, intent(in),optional :: u1,u2,u3
        !> new array configuration after reallocation
        integer, dimension(size(array,1)),intent(in),optional :: p1
        integer, dimension(size(array,2)),intent(in),optional :: p2
        integer, dimension(size(array,3)),intent(in),optional :: p3

        integer :: size1,size2,size3

        integer, dimension(size(array,1)) :: indices1
        integer, dimension(size(array,2)) :: indices2
        integer, dimension(size(array,3)) :: indices3

        integer :: i

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

        indices1 = [ (i,i=1,size(temp_array,1)) ]
        indices2 = [ (i,i=1,size(temp_array,2)) ]
        indices3 = [ (i,i=1,size(temp_array,3)) ]

        ! If the argument is present, define the new positions via the array
        if(present(p1)) indices1=p1
        if(present(p2)) indices2=p2
        if(present(p3)) indices3=p3

        ! Reassign the old values
        array(indices1,indices2,indices3) = temp_array

    end subroutine

    !> Reallocate a 1D array of integers
    subroutine reallocate_1_i(array,u1,p1)
        implicit none
        !> Array to be reallocated
        integer, dimension(:),allocatable, intent(inout) :: array
        !> New upper bound
        integer, intent(in),optional :: u1
        !> new array configuration after reallocation
        integer, dimension(size(array,1)),intent(in),optional :: p1

        integer :: size1

        integer, dimension(size(array,1)) :: indices1

        integer :: i

        integer, dimension(size(array,1)) :: temp_array

        ! Default reallocation size is the size of the original array
        size1 = size(array,1)

        ! If the argument is present, then reallocate to the new size
        if( present(u1) ) size1 = u1

        temp_array = array                  ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(size1))  ! Re-allocate with new size

        indices1 = [ (i,i=1,size(temp_array,1)) ]

        ! If the argument is present, define the new positions via the array
        if(present(p1)) indices1=p1

        ! Reassign the old values
        array(indices1) = temp_array

    end subroutine

    !> Reallocate a 2D array of integers
    subroutine reallocate_2_i(array,u1,u2,p1,p2)
        implicit none
        !> Array to be reallocated
        integer, dimension(:,:),allocatable, intent(inout) :: array
        !> New upper bound
        integer, intent(in),optional :: u1,u2
        !> new array configuration after reallocation
        integer, dimension(size(array,1)),intent(in),optional :: p1
        integer, dimension(size(array,2)),intent(in),optional :: p2

        integer :: size1,size2

        integer, dimension(size(array,1)) :: indices1
        integer, dimension(size(array,2)) :: indices2

        integer :: i

        integer, dimension(size(array,1),size(array,2)) :: temp_array

        ! Default reallocation size is the size of the original array
        size1 = size(array,1)
        size2 = size(array,2)

        ! If the argument is present, then reallocate to the new size
        if( present(u1) ) size1 = u1
        if( present(u2) ) size2 = u2

        temp_array = array                  ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(size1,size2))  ! Re-allocate with new size

        indices1 = [ (i,i=1,size(temp_array,1)) ]
        indices2 = [ (i,i=1,size(temp_array,2)) ]

        ! If the argument is present, define the new positions via the array
        if(present(p1)) indices1=p1
        if(present(p2)) indices2=p2

        ! Reassign the old values
        array(indices1,indices2) = temp_array

    end subroutine

    !> Reallocate a 3D array of integers
    subroutine reallocate_3_i(array,u1,u2,u3,p1,p2,p3)
        implicit none
        !> Array to be reallocated
        integer, dimension(:,:,:),allocatable, intent(inout) :: array
        !> New upper bound
        integer, intent(in),optional :: u1,u2,u3
        !> new array configuration after reallocation
        integer, dimension(size(array,1)),intent(in),optional :: p1
        integer, dimension(size(array,2)),intent(in),optional :: p2
        integer, dimension(size(array,3)),intent(in),optional :: p3

        integer :: size1,size2,size3

        integer, dimension(size(array,1)) :: indices1
        integer, dimension(size(array,2)) :: indices2
        integer, dimension(size(array,3)) :: indices3

        integer :: i

        integer, dimension(size(array,1),size(array,2),size(array,3)) :: temp_array

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

        indices1 = [ (i,i=1,size(temp_array,1)) ]
        indices2 = [ (i,i=1,size(temp_array,2)) ]
        indices3 = [ (i,i=1,size(temp_array,3)) ]

        ! If the argument is present, define the new positions via the array
        if(present(p1)) indices1=p1
        if(present(p2)) indices2=p2
        if(present(p3)) indices3=p3

        ! Reassign the old values
        array(indices1,indices2,indices3) = temp_array

    end subroutine






    subroutine add_point(point,array,narray,cluster_id)
        implicit none
        double precision, dimension(:), intent(in) :: point                     !> Point to be added to end of array
        double precision, dimension(:,:,:), allocatable, intent(inout) :: array !> Array to be added to
        integer,dimension(:), allocatable, intent(inout) :: narray              !> number of points in array (second index)
        integer, intent(in) :: cluster_id                                       !> cluster identity (third index)

        narray(cluster_id) = narray(cluster_id) + 1         ! Increase the number of points in the cluster

        ! If this takes us over the size of the array, then double it
        if(narray(cluster_id) > size(array,2) ) call reallocate_3_d(array, u2=size(array,2)*2 )

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
