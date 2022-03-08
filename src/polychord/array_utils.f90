module array_module
    use utils_module, only: dp

    interface reallocate
        module procedure reallocate_1_d, reallocate_2_d, reallocate_3_d, reallocate_1_i
    end interface reallocate
    contains

    !> Reallocate a 1D array of doubles
    !!
    !! The array is allocated with a new size defined by new_size1 (optional variable)
    !!
    !! The indices to be transferred across are defined by save_indices1 (optional variable)
    !!
    !! The positions of the data in the reallocated array are defined by target_indices1 (optional variable)
    subroutine reallocate_1_d(array,new_size1,save_indices1,target_indices1)
        use abort_module, only: halt_program
        implicit none
        real(dp), dimension(:),allocatable, intent(inout) :: array !> Array to be reallocated
        integer, intent(in),optional :: new_size1                          !> New size of the array
        integer, dimension(:),intent(in),optional :: save_indices1         !> indices of old array to be transferred
        integer, dimension(:),intent(in),optional :: target_indices1       !> new array indices after reallocation

        ! Temporary versions of the above variables
        real(dp), dimension(size(array,1)) :: a
        integer :: n1
        integer, allocatable, dimension(:) :: s1
        integer, allocatable, dimension(:) :: t1

        ! constructor variable
        integer :: i

        ! Check to see that it is already allocated
        if(.not.allocated(array)) call halt_program('reallocate_1_d error: array is not allocated') 

        ! Define the new size of the array
        if( present(new_size1) ) then
            n1 = new_size1     ! If the argument is present, then this is the reallocation size
        else
            n1 = size(array,1) ! Default reallocation size is the size of the original array
        end if

        ! Define the positions of where data is coming from and going to
        if(present(save_indices1) .and. present(target_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            if(size(s1)/=size(t1)) call halt_program('reallocate_1_d error: save and target indices must be equal in size') 
        else if(present(save_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(s1)))
            t1 = [ (i,i=1,size(s1)) ]
        else if(present(target_indices1)) then
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            allocate(s1(size(t1)))
            s1 = [ (i,i=1,size(t1)) ]
        else
            allocate(t1(size(array,1)),s1(size(array,1)))
            s1 = [ (i,i=1,size(array,1)) ]
            t1 = [ (i,i=1,size(array,1)) ]
        end if

        a = array                           ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(n1))                 ! Re-allocate with new size
        array = 0
        array(t1) = a(s1)                   ! Reassign the old values

    end subroutine reallocate_1_d


    !> Reallocate a 2D array of doubles
    !!
    !! The array is allocated with a new size defined by new_size1,new_size2 (optional variables)
    !!
    !! The indices to be transferred across are defined by save_indices1,save_indices2 (optional variables)
    !!
    !! The positions of the data in the reallocated array are defined by target_indices1,target_indices2, (optional variables)
    subroutine reallocate_2_d(array,new_size1,new_size2,save_indices1,save_indices2,target_indices1,target_indices2)
        use abort_module, only: halt_program
        implicit none
        !> Array to be reallocated
        real(dp), dimension(:,:),allocatable, intent(inout) :: array
        !> New size of the array 
        integer, intent(in),optional :: new_size1,new_size2
        !> indices of old array to be transferred
        integer, dimension(:),intent(in),optional :: save_indices1,save_indices2
        !> new array indices after reallocation 
        integer, dimension(:),intent(in),optional :: target_indices1,target_indices2

        ! Temporary versions of the above variables
        real(dp), dimension(size(array,1),size(array,2)) :: a
        integer :: n1,n2
        integer, allocatable, dimension(:) :: s1,s2
        integer, allocatable, dimension(:) :: t1,t2

        ! constructor variable
        integer :: i

        ! Check to see that it is already allocated
        if(.not.allocated(array)) call halt_program('reallocate_2_d error: array is not allocated') 

        ! Define the new size of the array
        if( present(new_size1) ) then
            n1 = new_size1     ! If the argument is present, then this is the reallocation size
        else
            n1 = size(array,1) ! Default reallocation size is the size of the original array
        end if
        if( present(new_size2) ) then
            n2 = new_size2     ! If the argument is present, then this is the reallocation size
        else
            n2 = size(array,2) ! Default reallocation size is the size of the original array
        end if

        ! Define the positions of where data is coming from and going to
        if(present(save_indices1) .and. present(target_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            if(size(s1)/=size(t1)) call halt_program('reallocate_2_d error: save and target indices must be equal in size') 
        else if(present(save_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(s1)))
            t1 = [ (i,i=1,size(s1)) ]
        else if(present(target_indices1)) then
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            allocate(s1(size(t1)))
            s1 = [ (i,i=1,size(t1)) ]
        else
            allocate(t1(size(array,1)),s1(size(array,1)))
            s1 = [ (i,i=1,size(array,1)) ]
            t1 = [ (i,i=1,size(array,1)) ]
        end if


        ! Define the positions of where data is coming from and going to
        if(present(save_indices2) .and. present(target_indices2)) then
            allocate(s2(size(save_indices2)))
            s2 = save_indices2
            allocate(t2(size(target_indices2)))
            t2 = target_indices2
            if(size(s2)/=size(t2)) call halt_program('reallocate_2_d error: save and target indices must be equal in size') 
        else if(present(save_indices2)) then
            allocate(s2(size(save_indices2)))
            s2 = save_indices2
            allocate(t2(size(s2)))
            t2 = [ (i,i=1,size(s2)) ]
        else if(present(target_indices2)) then
            allocate(t2(size(target_indices2)))
            t2 = target_indices2
            allocate(s2(size(t2)))
            s2 = [ (i,i=1,size(t2)) ]
        else
            allocate(t2(size(array,2)),s2(size(array,2)))
            s2 = [ (i,i=1,size(array,2)) ] 
            t2 = [ (i,i=1,size(array,2)) ] 
        end if

        a = array                           ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(n1,n2))              ! Re-allocate with new size
        array = 0
        array(t1,t2)= a(s1,s2)              ! Reassign the old values

    end subroutine reallocate_2_d


    !> Reallocate a 3D array of doubles
    !!
    !! The array is allocated with a new size defined by new_size1,new_size2,new_size3 (optional variables)
    !!
    !! The indices to be transferred across are defined by save_indices1,save_indices2,save_indices3 (optional variables)
    !!
    !! The positions of the data in the reallocated array are defined by target_indices1,target_indices2,target_indices3 (optional variables)
    subroutine reallocate_3_d(array,&
                                new_size1,new_size2,new_size3,&
                                save_indices1,save_indices2,save_indices3,&
                                target_indices1,target_indices2,target_indices3)
        use abort_module, only: halt_program
        implicit none
        !> Array to be reallocated
        real(dp), dimension(:,:,:),allocatable, intent(inout) :: array
        !> New size of the array 
        integer, intent(in),optional :: new_size1,new_size2,new_size3
        !> indices of old array to be transferred
        integer, dimension(:),intent(in),optional :: save_indices1,save_indices2,save_indices3
        !> new array indices after reallocation 
        integer, dimension(:),intent(in),optional :: target_indices1,target_indices2,target_indices3

        ! Temporary versions of the above variables
        real(dp), dimension(size(array,1),size(array,2),size(array,3)) :: a
        integer :: n1,n2,n3
        integer, allocatable, dimension(:) :: s1,s2,s3
        integer, allocatable, dimension(:) :: t1,t2,t3

        ! constructor variable
        integer :: i

        ! Check to see that it is already allocated
        if(.not.allocated(array)) call halt_program('reallocate_3_d error: array is not allocated') 

        ! Define the new size of the array
        if( present(new_size1) ) then
            n1 = new_size1     ! If the argument is present, then this is the reallocation size
        else
            n1 = size(array,1) ! Default reallocation size is the size of the original array
        end if
        if( present(new_size2) ) then
            n2 = new_size2     ! If the argument is present, then this is the reallocation size
        else
            n2 = size(array,2) ! Default reallocation size is the size of the original array
        end if
        if( present(new_size3) ) then
            n3 = new_size3     ! If the argument is present, then this is the reallocation size
        else
            n3 = size(array,3) ! Default reallocation size is the size of the original array
        end if

        ! Define the positions of where data is coming from and going to
        if(present(save_indices1) .and. present(target_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            if(size(s1)/=size(t1)) call halt_program('reallocate_3_d error: save and target indices must be equal in size') 
        else if(present(save_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(s1)))
            t1 = [ (i,i=1,size(s1)) ]
        else if(present(target_indices1)) then
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            allocate(s1(size(t1)))
            s1 = [ (i,i=1,size(t1)) ]
        else
            allocate(t1(size(array,1)),s1(size(array,1)))
            s1 = [ (i,i=1,size(array,1)) ]
            t1 = [ (i,i=1,size(array,1)) ]
        end if


        ! Define the positions of where data is coming from and going to
        if(present(save_indices2) .and. present(target_indices2)) then
            allocate(s2(size(save_indices2)))
            s2 = save_indices2
            allocate(t2(size(target_indices2)))
            t2 = target_indices2
            if(size(s2)/=size(t2)) call halt_program('reallocate_3_d error: save and target indices must be equal in size') 
        else if(present(save_indices2)) then
            allocate(s2(size(save_indices2)))
            s2 = save_indices2
            allocate(t2(size(s2)))
            t2 = [ (i,i=1,size(s2)) ]
        else if(present(target_indices2)) then
            allocate(t2(size(target_indices2)))
            t2 = target_indices2
            allocate(s2(size(t2)))
            s2 = [ (i,i=1,size(t2)) ]
        else
            allocate(t2(size(array,2)),s2(size(array,2)))
            s2 = [ (i,i=1,size(array,2)) ] 
            t2 = [ (i,i=1,size(array,2)) ] 
        end if


        ! Define the positions of where data is coming from and going to
        if(present(save_indices3) .and. present(target_indices3)) then
            allocate(s3(size(save_indices3)))
            s3 = save_indices3
            allocate(t3(size(target_indices3)))
            t3 = target_indices3
            if(size(s3)/=size(t3)) call halt_program('reallocate_3_d error: save and target indices must be equal in size') 
        else if(present(save_indices3)) then
            allocate(s3(size(save_indices3)))
            s3 = save_indices3
            allocate(t3(size(s3)))
            t3 = [ (i,i=1,size(s3)) ]
        else if(present(target_indices3)) then
            allocate(t3(size(target_indices3)))
            t3 = target_indices3
            allocate(s3(size(t3)))
            s3 = [ (i,i=1,size(t3)) ]

        else
            allocate(t3(size(array,3)),s3(size(array,3)))
            s3 = [ (i,i=1,size(array,3)) ] 
            t3 = [ (i,i=1,size(array,3)) ] 
        end if


        a = array                           ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(n1,n2,n3))           ! Re-allocate with new size
        array = 0
        array(t1,t2,t3)= a(s1,s2,s3)        ! Reassign the old values

    end subroutine reallocate_3_d




    !> Reallocate a 1D array of integers
    !!
    !! The array is allocated with a new size defined by new_size1 (optional variable)
    !!
    !! The indices to be transferred across are defined by save_indices1 (optional variable)
    !!
    !! The positions of the data in the reallocated array are defined by target_indices1 (optional variable)
    subroutine reallocate_1_i(array,new_size1,save_indices1,target_indices1)
        use abort_module, only: halt_program
        implicit none
        integer,          dimension(:),allocatable, intent(inout) :: array !> Array to be reallocated
        integer, intent(in),optional :: new_size1                          !> New size of the array
        integer, dimension(:),intent(in),optional :: save_indices1         !> indices of old array to be transferred
        integer, dimension(:),intent(in),optional :: target_indices1       !> new array indices after reallocation

        ! Temporary versions of the above variables
        integer, dimension(size(array,1)) :: a
        integer :: n1
        integer, allocatable, dimension(:) :: s1
        integer, allocatable, dimension(:) :: t1

        ! constructor variable
        integer :: i

        ! Check to see that it is already allocated
        if(.not.allocated(array)) call halt_program('reallocate_1_i error: array is not allocated') 

        ! Define the new size of the array
        if( present(new_size1) ) then
            n1 = new_size1     ! If the argument is present, then this is the reallocation size
        else
            n1 = size(array,1) ! Default reallocation size is the size of the original array
        end if

        ! Define the positions of where data is coming from and going to
        if(present(save_indices1) .and. present(target_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            if(size(s1)/=size(t1)) call halt_program('reallocate_1_i error: save and target indices must be equal in size') 
        else if(present(save_indices1)) then
            allocate(s1(size(save_indices1)))
            s1 = save_indices1
            allocate(t1(size(s1)))
            t1 = [ (i,i=1,size(s1)) ]
        else if(present(target_indices1)) then
            allocate(t1(size(target_indices1)))
            t1 = target_indices1
            allocate(s1(size(t1)))
            s1 = [ (i,i=1,size(t1)) ]
        else
            allocate(t1(size(array,1)),s1(size(array,1)))
            s1 = [ (i,i=1,size(array,1)) ]
            t1 = [ (i,i=1,size(array,1)) ]
        end if

        a = array                           ! Save the old array 
        deallocate(array)                   ! Deallocate it      
        allocate(array(n1))                 ! Re-allocate with new size
        array = 0
        array(t1) = a(s1)                   ! Reassign the old values

    end subroutine reallocate_1_i

    subroutine add_point(point,array,narray,nfixed)
        implicit none
        real(dp), dimension(:), intent(in) :: point                   !> Point to be added to end of array
        real(dp), dimension(:,:), allocatable, intent(inout) :: array !> Array to be added to
        integer, intent(inout) :: narray                              !> number of points in array 
        integer, intent(in), optional :: nfixed                       !> number of points in array (unincrementable -- useful for clusters)

        narray = narray + 1         ! Increase the number of points
        if (present(nfixed)) then
            ! If this takes us over the size of the array, then double it
            if(nfixed+1 > size(array,2) ) call reallocate(array, new_size2=max(1,size(array,2)*2) )
            array(:,nfixed+1) = point       ! Add the point to the end position
        else
            ! If this takes us over the size of the array, then double it
            if(narray > size(array,2) ) call reallocate(array, new_size2=max(1,size(array,2)*2) )
            array(:,narray) = point       ! Add the point to the end position
        end if

    end subroutine add_point


    function delete_point(i_point,array,narray,nfixed) result(point)
        implicit none
        integer, intent(in) :: i_point                                  !> Position of point to be deleted from array
        real(dp), dimension(:,:), allocatable, intent(inout) :: array   !> Array to be delete from
        integer, intent(inout) :: narray                                !> number of points in array (second index)
        integer, intent(in), optional :: nfixed                         !> number of points in array (undecrementable -- useful for clusters)
        real(dp), dimension(size(array,1)) :: point                     ! The point we have just deleted

        point = array(:,i_point)           ! Output the point to be deleted
        if (present(nfixed)) then
            array(:,i_point) = array(:,nfixed) ! delete the point by overwriting it with the point at the end 
            array(:,nfixed) = 0
        else
            array(:,i_point) = array(:,narray) ! delete the point by overwriting it with the point at the end 
            array(:,narray) = 0
        end if
        narray = narray - 1                ! reduce the number of points in the cluster

    end function delete_point


    function concat(a, b) result(c)
        implicit none
        real(dp), dimension(:,:), intent(in) :: a,b
        real(dp), dimension(size(a,1),size(a,2)+size(b,2)) :: c

        c(:,:size(a,2)) = a
        c(:,size(a,2)+1:) = b

    end function concat


    function sel(mask) 
        implicit none
        logical, dimension(:) :: mask
        integer i
        integer, dimension(count(mask)) :: sel
        sel = pack([(i,i=1,size(mask))],mask)
    end function


end module array_module
