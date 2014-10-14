module cluster_module

    implicit none
    contains


    subroutine SNN_clustering(settings,live_points,stack_size)
        use settings_module, only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer, intent(in) :: stack_size
        double precision, dimension(settings%nTotal,stack_size),intent(inout) :: live_points

        integer, dimension(settings%SNN_k+1,stack_size) :: knn

        integer, dimension(stack_size) :: cluster_list

        integer :: i_point,j_point

        integer :: cluster_orig_i,cluster_orig_j

        integer :: clusters(10)

        ! Begin by computing the k nearest neighbors for each point
        do i_point=1,stack_size
            knn(:,i_point) = compute_knn(settings%SNN_k,i_point,live_points(settings%h0:settings%h1,:stack_size))
        end do

        ! Set up the cluster list
        cluster_list = (/ (i_point,i_point=1,stack_size) /)

        ! Loop through all pairs of points
        ! But skip points that are already in the same cluster
        do i_point=1,stack_size
            do j_point=i_point,stack_size
                cluster_orig_i = cluster_list(i_point)
                cluster_orig_j = cluster_list(j_point)

                ! If they're not in the same cluster already...
                if(cluster_orig_i/=cluster_orig_j) then
                    ! ... check to see if they share kt nearest neighbors and each other
                    if(same_cluster(knn(:,i_point),knn(:,j_point),settings%SNN_kt) ) then
                        ! If they do, then set all of the points in the higher
                        ! cluster type to the points in the lower cluster type
                        if(cluster_orig_i>cluster_orig_j) then
                            where(cluster_list==cluster_orig_i) cluster_list=cluster_orig_j
                        else
                            where(cluster_list==cluster_orig_j) cluster_list=cluster_orig_i
                        end if
                    end if
                end if

            end do
        end do


        write(*,'(<size(cluster_list)>I5)') cluster_list

        ! Pass on the new cluster list
        live_points(settings%cluster,:) = cluster_list



    end subroutine SNN_clustering




    function compute_knn(k,i,data_array) result(knn)
        implicit none

        !> The length of vector to be computed
        integer,intent(in) :: k
        !> The index of the vector to be compared to
        integer,intent(in) :: i
        !> The data to compute on
        double precision, intent(in),dimension(:,:) :: data_array

        ! The indices of the k nearest neighbors
        integer, dimension(k+1) :: knn

        integer :: j
        double precision :: distance2
        double precision,dimension(k) :: distance2s

        integer :: insert_index(1)


        distance2s = huge(1d0)

        do j = 1,size(data_array,2)
            ! Find the squared distance between the two
            distance2 = dot_product(data_array(:,i)-data_array(:,j),data_array(:,i)-data_array(:,j))

            ! We want to insert this j into the the relevent arrays;
            ! Minloc finds the lowest index satisfying mask
            ! We want to find 
            insert_index = minloc(distance2s,mask=distance2s>distance2)


            ! If this point belongs in the array, then update the array
            if(insert_index(1)/=0) then
                distance2s(insert_index(1):) = (/ distance2, distance2s(insert_index(1):) /)
                knn(insert_index(1):) = (/ j, knn(insert_index(1):) /)
            end if

        end do

    end function compute_knn


    function same_cluster(knn1,knn2,threshold)
        implicit none
        integer,intent(in), dimension(:)          :: knn1
        integer,intent(in), dimension(size(knn1)) :: knn2
        integer,intent(in)                        :: threshold

        logical :: same_cluster

        integer i,j
        integer :: matches

        integer :: same_list(1)

        ! Check to see if they're in each others neigbors lists
        same_cluster=.false.

        same_list = minloc(knn1,mask=knn1==knn2(1))
        if(same_list(1)==0) return
        same_list = minloc(knn2,mask=knn2==knn1(1))
        if(same_list(1)==0) return

        matches=0
        do i=1,size(knn1)
            do j=1,size(knn1)
                if(knn1(i)==knn2(j)) matches = matches+1
                if(matches>=threshold) then
                    same_cluster=.true.
                    return
                end if
            end do
        end do



    end function same_cluster


end module cluster_module
