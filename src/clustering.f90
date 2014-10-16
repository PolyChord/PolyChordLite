module cluster_module

    implicit none
    contains


    subroutine SNN_clustering(settings,live_points,stack_size)
        use settings_module, only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        integer, intent(in) :: stack_size
        double precision, dimension(settings%nTotal,stack_size),intent(inout) :: live_points

        integer, dimension(settings%SNN_k+1) :: knn

        integer, dimension(stack_size) :: cluster_list
        integer, dimension(stack_size) :: cluster_list_final

        integer :: i_point,j_point

        integer :: cluster_orig_i,cluster_orig_j

        integer, dimension(settings%SNN_k+1) :: cluster_indices
        integer :: max_cluster_index
        integer :: new_cluster_index

        integer :: clusters(10)
        integer :: num_clusters

        integer :: i_cluster
        integer :: num_computes

        ! Initialise the cluster list
        cluster_list = 0 !(/ (i_point,i_point=1,stack_size) /)

        max_cluster_index=0

        num_computes = 0

        ! Loop through all pairs of points
        do i_point=1,stack_size

            ! If the point is unclustered
            if(cluster_list(i_point)==0) then

                ! Find the k nearest neighbors to this point
                knn = compute_knn(settings%SNN_k,i_point,live_points(settings%h0:settings%h1,:stack_size))
                num_computes=num_computes+1

                cluster_indices = cluster_list(knn)

                ! Find the lowest cluster index
                new_cluster_index = minval(cluster_indices)

                ! If none of them are assigned clusters...
                if(new_cluster_index==0) then
                    !... Create a new cluster index
                    max_cluster_index=max_cluster_index+1
                    new_cluster_index = max_cluster_index
                end if

                cluster_list(knn) = new_cluster_index

                ! Set the unassigned cluster indices to the new index
                where(cluster_indices==0) cluster_indices=new_cluster_index

                do j_point=1,stack_size
                    if(any(cluster_list(j_point)==cluster_indices)) cluster_list(j_point)=new_cluster_index
                end do

            end if

        end do

        ! Find the cluster
        num_clusters=1
        clusters(1) = cluster_list(1)

        do i_point=1,stack_size
            if(all(cluster_list(i_point)/=clusters(1:num_clusters))) then
                num_clusters=num_clusters+1
                clusters(num_clusters) = cluster_list(i_point)
            end if
        end do

        do i_cluster=1,num_clusters
            where(cluster_list==clusters(i_cluster)) cluster_list_final=i_cluster
        end do


        write(*,'(" num_clusters: ", I4)') num_clusters
        write(*,'(" num_computes: ", I4, "/", I6)') num_computes,stack_size

        ! Pass on the new cluster list
        live_points(settings%cluster,:) = cluster_list_final



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


    function same_cluster(knn1,knn2)
        implicit none
        integer,intent(in), dimension(:)          :: knn1
        integer,intent(in), dimension(size(knn1)) :: knn2

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

        same_cluster=.true.


    end function same_cluster


end module cluster_module
