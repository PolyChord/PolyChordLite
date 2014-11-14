module cluster_module

    implicit none
    contains


    function SNN_clustering(similarity_matrix,k,kt) result(cluster_list_final)
        implicit none
        double precision, intent(in), dimension(:,:) :: similarity_matrix
        integer, intent(in) :: k
        integer, intent(in) :: kt

        integer, dimension(size(similarity_matrix,1)) :: cluster_list
        integer, dimension(size(similarity_matrix,1)) :: cluster_list_final

        integer, dimension(k,size(similarity_matrix,1)) :: knn

        integer :: i_point,j_point,i_cluster

        integer :: cluster_orig_i,cluster_orig_j

        integer :: nlive

        integer :: num_clusters
        integer :: cluster_map(size(similarity_matrix,1))


        nlive=size(similarity_matrix,1)


        ! compute the k nearest neighbors for each point
        knn = compute_knn(similarity_matrix,k)

        ! Set up the cluster list
        cluster_list = [( i_point,i_point=1,nlive )]

        ! Loop through all pairs of points
        do i_point=1,nlive
            do j_point=i_point+1,nlive

                cluster_orig_i = cluster_list(i_point)
                cluster_orig_j = cluster_list(j_point)

                ! If they're not in the same cluster already...
                if(cluster_orig_i/=cluster_orig_j) then
                    ! ... check to see if they are each others nearest neihbors...
                    if( neighbors( knn(:,i_point),knn(:,j_point) ) ) then
                        ! ... and check to see if they share n nearest neighbors
                        if( matches( knn(:,i_point),knn(:,j_point) ) > kt ) then

                            if(cluster_orig_i>cluster_orig_j) then
                                where(cluster_list==cluster_orig_i) cluster_list=cluster_orig_j
                            else
                                where(cluster_list==cluster_orig_j) cluster_list=cluster_orig_i
                            end if
                        endif
                    end if
                end if

            end do
        end do

        ! Set clusters of the form 1,2,3
        num_clusters=1
        cluster_map(1) = cluster_list(1)

        do i_point=1,nlive
            if(all(cluster_list(i_point)/=cluster_map(1:num_clusters))) then
                num_clusters=num_clusters+1
                cluster_map(num_clusters) = cluster_list(i_point)
            end if
        end do

        do i_cluster=1,num_clusters
            where(cluster_list==cluster_map(i_cluster)) cluster_list_final=i_cluster
        end do

    end function SNN_clustering


    function NN_clustering(similarity_matrix,k,cluster_list_final) result(num_clusters)
        implicit none
        double precision, intent(in), dimension(:,:) :: similarity_matrix
        integer, intent(in) :: k
        integer, dimension(size(similarity_matrix,1)),intent(out) :: cluster_list_final

        integer :: num_clusters



        integer, dimension(size(similarity_matrix,1)) :: cluster_list

        integer, dimension(k,size(similarity_matrix,1)) :: knn

        integer :: i_point,j_point,i_cluster

        integer :: cluster_orig_i,cluster_orig_j

        integer :: nlive

        integer :: cluster_map(size(similarity_matrix,1))


        nlive=size(similarity_matrix,1)


        ! compute the k nearest neighbors for each point
        knn = compute_knn(similarity_matrix,k)

        ! Set up the cluster list
        cluster_list = [( i_point,i_point=1,nlive )]

        ! Loop through all pairs of points
        do i_point=1,nlive
            do j_point=i_point+1,nlive

                cluster_orig_i = cluster_list(i_point)
                cluster_orig_j = cluster_list(j_point)

                ! If they're not in the same cluster already...
                if(cluster_orig_i/=cluster_orig_j) then
                    ! ... check to see if they are within each others k nearest neihbors...
                    if( neighbors( knn(:k,i_point),knn(:k,j_point) ) ) then
                        if(cluster_orig_i>cluster_orig_j) then
                            where(cluster_list==cluster_orig_i) cluster_list=cluster_orig_j
                        else
                            where(cluster_list==cluster_orig_j) cluster_list=cluster_orig_i
                        end if
                    end if
                end if

            end do
        end do

        ! We now wish to relabel the cluster_list that we've got out, so that it
        ! uses integers 1,2,3,..., rather than a set of random integers between
        ! 1 and the number of points
        !
        ! The array cluster_map will detail this mapping between them

        ! Initialise the counter for the number of clusters
        num_clusters=1
        ! We will re-label the cluster type in cluster_list(1) with the integer 1
        cluster_map(num_clusters) = cluster_list(1)

        do i_point=1,nlive
            ! If the cluster type for i_point is not already included in the
            ! cluster_map, then add it
            if(all(cluster_list(i_point)/=cluster_map(1:num_clusters))) then
                num_clusters=num_clusters+1
                cluster_map(num_clusters) = cluster_list(i_point)
            end if
        end do

        ! cluster_map now contains the random integers that are found in cluster_list

        ! We now relabel according to 
        do i_cluster=1,num_clusters
            where(cluster_list==cluster_map(i_cluster)) cluster_list_final=i_cluster
        end do


    end function NN_clustering



    function compute_knn(similarity_matrix,n) result(knn)
        use utils_module, only: loginf
        implicit none

        !> The data to compute on
        double precision, intent(in),dimension(:,:) :: similarity_matrix
        !> The number of knn to compute
        integer,intent(in) :: n

        ! The indices of the k nearest neighbors
        integer, dimension(n,size(similarity_matrix,1)) :: knn

        integer :: nPoints,i,j

        integer :: insert_index(1)

        double precision, dimension(n) :: distance2s

        nPoints = size(similarity_matrix,1)
        knn=0


        ! Loop over each point
        do i=1,nPoints
            ! Find the k nearest neighbors for each point
            distance2s = loginf
            do j=1,nPoints
                ! If the distance between i and j is too large to be considered,
                ! this returns 0
                ! otherwise this returns the position to insert
                insert_index = minloc(distance2s, mask=distance2s>similarity_matrix(i,j))
                ! If it needs to be inserted, insert into both the integer
                ! array, and the local distance2s array
                if(insert_index(1)/=0) then
                    distance2s(insert_index(1):) =  eoshift( distance2s(insert_index(1):),  -1,dim=1,boundary=similarity_matrix(i,j))
                    knn(insert_index(1):,i) =  eoshift( knn(       insert_index(1):,i),-1,dim=1,boundary=j)
                end if
            end do
            if(knn(1,i)/=i) write(*,*) 'Catastrophic error'

        end do

    end function compute_knn


    ! Return whether they're each others n nearest neighbor list
    function neighbors(knn1,knn2) result(same_list)
        implicit none
        integer,intent(in), dimension(:) :: knn1
        integer,intent(in), dimension(:) :: knn2 

        logical :: same_list

        ! Check to see if they're in each others neighbors lists
        same_list= any(knn1==knn2(1)) .or. any(knn2==knn1(1))

    end function neighbors



    ! Return the number of matches in the n nearest neighbor list
    function matches(knn1,knn2)
        implicit none
        integer,intent(in), dimension(:) :: knn1
        integer,intent(in), dimension(:) :: knn2 

        integer :: matches

        integer i,j

        matches = count( [( [( knn1(i)==knn2(j), i=1,size(knn1) )], j=1,size(knn2) )] )

    end function matches


end module cluster_module
