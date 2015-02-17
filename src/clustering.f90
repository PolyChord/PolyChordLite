module KNN_clustering

    implicit none
    contains

    !> This function returns a clustering from a similarity matrix based on
    !! 'nearest neighbor' clustering.
    !!
    !! Points belong to the same cluster if they are in either of each others k
    !! nearest neighbor sets. 
    !!
    !! The algorithm computes the k nearest neihbor sets from the similarity
    !! matrix, and then tests
    recursive function NN_clustering(similarity_matrix,cluster_list_output) result(num_clusters_output)
        implicit none

        double precision, intent(in), dimension(:,:) :: similarity_matrix

        integer, dimension(size(similarity_matrix,1)),intent(out) :: cluster_list_output
        integer :: num_clusters_output


        integer, dimension(size(similarity_matrix,1),size(similarity_matrix,1)) :: cluster_list
        integer, dimension(size(similarity_matrix,1)) :: num_clusters

        integer, dimension(size(similarity_matrix,1),size(similarity_matrix,1)) :: knn

        integer :: nlive
        integer :: n


        ! Get the number of points to be clustered
        nlive=size(similarity_matrix,1)

        ! compute the k nearest neighbors for each point
        knn = compute_knn(similarity_matrix)

        ! Loop through all pairs of points
        do n=2,nlive

            ! Get a raw clustering 
            cluster_list(:,n) = do_clustering_k( knn(:,:n) )

            ! Re-label the cluster list using integers 1,2,3,....
            cluster_list(:,n) = relabel_clustering(cluster_list(:,n),1,num_clusters(n))

            if(num_clusters(n) == 1 ) then

                exit  ! If we're down to a single cluster, then just return

            else if(n>3) then

                ! Otherwise, check that the clustering hasn't changed for the past 3 passes.
                if( num_clusters(n) == num_clusters(n-1) .and. num_clusters(n-1) == num_clusters(n-2) ) then
                    if( all( cluster_list(:,n) == cluster_list(:,n-1) .and. cluster_list(:,n-1) == cluster_list(:,n-2) ) ) exit
                end if

            end if

        end do

        num_clusters_output = num_clusters_output + NN_clustering(similarity_matrix,cluster_list_output)

        num_clusters_output = num_clusters(n)
        cluster_list_output = cluster_list(:,n)


    end function NN_clustering


    function do_clustering_k(knn) result(cluster)
        implicit none

        integer, dimension(:,:) :: knn
        integer, dimension(size(knn,1)) :: cluster

        integer :: nlive

        integer :: cluster_i,cluster_j
        integer :: i_point,j_point

        nlive = size(knn,1)

        ! Set up the cluster list
        cluster = [( i_point,i_point=1,nlive )]

        do i_point=1,nlive
            do j_point=i_point+1,nlive

                cluster_i = cluster(i_point)
                cluster_j = cluster(j_point)

                ! If they're not in the same cluster already...
                if(cluster_i/=cluster_j) then
                    ! ... check to see if they are within each others k nearest neihbors...
                    if( neighbors( knn(:,i_point),knn(:,j_point) ) ) then

                        if(cluster_i>cluster_j) then
                            where(cluster==cluster_i) cluster=cluster_j
                        else
                            where(cluster==cluster_j) cluster=cluster_i
                        end if

                    end if
                end if

            end do
        end do

    end function do_clustering_k


    function relabel_clustering(cluster,i_start,i_finish) result(cluster_relabel)
        implicit none
        integer,intent(in),dimension(:)  :: cluster
        integer,intent(in)               :: i_start
        integer,intent(out)              :: i_finish

        integer,dimension(size(cluster)) :: cluster_relabel

        integer,dimension(size(cluster)) :: cluster_map

        integer :: npoints
        integer :: i_point
        integer :: i_cluster

        ! Find the number of points
        npoints = size(cluster)

        ! We will re-label the cluster type in cluster(1) with the integer 1
        cluster_map(1) = cluster(1)
        i_finish = 1

        do i_point=1,npoints
            ! If the cluster type for i_point is not already included in the
            ! cluster_map, then add it
            if( all(cluster(i_point)/=cluster_map(1:i_finish)) ) then
                i_finish=i_finish+1
                cluster_map(i_finish) = cluster(i_point)
            end if
        end do

        ! cluster_map now contains the random integers that are found in cluster

        ! We now relabel according to 
        do i_cluster=1,i_finish
            where(cluster==cluster_map(i_cluster)) cluster_relabel=i_cluster
        end do

        ! Add the offset of i_start
        cluster_relabel = cluster_relabel+i_start-1
        i_finish        = i_finish       +i_start-1

    end function

    function compute_knn(similarity_matrix) result(knn)
        use utils_module, only: loginf
        implicit none

        !> The data to compute on
        double precision, intent(in),dimension(:,:) :: similarity_matrix

        ! The indices of the k nearest neighbors
        integer, dimension(size(similarity_matrix,1),size(similarity_matrix,1)) :: knn

        integer :: nPoints,i,j

        integer :: insert_index(1)

        double precision, dimension(size(similarity_matrix,1)) :: distance2s

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


end module KNN_clustering











module cluster_module
    implicit none
    contains

    subroutine do_clustering(settings,RTI)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info
        use calculate_module,  only: calculate_similarity_matrix
        use KNN_clustering,    only: NN_clustering
        implicit none

        !> Program settings
        type(program_settings), intent(in) :: settings
        !> The evidence storage
        type(run_time_info), intent(in) :: RTI

        ! Similarity matrix
        double precision,dimension(settings%nlive,settings%nlive) :: similarity_matrix
        integer,dimension(settings%nlive) :: clusters

        ! number of live points
        integer :: nlive

        integer :: i_cluster


        i_cluster = 1
        do while(i_cluster<=RTI%ncluster)

            nlive = RTI%nlive(i_cluster) ! Get the number of live points in a temp variable

            ! Calculate the similarity matrix for this cluster
            similarity_matrix(:nlive,:nlive) = calculate_similarity_matrix(RTI%live(settings%h0:settings%h1,:,i_cluster))

            ! Do clustering on this 
            if ( NN_clustering(similarity_matrix(:nlive,:nlive),clusters(:nlive)) > 1 ) then
                ! If we find a cluster

                ! ... do re-organising ...
                ! we split this cluster into n, and move all the other
                ! clusters and files up by n-1
                
            end if


            
        end do


    end subroutine do_clustering



end module cluster_module












