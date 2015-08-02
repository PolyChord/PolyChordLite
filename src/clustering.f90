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
    recursive function NN_clustering(similarity_matrix,num_clusters) result(cluster_list)
        use utils_module, only: relabel
        use abort_module, only: halt_program
        implicit none

        double precision, intent(in), dimension(:,:) :: similarity_matrix

        integer, dimension(size(similarity_matrix,1)) :: cluster_list
        integer, intent(out):: num_clusters

        integer :: num_clusters_new

        integer, dimension(size(similarity_matrix,1)) :: cluster_list_old
        integer :: num_clusters_old

        integer, dimension(size(similarity_matrix,1),size(similarity_matrix,1)) :: knn

        integer :: k
        integer :: nlive
        integer :: n
        integer :: i_cluster
        integer :: i_point

        integer, allocatable, dimension(:) :: points

        ! Get the number of points to be clustered
        nlive=size(similarity_matrix,1)

        ! 10 degrees of separation are usually fine, we'll expand this if necessary
        k = min(nlive,10)

        ! compute the k nearest neighbors for each point
        knn(:k,:) = compute_knn(similarity_matrix,k)

        ! Set up the cluster list
        cluster_list_old = [( i_point,i_point=1,nlive )]
        num_clusters_old = nlive

        do n=2,k

            cluster_list = do_clustering_k( knn(:n,:) )     ! Get a raw clustering 

            ! Re-label the cluster list using integers 1,2,3,....
            cluster_list = relabel(cluster_list,num_clusters)

            if(num_clusters<=0) then
                call halt_program("Cluster error: cannot have fewer than 1 clusters")
            else if( num_clusters == 1 ) then
                return  ! If we're down to a single cluster, then just return
            else if( all( cluster_list == cluster_list_old ) ) then
                exit ! check that the clustering hasn't changed since the last pass
            else if(n==k) then
                ! If we need to cluster further, then expand the knn list
                k=min(k*2,nlive)
                knn(:k,:) = compute_knn(similarity_matrix,k)
            end if

            ! Save the old cluster list for later.
            cluster_list_old = cluster_list
            num_clusters_old = num_clusters

        end do


        ! If we've found clusters, then search within these clusters
        if(num_clusters>1) then
            i_cluster=1
            do while(i_cluster<=num_clusters)
                ! Get the indices of cluster i_cluster
                call get_indices_of_cluster(cluster_list,points,i_cluster)

                ! Call this function again on the similarity sub matrix, adding an offset
                cluster_list(points) = num_clusters + NN_clustering(similarity_matrix(points,points),num_clusters_new)

                ! If we didn't find any new clusters, then move on to the next one
                if(num_clusters_new==1) i_cluster=i_cluster+1

                ! re-label the clusters
                cluster_list = relabel(cluster_list,num_clusters)
            end do
        end if

    end function NN_clustering


    function do_clustering_k(knn) result(c)
        implicit none

        integer, dimension(:,:) :: knn
        integer, dimension(size(knn,2)) :: c

        integer :: i,j

        ! Set up the cluster list
        c = [( i,i=1,size(knn,2)  )]

        do i=1,size(knn,2)
            do j=i+1,size(knn,2)

                ! If they're not in the same cluster already...
                if(c(i)/=c(j)) then
                    ! ... check to see if they are within each others k nearest neihbors...
                    if( neighbors( knn(:,i),knn(:,j) ) ) then

                        ! If they are then relabel cluster_i and cluster_j to the smaller of the two
                        where(c==c(i).or.c==c(j)) 
                            c=min(c(i),c(j))
                        end where

                    end if
                end if

            end do
        end do

    end function do_clustering_k



    function compute_knn(similarity_matrix,k) result(knn)
        use utils_module, only: loginf
        implicit none

        !> The data to compute on
        double precision, intent(in),dimension(:,:) :: similarity_matrix

        !> The number of nearest neighbors to compute
        integer, intent(in) :: k

        ! The indices of the k nearest neighbors to output
        integer, dimension(k,size(similarity_matrix,1)) :: knn

        integer :: nPoints,i,j

        integer :: insert_index(1)

        double precision, dimension(k) :: distance2s

        nPoints = size(similarity_matrix,1)

        knn=0

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
                    distance2s(insert_index(1):) =  eoshift( distance2s(insert_index(1):), -1,dim=1,boundary=similarity_matrix(i,j))
                    knn(insert_index(1):,i) =  eoshift( knn(insert_index(1):,i) ,-1 ,dim=1, boundary=j)
                end if
            end do

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


    !> This subroutine returns points: an allocated array of indices of the points in cluster_list
    !! which belong to i_cluster
    subroutine get_indices_of_cluster(cluster_list,points,i_cluster)
        implicit none
        integer, dimension(:), intent(in)                 :: cluster_list
        integer, intent(out), allocatable, dimension(:)   :: points
        integer, intent(in)                               :: i_cluster

        integer :: nlive,npoints,j

        nlive = count(cluster_list==i_cluster)

        if(allocated(points)) deallocate(points)

        allocate( points(nlive) )

        npoints=0
        do j=1,size(cluster_list)
            if(cluster_list(j)==i_cluster) then
                npoints = npoints+1
                points(npoints) = j
            end if
        end do


    end subroutine get_indices_of_cluster



end module KNN_clustering











module cluster_module
    implicit none
    contains

    function do_clustering(settings,RTI,sub_dimensions)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info,add_cluster
        use calculate_module,  only: calculate_similarity_matrix
        use KNN_clustering,    only: NN_clustering
        implicit none

        !> Program settings
        type(program_settings), intent(in) :: settings
        !> The evidence storage
        type(run_time_info), intent(inout) :: RTI
        !> Dimensions to cluster on
        integer,dimension(:),optional,intent(in) :: sub_dimensions
        !> Whether or not a cluster has been found
        logical :: do_clustering


        ! Similarity matrix
        double precision,dimension(settings%nlive,settings%nlive) :: similarity_matrix
        integer,dimension(settings%nlive) :: clusters

        integer :: num_clusters
        integer :: num_old_clusters

        ! number of live points
        integer :: nlive

        integer :: i_cluster

        ! Initially no clusters found
        do_clustering=.false.

        ! Get the number of old clusters
        num_old_clusters=RTI%ncluster

        i_cluster = 1
        do while(i_cluster<=num_old_clusters)

            nlive = RTI%nlive(i_cluster) ! Get the number of live points in a temp variable

            if(nlive>2) then 
                ! Calculate the similarity matrix for this cluster
                if(present(sub_dimensions)) then
                    similarity_matrix(:nlive,:nlive) =&
                        calculate_similarity_matrix(RTI%live(sub_dimensions,:nlive,i_cluster))
                else 
                    similarity_matrix(:nlive,:nlive) =&
                        calculate_similarity_matrix(RTI%live(settings%h0:settings%h1,:nlive,i_cluster))
                end if

                clusters(:nlive) = NN_clustering(similarity_matrix(:nlive,:nlive),num_clusters)

                ! Do clustering on this 
                if ( num_clusters>1 ) then
                    ! If we find a cluster
                    do_clustering=.true.

                    ! we split this cluster into n, and put the others at the end
                    call add_cluster(settings,RTI,i_cluster,clusters(:nlive),num_clusters)

                else
                    i_cluster=i_cluster+1
                end if
            else
                i_cluster=i_cluster+1
            end if


            
        end do

    end function do_clustering



end module cluster_module












