module KNN_clustering
    use utils_module, only: dp

    implicit none
    contains

    !> This function returns a clustering from a distance^2 matrix based on
    !! 'nearest neighbour' clustering.
    !!
    !! Points belong to the same cluster if they are in either of each others k
    !! nearest neighbour sets. 
    !!
    !! The algorithm computes the k nearest neighbor sets from the distance^2
    !! matrix, and then tests
    recursive function NN_clustering(distance2_matrix) result(cluster_list)
        use utils_module, only: relabel
        use abort_module, only: halt_program
        implicit none

        real(dp), intent(in), dimension(:,:) :: distance2_matrix

        integer, dimension(size(distance2_matrix,1)) :: cluster_list

        integer :: num_clusters_new

        integer, dimension(size(distance2_matrix,1)) :: cluster_list_old
        integer :: num_clusters_old

        integer, dimension(size(distance2_matrix,1),size(distance2_matrix,1)) :: knn

        integer :: k
        integer :: nlive
        integer :: n
        integer :: i_cluster
        integer :: i_point
        integer :: num_clusters

        integer, allocatable, dimension(:) :: points

        ! Get the number of points to be clustered
        nlive=size(distance2_matrix,1)

        ! 10 degrees of separation are usually fine, we'll expand this if necessary
        k = min(nlive,10)

        ! compute the k nearest neighbours for each point
        knn(:k,:) = compute_knn(distance2_matrix,k)

        ! Set up the cluster list
        cluster_list_old = [( i_point,i_point=1,nlive )]

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
                knn(:k,:) = compute_knn(distance2_matrix,k)
            end if

            ! Save the old cluster list for later.
            cluster_list_old = cluster_list

        end do


        ! If we've found clusters, then search within these clusters
        if(num_clusters>1) then
            i_cluster=1
            do while(i_cluster<=num_clusters)
                ! Get the indices of cluster i_cluster
                call get_indices_of_cluster(cluster_list,points,i_cluster)

                ! Call this function again on the distance^2 sub matrix, adding an offset
                cluster_list(points) = num_clusters + NN_clustering(distance2_matrix(points,points))
                num_clusters_new = maxval(cluster_list(points))-num_clusters

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
                    ! ... check to see if they are within each others k nearest neighbours...
                    if( neighbours( knn(:,i),knn(:,j) ) ) then

                        ! If they are then relabel cluster_i and cluster_j to the smaller of the two
                        where(c==c(i).or.c==c(j)) 
                            c=min(c(i),c(j))
                        end where

                    end if
                end if

            end do
        end do

    end function do_clustering_k



    function compute_knn(distance2_matrix,k) result(knn)
        implicit none

        !> The data to compute on
        real(dp), intent(in),dimension(:,:) :: distance2_matrix

        !> The number of nearest neighbours to compute
        integer, intent(in) :: k

        ! The indices of the k nearest neighbours to output
        integer, dimension(k,size(distance2_matrix,1)) :: knn

        integer :: nPoints,i,j

        integer :: insert_index(1)

        real(dp), dimension(k) :: distance2s

        nPoints = size(distance2_matrix,1)

        knn=0

        do i=1,nPoints
            ! Find the k nearest neighbours for each point
            distance2s = huge(1d0)
            do j=1,nPoints
                ! If the distance between i and j is too large to be considered,
                ! this returns 0
                ! otherwise this returns the position to insert
                insert_index = minloc(distance2s, mask=distance2s>distance2_matrix(i,j))
                ! If it needs to be inserted, insert into both the integer
                ! array, and the local distance2s array
                if(insert_index(1)/=0) then
                    distance2s(insert_index(1):) =  eoshift( distance2s(insert_index(1):), -1,dim=1,boundary=distance2_matrix(i,j))
                    knn(insert_index(1):,i) =  eoshift( knn(insert_index(1):,i) ,-1 ,dim=1, boundary=j)
                end if
            end do

        end do

    end function compute_knn


    ! Return whether they're each others n nearest neighbour list
    function neighbours(knn1,knn2) result(same_list)
        implicit none
        integer,intent(in), dimension(:) :: knn1
        integer,intent(in), dimension(:) :: knn2 

        logical :: same_list

        ! Check to see if they're in each others neighbours lists
        same_list= any(knn1==knn2(1)) .or. any(knn2==knn1(1))

    end function neighbours



    ! Return the number of matches in the n nearest neighbour list
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
    use utils_module, only: dp
    implicit none
    contains

    function do_clustering(settings,RTI,cluster,sub_dimensions)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info,add_cluster
        use calculate_module,  only: calculate_distance2_matrix
        use KNN_clustering,    only: NN_clustering
        implicit none

        !> Program settings
        type(program_settings), intent(in) :: settings
        !> The evidence storage
        type(run_time_info), intent(inout) :: RTI
        !> Dimensions to cluster on
        integer,dimension(:),optional,intent(in) :: sub_dimensions

        interface
            function cluster(points) result(cluster_list)
                import :: dp
                real(dp), intent(in), dimension(:,:) :: points
                integer, dimension(size(points,2)) :: cluster_list
            end function
        end interface

        !> Whether or not a cluster has been found
        logical :: do_clustering


        ! distance^2 matrix
        real(dp),dimension(settings%nDims,sum(RTI%nlive)) :: points
        integer,dimension(sum(RTI%nlive)) :: clusters

        integer :: num_clusters
        integer :: num_old_clusters

        ! number of live points
        integer :: nlive
        ! number of dimensions
        integer :: nDims

        integer :: i_cluster

        ! Initially no clusters found
        do_clustering=.false.

        ! Get the number of old clusters
        num_old_clusters=RTI%ncluster

        ! get the number of dimensions
        nDims = settings%nDims

        i_cluster = 1
        do while(i_cluster<=num_old_clusters)

            nlive = RTI%nlive(i_cluster) ! Get the number of live points in a temp variable

            if(nlive>2) then 
                ! get the position matrix for this cluster
                if(present(sub_dimensions)) then
                    points(:nDims, :nlive) =&
                        RTI%live(sub_dimensions,:nlive,i_cluster)
                else 
                    points(:nDims, :nlive) =&
                        RTI%live(settings%h0:settings%h1,:nlive,i_cluster)
                end if

                clusters(:nlive) = cluster(points(:nDims, :nlive))

                ! default to KNN clustering
                if (any(clusters(:nlive)<=0)) then
                    clusters(:nlive) = NN_clustering(&
                        calculate_distance2_matrix(points(:nDims, :nlive)))
                end if
                num_clusters = maxval(clusters(:nlive))

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












