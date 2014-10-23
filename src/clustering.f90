module cluster_module

    implicit none
    contains




    subroutine Skilling_clustering(settings,live_points,cholesky)
        use settings_module, only: program_settings
        use utils_module, only: logzero,distance2
        use random_module, only: random_distinct_integers

        implicit none

        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nTotal,settings%nlive),intent(inout) :: live_points
        double precision, dimension(settings%nDims,settings%nDims), intent(in) :: cholesky

        double precision :: similarity_matrix(settings%nlive,settings%nlive)
        integer :: i,j

        double precision :: logsum_d
        double precision :: logsum_m
        double precision :: logsum_m_new
        double precision :: logprod_d
        double precision :: logprod_d_new
        double precision :: threshold

        integer :: numsub,numadd 
        double precision :: logprod_dadd,logprod_dsub
        logical :: keep_going

        double precision, parameter :: max_clusters=sqrt(1.5d0)

        logical, dimension(settings%nlive) :: cluster
        
        integer, dimension(settings%nlive) :: random_deck

        double precision, dimension(settings%nDims) :: displacement
        double precision, dimension(settings%nDims,settings%nDims) :: invcovmat
        !double precision, dimension(settings%nDims,settings%nDims) :: covmat
        integer :: info


        similarity_matrix = logzero

        invcovmat=cholesky
        call dpotri('L',settings%nDims,invcovmat,settings%nDims,info)
        do i=1,settings%nDims
            do j=i+1,settings%nDims
                invcovmat(i,j) = invcovmat(j,i)
            end do
        end do
        !covmat = matmul(transpose(cholesky),cholesky)
        !write(*,'(<settings%nDims>I3)') nint(matmul(covmat,invcovmat))


        ! Start by computing the similarity matrix
        do i=1,settings%nlive
            do j= i+1,settings%nlive
                displacement = live_points(settings%h0:settings%h1,i)-live_points(settings%h0:settings%h1,j)
                similarity_matrix(i,j) = 0.5d0 * log(dot_product(displacement,matmul(invcovmat,displacement)))
                similarity_matrix(j,i) = similarity_matrix(i,j)
            end do
        end do

        ! Compute our threshold
        logsum_d= log(sum(exp(similarity_matrix)))
        threshold = logsum_d-log(max_clusters+0d0)*2

        ! Assign them all to the initial cluster
        cluster = .true.

        ! Initially the sum_mij=0
        logprod_d = 1
        logsum_m = logzero

        keep_going=.true.

        do while(keep_going)

            keep_going=.false.

            ! Loop over the points in the cluster

            ! Generate a random ordering to go through
            random_deck = random_distinct_integers(settings%nlive,settings%nlive)


            do j=1,settings%nlive
                i = random_deck(j)
                numsub = count(cluster==cluster(i)) - 1
                numadd = count(cluster/=cluster(i)) + 1

                ! Find the contributions to subtract
                logprod_dsub = sum(similarity_matrix(:i-1,i),mask=cluster(:i-1)/=cluster(i)) + sum(similarity_matrix(i+1:,i),mask=cluster(i+1:)/=cluster(i)) 
                ! Find the contributions to add
                logprod_dadd = sum(similarity_matrix(:i-1,i),mask=cluster(:i-1)==cluster(i)) + sum(similarity_matrix(i+1:,i),mask=cluster(i+1:)==cluster(i)) 

                ! calculate the change in logprod_d
                logprod_d_new = logprod_d + logprod_dadd - logprod_dsub

                ! Calculate the new sum on the model values
                logsum_m_new = logprod_d_new/(numsub*numadd) + log(numsub+0d0) + log(numadd+0d0) + log(2d0)

                if(logsum_m_new>logsum_m) then
                    keep_going=.true.
                    cluster(i) = .not.cluster(i)

                    logsum_m = logsum_m_new
                    logprod_d = logprod_d_new
                end if

            end do
        end do

        if(.true.) then 
            !write(*,'("Cost: ", F10.4,"/",F10.4)') 2-1 + max_clusters**2  * (1 - exp(logsum_m-logsum_d) ), max_clusters**2
            write(*,'(F10.4)') exp(logsum_d-logsum_m) 
            if(2-1 + max_clusters**2  * (1 - exp(logsum_m-logsum_d) ) < max_clusters**2) write(*,*) 'cluster found'
            do i=1,settings%nlive
                if(cluster(i)) then
                    live_points(settings%cluster,i) = 1d0 
                else
                    live_points(settings%cluster,i) = 0d0 
                end if
            end do
        else
            live_points(settings%cluster,:) = 0
        end if
        !call sleep(1)


    end subroutine Skilling_clustering




    subroutine SNN_clustering(settings,live_points)
        use settings_module, only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nTotal,settings%nlive),intent(inout) :: live_points

        integer, dimension(settings%SNN_k+1) :: knn

        integer, dimension(settings%nlive) :: cluster_list
        integer, dimension(settings%nlive) :: cluster_list_final

        integer :: i_point,j_point

        integer, dimension(settings%SNN_k+1) :: cluster_indices
        integer :: max_cluster_index
        integer :: new_cluster_index

        integer :: clusters(100)
        integer :: num_clusters

        integer :: i_cluster

        ! Initialise the cluster list
        cluster_list = 0 

        max_cluster_index=0

        ! Loop through all pairs of points
        do i_point=1,settings%nlive

            ! Find the k nearest neighbors to this point
            knn = compute_knn(settings%SNN_k,i_point,live_points(settings%h0:settings%h1,:))

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

                do j_point=1,settings%nlive
                    if(any(cluster_list(j_point)==cluster_indices)) cluster_list(j_point)=new_cluster_index
                end do

        end do

        ! Find the cluster
        num_clusters=1
        clusters(1) = cluster_list(1)

        do i_point=1,settings%nlive
            if(all(cluster_list(i_point)/=clusters(1:num_clusters))) then
                num_clusters=num_clusters+1
                clusters(num_clusters) = cluster_list(i_point)
            end if
        end do

        do i_cluster=1,num_clusters
            where(cluster_list==clusters(i_cluster)) cluster_list_final=i_cluster
        end do

        write(*,'("ncluster  = ", I20                  )') num_clusters

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
