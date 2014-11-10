module cluster_module

    implicit none
    contains


    function SNN_clustering(settings,live_points) result(cluster_list)
        use settings_module, only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nTotal,settings%nlive),intent(in) :: live_points

        integer, dimension(settings%nlive) :: cluster_list

        integer, dimension(settings%SNN_k,settings%nlive) :: knn

        integer :: i_point,j_point

        integer :: cluster_orig_i,cluster_orig_j

        integer :: clusters(10)

        double precision, dimension(settings%nlive,settings%nlive) :: similarity_matrix


        ! Compute the similarity matrix
        similarity_matrix = spread( [( dot_product(live_points(settings%h0:settings%h1,i_point),live_points(settings%h0:settings%h1,i_point)),i_point=1,settings%nlive )], dim=2,ncopies=settings%nlive )
        similarity_matrix = similarity_matrix + transpose(similarity_matrix) - 2d0 * matmul( transpose(live_points(settings%h0:settings%h1,:)),live_points(settings%h0:settings%h1,:) )

        ! computing the k nearest neighbors for each point
        knn = compute_knn(similarity_matrix,settings%SNN_k)

        ! Set up the cluster list
        cluster_list = [( i_point,i_point=1,settings%nlive )]

        ! Loop through all pairs of points
        ! But skip points that are already in the same cluster

        do i_point=1,settings%nlive
            do j_point=i_point+1,settings%nlive

                cluster_orig_i = cluster_list(i_point)
                cluster_orig_j = cluster_list(j_point)

                ! If they're not in the same cluster already...
                if(cluster_orig_i/=cluster_orig_j) then
                    ! ... check to see if they are each others nearest neihbors...
                    if( neighbors( knn(:,i_point),knn(:,j_point) ) ) then
                        ! ... and check to see if they share n nearest neighbors
                        if( matches( knn(:,i_point),knn(:,j_point) ) > settings%SNN_kt ) then

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

    end function SNN_clustering




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
        integer :: same_list_arr(1)

        ! Check to see if they're in each others neighbors lists
        same_list=.false.

        same_list_arr = minloc(knn1,mask=knn1==knn2(1))
        if(same_list_arr(1)==0) return

        same_list_arr = minloc(knn2,mask=knn2==knn1(1))
        if(same_list_arr(1)==0) return

        same_list=.true.

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
