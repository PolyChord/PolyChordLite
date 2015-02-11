module nested_sampling_module

#ifdef MPI
    use mpi_module
#endif

    implicit none


    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood,priors,settings,mpi_communicator) result(output_info)
        use priors_module,     only: prior,prior_log_volume
        use utils_module,      only: stdout_unit
        use settings_module
        use utils_module,      only: logsumexp,calc_similarity_matrix,write_untxt_unit
        use read_write_module, only: write_resume_file,write_posterior_file,write_phys_live_points,read_resume_file,resume_file,write_stats_file,posterior_file
        use feedback_module
        use evidence_module,   only: run_time_info,allocate_run_time_info
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: NN_clustering
        use generate_module

        implicit none

        interface
            function loglikelihood(theta,phi,context)
                double precision, intent(in),  dimension(:) :: theta
                double precision, intent(out),  dimension(:) :: phi
                integer,          intent(in)                 :: context
                double precision :: loglikelihood
            end function
        end interface

        type(prior), dimension(:), intent(in) :: priors
        type(program_settings), intent(in) :: settings

        integer, intent(in) :: mpi_communicator

        ! Output of the program
        ! 1) log(evidence)
        ! 2) error(log(evidence))
        ! 3) ndead
        ! 4) number of likelihood calls
        ! 5) log(evidence) + log(prior volume)
        double precision, dimension(5) :: output_info

        ! The run time info (very important, see src/run_time_info.f90)
        type(run_time_info) :: RTI

        logical :: more_samples_needed

        ! The new-born baby points
        double precision, dimension(settings%nTotal,settings%num_babies) :: baby_points

        ! Point to seed a new one from
        double precision, dimension(settings%nTotal) :: seed_point

        ! Cholesky matrix to send
        double precision, dimension(settings%nDims,settings%nDims)   :: seed_point


        logical :: linear_mode ! Whether to run in linear mode if nprocs==1
        logical :: resume      ! Whether to resume from file


        integer :: nprocs
        integer :: myrank
        integer :: root

        integer :: i_cluster
        integer :: clusters(settings%nlive)
        integer :: num_new_clusters
        double precision, dimension(settings%nlive,settings%nlive) :: similarity_matrix

#ifdef MPI
        integer, dimension(MPI_STATUS_SIZE) :: mpi_status
        integer :: i_slave

        ! MPI initialisation
        nprocs = get_nprocs(mpi_communicator)         ! Get the number of MPI processes
        myrank = get_rank(mpi_communicator)           ! Get the identity of this process
        root   = assign_root(myrank,mpi_communicator) ! Assign the root process as the minimum integer

#else 

        ! non-MPI initialisation
        root=0      ! Define root node
        myrank=root ! set the only node's rank to be the root
        nprocs=1    ! there's only one process
#endif

        linear_mode = nprocs==1 ! Run in linear mode if only one process

        if(myrank==root) call write_opening_statement(settings) 




        !======= 1) Initialisation =====================================

        ! Check if we actually want to resume
        resume = settings%read_resume .and. resume_file_exists(settings)

        if (resume) then 

            ! Read the resume file on root
            if(myrank==root) call read_resume_file(settings,RTI) 

        else 
            
            ! Intialise the run by setting all of the relevant run time info, and generating live points
            call GenerateLivePoints(loglikelihood,priors,settings,RTI,mpi_communicator,nprocs,myrank,root)

            ! Write a resume file (as the generation of live points can be intensive)
            if(myrank==root.and.settings%write_resume) then
                call delete_files(settings)     ! Delete any existing files
                write_resume_file(settings,RTI) ! Write a new resume file
            end if

        end if !(.not.resume)


        if(myrank==root) then


            ! --- consider revising
            call calc_covmats(settings,RTI)

            !======= 2) Main loop body =====================================

            ! -------------------------------------------- !
            call write_started_sampling(settings%feedback)
            ! -------------------------------------------- !

            ! definitely more samples needed than this
            more_samples_needed = .true.

            ! Delete the first outer point
            !---consider revising
            call delete_outer_point(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior)


            do while ( more_samples_needed )

                ! (1) Generate a fresh incubation stack of baby_points
                if(settings%feedback>=2) write(stdout_unit,'(" Generating incubation stack ")')

                if(linear_mode) then

                    ! Generate a seed point --- update this
                    seed_point = GenerateSeed(settings,info,live_points,i_cluster)

                    ! Choose the cholesky decomposition for the cluster
                    cholesky = choleskys(:,:,i_cluster)

                    ! Generate a new set of points within the likelihood bound of the late point
                    baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)

#ifdef MPI
                else !(.not.linear_mode)

                    ! Recieve any new baby points from any slave currently sending
                    i_slave = catch_babies(baby_points,mpi_communicator)

                    seed_point = GenerateSeed(settings,info,live_points,i_cluster)  ! generate a new seed point --- update this
                    call throw_seed(seed_point,mpi_communicator,i_slave,.true.)     ! and throw back to slave (true => keep going)
                    
                    cholesky = RTI%cholesky(:,:,i_cluster)                 ! Find the cholesky decomp for seed
                    call throw_cholesky(cholesky,mpi_communicator,i_slave) ! Throw this back to slave

#endif
                end if !(linear_mode / .not.linear_mode)

                
                !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

                ! (2) Add the babies to the array, testing to see if these
                ! constitute a valid point. 
                if( add_babies(settings,info,live_points,phantom_points,nphantom,baby_points) ) then

                    ! Record that we have a new dead point
                    ndead = ndead + 1

                    ! (5) Update the covariance matrix of the distribution of live points
                    if(mod(ndead,settings%nlive) .eq.0) then

                        if(settings%do_clustering) then 

                            if(settings%feedback>=2) write(stdout_unit,'(" Doing Clustering ")')
                            i_cluster=1
                            do while(i_cluster<=info%ncluster_A)
                                ! For each active cluster, see if it is further sub-clustered

                                ! Calculate a similarity matrix
                                similarity_matrix(:info%n(i_cluster),:info%n(i_cluster)) &
                                        = calc_similarity_matrix(live_points(settings%h0:settings%h1,:info%n(i_cluster),i_cluster))

                                ! Do clustering on this matrix
                                num_new_clusters = NN_clustering( &
                                        similarity_matrix(:info%n(i_cluster),:info%n(i_cluster)), &
                                        settings%SNN_k,clusters(:info%n(i_cluster)))

                                ! If we've found new clusters, then we should bifurcate the algorithm at this point
                                if(num_new_clusters>1) then

                                    if( num_new_clusters+info%ncluster_A>settings%ncluster ) then
                                        call halt_program(" Too many clusters. Consider increasing settings%ncluster")
                                    else if (num_new_clusters + info%ncluster_T > settings%ncluster*2 ) then
                                        call halt_program(" Too many clusters. Consider increasing settings%nclustertot")
                                    else
                                        call create_new_clusters(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior,i_cluster,clusters(:info%n(i_cluster)),num_new_clusters)

                                        write(stdout_unit,'( I8, " clusters found at iteration ", I8)') info%ncluster_A, ndead
                                    end if
                                else
                                    ! Otherwise move on to the next cluster
                                    i_cluster=i_cluster+1
                                end if


                            end do

                        end if
                        ! Calculate the covariance matrices
                        covmats = calc_covmats(settings,info,live_points,phantom_points,nphantom)
                        ! Calculate the cholesky decomposition
                        choleskys = calc_choleskys(covmats)

                    end if






                    ! (3) Feedback to command line every nlive iterations
                    ! Test to see if we need to finish
                    more_samples_needed =  (live_logZ(settings,info,live_points) > log(settings%precision_criterion) + info%logevidence ) 

                    ! (4) Update the resume and posterior files every update_resume iterations, or at program termination
                    if ( (mod(ndead,settings%update_resume) == 0) .or.  (more_samples_needed.eqv..false.) )  then


                        ! ---------------------------------------------------------------------- !
                        call write_intermediate_results(settings,info,ndead,nphantom,nposterior,&
                                mean_likelihood_calls(settings,info,live_points) ) 
                        ! ---------------------------------------------------------------------- !

                        if(settings%feedback>=2) write(stdout_unit,'(" Writing resume files ")')
                        if(settings%calculate_posterior) call write_posterior_file(settings,info,posterior_points,nposterior)  
                        if(settings%write_resume)        call write_resume_file(settings,info,live_points,nphantom,phantom_points,&
                                ndead,total_likelihood_calls)
                        if(settings%write_live)          call write_phys_live_points(settings,info,live_points)
                        call write_stats_file(settings,info,ndead) 

                    end if

                    ! If we've put a limit on the maximum number of iterations, then
                    ! check to see if we've reached this
                    if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.



                    ! (6] delete the next outer point.
                    if(settings%feedback>=2) write(stdout_unit,'(" Deleting outer point ")')
                    call delete_outer_point(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior)


                end if

            end do ! End main loop



#ifdef MPI
            if(.not.linear_mode) then

                ! Kill off the final slaves
                ! If we're done, then clean up by receiving the last piece of
                ! data from each node (and throw it away) and then send a kill signal back to it
                do active_slaves=nprocs-1,1,-1

                    ! Recieve baby point from slave i_slave
                    i_slave = catch_babies(baby_points,mpi_communicator)

                    ! Send kill signal to slave i_slave (note that we no longer care about seed_point, so we'll just use the last one
                    call throw_seed(seed_point,mpi_communicator,i_slave,.false.) 

                end do

            end if !(.not.linear_mode / linear_mode )
#endif


            ! Create the output array
            ! (1) log evidence
            ! (2) Error in the log evidence
            ! (3) Number of dead points
            ! (4) Number of likelihood calls
            ! (5) log(evidence * prior volume)
            output_info(1) = 2*info%logevidence - 0.5*info%logevidence2
            output_info(2) = sqrt(info%logevidence2 - 2*info%logevidence)
            output_info(3) = ndead
            output_info(4) = total_likelihood_calls
            output_info(5) = output_info(1)+prior_log_volume(priors)

            ! ------------------------------------------------------------ !
            call write_final_results(output_info,settings%feedback,priors)
            ! ------------------------------------------------------------ !





#ifdef MPI
        else !(myrank/=root)

            ! These are the slave tasks
            ! -------------------------
            !
            ! This is considerably simpler than that of the master.
            ! All slaves do is:
            ! 1) recieve a seed point from the master
            ! 2) recieve a cholesky decomposition from the master
            ! 3) using the seed and cholesky, generate a new set of baby points
            ! 4) Send the baby points back to the master


            ! BEGIN) On the first loop, send a nonsense set of baby_points
            ! to indicate that we're ready to start receiving

            baby_points = 0d0
            baby_points(settings%l0) = logzero

            call throw_babies(baby_points,mpi_communicator,root)

            ! 1) Listen for a seed point being sent by the master
            !    Note that this also tests for a kill signal sent by the master
            do while(catch_seed(seed_point,mpi_communicator,root))


                ! 2) Listen for the cholesky decomposition sent by the master
                call catch_cholesky(cholesky,mpi_communicator,root) 

                ! 3) Generate a new set of baby points
                baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,seed_point)

                ! 4) Send the baby points back
                call throw_babies(baby_points,mpi_communicator,root)

            end do

#endif
        end if !(myrank==root / myrank/=root) 




    end function NestedSampling







    subroutine create_new_clusters(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior,i_cluster,clusters,num_new_clusters)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info,bifurcate_evidence
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: info

        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(inout)   :: live_points
        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster),intent(inout)  :: phantom_points
        integer, dimension(settings%ncluster),intent(inout) :: nphantom
        double precision, dimension(settings%nposterior,settings%nstack,0:settings%ncluster*2),intent(inout) :: posterior_points
        integer, dimension(0:settings%ncluster*2),intent(inout) :: nposterior

        integer, intent(in)                             :: i_cluster
        integer,dimension(info%n(i_cluster)),intent(in) :: clusters
        integer, intent(in)                             :: num_new_clusters
        

        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster) :: temp_phantoms
        integer, dimension(settings%ncluster) :: temp_nphantom

        ! Temporary live points
        double precision, dimension(settings%nTotal,info%n(i_cluster)) :: lives
        integer :: nlives
        integer :: i_live

        integer :: i_phantom

        integer :: j_cluster,k_cluster

        integer, dimension(num_new_clusters) :: ni
        integer :: ncluster_A_old

        ! Get the lives to split
        nlives = info%n(i_cluster)
        lives(:,:)  = live_points(:,:info%n(i_cluster),i_cluster)

        ! Create a temporary set of phantom points for re-assigning to clusters
        ncluster_A_old= info%ncluster_A
        temp_phantoms = phantom_points
        temp_nphantom = nphantom


        ! Set the number of points in the new clusters to 0 initially
        ni = 0

        ! Split the live points stored in lives into their new clusters
        do i_live=1,nlives

            ! Iterate the number of points in cluster 'clusters(i_live)'
            ni(clusters(i_live))= ni(clusters(i_live)) + 1

            ! Add the 'i_live' point to the 'clusters(i_live)' cluster from the end
            live_points(:,ni(clusters(i_live)),clusters(i_live) + info%ncluster_A) = lives(:,i_live)
        end do

        ! Move up the posterior points
        posterior_points(:,:,info%ncluster_A+1+size(ni):info%ncluster_T+size(ni)) = posterior_points(:,:,info%ncluster_A+1:info%ncluster_T)
        nposterior(info%ncluster_A+1+size(ni):info%ncluster_T+size(ni)) = nposterior(info%ncluster_A+1:info%ncluster_T)
        nposterior(info%ncluster_A+1:info%ncluster_A+size(ni))=0

        ! Update the info variable for these new clusters
        call bifurcate_evidence(info,i_cluster,ni)

        ! Delete the old cluster
        call reorganise_clusters(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior,i_cluster,.false.) 

        ! Reassign all of the phantom points
        ! set the number of phantoms to zero
        nphantom=0

        ! Reassign the phantom points
        do j_cluster=1,ncluster_A_old
            do i_phantom=1,temp_nphantom(j_cluster)
                ! Identify the cluster this phantom point belongs to
                k_cluster = identify_cluster(settings,info,temp_phantoms(:,i_phantom,j_cluster),live_points)
                nphantom(k_cluster) = nphantom(k_cluster) + 1

                phantom_points(:,nphantom(k_cluster),k_cluster) = temp_phantoms(:,i_phantom,j_cluster)
            end do
        end do


    end subroutine create_new_clusters



    function identify_cluster(settings,info,point,live_points) result(cluster)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use utils_module,      only: loginf,distance2
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: info

        double precision, dimension(settings%nTotal),intent(in)   :: point
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(in)   :: live_points

        integer :: cluster

        integer :: i_cluster
        integer :: i_live

        double precision :: temp_distance2
        double precision :: closest_distance2

        closest_distance2=loginf

        ! Find the cluster this point is nearest to
        do i_cluster=1,info%ncluster_A
            do i_live=1,info%n(i_cluster)
                temp_distance2 = distance2(point(settings%h0:settings%h1),live_points(settings%h0:settings%h1,i_live,i_cluster) )
                if(temp_distance2 < closest_distance2) then
                    cluster = i_cluster
                    closest_distance2 = temp_distance2
                end if
            end do
        end do

    end function identify_cluster














    function add_babies(settings,info,live_points,phantom_points,nphantom,baby_points) result(babies_added)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: info

        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(inout)   :: live_points
        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster),intent(inout)  :: phantom_points
        integer, dimension(settings%ncluster),intent(inout) :: nphantom
        double precision, dimension(settings%nTotal,settings%num_babies),intent(in)                   :: baby_points

        logical :: babies_added

        integer,dimension(settings%num_babies) :: baby_cluster
        integer :: i_cluster
        integer :: i_baby

        ! Identify the clusters
        if(settings%do_clustering .and. info%ncluster_A>1 ) then
            do i_baby=1,settings%num_babies
                baby_cluster(i_baby) = identify_cluster(settings,info,baby_points(:,i_baby),live_points)
            end do
        else
            ! No clustering needed
            baby_cluster = 1
        end if

        ! If the baby has a likelihood >= contour, then add it to the end of the
        ! live points (this is where the lowest likelihood point is found).
        i_baby = settings%num_babies
        i_cluster = baby_cluster(i_baby)
        if( baby_points(settings%l0,i_baby) >= info%logL(i_cluster) ) then

            ! record the increase in the number of live points for this cluster
            info%n(i_cluster) = info%n(i_cluster) + 1

            ! Add it to the end of the live points in cluster 'i_cluster'
            live_points(:,info%n(i_cluster),i_cluster) = baby_points(:,settings%num_babies)

            ! Record that we've added a baby point
            babies_added=.true.
        else
            babies_added=.false.
        end if

        ! Now add all of the rest to the end of the phantom array (these will be pruned later if necessary)
        do i_baby=1,settings%num_babies-1
            ! If the baby is within the likelihood contour
            if( baby_points(settings%l0,i_baby) >= info%logL(i_cluster) ) then
                ! Note the cluster number
                i_cluster = baby_cluster(i_baby)

                ! record the increase in number of phantom points in that cluster
                nphantom(i_cluster) = nphantom(i_cluster) + 1

                ! Add the baby point to the correct part of the phantom array
                phantom_points(:,nphantom(i_cluster),i_cluster) = baby_points(:,i_baby)
            end if
        end do

    end function add_babies



    subroutine delete_outer_point(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior) 
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info,update_evidence
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: info

        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(inout)   :: live_points
        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster),intent(inout)  :: phantom_points
        integer, dimension(settings%ncluster),intent(inout) :: nphantom

        double precision, dimension(settings%nposterior,settings%nstack,0:settings%ncluster*2),intent(inout) :: posterior_points
        integer, dimension(0:settings%ncluster*2),intent(inout) :: nposterior

        double precision :: min_loglike
        integer :: min_cluster

        double precision, dimension(settings%nTotal) :: dead_point

        double precision :: logweight


        ! Find the outermost point, and rearrange the live points to 'delete'
        ! the outermost point
        dead_point = find_outer_point(settings,info,live_points,min_loglike,min_cluster)

        ! Get the logweighting for use in posterior calculation
        logweight = info%logX(min_cluster) - log( info%n(min_cluster) + 1d0 )

        ! Update the evidence and other stored information 
        call update_evidence(info,min_cluster,min_loglike)

        evidence = info%logevidence
        ! Update the posterior and phantom_arrays
        call update_posterior_and_phantom(settings,info,posterior_points,nposterior,dead_point,phantom_points,nphantom,min_cluster,min_loglike,logweight) 
        ! Decrease the chord estimate for this cluster
        !live_points(settings%last_chord,:,min_cluster) = live_points(settings%last_chord,:,min_cluster) *info%n(min_cluster) / (info%n(min_cluster) + 1d0) 

        ! If we've deleted a cluster, we should re-organise live and phantom points
        if(info%n(min_cluster)==0) then
            write(*,'(" Deleting cluster", I4)') min_cluster
            call reorganise_clusters(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior,min_cluster,.true.)
        end if

    end subroutine delete_outer_point


    subroutine reorganise_clusters(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior,min_cluster,deleted)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info,delete_evidence
        use read_write_module, only: posterior_file
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: info

        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(inout)   :: live_points
        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster),intent(inout)  :: phantom_points
        integer, dimension(settings%ncluster),intent(inout) :: nphantom

        double precision, dimension(settings%nposterior,settings%nstack,0:settings%ncluster*2),intent(inout) :: posterior_points
        integer, dimension(0:settings%ncluster*2),intent(inout) :: nposterior

        integer,intent(in) :: min_cluster
        logical,intent(in) :: deleted

        integer :: top_bound

        integer :: i_cluster
        
        top_bound = info%ncluster_A

        ! Cyclically shift the active points down one
        live_points(     :,:,min_cluster:top_bound) = cshift(live_points(     :,:,min_cluster:top_bound),shift=1,dim=3)

        phantom_points(  :,:,min_cluster:top_bound) = cshift(phantom_points(  :,:,min_cluster:top_bound),shift=1,dim=3)
        nphantom(            min_cluster:top_bound) = cshift(nphantom(            min_cluster:top_bound),shift=1,dim=1) 

        if(.not. deleted) top_bound = top_bound + info%ncluster_P

        posterior_points(:,:,min_cluster:top_bound) = cshift(posterior_points(:,:,min_cluster:top_bound),shift=1,dim=3)
        nposterior(          min_cluster:top_bound) = cshift(nposterior(          min_cluster:top_bound),shift=1,dim=1) 

        ! Rename the posterior files
        call rename(trim(posterior_file(settings,.true.,min_cluster)), 'temp_cluster.dat')
        do i_cluster=min_cluster,top_bound-1
            call rename(trim(posterior_file(settings,.true.,i_cluster+1)),trim(posterior_file(settings,.true.,i_cluster)))
        end do
        call rename( 'temp_cluster.dat',trim(posterior_file(settings,.true.,top_bound)))

        ! Update the evidence
        call delete_evidence(info,min_cluster,deleted)

    end subroutine reorganise_clusters





    function find_outer_point(settings,info,live_points,min_loglike,min_cluster) result(dead_point)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use utils_module,      only: loginf
        implicit none

        ! Inputs
        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(in) :: info
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(inout)   :: live_points
        double precision, dimension(settings%nTotal) :: dead_point

        double precision, intent(out)   :: min_loglike
        integer, intent(out)            :: min_cluster

        ! Local variables
        integer :: i_cluster
        double precision :: loglike
        integer          :: pos(1)
        integer          :: min_pos


        ! Initialise min_loglike at the largest possible value
        min_loglike = loginf

        ! Iterate through the clusters to find the cluster with the lowest
        ! outermost loglikelihood
        do i_cluster=1,info%ncluster_A

            ! Find the position of the lowest loglikelihood point in cluster 'i_cluster'
            pos = minloc(live_points(settings%l0,:info%n(i_cluster),i_cluster))
            ! Find that lowest loglikelihood
            loglike = live_points(settings%l0,pos(1),i_cluster)

            ! If this loglikelihood is lower than the current lowest, then update the outputs
            if(loglike < min_loglike) then
                min_loglike = loglike
                min_cluster = i_cluster
                min_pos     = pos(1)
            end if

        end do

        ! Record the point that has just died
        dead_point = live_points(:,min_pos,min_cluster)

        ! move the end point of min_cluster to replace the outermost point
        live_points(:,min_pos,min_cluster) = live_points(:,info%n(min_cluster),min_cluster)

    end function find_outer_point


    !> This function runs through the phantom array, removing any points from
    !! the relevant cluster that are now below the loglikelihood contour. It
    !! adds a fraction of the discarded phantoms to the posterior array.
    subroutine update_posterior_and_phantom(settings,info,posterior_points,nposterior,dead_point,phantom_points,nphantom,min_cluster,min_loglike,logweight)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use utils_module,      only: logincexp,logsumexp
        use random_module,     only: random_real,random_logicals
        implicit none

        ! Inputs
        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(in) :: info
        double precision, dimension(settings%nTotal),intent(in) :: dead_point
        integer, intent(in)            :: min_cluster
        double precision, intent(in)            :: min_loglike
        double precision, intent(in)            :: logweight

        ! Outputs
        double precision, dimension(settings%nposterior,settings%nstack,0:settings%ncluster*2),intent(inout) :: posterior_points
        integer, dimension(0:settings%ncluster*2),intent(inout) :: nposterior
        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster),intent(inout)  :: phantom_points
        integer, dimension(settings%ncluster),intent(inout) :: nphantom


        ! Local variables
        double precision, dimension(settings%nposterior) :: posterior_point

        integer :: i_phantom

        integer :: num_new
        double precision :: lognum_new
        integer :: i_new

        logical,dimension(settings%nstack) :: keeps






        if(settings%calculate_posterior) then
            ! Count the number of new posterior points
            num_new = count( phantom_points(settings%l0,:nphantom(min_cluster),min_cluster) < min_loglike) + 1
            i_new=1
            keeps = (/.true. , random_logicals(num_new-1,settings%thin_posterior) /)
            num_new = count(keeps(:num_new))
            lognum_new =log( num_new  + 0d0 )


            nposterior(0) = nposterior(0) + 1
            ! Now update the posterior information for the dead point
            !   - calculate the new posterior point
            posterior_point = calc_posterior_point(settings,dead_point,logweight-lognum_new,info%logevidence)
            !   - add this point to the end of the arrays
            posterior_points(:,nposterior(0)          ,0          ) = posterior_point

            if(settings%do_clustering) then
                nposterior(min_cluster) = nposterior(min_cluster) + 1
                posterior_point = calc_posterior_point(settings,dead_point,logweight-lognum_new,info%logZ(min_cluster))
                posterior_points(:,nposterior(min_cluster),min_cluster) = posterior_point
            end if


        end if





        ! Now run through the stack and strip out any points that are less
        ! than the min_loglike, replacing them with points drawn from
        ! the end 
        i_phantom=1
        do while(i_phantom<=nphantom(min_cluster))
            if( phantom_points(settings%l0,i_phantom,min_cluster) < min_loglike ) then

                if(settings%calculate_posterior) then
                    i_new=i_new+1
                    if(keeps(i_new)) then
                        nposterior(0) = nposterior(0) + 1
                        ! Now update the posterior information
                        !   - calculate the new posterior point
                        posterior_point = calc_posterior_point(settings,phantom_points(:,i_phantom,min_cluster),logweight-lognum_new,info%logevidence)
                        !   - add this point to the end of the arrays
                        posterior_points(:,nposterior(0)          ,0          ) = posterior_point

                        if(settings%do_clustering) then
                            nposterior(min_cluster) = nposterior(min_cluster) + 1
                            posterior_point = calc_posterior_point(settings,phantom_points(:,i_phantom,min_cluster),logweight-lognum_new,info%logZ(min_cluster))
                            posterior_points(:,nposterior(min_cluster),min_cluster) = posterior_point
                        end if


                    end if
                end if

                ! Overwrite the discarded point with a point from the end...
                phantom_points(:,i_phantom,min_cluster) = phantom_points(:,nphantom(min_cluster),min_cluster)
                ! ...and reduce the number of phantom points
                nphantom(min_cluster)=nphantom(min_cluster)-1
            else
                i_phantom=i_phantom+1
            end if
        end do


    end subroutine update_posterior_and_phantom



    function insert_point_at_min(vector,array,column,n) result(discarded_value)
        use utils_module, only: logzero
        implicit none
        double precision, intent(in),    dimension(:)   :: vector
        double precision, intent(inout), dimension(:,:) :: array
        integer, intent(in)    :: column
        integer, intent(inout) :: n

        double precision :: discarded_value

        integer :: row(1)

        discarded_value = logzero

        if(n<size(array,2)) then
            ! If there's space left in the array, then tag this point on the end
            n=n+1
            array(:,n) = vector
        else
            ! Otherwise find the minimum position
            row = minloc(array(column,:),mask=array(column,:)<vector(column))

            if(row(1)/=0) then
                ! Save the value in this discarded point
                discarded_value = array(column,row(1))
                ! Overwrite the array
                array(:,row(1)) = vector
            else
                discarded_value = vector(column)
            end if

        end if

    end function insert_point_at_min



    !> Calculate a posterior point from a live/phantom point, suitable for
    !! adding to a .txt file
    function calc_posterior_point(settings,point,logweight,evidence) result(posterior_point)
        use settings_module,   only: program_settings
        use utils_module, only: logincexp
        implicit none

        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nTotal),intent(in) :: point
        double precision,intent(in) :: logweight
        double precision,intent(in) :: evidence
        double precision, dimension(settings%nposterior) :: posterior_point


        ! Un-normalised weighting (needs to be unnormalised since the evidence is only correct at the end)
        posterior_point(settings%pos_w)  = point(settings%l0) + logweight
        ! un-normalise cumulative weighting
        posterior_point(settings%pos_Z)  = evidence
        ! Likelihood
        posterior_point(settings%pos_l)  = point(settings%l0)
        ! Physical parameters
        posterior_point(settings%pos_p0:settings%pos_p1) = point(settings%p0:settings%p1)
        ! Derived parameters
        posterior_point(settings%pos_d0:settings%pos_d1) = point(settings%d0:settings%d1)

    end function calc_posterior_point














    function calc_covmats(settings,info,live_points,phantom_points,nphantom) result(covmats)
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use utils_module,      only: calc_covmat
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),intent(in) :: info
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(in) :: live_points
        double precision, dimension(settings%nTotal,settings%nstack,settings%ncluster),intent(in) :: phantom_points
        integer, dimension(settings%ncluster),intent(in) :: nphantom

        double precision, dimension(settings%nDims,settings%nDims,settings%ncluster) :: covmats

        integer :: i,j

        covmats = 0d0
        do i=1,settings%nDims
            covmats(i,i,:) = 1d0
        end do

        do i=1,info%ncluster_A
            if(info%n(i)+nphantom(i) >= settings%nDims*(settings%nDims+1)/2) then
                covmats(:,:,i) = calc_covmat(&
                    live_points(settings%h0:settings%h1,:info%n(i),i),&
                    phantom_points(settings%h0:settings%h1,:nphantom(i),i)&
                    )
            else
                do j=1,settings%nDims
                    covmats(j,j,i) = 1d0
                end do
            end if
        end do

    end function calc_covmats

    function calc_choleskys(covmats) result(choleskys)
        use utils_module, only: calc_cholesky
        implicit none
        double precision, intent(in), dimension(:,:,:) :: covmats
        double precision, dimension(size(covmats,1),size(covmats,2),size(covmats,3)) :: choleskys

        integer :: i

        choleskys=0d0

        do i=1,size(covmats,3)
            choleskys(:,:,i) = calc_cholesky( covmats(:,:,i) )
        end do

    end function calc_choleskys


    function mean_likelihood_calls(settings,info,live_points) 
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),intent(in) :: info
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(in) :: live_points
        
        double precision mean_likelihood_calls
        integer i_cluster

        mean_likelihood_calls=0d0

        do i_cluster=1,info%ncluster_A
            mean_likelihood_calls = mean_likelihood_calls + sum(live_points(settings%nlike,:info%n(i_cluster),i_cluster))
        end do

        mean_likelihood_calls = mean_likelihood_calls/(settings%nlive + 0d0)


    end function mean_likelihood_calls 

    function live_logZ(settings,info,live_points) 
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use utils_module,      only: logzero,logsumexp,logincexp
        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),intent(in) :: info
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(in) :: live_points
        
        double precision live_logZ
        integer i_cluster

        ! Sum up the live log evidence
        live_logZ=logzero

        do i_cluster=1,info%ncluster_A

            call logincexp( live_logZ, logsumexp(live_points(settings%l0,:info%n(i_cluster),i_cluster)) &
                                       - log( info%n(i_cluster) + 0d0) &
                                       + info%logX(i_cluster)   )


        end do


    end function live_logZ




    function mean_w(settings,info,live_points) 
        use settings_module,   only: program_settings
        use evidence_module,   only: run_time_info
        use utils_module,      only: loginf,logzero
        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),intent(in) :: info
        double precision, dimension(settings%nTotal,settings%nlive,settings%ncluster),intent(in) :: live_points
        
        double precision mean_w(3)
        integer i_cluster

        mean_w(1)=0
        mean_w(2)=logzero
        mean_w(3)=loginf

        do i_cluster=1,info%ncluster_A
            mean_w(1) = mean_w(1) + sum(live_points(settings%last_chord,:info%n(i_cluster),i_cluster))
            mean_w(2) = max(mean_w(2),maxval(live_points(settings%last_chord,:info%n(i_cluster),i_cluster)))
            mean_w(3) = min(mean_w(3),minval(live_points(settings%last_chord,:info%n(i_cluster),i_cluster)))
        end do



    end function mean_w 





end module nested_sampling_module
