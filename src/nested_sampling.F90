module nested_sampling_module

#ifdef MPI
    use mpi_module
#endif

    implicit none


    contains

    !> Main subroutine for computing a generic nested sampling algorithm
    function NestedSampling(loglikelihood,priors,settings,mpi_communicator) result(output_info)
        use priors_module,     only: prior,prior_log_volume
        use utils_module,      only: stdout_unit,swap_int
        use settings_module,   only: program_settings
        use utils_module,      only: logsumexp,calc_similarity_matrix,write_untxt_unit
        use read_write_module, only: write_resume_file
        use feedback_module
        use run_time_module,   only: run_time_info
        use chordal_module,    only: SliceSampling
        use random_module,     only: random_integer,random_direction
        use cluster_module,    only: NN_clustering
        use generate_module,   only: GenerateSeed,GenerateLivePoints

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
        double precision, dimension(settings%nDims,settings%nDims)   :: cholesky


        logical :: linear_mode ! Whether to run in linear mode if nprocs==1
        logical :: resume      ! Whether to resume from file

        integer :: nlike ! Temporary storage for number of likelihood calls


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
        integer, dimension(:), allocatable :: slave_cluster       !> The cluster the slave is currently working on

        ! MPI initialisation
        nprocs = get_nprocs(mpi_communicator)         ! Get the number of MPI processes
        myrank = get_rank(mpi_communicator)           ! Get the identity of this process
        root   = assign_root(myrank,mpi_communicator) ! Assign the root process as the minimum integer

        allocate(slave_cluster(nprocs-1)) ! Allocate the slave arrays

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














        !======= 2) Main loop body =====================================

        if(myrank==root) then


            ! -------------------------------------------- !
            call write_started_sampling(settings%feedback)
            ! -------------------------------------------- !

            ! definitely more samples needed than this
            more_samples_needed = .true.


            do while ( more_samples_needed )

                ! Generate a seed point --- update this
                seed_point = GenerateSeed(settings,RTI,i_cluster)

                ! Choose the cholesky decomposition for the cluster
                cholesky = RTI%cholesky(:,:,i_cluster)

                ! Get the loglikelihood contour we're generating from
                logL = RTI%logLp(i_cluster)


                if(linear_mode) then
                    ! Generate a new set of points within the likelihood bound of the late point
                    baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,seed_point,logL,nlike)
#ifdef MPI
                else !(.not.linear_mode)

                    ! Recieve any new baby points from any slave currently sending
                    i_slave = catch_babies(baby_points,nlike,mpi_communicator)

                    ! and throw seeding information back to slave (true => keep going)
                    call throw_seed(seed_point,cholesky,logL,mpi_communicator,i_slave,.true.)

                    ! set logL to be the bound of the babies just recieved (saved in slave_logL)
                    ! and set slave_logL to be the bound just sent off
                    call swap_int(i_cluster,slave_cluster(i_slave))

#endif
                end if !(linear_mode / .not.linear_mode)

                ! Add the likelihood calls to our counter
                RTI%nlike = RTI%nlike + nlike

                

                if( replace_point(settings,RTI,baby_points,i_cluster) ) then

                    call check_end

                    call do_clustering

                    call write_info

                end if
            end do
                !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv







!                ! (2) Add the babies to the array, testing to see if these
!                ! constitute a valid point. 
!                if( add_babies(settings,info,live_points,phantom_points,nphantom,baby_points) ) then
!
!                    ! Record that we have a new dead point
!                    ndead = ndead + 1
!
!                    ! (5) Update the covariance matrix of the distribution of live points
!                    if(mod(ndead,settings%nlive) .eq.0) then
!
!                        if(settings%do_clustering) then 
!
!                            if(settings%feedback>=2) write(stdout_unit,'(" Doing Clustering ")')
!                            i_cluster=1
!                            do while(i_cluster<=info%ncluster_A)
!                                ! For each active cluster, see if it is further sub-clustered
!
!                                ! Calculate a similarity matrix
!                                similarity_matrix(:info%n(i_cluster),:info%n(i_cluster)) &
!                                        = calc_similarity_matrix(live_points(settings%h0:settings%h1,:info%n(i_cluster),i_cluster))
!
!                                ! Do clustering on this matrix
!                                num_new_clusters = NN_clustering( &
!                                        similarity_matrix(:info%n(i_cluster),:info%n(i_cluster)), &
!                                        settings%SNN_k,clusters(:info%n(i_cluster)))
!
!                                ! If we've found new clusters, then we should bifurcate the algorithm at this point
!                                if(num_new_clusters>1) then
!
!                                    if( num_new_clusters+info%ncluster_A>settings%ncluster ) then
!                                        call halt_program(" Too many clusters. Consider increasing settings%ncluster")
!                                    else if (num_new_clusters + info%ncluster_T > settings%ncluster*2 ) then
!                                        call halt_program(" Too many clusters. Consider increasing settings%nclustertot")
!                                    else
!                                        call create_new_clusters(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior,i_cluster,clusters(:info%n(i_cluster)),num_new_clusters)
!
!                                        write(stdout_unit,'( I8, " clusters found at iteration ", I8)') info%ncluster_A, ndead
!                                    end if
!                                else
!                                    ! Otherwise move on to the next cluster
!                                    i_cluster=i_cluster+1
!                                end if
!
!
!                            end do
!
!                        end if
!                        ! Calculate the covariance matrices
!                        covmats = calc_covmats(settings,info,live_points,phantom_points,nphantom)
!                        ! Calculate the cholesky decomposition
!                        choleskys = calc_choleskys(covmats)
!
!                    end if
!
!
!
!
!
!
!                    ! (3) Feedback to command line every nlive iterations
!                    ! Test to see if we need to finish
!                    more_samples_needed =  (live_logZ(settings,info,live_points) > log(settings%precision_criterion) + info%logevidence ) 
!
!                    ! (4) Update the resume and posterior files every update_resume iterations, or at program termination
!                    if ( (mod(ndead,settings%update_resume) == 0) .or.  (more_samples_needed.eqv..false.) )  then
!
!
!                        ! ---------------------------------------------------------------------- !
!                        call write_intermediate_results(settings,info,ndead,nphantom,nposterior,&
!                                mean_likelihood_calls(settings,info,live_points) ) 
!                        ! ---------------------------------------------------------------------- !
!
!                        if(settings%feedback>=2) write(stdout_unit,'(" Writing resume files ")')
!                        if(settings%calculate_posterior) call write_posterior_file(settings,info,posterior_points,nposterior)  
!                        if(settings%write_resume)        call write_resume_file(settings,info,live_points,nphantom,phantom_points,&
!                                ndead,total_likelihood_calls)
!                        if(settings%write_live)          call write_phys_live_points(settings,info,live_points)
!                        call write_stats_file(settings,info,ndead) 
!
!                    end if
!
!                    ! If we've put a limit on the maximum number of iterations, then
!                    ! check to see if we've reached this
!                    if (settings%max_ndead >0 .and. ndead .ge. settings%max_ndead) more_samples_needed = .false.
!
!
!
!                    ! (6] delete the next outer point.
!                    if(settings%feedback>=2) write(stdout_unit,'(" Deleting outer point ")')
!                    call delete_outer_point(settings,info,live_points,phantom_points,nphantom,posterior_points,nposterior)
!
!
!                end if
!
!            end do ! End main loop



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
            call calculate_logZ_estimate(RTI,output_info(1),output_info(2))
            output_info(3) = RTI%ndead
            output_info(4) = RTI%nlike
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
            nlike = 0

            call throw_babies(baby_points,nlike,mpi_communicator,root)

            ! 1) Listen for a seed point being sent by the master
            !    Note that this also tests for a kill signal sent by the master
            do while(catch_seed(seed_point,cholesky,logL,mpi_communicator,root))

                ! 2) Generate a new set of baby points
                baby_points = SliceSampling(loglikelihood,priors,settings,cholesky,logL,seed_point,nlike)

                ! 3) Send the baby points back
                call throw_babies(baby_points,nlike,mpi_communicator,root)

            end do

#endif
        end if !(myrank==root / myrank/=root) 




    end function NestedSampling







    function identify_cluster(settings,RTI,point) result(cluster)
        use settings_module,   only: program_settings
        use run_time_module,   only: run_time_info
        use utils_module,      only: loginf,distance2
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info), intent(inout) :: RTI

        double precision, dimension(settings%nTotal),intent(in)   :: point

        integer :: cluster

        integer :: i_cluster
        integer :: i_live

        double precision :: temp_distance2
        double precision :: closest_distance2

        closest_distance2=loginf

        ! Find the cluster this point is nearest to
        do i_cluster=1,RTI%ncluster
            do i_live=1,RTI%nlive(i_cluster)
                temp_distance2 = distance2(point(settings%h0:settings%h1),RTI%live(settings%h0:settings%h1,i_live,i_cluster) )
                if(temp_distance2 < closest_distance2) then
                    cluster = i_cluster
                    closest_distance2 = temp_distance2
                end if
            end do
        end do

    end function identify_cluster






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

















end module nested_sampling_module
