module mpi_module

#ifdef MPI
    use mpi
#endif
    implicit none

    integer :: mpierror

    integer, parameter :: tag_gen_point=1
    integer, parameter :: tag_gen_request=2
    integer, parameter :: tag_gen_stop=3

    integer, parameter :: tag_run_baby=4
    integer, parameter :: tag_run_seed=5
    integer, parameter :: tag_run_cholesky=6
    integer, parameter :: tag_run_logL=7
    integer, parameter :: tag_run_nlike=8
    integer, parameter :: tag_run_epoch_seed=9
    integer, parameter :: tag_run_epoch_babies=10
    integer, parameter :: tag_run_stop=11

    type mpi_bundle
        integer :: rank
        integer :: nprocs
        integer :: root
        integer :: colour
        integer :: communicator

    end type mpi_bundle



    contains

    !> Subroutine to get all mpi information
    function get_mpi_information(mpi_communicator,colour) result(mpi_information)
        implicit none
        integer, intent(in)           :: mpi_communicator
        integer, intent(in), optional :: colour
        type(mpi_bundle)                :: mpi_information

        mpi_information%rank         = get_rank(mpi_communicator)
        mpi_information%nprocs       = get_nprocs(mpi_communicator)
        mpi_information%root         = get_root(mpi_communicator)
        mpi_information%communicator = mpi_communicator

        if(present(colour)) then
            mpi_information%colour = colour
        else
            mpi_information%colour = 0
        end if

    end function

    !> Returns whether this is the root node
    function is_root(mpi_information) 
        implicit none
        type(mpi_bundle),intent(in) :: mpi_information
        logical :: is_root

        is_root = mpi_information%rank==mpi_information%root

    end function is_root

    !> Returns whether there is a single processor
    function linear_mode(mpi_information) 
        implicit none
        type(mpi_bundle),intent(in) :: mpi_information
        logical :: linear_mode

        linear_mode = mpi_information%nprocs==1

    end function linear_mode

    !> Procedure to get the number of mpi processors
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Comm_size.html
    function get_nprocs(mpi_communicator) result(nprocs)
        implicit none
        integer, intent(in) :: mpi_communicator
        integer :: nprocs

#ifdef MPI
        call MPI_COMM_SIZE(   & 
            mpi_communicator, &!handle
            nprocs,           &!return number of processors
            mpierror          &!error flag
            )
#else
        nprocs = 1
#endif

    end function get_nprocs

    !> Procedure to get the number of mpi processors
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Comm_rank.html
    function get_rank(mpi_communicator) result(myrank)
        implicit none
        integer, intent(in) :: mpi_communicator
        integer :: myrank

#ifdef MPI
        call MPI_COMM_RANK(   & 
            mpi_communicator, &!handle
            myrank,           &!return rank of calling processor 
            mpierror          &!error flag
            )
#else
        myrank = 0
#endif

    end function get_rank


    !> This defines the root node for this mpi communicator
    !!
    !! All processes call this simultaneously, and the root is defined
    !! as the process with the lowest rank. This value is returned
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Allreduce.html
    !!
    function get_root(mpi_communicator) result(root)
        implicit none
        integer, intent(in) :: mpi_communicator
        integer :: root
        integer :: myrank

        ! Get the rank of the process
        myrank = get_rank(mpi_communicator)

#ifdef MPI
        call MPI_ALLREDUCE(    &
            myrank,            &!send buffer 
            root,              &!recieve buffer
            1,                 &!number of elements sent
            MPI_INTEGER,       &!type of element sent
            MPI_MIN,           &!reduce by finding the minimum
            mpi_communicator,  &!handle
            mpierror           &!error flag
            )
#else
        root = myrank
#endif

    end function get_root


#ifdef MPI
    !> Procedure to initialise mpi
    subroutine initialise_mpi()
        implicit none

        call MPI_INIT(mpierror)

    end subroutine initialise_mpi

    !> Procedure to finalise mpi
    subroutine finalise_mpi()
        implicit none

        call MPI_FINALIZE(mpierror)

    end subroutine finalise_mpi




    !> Split a communicator into n even groups
    function mpi_split(n,mpi_communicator) result(mpi_information)
        implicit none
        integer,intent(in)  :: n
        integer,intent(in)  :: mpi_communicator
        type(mpi_bundle)      :: mpi_information

        integer :: new_mpi_communicator
        integer :: numprocs
        integer :: colour
        integer :: key

        ! Get the original mpi info
        mpi_information = get_mpi_information(mpi_communicator)

        ! Define the new number of processors
        numprocs= ceiling(dble(mpi_information%nprocs)/dble(n))

        ! Define the 'colour' of this process by dividing it into n adjacent processes
        colour = mpi_information%rank / numprocs

        ! The new rank is just modulo this
        key  = mod(mpi_information%rank, numprocs)

        ! Split up the communicator
        call MPI_COMM_SPLIT(mpi_communicator,colour,key,new_mpi_communicator,mpierror)

        ! Assign the new mpi info
        mpi_information = get_mpi_information(new_mpi_communicator,colour)

    end function mpi_split

    subroutine mpi_synchronise(mpi_information)
        implicit none
        type(mpi_bundle), intent(in) :: mpi_information

#ifdef MPI
        call MPI_BARRIER(mpi_information%communicator,mpierror)
#endif 

    end subroutine mpi_synchronise


    !> This sums a whole set of integers across all processes
    !!
    !! All processes call this simultaneously
    !!
    !! http://www.mpich.org/static/docs/v3.1/www3/MPI_Allreduce.html
    !!
    function sum_integers(intgr_local,mpi_information) result(intgr)
        implicit none
        integer, intent(in) :: intgr_local
        type(mpi_bundle), intent(in) :: mpi_information
        integer :: intgr

        call MPI_ALLREDUCE(    &
            intgr_local,       &!send buffer 
            intgr,             &!recieve buffer
            1,                 &!number of elements sent
            MPI_INTEGER,       &!type of element sent
            MPI_SUM,           &!reduce by finding the minimum
            mpi_information%communicator,  &!handle
            mpierror           &!error flag
            )

    end function sum_integers

    function sum_doubles(db_local,mpi_information) result(db)
        implicit none
        double precision, intent(in) :: db_local
        type(mpi_bundle), intent(in) :: mpi_information
        double precision :: db

        call MPI_ALLREDUCE(       &
            db_local,             &!send buffer 
            db,                   &!recieve buffer
            1,                    &!number of elements sent
            MPI_DOUBLE_PRECISION, &!type of element sent
            MPI_SUM,              &!reduce by finding the minimum
            mpi_information%communicator,     &!handle
            mpierror              &!error flag
            )

    end function sum_doubles

    subroutine broadcast_doubles(doubles,mpi_information)
        implicit none
        double precision, dimension(:), intent(inout) :: doubles
        type(mpi_bundle), intent(in) :: mpi_information

        call MPI_BCAST(            & 
            doubles,               &!broadcast buffer
            size(doubles),         &!size of buffer
            MPI_DOUBLE_PRECISION,  &!type of element sent
            mpi_information%root,         &!root doing the sending
            mpi_information%communicator, &!handle
            mpierror               &!error flag
            )

    end subroutine broadcast_doubles

    subroutine broadcast_integers(integers,mpi_information)
        implicit none
        integer, dimension(:), intent(inout) :: integers
        type(mpi_bundle), intent(in) :: mpi_information

        call MPI_BCAST(            & 
            integers,              &!broadcast buffer
            size(integers),        &!size of buffer
            MPI_INT,               &!type of element sent
            mpi_information%root,         &!root doing the sending
            mpi_information%communicator, &!handle
            mpierror               &!error flag
            )

    end subroutine broadcast_integers

    !============== Throwing and catching routines ====================
    ! There are four points in the algorithm where arrays need to be transferred 
    ! between master and slaves
    !
    ! 1) live point during initial generation
    !     any slave    ---->  root
    !     throw_point         catch_point
    !
    ! 2) baby points for use in nested sampling
    !     any slave   ----> root
    !     throw_babies      catch_babies
    !
    ! 3) seed information for generation of babies
    !     root      ----> specific slave
    !     throw_seed      catch_seed


    !> Root catch
    !!
    !! This a process by which the root node 'catches' a thrown point from 
    !! any slave, and returns the slave identifier for use to throw back
    function catch_point(live_point,mpi_information) result(slave_id)
        implicit none

        double precision,intent(out),dimension(:) :: live_point !> The caught live point
        type(mpi_bundle), intent(in) :: mpi_information

        integer :: slave_id ! slave identifier

        integer, dimension(MPI_STATUS_SIZE) :: mpistatus ! status identifier

        call MPI_RECV(             &!
            live_point,            &!
            size(live_point),      &!
            MPI_DOUBLE_PRECISION,  &!
            MPI_ANY_SOURCE,        &!
            tag_gen_point,         &!
            mpi_information%communicator,      &!
            mpistatus,            &!
            mpierror               &!
            )

        slave_id = mpistatus(MPI_SOURCE) ! Pass on the slave id

    end function catch_point


    !> Slave throws a point to the master
    !!
    !! This a process by which a slave node 'throws' a point to the root
    subroutine throw_point(live_point,mpi_information)
        implicit none

        double precision,intent(in),dimension(:) :: live_point !> live point to throw
        type(mpi_bundle), intent(in) :: mpi_information

        call MPI_SEND(             &!
            live_point,            &!
            size(live_point),      &!
            MPI_DOUBLE_PRECISION,  &!
            mpi_information%root,                  &!
            tag_gen_point,         &!
            mpi_information%communicator,      &!
            mpierror               &!
            )

    end subroutine throw_point







    !> Master catches babies thrown by any slave, and returns the slave identity that did the throwing
    function catch_babies(baby_points,nlike,epoch,mpi_information) result(slave_id)
        implicit none

        double precision,intent(out),dimension(:,:) :: baby_points !> The babies to be caught
        integer, dimension(:), intent(out)          :: nlike       !> The number of likelihood evaluations to be caught
        integer,               intent(out)          :: epoch       !> The epoch the points were generated in
        type(mpi_bundle), intent(in)                :: mpi_information    !> The mpi communicator

        integer :: slave_id ! slave identifier

        integer, dimension(MPI_STATUS_SIZE) :: mpistatus ! status identifier

        call MPI_RECV(                               &!
            baby_points,                             &!
            size(baby_points,1)*size(baby_points,2), &!
            MPI_DOUBLE_PRECISION,                    &!
            MPI_ANY_SOURCE,                          &!
            tag_run_baby,                            &!
            mpi_information%communicator,            &!
            mpistatus,                               &!
            mpierror                                 &!
            )

        ! Pass on the slave id
        slave_id = mpistatus(MPI_SOURCE)

        call MPI_RECV(            &! 
            nlike,                &! 
            size(nlike),          &! 
            MPI_INT,              &! 
            slave_id,             &! 
            tag_run_nlike,        &! 
            mpi_information%communicator,&! 
            mpistatus,            &! 
            mpierror              &! 
            )

        call MPI_RECV(            &! 
            epoch,                &! 
            1,                    &! 
            MPI_INT,              &! 
            slave_id,             &! 
            tag_run_epoch_babies, &! 
            mpi_information%communicator,&! 
            mpistatus,            &! 
            mpierror              &! 
            )

    end function catch_babies

    !> Slave throws babies to the master
    subroutine throw_babies(baby_points,nlike,epoch,mpi_information)
        implicit none

        double precision,intent(in),dimension(:,:) :: baby_points !> The babies to be thrown
        integer, dimension(:), intent(in) :: nlike                !> The number of likelihood evaluations to be caught
        integer,               intent(in) :: epoch                !> The epoch the babies were generated in
        type(mpi_bundle), intent(in) :: mpi_information

        call MPI_SEND(                               &! 
            baby_points,                             &! 
            size(baby_points,1)*size(baby_points,2), &! 
            MPI_DOUBLE_PRECISION,                    &! 
            mpi_information%root,                                    &! 
            tag_run_baby,                            &! 
            mpi_information%communicator,                        &! 
            mpierror                                 &! 
            )                                        
        call MPI_SEND(         &!  
            nlike,             &!  
            size(nlike),       &!  
            MPI_INT,           &!  
            mpi_information%root,              &!  
            tag_run_nlike,     &!  
            mpi_information%communicator,  &!  
            mpierror           &!  
            )
        call MPI_SEND(         &!  
            epoch,             &!  
            1,                 &!  
            MPI_INT,           &!  
            mpi_information%root,              &!  
            tag_run_epoch_babies,&!  
            mpi_information%communicator,  &!  
            mpierror           &!  
            )

    end subroutine throw_babies






    !> slave catches seed thrown by master
    function catch_seed(seed_point,cholesky,logL,epoch,mpi_information) result(more_points_needed)
        use abort_module, only: halt_program
        implicit none


        double precision,intent(out),dimension(:) :: seed_point  !> The seed point to be caught
        double precision,intent(out),dimension(:,:) :: cholesky  !> Cholesky matrix to be caught
        double precision,intent(out)               :: logL       !> loglikelihood contour to be caught
        integer,         intent(out)               :: epoch
        type(mpi_bundle), intent(in)               :: mpi_information

        logical :: more_points_needed ! whether or not we need more points

        integer, dimension(MPI_STATUS_SIZE) :: mpistatus ! status identifier


        call MPI_RECV(                      &!
            seed_point,                     &!
            size(seed_point),               &!
            MPI_DOUBLE_PRECISION,           &!
            mpi_information%root,                           &!
            MPI_ANY_TAG,                    &!
            mpi_information%communicator,               &!
            mpistatus,                     &!
            mpierror                        &!
            )
        if(mpistatus(MPI_TAG) == tag_run_stop ) then
            more_points_needed = .false.
            return
        else if(mpistatus(MPI_TAG) == tag_run_seed) then
            more_points_needed = .true.
        else
            call halt_program('slave error: unrecognised tag')
        end if

        call MPI_RECV(                        &!
            cholesky,                         &!
            size(cholesky,1)*size(cholesky,1),&!
            MPI_DOUBLE_PRECISION,             &!
            mpi_information%root,                             &!
            tag_run_cholesky,                 &!
            mpi_information%communicator,                 &!
            mpistatus,                       &!
            mpierror                          &!
            )
        call MPI_RECV(                 &!
            logL,                      &!
            1,                         &!
            MPI_DOUBLE_PRECISION,      &!
            mpi_information%root,      &!
            tag_run_logL,              &!
            mpi_information%communicator,          &!
            mpistatus,                 &!
            mpierror                   &!
            )
        call MPI_RECV(            &! 
            epoch,                &! 
            1,                    &! 
            MPI_INT,              &! 
            mpi_information%root, &! 
            tag_run_epoch_seed,   &! 
            mpi_information%communicator,&! 
            mpistatus,            &! 
            mpierror              &! 
            )

    end function catch_seed



    !> root throws seed to slave
    subroutine throw_seed(seed_point,cholesky,logL,mpi_information,slave_id,epoch,keep_going)
        implicit none

        double precision,intent(in),dimension(:) :: seed_point   !> seed to be thrown
        double precision,intent(in),dimension(:,:) :: cholesky   !> cholesky to be thrown
        double precision,intent(in)                :: logL       !> loglikelihood contour to be thrown
        type(mpi_bundle),intent(in) :: mpi_information           !> mpi handle
        integer, intent(in) :: slave_id                          !> identity of target slave
        integer, intent(in) :: epoch                             !> epoch of seed
        logical, intent(in) :: keep_going                        !> Further signal whether to keep going 

        integer :: tag ! tag variable to

        tag = tag_run_stop                 ! Default tag is stop tag
        if(keep_going) tag = tag_run_seed  ! If we want to keep going then change this to the seed tag


        call MPI_SEND(             &!  
            seed_point,            &!  
            size(seed_point),      &!  
            MPI_DOUBLE_PRECISION,  &!  
            slave_id,              &!  
            tag,                   &!  
            mpi_information%communicator,      &!  
            mpierror               &!  
            )

        if(.not. keep_going) return ! Stop here if we're wrapping up

        call MPI_SEND(                          &!  
            cholesky,                           &!  
            size(cholesky,1)*size(cholesky,2),  &!  
            MPI_DOUBLE_PRECISION,               &!  
            slave_id,                           &!  
            tag_run_cholesky,                   &!  
            mpi_information%communicator,                   &!  
            mpierror                            &!  
            )
        call MPI_SEND(              &!  
            logL,                   &!  
            1,                      &!  
            MPI_DOUBLE_PRECISION,   &!  
            slave_id,               &!  
            tag_run_logL,           &!  
            mpi_information%communicator,       &!  
            mpierror                &!  
            )
        call MPI_SEND(         &!  
            epoch,             &!  
            1,                 &!  
            MPI_INT,           &!  
            slave_id,          &!  
            tag_run_epoch_seed,&!  
            mpi_information%communicator,  &!  
            mpierror           &!  
            )


    end subroutine throw_seed




    !============== Pure messaging routines ===========================
    ! During initial live point generation, the master needs to signal to the slaves
    ! whether or not to keep generating live points, or whether to stop
    !
    ! master         ----> slave
    ! request_point        more_points_needed -> true
    ! no_more_points       more_points_needed -> false
    ! 

    !> Request point
    !!
    !! This subroutine is used by the root node to request a new live point
    subroutine request_point(mpi_information,slave_id)
        implicit none
        type(mpi_bundle), intent(in) :: mpi_information
        integer, intent(in) :: slave_id         !> Slave to request a new point from


        integer :: empty_buffer(0) ! empty buffer to send

        call MPI_SEND(         &
            empty_buffer,      &! not sending anything
            0,                 &! size of nothing
            MPI_INT,           &! sending no integers
            slave_id,          &! process id to send to
            tag_gen_request,   &! continuation tag
            mpi_information%communicator,  &! mpi handle
            mpierror           &! error flag
            )

    end subroutine request_point


    !> No more points please
    !!
    !! This subroutine is used by the root node to signal that no more points are required
    subroutine no_more_points(mpi_information,slave_id)
        implicit none
        type(mpi_bundle), intent(in) :: mpi_information
        integer, intent(in) :: slave_id         !> Slave to request a new point from


        integer :: empty_buffer(0) ! empty buffer to send

        call MPI_SEND(         &
            empty_buffer,      &! not sending anything
            0,                 &! size of nothing
            MPI_INT,           &! sending no integers
            slave_id,          &! process id to send to
            tag_gen_stop,      &! continuation tag
            mpi_information%communicator,  &! mpi handle
            mpierror           &! error flag
            )

    end subroutine no_more_points

    !> See if more points are needed
    !!
    !! This subroutine is used by the root node to request a new live point
    function more_points_needed(mpi_information)
        use abort_module
        implicit none
        type(mpi_bundle), intent(in) :: mpi_information

        integer :: empty_buffer(0) ! empty buffer to send
        integer, dimension(MPI_STATUS_SIZE) :: mpistatus  ! status identifier

        logical :: more_points_needed !> Whether we need more points or not

        call MPI_RECV(       &!
            empty_buffer,    &!
            0,               &!
            MPI_INT,         &!
            mpi_information%root,            &!
            MPI_ANY_TAG,     &!
            mpi_information%communicator,&!
            mpistatus,      &!
            mpierror         &!
            )

        ! If we've recieved a kill signal, then exit this loop
        if(mpistatus(MPI_TAG) == tag_gen_stop ) then
            more_points_needed = .false.
        else if(mpistatus(MPI_TAG) == tag_gen_request) then
            more_points_needed = .true.
        else
            call halt_program('generate error: unrecognised tag')
        end if

    end function more_points_needed

#endif








end module mpi_module
