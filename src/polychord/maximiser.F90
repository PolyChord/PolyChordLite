module maximise_module
    use utils_module, only: dp

    implicit none

    contains

    subroutine maximise(loglikelihood,prior,settings,RTI,num_repeats,mpi_information)
        use utils_module, only: stdout_unit
        use settings_module,  only: program_settings
        use run_time_module,   only: run_time_info
        use read_write_module,   only: write_max_file, mean
        use chordal_module, only: generate_nhats, slice_sample
#ifdef MPI
        use mpi_module, only: mpi_bundle,is_root, throw_point, catch_point, mpi_synchronise, throw_seed, catch_seed
#else
        use mpi_module, only: mpi_bundle,is_root
#endif
        interface
            function loglikelihood(theta,phi)
                import :: dp
                real(dp), intent(in), dimension(:)  :: theta
                real(dp), intent(out), dimension(:) :: phi
                real(dp) :: loglikelihood
            end function
        end interface
        interface
            function prior(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function
        end interface

        !> Program settings
        type(program_settings), intent(in) :: settings

        ! The run time info (very important, see src/run_time_info.f90)
        type(run_time_info) :: RTI
        integer, dimension(size(settings%grade_dims)), intent(in) :: num_repeats

        type(mpi_bundle),intent(in) :: mpi_information

        real(dp) :: max_loglike, max_loglike0
        integer :: imax(2)
        integer :: cluster_id, i, j, epoch
        real(dp), dimension(settings%nTotal) :: max_point, max_point0, max_point1, mean_point
        real(dp), dimension(settings%nDims,settings%nDims)   :: cholesky
        integer,   dimension(sum(num_repeats))   :: speeds ! The speed of each nhat
        real(dp),    dimension(settings%nDims,sum(num_repeats))   :: nhats
        integer, dimension(size(settings%grade_dims)) :: nlike      ! Temporary storage for number of likelihood calls
        real(dp) :: w, dx
        real(dp),    dimension(settings%nDims)   :: nhat, x0, x1
#ifdef MPI
            call mpi_synchronise(mpi_information)
#endif

        ! Get highest likelihood point
        if(is_root(mpi_information)) then
            write(stdout_unit,'("Maximising")') 
            imax = maxloc(RTI%live(settings%l0,:,:))
            cluster_id = imax(2)
            max_point = RTI%live(:,imax(1),cluster_id)
            max_loglike = max_point(settings%l0)
            cholesky = RTI%cholesky(:,:,cluster_id)
        end if

        do while(.true.)

#ifdef MPI
            call mpi_synchronise(mpi_information)
            ! Synchronise details for maximising
            if(is_root(mpi_information)) then
                do i=1,mpi_information%nprocs-1
                    call throw_seed(max_point,cholesky,max_loglike,mpi_information,i,j,.true.)
                end do
            else
                if(.not. catch_seed(max_point,cholesky,max_loglike,j,mpi_information)) exit
            end if
#endif

            ! Do maximisation
            nhats = generate_nhats(settings,num_repeats,speeds)
            nhats = matmul(cholesky,nhats)
            max_point0 = max_point
            max_loglike0 = max_loglike
            nlike = 0
            do i=1,size(nhats,2)
                nhat = nhats(:,i)
                w = sqrt(dot_product(nhat,nhat))
                w = w * 3d0 
                max_point1 = slice_sample(loglikelihood,prior,max_loglike,nhat,max_point,1d0,settings,nlike(speeds(i))) 
                if (max_point1(settings%l0) > max_loglike) then
                    max_loglike = max_point1(settings%l0)
                    max_point = max_point1
                end if
            end do

#ifdef MPI
            ! Catch maxima from all cores
            if(is_root(mpi_information)) then
                do i=1,mpi_information%nprocs-1
                    j = catch_point(max_point1,mpi_information)
                    if (max_point1(settings%l0) > max_loglike) then
                        max_loglike = max_point1(settings%l0)
                        max_point = max_point1
                    end if
                end do
            else
                call throw_point(max_point,mpi_information)
            end if
            call mpi_synchronise(mpi_information)
#endif

            if(is_root(mpi_information)) then
                x0 = max_point0(settings%h0:settings%h1) 
                x1 = max_point(settings%h0:settings%h1) 
                dx = maxval(abs(x0-x1))
                if (dx < 1e-5 .or. max_loglike - max_loglike0 < 1e-5) then
#ifdef MPI
                call mpi_synchronise(mpi_information)
                    do i=1,mpi_information%nprocs-1
                        call throw_seed(max_point,cholesky,max_loglike,mpi_information,i,j,.false.)
                    end do
#endif
                    exit
                endif
                write(stdout_unit,'("Loglike: ", F15.5, " change: ", F15.5 )') max_loglike, max_loglike - max_loglike0 
                write(stdout_unit,*) nlike
            end if
        end do
        if(is_root(mpi_information)) then
            if (settings%posteriors) then
                mean_point(settings%p0:settings%d1) = mean(RTI, settings)
                mean_point(settings%l0) = loglikelihood(mean_point(settings%p0:settings%p1),mean_point(settings%d0:settings%d1)) 
                call write_max_file(settings, max_point, mean_point)
            else
                call write_max_file(settings, max_point)
            end if
        end if


    end subroutine maximise

end module maximise_module
