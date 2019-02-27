module maximise_module
    use utils_module, only: dp

    implicit none

    contains

    subroutine maximise(loglikelihood,prior,settings,RTI,mpi_information)
        use utils_module, only: stdout_unit
        use settings_module,  only: program_settings
        use run_time_module,   only: run_time_info
        use read_write_module,   only: write_max_file
        use chordal_module, only: generate_nhats, slice_sample
#ifdef MPI
        use mpi_module, only: mpi_bundle,is_root,linear_mode
#else
        use mpi_module, only: mpi_bundle,is_root,linear_mode
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

        type(mpi_bundle),intent(in) :: mpi_information

        real(dp) :: max_loglike, max_loglike0
        integer :: imax(2)
        integer :: cluster_id, i
        real(dp), dimension(settings%nTotal) :: max_point, max_point0
        real(dp), dimension(settings%nDims,settings%nDims)   :: cholesky
        integer,   dimension(sum(RTI%num_repeats))   :: speeds ! The speed of each nhat
        real(dp),    dimension(settings%nDims,sum(RTI%num_repeats))   :: nhats
        integer, dimension(size(settings%grade_dims)) :: nlike      ! Temporary storage for number of likelihood calls
        real(dp) :: w, dx
        real(dp),    dimension(settings%nDims)   :: nhat, x0, x1

        ! Get highest likelihood point
        imax = maxloc(RTI%live(settings%l0,:,:))
        cluster_id = imax(2)
        max_point = RTI%live(:,imax(1),cluster_id)
        max_loglike = max_point(settings%l0)
        cholesky = RTI%cholesky(:,:,cluster_id)

        write(stdout_unit,'("Maximising")') 
        do while(.true.)

            nhats = generate_nhats(settings,RTI%num_repeats,speeds)
            nhats = matmul(cholesky,nhats)
            max_point0 = max_point
            max_loglike0 = max_loglike
            do i=1,size(nhats,2)
                nhat = nhats(:,i)
                w = sqrt(dot_product(nhat,nhat))
                w = w * 3d0 
                max_point = slice_sample(loglikelihood,prior,max_loglike,nhat,max_point,1d0,settings,nlike(speeds(i))) 
                max_loglike = max(max_point(settings%l0),max_loglike)
            end do

            x0 = max_point0(settings%h0:settings%h1) 
            x1 = max_point(settings%h0:settings%h1) 
            dx = sqrt(dot_product(x0-x1,x0-x1))
            if (dx < 1e-5) exit
            if (max_loglike - max_loglike0 < 1e-4) exit
            write(stdout_unit,'("Loglike: ", F15.5, " change: ", F15.5 )') max_loglike, max_loglike - max_loglike0 
        end do
        call write_max_file(settings, max_point)


    end subroutine maximise

end module maximise_module
