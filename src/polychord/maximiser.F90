module maximise_module
    use utils_module, only: dp

    implicit none

    contains

    subroutine maximise(loglikelihood,prior,settings,RTI,num_repeats)
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


        real(dp) :: max_loglike, max_loglike0
        integer :: imax(2)
        integer :: cluster_id, i, j, epoch
        real(dp), dimension(settings%nTotal) :: max_point, max_posterior_point, mean_point
        real(dp), parameter :: dl=1d-8

        ! Get highest likelihood point
        write(stdout_unit,'("-------------------------------------")') 
        write(stdout_unit,'("Maximising Likelihood")') 
        max_point = do_maximisation(loglikelihood,prior,settings,RTI,num_repeats, dl, .false.) 

        write(stdout_unit,'("-------------------------------------")') 
        write(stdout_unit,'("Maximising Posterior")') 
        max_posterior_point = do_maximisation(loglikelihood,prior,settings,RTI,num_repeats, dl, .true.) 

        if (settings%posteriors) then
            mean_point(settings%p0:settings%d1) = mean(RTI, settings)
            mean_point(settings%l0) = loglikelihood(mean_point(settings%p0:settings%p1),mean_point(settings%d0:settings%d1)) 
            call write_max_file(settings, max_point, max_posterior_point, dXdtheta(prior,max_posterior_point(settings%h0:settings%h1)),  mean_point)
        else
            call write_max_file(settings, max_point, max_posterior_point, dXdtheta(prior,max_posterior_point(settings%h0:settings%h1)))
        end if


    end subroutine maximise


    function do_maximisation(loglikelihood,prior,settings,RTI,num_repeats, dl, posterior) result(max_point)
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
        type(run_time_info), intent(in) :: RTI
        integer, dimension(size(settings%grade_dims)), intent(in) :: num_repeats
        logical, intent(in) :: posterior
        real(dp), intent(in) :: dl


        real(dp) :: max_loglike, max_loglike0
        integer :: imax(2)
        integer :: cluster_id, i, j, epoch
        real(dp), dimension(settings%nTotal) :: max_point, max_point0, max_point1, mean_point
        real(dp), dimension(settings%nDims,settings%nDims)   :: cholesky
        integer,   dimension(sum(num_repeats))   :: speeds ! The speed of each nhat
        real(dp),    dimension(settings%nDims,sum(num_repeats))   :: nhats
        integer, dimension(size(settings%grade_dims)) :: nlike      ! Temporary storage for number of likelihood calls
        real(dp) :: w
        real(dp),    dimension(settings%nDims)   :: nhat, x0, x1
        logical temp

        ! Get highest likelihood point
        imax = maxloc(RTI%live(settings%l0,:,:))
        cluster_id = imax(2)
        max_point = RTI%live(:,imax(1),cluster_id)
        max_loglike = max_point(settings%l0)
        if (posterior) max_loglike = max_loglike + dXdtheta(prior, max_point(settings%h0:settings%h1)) 
        cholesky = RTI%cholesky(:,:,cluster_id)

        do while(.true.)

            ! Do maximisation
            nhats = generate_nhats(settings,num_repeats,speeds)
            nhats = matmul(cholesky,nhats)
            max_point0 = max_point
            max_loglike0 = max_loglike
            nlike = 0
            do i=1,size(nhats,2)
                nhat = nhats(:,i)
                max_point = maximise_direction(loglikelihood,prior,settings, max_point, nhat, dl, posterior)
            end do

            max_loglike = max_point(settings%l0)
            if (max_loglike - max_loglike0 < 1e-5) exit
            write(stdout_unit,'("Loglike: ", F15.8, " change: ", F15.8 )') max_loglike, max_loglike - max_loglike0 
        end do


    end function do_maximisation

    function maximise_direction(loglikelihood,prior,settings, point, n, dl, posterior) result(max_point)
        use utils_module, only: stdout_unit
        use settings_module,  only: program_settings
        use run_time_module,   only: run_time_info
        use read_write_module,   only: write_max_file, mean
        use chordal_module, only: generate_nhats, slice_sample
        use calculate_module, only: calculate_point
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
        logical, intent(in) :: posterior
        real(dp), intent(in) :: dl

        real(dp), dimension(settings%nTotal) :: max_point, point, point_a, point_b, point_c, point_d
        real(dp),    dimension(settings%nDims)   :: n, x, a, b, c, d
        integer nlike
        real(dp),parameter :: phi = (1+sqrt(5.))/2


        ! Construct initial bracket
        point_a = point
        point_b = point
        if (posterior) point_b(settings%l0) = point_b(settings%l0)  + dXdtheta(prior, point_b(settings%h0:settings%h1))

        do while (.true.)
            point_a(settings%h0:settings%h1) = point_a(settings%h0:settings%h1) - n
            call calculate_point(loglikelihood,prior,point_a,settings,nlike) 
            if (posterior) point_a(settings%l0) = point_a(settings%l0)  + dXdtheta(prior, point_a(settings%h0:settings%h1))
            if (point_a(settings%l0) <= point(settings%l0)) exit
            point_b = point
            point = point_a
        end do

        do while (point_b(settings%l0) > point(settings%l0))
            point = point_b
            point_b(settings%h0:settings%h1) = point_b(settings%h0:settings%h1) + n
            call calculate_point(loglikelihood,prior,point_b,settings,nlike) 
            if (posterior) point_b(settings%l0) = point_b(settings%l0)  + dXdtheta(prior, point_b(settings%h0:settings%h1))
        end do

        ! Now do golden section search
        a = point_a(settings%h0:settings%h1)
        b = point_b(settings%h0:settings%h1)

        point_c(settings%h0:settings%h1) = b - (b-a)/phi
        call calculate_point(loglikelihood,prior,point_c,settings,nlike) 
        if (posterior) point_c(settings%l0) = point_c(settings%l0)  + dXdtheta(prior, point_c(settings%h0:settings%h1))

        point_d(settings%h0:settings%h1) = a + (b-a)/phi
        call calculate_point(loglikelihood,prior,point_d,settings,nlike) 
        if (posterior) point_d(settings%l0) = point_d(settings%l0)  + dXdtheta(prior, point_d(settings%h0:settings%h1))

        do while ( abs(point_a(settings%l0)-point_b(settings%l0)) > dl ) 
            if (point_c(settings%l0) > point_d(settings%l0)) then
                point_b = point_d
                point_d = point_c
                a = point_a(settings%h0:settings%h1)
                b = point_b(settings%h0:settings%h1)
                point_c(settings%h0:settings%h1) = b - (b-a)/phi
                call calculate_point(loglikelihood,prior,point_c,settings,nlike) 
                if (posterior) point_c(settings%l0) = point_c(settings%l0)  + dXdtheta(prior, point_c(settings%h0:settings%h1))
            else
                point_a = point_c
                point_c = point_d
                a = point_a(settings%h0:settings%h1)
                b = point_b(settings%h0:settings%h1)
                point_d(settings%h0:settings%h1) = a + (b-a)/phi
                call calculate_point(loglikelihood,prior,point_d,settings,nlike) 
                if (posterior) point_d(settings%l0) = point_d(settings%l0)  + dXdtheta(prior, point_d(settings%h0:settings%h1))
            end if
        end do

        if (point_c(settings%l0) > point_d(settings%l0)) then
            max_point = point_c
        else
            max_point = point_d
        end if

    end function maximise_direction

    function dXdtheta(prior, cube, dx_)
        implicit none
        interface
            function prior(cube) result(theta)
                import :: dp
                real(dp), intent(in), dimension(:) :: cube
                real(dp), dimension(size(cube))    :: theta
            end function
        end interface
        real(dp), intent(in), dimension(:) :: cube
        real(dp), optional, intent(in) ::  dx_
        real(dp) dXdtheta

        real(dp), dimension(size(cube)) :: cube0
        real(dp), dimension(size(cube),size(cube)) :: dtheta
        integer :: i, s
        real(dp) :: dx

        dx = 1e-8
        if (present(dx_)) dx=dx_

        s=1
        do i=1,size(cube)
            cube0 = cube
            if (cube0(i) + dx >= 1) then
                cube0(i) = cube0(i) - dx
                s = -s 
            else
                cube0(i) = cube0(i) + dx
            end if
            dtheta(:,i) = prior(cube0) - prior(cube)
        end do
        dXdtheta =  size(cube)*log(dx) - log(s*det(dtheta))

    end function dXdtheta

    function det(matrix)
        implicit none
        real(dp), dimension(:,:) :: matrix
        real(dp) det
        real(dp) :: m, temp
        integer :: n, i, j, k, l
        logical :: detexists = .true.
        n = size(matrix,1)
        l = 1
        !convert to upper triangular form
        do k = 1, n-1
            if (matrix(k,k) == 0) then
                detexists = .false.
                do i = k+1, n
                    if (matrix(i,k) /= 0) then
                        do j = 1, n
                            temp = matrix(i,j)
                            matrix(i,j)= matrix(k,j)
                            matrix(k,j) = temp
                        end do
                        detexists = .true.
                        l=-l
                        exit
                    endif
                end do
                if (detexists .eqv. .false.) then
                    det = 0
                    return
                end if
            endif
            do j = k+1, n
                m = matrix(j,k)/matrix(k,k)
                do i = k+1, n
                    matrix(j,i) = matrix(j,i) - m*matrix(k,i)
                end do
            end do
        end do
        
        !calculate determinant by finding product of diagonal elements
        det = l
        do i = 1, n
            det = det * matrix(i,i)
        end do
        
    end function det



end module maximise_module
