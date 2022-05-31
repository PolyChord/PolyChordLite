module maximise_module
    use utils_module, only: dp

    implicit none

    contains

    subroutine maximise(loglikelihood,prior,settings,RTI)
        use utils_module, only: stdout_unit
        use settings_module,  only: program_settings
        use run_time_module,   only: run_time_info
        use read_write_module,   only: write_max_file, mean
        use chordal_module, only: generate_nhats, slice_sample
#ifdef USE_MPI
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
        real(dp), dimension(settings%nTotal) :: max_point, max_posterior_point, mean_point

        ! Get highest likelihood point
        write(stdout_unit,'("-------------------------------------")') 
        write(stdout_unit,'("Maximising Likelihood")') 
        max_point = do_maximisation(loglikelihood,prior,settings,RTI, .false.) 

        write(stdout_unit,'("-------------------------------------")') 
        write(stdout_unit,'("Maximising Posterior")') 
        max_posterior_point = do_maximisation(loglikelihood,prior,settings,RTI, .true.) 

        if (settings%posteriors) then
            mean_point(settings%p0:settings%d1) = mean(RTI, settings)
            mean_point(settings%l0) = loglikelihood(mean_point(settings%p0:settings%p1),mean_point(settings%d0:settings%d1)) 
            call write_max_file(settings, max_point, max_posterior_point,&
                                dXdtheta(prior,max_posterior_point(settings%h0:settings%h1)),  mean_point)
        else
            call write_max_file(settings, max_point, max_posterior_point,&
                                dXdtheta(prior,max_posterior_point(settings%h0:settings%h1)))
        end if


    end subroutine maximise


    function do_maximisation(loglikelihood,prior,settings,RTI, posterior) result(max_point)
        use calculate_module, only: calculate_point
        use utils_module, only: sort_doubles
        use settings_module,  only: program_settings
        use run_time_module,   only: run_time_info
        use read_write_module,   only: write_max_file, mean
        use chordal_module, only: generate_nhats, slice_sample
        use nelder_mead_module, only: nelder_mead
#ifdef USE_MPI
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
        logical, intent(in) :: posterior


        real(dp) :: max_loglike
        integer :: cluster_id
        real(dp), dimension(settings%nTotal) :: max_point
        integer :: nlike, j
        real(dp), dimension(settings%nDims)   :: x
        real(dp), dimension(settings%nDims,settings%nDims+1) :: simplex
        real(dp), dimension(settings%nDims+1) :: f
        real(dp), dimension(settings%nlive) :: l
        integer, dimension(settings%nlive) :: i
        ! f needs to be the posterior, not the likelihood

        ! Get highest likelihood points
        max_loglike = settings%logzero
        do cluster_id=1,RTI%ncluster
            if (RTI%nlive(cluster_id) >= settings%nDims +1) then
                l = RTI%live(settings%l0,:RTI%nlive(cluster_id), cluster_id)
                if (posterior) then
                    do j=1,RTI%nlive(cluster_id)
                        l(j) = l(j) + dXdtheta(prior, RTI%live(settings%h0:settings%h1,j, cluster_id))
                    end do
                end if
                i(:RTI%nlive(cluster_id)) = sort_doubles(l(:RTI%nlive(cluster_id)))

                if (l(i(RTI%nlive(cluster_id))) > max_loglike) then
                    max_loglike = l(i(RTI%nlive(cluster_id)))
                    simplex = RTI%live(settings%h0:settings%h1,i(RTI%nlive(cluster_id)-settings%nDims:RTI%nlive(cluster_id)),cluster_id)
                    f = l(i(RTI%nlive(cluster_id)-settings%nDims:RTI%nlive(cluster_id)))
                end if
            end if
        end do
        if (max_loglike > settings%logzero) then
            x = nelder_mead(func, simplex, f, 1d-5)
            max_point(settings%h0:settings%h1) = x
            call calculate_point(loglikelihood,prior,max_point,settings,nlike) 
        else
            write(*,*) "Could not construct simplex"
        end if


        contains
         function func(x)
            use calculate_module, only: calculate_point
            implicit none
            double precision, intent(in), dimension(:) :: x
            double precision :: func
            real(dp), dimension(settings%nTotal) :: point
            integer nlike

            point(settings%h0:settings%h1) = x
            call calculate_point(loglikelihood,prior,point,settings,nlike) 
            func  = point(settings%l0)
            if (posterior .and. func>settings%logzero) then
                func = func + dXdtheta(prior, point(settings%h0:settings%h1))
            end if

         end function


    end function do_maximisation

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

        dx = 1e-5
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
