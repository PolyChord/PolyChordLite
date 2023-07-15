!> Module for reading and writing to files
module read_write_module
    use utils_module, only: dp
    implicit none
    interface read_doubles
        module procedure read_doubles_1, read_doubles_2, read_doubles_3
    end interface read_doubles
    interface write_doubles
        module procedure write_doubles_1, write_doubles_2, write_doubles_3
    end interface write_doubles

    contains

    !> Check if a directory exists
    function directory_exists(directory)
        implicit none
        character(len=*), intent(in) :: directory
        logical :: directory_exists

#ifdef __INTEL_COMPILER
        inquire( directory=directory, exist=directory_exists )         ! Works with ifort, but not gfortran
#else
        inquire( file=trim(directory)//'/.', exist=directory_exists )  ! Works with gfortran, but not ifort
#endif

    end function directory_exists

    subroutine check_directories(settings)
        use settings_module, only: program_settings
        use abort_module, only: halt_program
        implicit none
        type(program_settings), intent(in) :: settings
        if (.not. directory_exists(settings%base_dir) ) then
            call halt_program('directory error: ' // trim(settings%base_dir) // ' does not exist')
        else if (settings%do_clustering .and. .not. directory_exists(cluster_dir(settings))) then
            call halt_program('directory error: ' // trim(cluster_dir(settings)) // ' does not exist')
        end if
    end subroutine

    !> Check to see whether resume file exists
    function resume_file_exists(settings) result(resume)
        use settings_module, only: program_settings
        implicit none

        !> settings variable to get the base_dir and root_name out of
        type(program_settings), intent(in) :: settings
        !> Whether or not we should resume from file
        logical :: resume

        ! See if the resume file with the correct name exists
        inquire(                                     &! Inquire function
            file=trim(resume_file(settings,.false.)),&! file name defined by subroutine resume_file
            exist=resume                             &! return whether this exists
            )

    end function resume_file_exists


    !> Remove resume files
    subroutine delete_files(settings)
        use settings_module, only: program_settings
        use utils_module, only: delete_file,verbose_fb
        implicit none

        !> settings variable to get the base_dir and root_name out of
        type(program_settings), intent(in) :: settings

        integer :: i_cluster ! cluster iterator

        logical :: deleted ! Whether a file has been deleted

        logical :: fb ! Temporary feedback variable

        fb = settings%feedback>=verbose_fb

        deleted = delete_file( stats_file(settings), fb )          ! Delete stats file
        deleted = delete_file( phys_live_file(settings), fb )      ! Delete phys_live file
        deleted = delete_file( resume_file(settings,.false.), fb ) ! Delete temp resume file
        deleted = delete_file( resume_file(settings,.true.), fb )  ! Delete resume file


        ! Delete posterior files
        deleted = delete_file( posterior_file(settings,.true.), fb )
        deleted = delete_file( posterior_file(settings,.false.), fb )

        i_cluster = 1
        do while ( &
                delete_file( posterior_file(settings,.true.,.false.,i_cluster), fb ) .or.  &
                delete_file( posterior_file(settings,.false.,.false.,i_cluster), fb )   &
                )
            i_cluster = i_cluster + 1
        end do


    end subroutine delete_files

    subroutine rename_files(settings,RTI)
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info
        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        integer i_cluster

        ! Move the resume file if resume
        if(settings%write_resume) call rename(trim(resume_file(settings,.true.)),trim(resume_file(settings,.false.)))

        ! Rename the files
        if(settings%equals) then
            call rename(trim(posterior_file(settings,.false.,.true.)),trim(posterior_file(settings,.false.,.false.)))
            do i_cluster = 1,RTI%ncluster+RTI%ncluster_dead
                call rename(trim(posterior_file(settings,.false.,.true.,i_cluster)),trim(posterior_file(settings,.false.,.false.,i_cluster)))
            end do
        end if

        if(settings%posteriors) then
            call rename(trim(posterior_file(settings,.true. ,.true.)),trim(posterior_file(settings,.true. ,.false.)))
            do i_cluster = 1,RTI%ncluster+RTI%ncluster_dead
                call rename(trim(posterior_file(settings,.true.,.true.,i_cluster)),trim(posterior_file(settings,.true.,.false.,i_cluster)))
            end do
        end if
    end subroutine rename_files
    








    subroutine write_integer(a,str)
        use utils_module, only: write_resume_unit,integer_format
        implicit none
        integer, intent(in) :: a
        character(len=*), intent(in),optional :: str

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        write(write_resume_unit,integer_format(1)) a

    end subroutine write_integer

    subroutine write_integers(arr,str)
        use utils_module, only: write_resume_unit,integer_format
        implicit none
        integer,dimension(:), intent(in) :: arr
        character(len=*), intent(in),optional :: str

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        if(size(arr)>0) write(write_resume_unit,integer_format(size(arr))) arr

    end subroutine write_integers
    
    subroutine write_double(a,str)
        use utils_module, only: write_resume_unit,double_format
        implicit none
        real(dp), intent(in) :: a
        character(len=*), intent(in),optional :: str

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        write(write_resume_unit,double_format(1)) a

    end subroutine write_double

    subroutine write_doubles_1(arr,str)
        use utils_module, only: write_resume_unit,double_format
        implicit none
        real(dp),dimension(:), intent(in) :: arr
        character(len=*), intent(in),optional :: str

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        if(size(arr)>0) write(write_resume_unit,double_format(size(arr))) arr

    end subroutine write_doubles_1

    subroutine write_doubles_2(arr,str,n)
        use utils_module, only: write_resume_unit,double_format
        implicit none
        real(dp),dimension(:,:), intent(in) :: arr
        integer, intent(in), optional :: n
        character(len=*), intent(in),optional :: str
        integer :: i

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        if(present(n)) then
            do i=1,n
                call write_doubles( arr(:,i) )
            end do
        else
            do i=1,size(arr,2)
                call write_doubles( arr(:,i) )
            end do
        end if

    end subroutine write_doubles_2

    subroutine write_doubles_3(arr,str,n)
        use utils_module, only: write_resume_unit,double_format
        implicit none
        real(dp),dimension(:,:,:), intent(in) :: arr
        integer, intent(in),dimension(size(arr,3)),optional :: n
        character(len=*), intent(in),optional :: str
        integer :: i

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        do i=1,size(arr,3)
            if(present(n)) then
                call write_doubles( arr(:,:,i), '---------------------------------------',n(i) )
            else
                call write_doubles( arr(:,:,i), '---------------------------------------' )
            end if
        end do

    end subroutine write_doubles_3



    subroutine write_resume_file(settings,RTI)
        use utils_module,    only: write_resume_unit
        use run_time_module, only: run_time_info
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        call check_directories(settings)

        ! Open the .resume file
        open(write_resume_unit,file=trim(resume_file(settings,.true.)), action='write') 

        ! write integers 
        call write_integer(settings%nDims,           "=== Number of dimensions ===")
        call write_integer(settings%nDerived,        "=== Number of derived parameters ===")                    
        call write_integer(RTI%ndead,                "=== Number of dead points/iterations ===")                    
        call write_integer(RTI%ncluster,             "=== Number of clusters ===")
        call write_integer(RTI%ncluster_dead,        "=== Number of dead clusters ===")
        call write_integer(RTI%nposterior_global,    "=== Number of global weighted posterior points ===")
        call write_integer(RTI%nequals_global,       "=== Number of global equally weighted posterior points ===")
        call write_integer(size(settings%grade_dims),"=== Number of grades ===")
        call write_integers(settings%grade_dims,     "=== positions of grades ===")
        call write_integers(RTI%num_repeats,         "=== Number of repeats ===")
        call write_integers(RTI%nlike,               "=== Number of likelihood calls ===")                    
        call write_integers(RTI%nlive,               "=== Number of live points in each cluster ===")
        call write_integers(RTI%nphantom,            "=== Number of phantom points in each cluster ===")                    
        call write_integers(RTI%nposterior,          "=== Number of weighted posterior points in each cluster ===")                    
        call write_integers(RTI%nequals,             "=== Number of equally weighted posterior points in each cluster ===")                    
        call write_integers(RTI%i,                   "=== Minimum loglikelihood positions ===")                    
        call write_integers(RTI%nposterior_dead,     "=== Number of weighted posterior points in each dead cluster ===")
        call write_integers(RTI%nequals_dead,        "=== Number of equally weighted posterior points in each dead cluster ===")
  
       ! write evidences
        call write_double(RTI%logZ,                  "=== global evidence -- log(<Z>) ===")                    
        call write_double(RTI%logZ2,                 "=== global evidence^2 -- log(<Z^2>) ===")                    
        call write_double(RTI%thin_posterior,        "=== posterior thin factor ===")                    
        call write_doubles(RTI%logLp,                "=== local loglikelihood bounds ===")                    
        call write_doubles(RTI%logXp,                "=== local volume -- log(<X_p>) ===")                    
        call write_double(RTI%logX_last_update,      "=== last update volume ===")                    
        call write_doubles(RTI%logZXp,               "=== global evidence volume cross correlation -- log(<ZX_p>) ===")                    
        call write_doubles(RTI%logZp,                "=== local evidence -- log(<Z_p>) ===")                    
        call write_doubles(RTI%logZp2,               "=== local evidence^2 -- log(<Z_p^2>) ===")                    
        call write_doubles(RTI%logZpXp,              "=== local evidence volume cross correlation -- log(<Z_pX_p>) ===")                    
        call write_doubles_2(RTI%logXpXq,            "=== local volume cross correlation -- log(<X_pX_q>) ===")                    
        call write_doubles(RTI%maxlogweight,         "=== maximum log weights -- log(w_p) ===")                    
  
        call write_doubles(RTI%logZp_dead,           "=== local dead evidence -- log(<Z_p>) ===")
        call write_doubles(RTI%logZp2_dead,          "=== local dead evidence^2 -- log(<Z_p^2>) ===")
        call write_doubles(RTI%maxlogweight_dead,    "=== maximum dead log weights -- log(w_p) ===")                    
  
        call write_doubles(RTI%covmat,             "=== covariance matrices ===")
        call write_doubles(RTI%cholesky,           "=== cholesky decompositions ===")
  
        call write_doubles(RTI%live,                  "=== live points ===",                               RTI%nlive )
        call write_doubles(RTI%dead,                  "=== dead points ===",                               RTI%ndead )
        call write_doubles(RTI%logweights(:RTI%ndead),"=== logweights of dead points ===")
        call write_doubles(RTI%phantom,               "=== phantom points ===",                            RTI%nphantom )
        call write_doubles(RTI%posterior,             "=== weighted posterior points ===",                 RTI%nposterior )
        call write_doubles(RTI%posterior_dead,        "=== dead weighted posterior points ===",            RTI%nposterior_dead )
        call write_doubles(RTI%posterior_global,      "=== global weighted posterior points ===",          RTI%nposterior_global )
        call write_doubles(RTI%equals,                "=== equally weighted posterior points ===",         RTI%nequals )
        call write_doubles(RTI%equals_dead,           "=== dead equally weighted posterior points ===",    RTI%nequals_dead )
        call write_doubles(RTI%equals_global,         "=== global equally weighted posterior points ===",  RTI%nequals_global )

        ! Close the writing file
        close(write_resume_unit)

    end subroutine write_resume_file






    subroutine read_integer(a,str)
        use utils_module, only: read_resume_unit,integer_format
        implicit none
        integer, intent(out) :: a
        character(len=*), intent(in),optional :: str

        if(present(str)) read(read_resume_unit,*)
        read(read_resume_unit,integer_format(1)) a

    end subroutine read_integer

    subroutine read_integers(arr,str,n)
        use utils_module, only: read_resume_unit,integer_format
        implicit none
        integer, allocatable, dimension(:), intent(out) :: arr
        integer, intent(in) :: n
        character(len=*), intent(in),optional :: str

        if(present(str)) read(read_resume_unit,*)
        allocate(arr(n))
        if(n>0) read(read_resume_unit,integer_format(n)) arr

    end subroutine read_integers

    subroutine read_double(a,str)
        use utils_module, only: read_resume_unit,double_format
        implicit none
        real(dp), intent(out) :: a
        character(len=*), intent(in),optional :: str

        if(present(str)) read(read_resume_unit,*)
        read(read_resume_unit,double_format(1)) a

    end subroutine read_double

    subroutine read_doubles_1(arr,str,n)
        use utils_module, only: read_resume_unit,double_format
        implicit none
        real(dp),allocatable,dimension(:), intent(out) :: arr
        integer,intent(in) :: n
        character(len=*), intent(in),optional :: str

        if(present(str)) read(read_resume_unit,*)
        allocate(arr(n))
        if(n>0) read(read_resume_unit,double_format(n)) arr

    end subroutine read_doubles_1

    subroutine read_doubles_2(arr,str,n1,n2)
        use utils_module, only: read_resume_unit,double_format
        implicit none
        real(dp),allocatable,dimension(:,:), intent(out) :: arr
        integer,intent(in) :: n1,n2
        character(len=*), intent(in),optional :: str
        integer :: i2

        if(present(str)) read(read_resume_unit,*)
        allocate(arr(n1,n2))
        do i2=1,n2
            read(read_resume_unit,double_format(n1)) arr(:,i2)
        end do

    end subroutine read_doubles_2

    subroutine read_doubles_3(arr,str,n1,n2,n3,n)
        use utils_module, only: read_resume_unit,double_format
        implicit none
        real(dp),allocatable,dimension(:,:,:), intent(out) :: arr
        integer,intent(in) :: n1,n2,n3
        integer,optional,intent(in),dimension(n3) :: n
        character(len=*), intent(in),optional :: str
        integer :: i2,i3,m

        if(present(str)) read(read_resume_unit,*)
        allocate(arr(n1,n2,n3))
        do i3=1,n3
            read(read_resume_unit,*)
            if(present(n)) then
                m=n(i3)
            else
                m=n2
            end if
            do i2=1,m
                read(read_resume_unit,double_format(n1)) arr(:,i2,i3)
            end do
        end do

    end subroutine read_doubles_3

    subroutine read_resume_file(settings,RTI)
        use utils_module,    only: read_resume_unit
        use run_time_module, only: run_time_info
        use settings_module, only: program_settings
        use abort_module,    only: halt_program
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(out) :: RTI

        integer :: i_temp
        integer,allocatable,dimension(:) :: i_temps
        integer :: ngrades

        ! Open the .resume file
        open(read_resume_unit,file=trim(resume_file(settings,.false.)), action='read') 

        call read_integer(i_temp,'-') ! number of dimensions
        if(settings%nDims/=i_temp) call halt_program('resume error: nDims does not match')
  
        call read_integer(i_temp,'-') ! number of derived parameters
        if(settings%nDerived/=i_temp) call halt_program('resume error: nDerived does not match')
  
        call read_integer(RTI%ndead,'-')             ! number of dead points
        call read_integer(RTI%ncluster,'-')          ! number of clusters
        call read_integer(RTI%ncluster_dead,'-')     ! number of dead clusters
        call read_integer(RTI%nposterior_global,'-') ! number of weighted posteriors 
        call read_integer(RTI%nequals_global,'-')    ! number of equally weighted posteriors 
  
        call read_integer(ngrades,'-') ! check the number of grades
        if(size(settings%grade_dims)/=ngrades) call halt_program('resume error: number of grades does not match')
  
        call read_integers(i_temps,'-',ngrades)    ! check the grades themselves
        if(any(settings%grade_dims/=i_temps)) call halt_program('resume error: Grades do not match') 
  
        call read_integers(RTI%num_repeats,'-',ngrades)  ! number of repeats per grade
        call read_integers(RTI%nlike,'-',ngrades)        ! number of likelihood calls per grade
  
        call read_integers(RTI%nlive,'-',RTI%ncluster)
        call read_integers(RTI%nphantom,'-',RTI%ncluster)
        call read_integers(RTI%nposterior,'-',RTI%ncluster)
        call read_integers(RTI%nequals,'-',RTI%ncluster)
        call read_integers(RTI%i,'-',RTI%ncluster)
  
        ! Check to see if this is consistent with settings
        !if(settings%nlive/=sum(RTI%nlive)) call halt_program('resume error: nlive does not match')
  
        call read_integers(RTI%nposterior_dead,'-',RTI%ncluster_dead)
        call read_integers(RTI%nequals_dead,'-',RTI%ncluster_dead)
        call read_double(RTI%logZ,'-')
        call read_double(RTI%logZ2,'-')
        call read_double(RTI%thin_posterior,'-')
        
        call read_doubles(RTI%logLp,'-',RTI%ncluster)
        call read_doubles(RTI%logXp,'-',RTI%ncluster)
        call read_double(RTI%logX_last_update,'-')
        call read_doubles(RTI%logZXp,'-',RTI%ncluster)
        call read_doubles(RTI%logZp,'-',RTI%ncluster)
        call read_doubles(RTI%logZp2,'-',RTI%ncluster)
        call read_doubles(RTI%logZpXp,'-',RTI%ncluster)
        call read_doubles(RTI%logXpXq,'-',RTI%ncluster,RTI%ncluster)
        call read_doubles(RTI%maxlogweight,'-',RTI%ncluster)
  
        call read_doubles(RTI%logZp_dead,'-',RTI%ncluster_dead)
        call read_doubles(RTI%logZp2_dead,'-',RTI%ncluster_dead)
        call read_doubles(RTI%maxlogweight_dead,'-',RTI%ncluster_dead)
  
        call read_doubles(RTI%covmat,'-',settings%nDims,settings%nDims,RTI%ncluster)
        call read_doubles(RTI%cholesky,'-',settings%nDims,settings%nDims,RTI%ncluster)
  
        call read_doubles(RTI%live,'-',settings%nTotal,maxval(RTI%nlive),RTI%ncluster,RTI%nlive)
        call read_doubles(RTI%dead,'-',settings%nTotal,RTI%ndead)
        call read_doubles(RTI%logweights,'-',RTI%ndead)
        call read_doubles(RTI%phantom,'-',settings%nTotal,maxval(RTI%nphantom),RTI%ncluster,RTI%nphantom)
  
        call read_doubles(RTI%posterior,'-',settings%nposterior,maxval(RTI%nposterior),RTI%ncluster,RTI%nposterior)
        call read_doubles(RTI%posterior_dead,'-',settings%nposterior,maxval(RTI%nposterior_dead),RTI%ncluster_dead,RTI%nposterior_dead)
        call read_doubles(RTI%posterior_global,'-',settings%nposterior,RTI%nposterior_global)
  
        call read_doubles(RTI%equals,'-',settings%np,maxval(RTI%nequals),RTI%ncluster,RTI%nequals)
        call read_doubles(RTI%equals_dead,'-',settings%np,maxval(RTI%nequals_dead),RTI%ncluster_dead,RTI%nequals_dead)
        call read_doubles(RTI%equals_global,'-',settings%np,RTI%nequals_global)

        ! Close the reading unit
        close(read_resume_unit)

        ! Allocate the posterior stack if we're calculating this
        allocate(RTI%posterior_stack(settings%nposterior,settings%nlive,RTI%ncluster),RTI%nposterior_stack(RTI%ncluster))
        RTI%nposterior_stack = 0 ! Initialise number of posterior points at 0

        RTI%maxlogweight_global = maxval(RTI%maxlogweight)

    end subroutine read_resume_file


    subroutine write_posterior_file(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_posterior_unit,write_equals_unit,sort_doubles
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in)    :: settings
        type(run_time_info),    intent(inout) :: RTI

        integer :: i_cluster
        integer :: i_post
        real(dp) :: weight

        integer, dimension(RTI%ncluster+RTI%ncluster_dead) :: ordering

        character(len=fmt_len) :: fmt_dbl

        call check_directories(settings)

        ! Assign the writing format
        write(fmt_dbl,'("(",I0,A,")")') settings%np,DB_FMT 

        ! Sort the indices into order
        ordering = sort_doubles([-RTI%logZp,-RTI%logZp_dead])

        ! ============= Equally weighted posteriors ================
        if(settings%equals) then
            ! ------------- global equally weighted posteriors ----------------

            ! Open the equally weighted global posterior file
            open(write_equals_unit,file=trim(posterior_file(settings,.false.,.true.)))

            ! Print out the posteriors
            do i_post=1,RTI%nequals_global
                write(write_equals_unit,fmt_dbl) 1d0,RTI%equals_global(settings%p_2l:,i_post)
            end do

            ! Close the equally weighted global posterior file
            close(write_equals_unit)

            ! ------------- cluster equally weighted posteriors ----------------

            if(settings%cluster_posteriors) then

                do i_cluster = 1,RTI%ncluster+RTI%ncluster_dead

                    if(ordering(i_cluster)<=RTI%ncluster) then
                        ! Open the equally weighted cluster posterior file
                        open(write_equals_unit,file=trim(posterior_file(settings,.false.,.true.,i_cluster)))

                        ! Print out the posterior for the active clusters
                        do i_post = 1,RTI%nequals(ordering(i_cluster))
                            write(write_equals_unit,fmt_dbl) exp(RTI%logZp(ordering(i_cluster))-RTI%logZ),RTI%equals(settings%p_2l:,i_post,ordering(i_cluster))
                        end do

                    else
                        ! Open the equally weighted cluster posterior file
                        open(write_equals_unit,file=trim(posterior_file(settings,.false.,.true.,i_cluster)))

                        ! Print out the posterior for the dead clusters
                        do i_post = 1,RTI%nequals_dead(ordering(i_cluster)-RTI%ncluster)
                            write(write_equals_unit,fmt_dbl) exp(RTI%logZp_dead(ordering(i_cluster)-RTI%ncluster)-RTI%logZ),RTI%equals_dead(settings%p_2l:,i_post,ordering(i_cluster)-RTI%ncluster)
                        end do
                    end if

                    ! Close the equally weighted cluster posterior file
                    close(write_equals_unit)

                end do

            end if
        end if

        ! ============= weighted posteriors ================
        if(settings%posteriors) then
            ! ------------- global weighted posteriors ----------------

            ! Open the weighted global posterior file
            open(write_posterior_unit,file=trim(posterior_file(settings,.true.,.true.)))

            ! Print out the posteriors
            do i_post=1,RTI%nposterior_global
                weight = exp(RTI%posterior_global(settings%pos_w,i_post) + RTI%posterior_global(settings%pos_l,i_post) - RTI%maxlogweight_global) 
                if( weight>0d0  ) write(write_posterior_unit,fmt_dbl) weight,-2*RTI%posterior_global(settings%pos_l,i_post),RTI%posterior_global(settings%pos_p0:,i_post) 
            end do

            ! Close the weighted global posterior file
            close(write_posterior_unit)

            ! ------------- cluster weighted posteriors ----------------
            if(settings%cluster_posteriors) then
                do i_cluster = 1,RTI%ncluster+RTI%ncluster_dead

                    if(ordering(i_cluster)<=RTI%ncluster) then

                        ! Open the weighted cluster posterior file
                        open(write_posterior_unit,file=trim(posterior_file(settings,.true.,.true.,i_cluster)))

                        ! Print out the posterior for the active clusters
                        do i_post = 1,RTI%nposterior(ordering(i_cluster))
                            weight = exp(RTI%posterior(settings%pos_w,i_post,ordering(i_cluster)) + RTI%posterior(settings%pos_l,i_post,ordering(i_cluster)) - RTI%maxlogweight(ordering(i_cluster)) +RTI%logZp(ordering(i_cluster))-RTI%logZ)
                            if( weight>0d0  ) write(write_posterior_unit,fmt_dbl) weight,-2*RTI%posterior(settings%pos_l,i_post,ordering(i_cluster)),RTI%posterior(settings%pos_p0:,i_post,ordering(i_cluster)) 
                        end do
                    else

                        ! Open the weighted cluster posterior file
                        open(write_posterior_unit,file=trim(posterior_file(settings,.true.,.true.,i_cluster+RTI%ncluster)))

                        ! Print out the posterior for the dead clusters
                        do i_post = 1,RTI%nposterior_dead(ordering(i_cluster)-RTI%ncluster)
                            weight = exp(RTI%posterior_dead(settings%pos_w,i_post,ordering(i_cluster)-RTI%ncluster) + RTI%posterior_dead(settings%pos_l,i_post,ordering(i_cluster)-RTI%ncluster) - RTI%maxlogweight_dead(ordering(i_cluster)-RTI%ncluster) +RTI%logZp_dead(ordering(i_cluster)-RTI%ncluster)-RTI%logZ)
                            if( weight>0d0  ) write(write_posterior_unit,fmt_dbl) weight,-2*RTI%posterior_dead(settings%pos_l,i_post,ordering(i_cluster)-RTI%ncluster),RTI%posterior_dead(settings%pos_p0:,i_post,ordering(i_cluster)-RTI%ncluster) 
                        end do

                    end if

                    ! Close the weighted cluster posterior file
                    close(write_posterior_unit)

                end do
            end if
        end if

        ! Rename the files
        if(settings%equals) then
            call rename(trim(posterior_file(settings,.false.,.true.)),trim(posterior_file(settings,.false.,.false.)))
            do i_cluster = 1,RTI%ncluster+RTI%ncluster_dead
                call rename(trim(posterior_file(settings,.false.,.true.,i_cluster)),trim(posterior_file(settings,.false.,.false.,i_cluster)))
            end do
        end if

        if(settings%posteriors) then
            call rename(trim(posterior_file(settings,.true. ,.true.)),trim(posterior_file(settings,.true. ,.false.)))
            do i_cluster = 1,RTI%ncluster+RTI%ncluster_dead
                call rename(trim(posterior_file(settings,.true.,.true.,i_cluster)),trim(posterior_file(settings,.true.,.false.,i_cluster)))
            end do
        end if

    end subroutine write_posterior_file



    subroutine write_phys_live_points(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_phys_unit,write_phys_cluster_unit, write_live_birth_unit
        use settings_module, only: program_settings 
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        integer i_live
        integer i_cluster

        character(len=fmt_len) :: fmt_dbl, fmt_dbl_1

        call check_directories(settings)

        ! Initialise the formats
        write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT
        write(fmt_dbl_1,'("(",I0,A,")")') settings%nDims+settings%nDerived+2, DB_FMT


        ! Open a new file for appending to
        open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')
        open(write_live_birth_unit,file=trim(phys_live_birth_file(settings)), action='write')

        do i_cluster = 1,RTI%ncluster

            if(settings%do_clustering) open(write_phys_cluster_unit,file=trim(phys_live_file(settings,i_cluster)), action='write') 

            do i_live=1,RTI%nlive(i_cluster)
                write(write_phys_unit,fmt_dbl) &
                    RTI%live(settings%p0:settings%d1,i_live,i_cluster), &
                    RTI%live(settings%l0,i_live,i_cluster)

                write(write_live_birth_unit,fmt_dbl_1) &
                    RTI%live(settings%p0:settings%d1,i_live,i_cluster), &
                    RTI%live(settings%l0,i_live,i_cluster), &
                    RTI%live(settings%b0,i_live,i_cluster)

                if(settings%do_clustering) then
                    write(write_phys_cluster_unit,fmt_dbl) &
                        RTI%live(settings%p0:settings%d1,i_live,i_cluster), &
                        RTI%live(settings%l0,i_live,i_cluster)
                end if

            end do

            if(settings%do_clustering) close(write_phys_cluster_unit)

        end do

        close(write_phys_unit)
        close(write_live_birth_unit)


    end subroutine write_phys_live_points


    subroutine write_dead_points(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_dead_unit, write_dead_birth_unit
        use settings_module, only: program_settings 
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        integer i_dead

        character(len=fmt_len) :: fmt_dbl

        call check_directories(settings)

        ! Initialise the formats
        write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT

        ! Open a new file for appending to
        open(write_dead_unit,file=trim(dead_file(settings)), action='write')

        do i_dead=1,RTI%ndead
            write(write_dead_unit,fmt_dbl) &
                RTI%dead(settings%l0,i_dead), &
                RTI%dead(settings%p0:settings%d1,i_dead)
        end do
        close(write_dead_unit)

        write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+2, DB_FMT

        open(write_dead_birth_unit,file=trim(dead_birth_file(settings)), action='write')
        do i_dead=1,RTI%ndead
            write(write_dead_birth_unit,fmt_dbl) &
                RTI%dead(settings%p0:settings%d1,i_dead), &
                RTI%dead(settings%l0,i_dead), &
                RTI%dead(settings%b0,i_dead)
        end do
        close(write_dead_birth_unit)


    end subroutine write_dead_points

    subroutine write_prior_file(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_prior_unit
        use settings_module, only: program_settings 
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        integer i_prior

        character(len=fmt_len) :: fmt_dbl

        call check_directories(settings)

        ! Initialise the formats
        write(fmt_dbl,'("(",I0,A,")")') settings%np,DB_FMT 

        ! Open a new file for appending to
        open(write_prior_unit,file=trim(prior_file(settings)), action='write')

        do i_prior=1,RTI%nlive(1)
            write(write_prior_unit,fmt_dbl) &
                1d0,&
                -2*RTI%live(settings%l0,i_prior,1), &
                RTI%live(settings%p0:settings%d1,i_prior,1)
        end do

        close(write_prior_unit)


    end subroutine write_prior_file

    subroutine write_max_file(settings,max_point, max_posterior_point, dXdtheta, mean_point)
        use utils_module, only: DB_FMT,fmt_len,write_max_unit
        use settings_module, only: program_settings 
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        real(dp), dimension(settings%nTotal), intent(in) :: max_point, max_posterior_point
        real(dp), dimension(settings%nTotal), intent(in), optional :: mean_point
        real(dp), intent(in) :: dXdtheta

        character(len=fmt_len) :: fmt_dbl

        call check_directories(settings)

        ! Initialise the formats
        write(fmt_dbl,'("(",I0,A,")")') settings%nDims + settings%nDerived + 1,DB_FMT 

        ! Open a new file for appending to
        open(write_max_unit,file=trim(max_file(settings)), action='write')
        write(write_max_unit, '("Maximum LogLikelihood:" )')
        write(fmt_dbl,'("(",I0,A,")")') 1,DB_FMT 
        write(write_max_unit, fmt_dbl) max_point(settings%l0)

        write(write_max_unit, '("Maximum Likelihood point:" )')
        write(fmt_dbl,'("(",I0,A,")")') settings%nDims + settings%nDerived,DB_FMT 
        write(write_max_unit, fmt_dbl) max_point(settings%p0:settings%d1)
        write(write_max_unit, '("")')


        write(write_max_unit, '("Maximum Posterior:" )')
        write(fmt_dbl,'("(",I0,A,")")') 1,DB_FMT 
        write(write_max_unit, fmt_dbl) max_posterior_point(settings%l0) + dXdtheta
        write(write_max_unit, '("Maximum Likelihood at posterior:" )')
        write(write_max_unit, fmt_dbl) max_posterior_point(settings%l0)

        write(write_max_unit, '("Maximum Posterior point:" )')
        write(fmt_dbl,'("(",I0,A,")")') settings%nDims + settings%nDerived,DB_FMT 
        write(write_max_unit, fmt_dbl) max_posterior_point(settings%p0:settings%d1)
        write(write_max_unit, '("")')

        if (present(mean_point)) then
            write(write_max_unit, '("LogLikelihood(mean):" )')
            write(fmt_dbl,'("(",I0,A,")")') 1,DB_FMT 
            write(write_max_unit, fmt_dbl) mean_point(settings%l0)

            write(write_max_unit, '("mean point:" )')
            write(fmt_dbl,'("(",I0,A,")")') settings%nDims + settings%nDerived,DB_FMT 
            write(write_max_unit, fmt_dbl) mean_point(settings%p0:settings%d1)
        end if

        close(write_max_unit)

    end subroutine write_max_file

    subroutine write_stats_file(settings,RTI,nlikesum)
        use utils_module, only: DB_FMT,fmt_len,write_stats_unit,logsubexp
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info,calculate_logZ_estimate
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI
        integer,dimension(:),   intent(in) :: nlikesum    !> number of likelihood calls since last call

        real(dp)                           :: logZ       
        real(dp)                           :: varlogZ  
        real(dp), dimension(RTI%ncluster)  :: logZp      
        real(dp), dimension(RTI%ncluster)  :: varlogZp 
        real(dp), dimension(RTI%ncluster_dead) :: logZp_dead      
        real(dp), dimension(RTI%ncluster_dead) :: varlogZp_dead 
        real(dp), dimension(settings%nDims+settings%nDerived) :: mu,sig
        real(dp)                           :: update_files

        integer :: p, i_dims, i_char

        character(len=fmt_len) :: fmt_Z,fmt_nlike

        call check_directories(settings)

        open(write_stats_unit,file=trim(stats_file(settings)), action='write') 

        call calculate_logZ_estimate(RTI,logZ,varlogZ,logZp,varlogZp,logZp_dead,varlogZp_dead)            


        write(fmt_Z,'("(""log(Z)       = "",", A, ","" +/- "",", A, ")")') DB_FMT,DB_FMT


        write(write_stats_unit, '("Evidence estimates:")')
        write(write_stats_unit, '("===================")')
        write(write_stats_unit, '("  - The evidence Z is a log-normally distributed, with location and scale parameters mu and sigma.")')
        write(write_stats_unit, '("  - We denote this as log(Z) = mu +/- sigma.")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Global evidence:")')
        write(write_stats_unit, '("----------------")')
        write(write_stats_unit,'("")')
        write(write_stats_unit,fmt_Z) logZ,sqrt(abs(varlogZ))
        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Local evidences:")')
        write(write_stats_unit, '("----------------")')
        write(write_stats_unit,'("")')


        do p=1,RTI%ncluster
            i_char = ceiling(log10(float(p+1)))
            write(fmt_Z,'("(""log(Z_"",I",I1,","")"",A",I1,",""= "",", A, ","" +/- "",", A, """ (Still Active)"")")') i_char,6-i_char,DB_FMT,DB_FMT
            write(write_stats_unit,fmt_Z) p, '', logZp(p), sqrt(abs(varlogZp(p)))
        end do
        do p=1,RTI%ncluster_dead
            i_char = ceiling(log10(float(p+RTI%ncluster+1)))
            write(fmt_Z,'("(""log(Z_"",I",I1,","")"",A",I1,",""= "",", A, ","" +/- "",", A, ")")') i_char,6-i_char,DB_FMT,DB_FMT
            write(write_stats_unit,fmt_Z) p+RTI%ncluster, '', logZp_dead(p), sqrt(abs(varlogZp_dead(p)))
        end do

        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Run-time information:")')
        write(write_stats_unit, '("---------------------")')
        write(write_stats_unit,'("")')
        write(write_stats_unit,'(" ncluster:   ", I8," /",I8           )') RTI%ncluster, RTI%ncluster+RTI%ncluster_dead
        write(write_stats_unit,'(" nposterior: ", I8                   )') RTI%nposterior_global
        write(write_stats_unit,'(" nequals:    ", I8                   )') RTI%nequals_global
        write(write_stats_unit,'(" ndead:      ", I8)') RTI%ndead
        write(write_stats_unit,'(" nlive:      ", I8)') sum(RTI%nlive)

        write(fmt_nlike,'("("" nlike:      "",",I0,"I8)")') size(nlikesum)
        write(write_stats_unit,fmt_nlike) RTI%nlike

        write(fmt_nlike,'(  "("" <nlike>:    "","  ,I0,   "F8.2,""   (""",I0,"F8.2 "" per slice )"")")') size(nlikesum), size(nlikesum)
        if (sum(RTI%nlive)>0) then
            update_files = -sum(RTI%nlive)*log(settings%compression_factor)
            write(write_stats_unit,fmt_nlike) dble(nlikesum)/update_files,dble(nlikesum)/dble(RTI%num_repeats*update_files)
        else
            write(write_stats_unit,fmt_nlike) [(0d0, p=1,size(nlikesum))], [(0d0, p=1,size(nlikesum))]
        end if

        if (settings%posteriors) then
            write(fmt_Z,'("(I3,", A, ","" +/- "",", A, ")")') DB_FMT,DB_FMT
            write(write_stats_unit,'("")')
            write(write_stats_unit,'("")')
            write(write_stats_unit, '("Dim No.       Mean        Sigma")')
            mu = mean(RTI,settings)
            sig = sqrt(variance(RTI,settings))
            do i_dims=1,settings%nDims
                write(write_stats_unit,fmt_Z) i_dims, mu(i_dims), sig(i_dims)
            end do
            write(write_stats_unit, '("-------------------------------")')
            do i_dims=settings%nDims+1,settings%nDims+settings%nDerived
                write(write_stats_unit,fmt_Z) i_dims, mu(i_dims), sig(i_dims)
            end do
        end if


        close(write_stats_unit)

    end subroutine write_stats_file

    function mean(RTI,settings) result(mu)
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info
        use utils_module,    only: logaddexp
        implicit none
        type(run_time_info),    intent(in) :: RTI
        type(program_settings), intent(in) :: settings

        integer i
        real(dp), dimension(settings%nDims+settings%nDerived) :: mu,x
        real(dp) :: logwsum, logw
        mu = 0
        logwsum = settings%logzero

        do i=1,RTI%nposterior_global
            x = RTI%posterior_global(settings%pos_p0:settings%pos_d1,i)
            logw = RTI%posterior_global(settings%pos_w,i) + RTI%posterior_global(settings%pos_l,i) 
            logwsum = logaddexp(logwsum,logw)
            mu = mu + exp(logw-logwsum)*(x-mu)
        end do
        

    end function mean

    function variance(RTI,settings) result(var)
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info
        use utils_module, only: logaddexp,logaddexp_
        implicit none
        type(run_time_info),    intent(in) :: RTI
        type(program_settings), intent(in) :: settings

        integer i
        real(dp), dimension(settings%nDims+settings%nDerived) :: mu,mu_old,x,var,logS
        real(dp) :: logw, logwsum
        mu = 0
        logwsum = settings%logzero
        logS = settings%logzero 

        do i=1,RTI%nposterior_global
            x = RTI%posterior_global(settings%pos_p0:settings%pos_d1,i)
            logw = RTI%posterior_global(settings%pos_w,i) + RTI%posterior_global(settings%pos_l,i) 
            mu_old = mu
            logwsum = logaddexp(logwsum,logw)
            mu = mu_old + exp(logw-logwsum)*(x-mu_old)
            if (all((x-mu_old)*(x-mu)>0)) logS = logaddexp_(logS,logw + log((x-mu_old)*(x-mu)))
        end do
        var = exp(logS - logwsum)

    end function variance


    subroutine write_paramnames_file(settings,params,derived_params)
        use priors_module, only: prior
        use settings_module,   only: program_settings
        use utils_module,  only: paramnames_unit
        use params_module, only: param_type
        implicit none
        
        type(program_settings),intent(in)                    :: settings       !> Program settings
        type(param_type),dimension(:),allocatable,intent(in) :: params         !> Parameter array
        type(param_type),dimension(:),allocatable,intent(in) :: derived_params !> Derived parameter array

        integer :: i

        call check_directories(settings)

        open(unit=paramnames_unit,file=trim(paramnames_file(settings)))

        do i=1,size(params)
            write(paramnames_unit,'(A,"      ",A)') trim(params(i)%paramname),trim(params(i)%latex) 
        end do

        do i=1,size(derived_params)
            write(paramnames_unit,'(A,"      ",A)') trim(derived_params(i)%paramname),trim(derived_params(i)%latex) 
        end do


        close(paramnames_unit)
        call write_properties_file(settings)
 
    end subroutine write_paramnames_file


    subroutine write_properties_file(settings)
        use priors_module, only: prior
        use settings_module,   only: program_settings
        use utils_module,  only: properties_unit
        use params_module, only: param_type
        implicit none
        
        type(program_settings),intent(in)                    :: settings       !> Program settings

        call check_directories(settings)

        open(unit=properties_unit,file=trim(properties_file(settings)))

        write(properties_unit,'("sampler=nested")')
        write(properties_unit,'("label=",A)') settings%file_root

        close(properties_unit)
 
    end subroutine write_properties_file




    ! File namers

    !> Name of the resume file
    function resume_file(settings,temp) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        !> Program settings (for base_dir and file_root)
        type(program_settings), intent(in) :: settings
        !> whether or not to create a temp file
        logical, intent(in) :: temp

        character(STR_LENGTH) :: file_name

        if(temp) then
            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_temp.resume'
        else
            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.resume'
        end if

    end function resume_file

    function stats_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.stats'

    end function stats_file

    function cluster_dir(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/clusters'

    end function cluster_dir

    function posterior_file(settings,weighted,temp,i) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        logical,intent(in) :: weighted
        logical,intent(in),optional :: temp !> whether to put a temporary prefix on it
        integer,intent(in),optional :: i

        character(STR_LENGTH) :: file_name

        character(STR_LENGTH) :: cluster_num

        if(present(i)) then
            write(cluster_num,'(I5)') i
            file_name = trim(cluster_dir(settings)) // '/' // trim(settings%file_root) // '_' // trim(adjustl(cluster_num)) 
        else 
            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root)
        end if

        if(present(temp)) then
            if(temp) file_name = trim(file_name) // '_temp'
        end if

        if(weighted) then
            file_name = trim(file_name) // '.txt'
        else
            file_name = trim(file_name) // '_equal_weights.txt'
        end if

    end function posterior_file







    function phys_live_file(settings,i) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in),optional :: i

        character(STR_LENGTH) :: file_name

        character(STR_LENGTH) :: cluster_num

        if(present(i)) then
            write(cluster_num,'(I5)') i
            file_name = trim(cluster_dir(settings)) // '/' // trim(settings%file_root) &
                // '_phys_live_' // trim(adjustl(cluster_num)) //'.txt'
        else 
            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_phys_live.txt'
        end if

    end function phys_live_file


    function phys_live_birth_file(settings,i) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in),optional :: i

        character(STR_LENGTH) :: file_name

        character(STR_LENGTH) :: cluster_num

        if(present(i)) then
            write(cluster_num,'(I5)') i
            file_name = trim(cluster_dir(settings)) // '/' // trim(settings%file_root) &
                // '_phys_live_' // trim(adjustl(cluster_num)) //'.txt'
        else 
            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_phys_live-birth.txt'
        end if

    end function phys_live_birth_file

    function prior_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_prior.txt'

    end function prior_file

    function max_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.maximum'

    end function max_file

    function dead_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_dead.txt'

    end function dead_file

    function dead_birth_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_dead-birth.txt'

    end function dead_birth_file

    function paramnames_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.paramnames'

    end function paramnames_file

    function properties_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.properties.ini'

    end function properties_file

    function prior_info_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.prior_info'

    end function prior_info_file

end module
