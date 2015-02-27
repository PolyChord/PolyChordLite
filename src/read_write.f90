!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

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
        deleted = delete_file( paramnames_file(settings), fb )     ! Delete paramnames file


        ! Delete normalised .txt files
        deleted = delete_file( posterior_file(settings), fb )
        deleted = delete_file( posterior_file(settings,.true.), fb )
        deleted = delete_file( posterior_file(settings,.false.), fb )

        i_cluster = 1
        do while ( &
                delete_file( posterior_file(settings,i=i_cluster), fb ) .or.  &
                delete_file( posterior_file(settings,.true.,i_cluster), fb ) .or.  &
                delete_file( posterior_file(settings,.false.,i_cluster), fb )   &
                )
            i_cluster = i_cluster + 1
        end do


    end subroutine delete_files


    subroutine write_resume_file(settings,RTI)
        use utils_module, only: DB_FMT,INT_FMT,fmt_len,write_resume_unit
        use run_time_module, only: run_time_info
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        integer :: i_cluster
        integer :: i_dims
        integer :: i_point


        character(len=fmt_len) :: fmt_int
        character(len=fmt_len) :: fmt_dbl



        ! Open the .resume file
        open(write_resume_unit,file=trim(resume_file(settings,.true.)), action='write') 

        ! write integers 
        write(fmt_int,'("(",I0,A,")")') 1,INT_FMT   ! define the integer format

        write(write_resume_unit,'("=== Number of dimensions ===")')                    
        write(write_resume_unit,fmt_int) settings%nDims

        write(write_resume_unit,'("=== Number of derived parameters ===")')                    
        write(write_resume_unit,fmt_int) settings%nDerived

        write(write_resume_unit,'("=== Number of dead points/iterations ===")')                    
        write(write_resume_unit,fmt_int) RTI%ndead
        write(write_resume_unit,'("=== Number of clusters ===")')                    
        write(write_resume_unit,fmt_int) RTI%ncluster

        ! Write out the grade information
        write(write_resume_unit,'("=== Number of grades ===")')                    
        write(write_resume_unit,fmt_int) size(settings%grade_dims)
        write(fmt_int,'("(",I0,A,")")') size(settings%grade_dims),INT_FMT   ! define the integer array format
        write(write_resume_unit,'("=== Total number of likelihood calls ===")')                    
        write(write_resume_unit,fmt_int) RTI%nlike
        write(write_resume_unit,'("=== positions of grades ===")')                    
        write(write_resume_unit,fmt_int) settings%grade_dims
        write(write_resume_unit,'("=== number of repeats ===")')                    
        write(write_resume_unit,fmt_int) RTI%num_repeats



        ! write number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%ncluster,INT_FMT   ! define the integer array format

        write(write_resume_unit,'("=== Number of live points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nlive           
        write(write_resume_unit,'("=== Number of phantom points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nphantom         
        write(write_resume_unit,'("=== Minimum loglikelihood positions ===")')                    
        write(write_resume_unit,fmt_int) RTI%i

        ! write evidences
        write(fmt_dbl,'("(",I0,A,")")') 1, DB_FMT              ! Initialise the double format
        write(write_resume_unit,'("=== global evidence -- log(<Z>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZ                
        write(write_resume_unit,'("=== global evidence^2 -- log(<Z^2>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZ2              



        write(fmt_dbl,'("(",I0,A,")")') RTI%ncluster, DB_FMT   ! Initialise the double array format
        write(write_resume_unit,'("=== local loglikelihood bounds ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logLp
        write(write_resume_unit,'("=== local volume -- log(<X_p>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logXp
        write(write_resume_unit,'("=== global evidence volume cross correlation -- log(<ZX_p>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZXp
        write(write_resume_unit,'("=== local evidence -- log(<Z_p>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZp
        write(write_resume_unit,'("=== local evidence^2 -- log(<Z_p^2>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZp2
        write(write_resume_unit,'("=== local evidence volume cross correlation -- log(<Z_pX_p>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZpXp
        write(write_resume_unit,'("=== local volume cross correlation -- log(<X_pX_q>) ===")')                    
        do i_cluster=1,RTI%ncluster
            write(write_resume_unit,fmt_dbl) RTI%logXpXq       
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nDims, DB_FMT   ! Initialise the double array format for matrices
        write(write_resume_unit,'("=== covariance matrices ===")') 
        do i_cluster=1,RTI%ncluster
            write(write_resume_unit,'("--- covariance matrix ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_dims=1,settings%nDims
                write(write_resume_unit,fmt_dbl) RTI%covmat(:,i_dims,i_cluster)
            end do
        end do
        write(write_resume_unit,'("=== cholesky decompositions ===")') 
        do i_cluster=1,RTI%ncluster
            write(write_resume_unit,'("--- cholesky decomposition ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_dims=1,settings%nDims
                write(write_resume_unit,fmt_dbl) RTI%cholesky(:,i_dims,i_cluster)
            end do
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nTotal, DB_FMT   ! Initialise the double array format for live points

        ! write live points
        write(write_resume_unit,'("=== live points ===")')                    
        do i_cluster=1,RTI%ncluster
            do i_point=1,RTI%nlive(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%live(:,i_point,i_cluster)
            end do
        end do

        ! write phantom points
        write(write_resume_unit,'("=== phantom points ===")')                    
        do i_cluster=1,RTI%ncluster
            do i_point=1,RTI%nphantom(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%phantom(:,i_point,i_cluster)
            end do
        end do

        ! Close the writing file
        close(write_resume_unit)

        ! Move the temp file onto the new file
        call rename(trim(resume_file(settings,.true.)),trim(resume_file(settings,.false.)))

    end subroutine write_resume_file


    subroutine read_resume_file(settings,RTI)
        use utils_module, only: DB_FMT,INT_FMT,fmt_len,read_resume_unit
        use run_time_module, only: run_time_info
        use settings_module, only: program_settings
        use feedback_module, only: write_resuming
        use abort_module, only: halt_program
        use array_module, only: add_point
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(out) :: RTI

        integer :: i_cluster
        integer :: i_dims
        integer :: i_point
        integer :: i_temp


        character(len=fmt_len) :: fmt_int
        character(len=fmt_len) :: fmt_dbl

        double precision, dimension(settings%nTotal) :: temp_point

        integer, allocatable,dimension(:) :: nlive
        integer, allocatable,dimension(:) :: nphantom

        ! -------------------------------------------- !
        call write_resuming(settings%feedback)
        ! -------------------------------------------- !


        ! Open the .resume file
        open(read_resume_unit,file=trim(resume_file(settings,.false.)), action='read') 

        ! Read in integers 
        write(fmt_int,'("(",I0,A,")")') 1,INT_FMT   ! define the integer format

        ! number of dimensions
        read(read_resume_unit,*)                    
        read(read_resume_unit,fmt_int) i_temp       
        if(settings%nDims/=i_temp) call halt_program('resume error: nDims does not match')

        ! number of derived parameters
        read(read_resume_unit,*)                    
        read(read_resume_unit,fmt_int) i_temp       
        if(settings%nDerived/=i_temp) call halt_program('resume error: nDerived does not match')

        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) RTI%ndead    ! number of dead points
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) RTI%ncluster ! number of clusters

        ! Allocate nlive and nphantom arrays based on these
        allocate(RTI%nlive(RTI%ncluster),RTI%nphantom(RTI%ncluster),RTI%nposterior(RTI%ncluster),RTI%i(RTI%ncluster),nphantom(RTI%ncluster),nlive(RTI%ncluster))
        RTI%nphantom=0
        RTI%nlive=0


        ! read in out the grade information
        read(read_resume_unit,*)                    
        read(read_resume_unit,fmt_int) i_temp
        if(size(settings%grade_dims)/=i_temp) call halt_program('resume error: Grades do not match')
        allocate(RTI%num_repeats(i_temp))
        write(fmt_int,'("(",I0,A,")")') size(settings%grade_dims),INT_FMT   ! define the integer array format
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) RTI%nlike    ! number of likelihood calls
        read(read_resume_unit,*)                    
        read(read_resume_unit,*)
        read(read_resume_unit,*)                    
        read(read_resume_unit,fmt_int) RTI%num_repeats

        ! Read in number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%ncluster,INT_FMT  ! define the integer array format

        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) nlive        ! temporary number of live points
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) nphantom     ! temporary number of phantom points
        read(read_resume_unit,*)                    !
        read(read_resume_unit,fmt_int) RTI%i        ! minimum loglikelihood positions


        ! Check to see if this is consistent with settings
        if(settings%nlive/=sum(nlive)) call halt_program('resume error: nlive does not match')


        ! Allocate the rest of the arrays
        allocate(                                                       &
            RTI%live(settings%nTotal,settings%nlive,RTI%ncluster),      &
            RTI%phantom(settings%nTotal,settings%nlive,RTI%ncluster),   &
            RTI%logLp(RTI%ncluster),                                    &
            RTI%logXp(RTI%ncluster),                                    &
            RTI%logZp(RTI%ncluster),                                    &
            RTI%logZXp(RTI%ncluster),                                   &
            RTI%logZp2(RTI%ncluster),                                   &
            RTI%logZpXp(RTI%ncluster),                                  &
            RTI%logXpXq(RTI%ncluster,RTI%ncluster),                     &
            RTI%covmat(settings%nDims,settings%nDims,RTI%ncluster),     &
            RTI%cholesky(settings%nDims,settings%nDims,RTI%ncluster)    &
            )


        ! Read in evidences
        write(fmt_dbl,'("(",I0,A,")")') 1, DB_FMT              ! Initialise the double format
        read(read_resume_unit,*)                               !
        read(read_resume_unit,fmt_dbl) RTI%logZ                ! global evidence estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZ2               ! global evidence^2 estimate



        write(fmt_dbl,'("(",I0,A,")")') RTI%ncluster, DB_FMT   ! Initialise the double array format
        read(read_resume_unit,*)                               !
        read(read_resume_unit,fmt_dbl) RTI%logLp               ! local loglikehood bound
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logXp               ! local volume estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZXp              ! global evidence volume cross correlation
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZp               ! local evidence estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZp2              ! local evidence^2 estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZpXp             ! local evidence volume cross correlation
        read(read_resume_unit,*)                               ! 
        do i_cluster=1,RTI%ncluster
            read(read_resume_unit,fmt_dbl) RTI%logXpXq         ! local volume cross correlation
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nDims, DB_FMT   ! Initialise the double array format for matrices
        read(read_resume_unit,*) 
        do i_cluster=1,RTI%ncluster
            read(read_resume_unit,*)
            do i_dims=1,settings%nDims
                read(read_resume_unit,fmt_dbl) RTI%covmat(:,i_dims,i_cluster)
            end do
        end do
        read(read_resume_unit,*) 
        do i_cluster=1,RTI%ncluster
            read(read_resume_unit,*)
            do i_dims=1,settings%nDims
                read(read_resume_unit,fmt_dbl) RTI%cholesky(:,i_dims,i_cluster)
            end do
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nTotal, DB_FMT   ! Initialise the double array format for live points

        ! Read in live points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster
            do i_point=1,nlive(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_point
                call add_point(temp_point,RTI%live,RTI%nlive,i_cluster)
            end do
        end do

        ! Read in phantom points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster
            do i_point=1,nphantom(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_point
                call add_point(temp_point,RTI%phantom,RTI%nphantom,i_cluster)
            end do
        end do

        ! Close the reading unit
        close(read_resume_unit)


        ! Allocate the posterior array if we're calculating this
        if(settings%calculate_posterior) allocate(RTI%posterior(settings%nposterior,settings%nlive,RTI%ncluster))
        RTI%nposterior = 0 ! Initialise number of posterior points at 0

    end subroutine read_resume_file


    subroutine write_posterior_file(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_untxt_unit,write_untxt_cluster_unit,write_equals_unit
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info 
        use random_module, only: bernoulli_trial
        use array_module, only: delete_point,add_point
        implicit none

        type(program_settings), intent(in)    :: settings
        type(run_time_info),    intent(inout) :: RTI

        integer :: i_cluster
        integer :: i_post

        character(len=fmt_len) :: fmt_dbl

        double precision, dimension(settings%np) :: posterior_point

        ! ============= Equally weighted posteriors ================
        open(write_equals_unit,file=trim(equals_file(settings)))

        write(fmt_dbl,'("(",I0,A,")")') settings%np,DB_FMT 

        ! strip out old posterior points
        i_post=1
        do while(i_post<=RTI%nequals(1)) 
            if(bernoulli_trial( exp(RTI%equals(settings%p_w,i_post,1)-RTI%maxlogweight) )) then
                RTI%equals(settings%p_w,i_post,1) = RTI%maxlogweight 
                write(write_equals_unit,fmt_dbl) 1d0,RTI%equals(settings%p_2l:,i_post,1)
                i_post=i_post+1
            else
                posterior_point = delete_point(i_post,RTI%equals,RTI%nequals,1)
            end if
        end do

        do i_cluster=1,RTI%ncluster
            do i_post=1,RTI%nposterior(i_cluster)
                if(bernoulli_trial( exp( RTI%posterior(settings%pos_w,i_post,i_cluster) + RTI%posterior(settings%pos_l,i_post,i_cluster) - RTI%maxlogweight) )) then
                    posterior_point(settings%p_w) = RTI%maxlogweight
                    posterior_point(settings%p_2l) = -2*RTI%posterior(settings%pos_l,i_post,i_cluster)
                    posterior_point(settings%p_p0:settings%p_d1) = RTI%posterior(settings%pos_p0:settings%pos_d1,i_post,i_cluster)
                    call add_point(posterior_point,RTI%equals,RTI%nequals,1)
                    write(write_equals_unit,fmt_dbl) 1d0,posterior_point(settings%p_2l:)
                end if
            end do
        end do

        close(write_equals_unit)
        write(*,*) RTI%maxlogweight
        


        ! ============= Unnormalised posteriors ====================

        if(settings%write_unnormalised_posterior) then

            ! Open the new unormalised .txt file for writing
            open(write_untxt_unit,file=trim(posterior_file(settings,.true.)), action='write',position='append') 

            ! Define the printing format
            write(fmt_dbl,'("(",I0,A,")")') settings%nposterior, DB_FMT

            do i_cluster=1,RTI%ncluster

                if(settings%do_clustering) open(write_untxt_cluster_unit,file=trim(posterior_file(settings,.true.,i_cluster)),action='write',position='append') 

                do i_post=1,RTI%nposterior(i_cluster)

                    ! Print each cluster sequentially to the main unnormalised chains file
                    write(write_untxt_unit,fmt_dbl) RTI%posterior(:,i_post,i_cluster)

                    ! If we're clustering, then print out separate cluster files
                    if(settings%do_clustering) write(write_untxt_cluster_unit,fmt_dbl) RTI%posterior(:,i_post,i_cluster) 

                end do

                if(settings%do_clustering) close(write_untxt_cluster_unit)

            end do
        end if

        ! Delete all of the posterior points
        RTI%nposterior = 0

        close(write_untxt_unit)

    end subroutine write_posterior_file



    subroutine write_phys_live_points(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_phys_unit,write_phys_cluster_unit
        use settings_module, only: program_settings 
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        integer i_live
        integer i_cluster

        character(len=fmt_len) :: fmt_dbl

        ! Initialise the formats
        write(fmt_dbl,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT


        ! Open a new file for appending to
        open(write_phys_unit,file=trim(phys_live_file(settings)), action='write')

        do i_cluster = 1,RTI%ncluster

            if(settings%do_clustering) open(write_phys_cluster_unit,file=trim(phys_live_file(settings,i_cluster)), action='write') 

            do i_live=1,RTI%nlive(i_cluster)
                write(write_phys_unit,fmt_dbl) &
                    RTI%live(settings%p0:settings%d1,i_live,i_cluster), &
                    RTI%live(settings%l0,i_live,i_cluster)

                if(settings%do_clustering) then
                    write(write_phys_cluster_unit,fmt_dbl) &
                        RTI%live(settings%p0:settings%d1,i_live,i_cluster), &
                        RTI%live(settings%l0,i_live,i_cluster)
                end if

            end do

            if(settings%do_clustering) close(write_phys_cluster_unit)

        end do

        close(write_phys_unit)


    end subroutine write_phys_live_points


    subroutine write_stats_file(settings,RTI)
        use utils_module, only: DB_FMT,fmt_len,write_stats_unit,logsubexp
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info,calculate_logZ_estimate
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI

        double precision                          :: logZ       
        double precision                          :: sigmalogZ  
        double precision, dimension(RTI%ncluster) :: logZp      
        double precision, dimension(RTI%ncluster) :: sigmalogZp 

        integer :: p

        character(len=fmt_len) :: fmt_Z

        open(write_stats_unit,file=trim(stats_file(settings)), action='write') 

        call calculate_logZ_estimate(RTI,logZ,sigmalogZ,logZp,sigmalogZp)            


        write(fmt_Z,'("(""log(Z)       = "",", A, ","" +/- "",", A, ")")') DB_FMT,DB_FMT


        write(write_stats_unit, '("Evidence estimates:")')
        write(write_stats_unit, '("===================")')
        write(write_stats_unit, '("  - The evidence Z is a log-normally distributed, with location and scale parameters mu and sigma.")')
        write(write_stats_unit, '("  - We denote this as log(Z) = mu +/- sigma.")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Global evidence:")')
        write(write_stats_unit, '("----------------")')
        write(write_stats_unit,'("")')
        write(write_stats_unit,fmt_Z) logZ,sigmalogZ
        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Local evidences:")')
        write(write_stats_unit, '("----------------")')
        write(write_stats_unit, '(" ", I2, " cluster still active")') RTI%ncluster
        write(write_stats_unit,'("")')


        write(fmt_Z,'("(""log(Z_"",",A,","")  = "",", A, ","" +/- "",", A, ")")') 'I2',DB_FMT,DB_FMT

        if(RTI%ncluster>1) then
            do p=1,RTI%ncluster
                write(write_stats_unit,fmt_Z) p, logZp(p), sigmalogZp(p)
            end do
        end if

        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Run-time information:")')
        write(write_stats_unit, '("---------------------")')
        write(write_stats_unit,'("")')
        write(write_stats_unit,'(" ndead:    ", I8)') RTI%ndead
        write(write_stats_unit,'(" nlive:    ", I8)') settings%nlive
        write(write_stats_unit,'(" active clusters: ", I8)') RTI%ncluster



        close(write_stats_unit)

    end subroutine write_stats_file


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

        open(unit=paramnames_unit,file=trim(paramnames_file(settings)))

        do i=1,size(params)
            write(paramnames_unit,'(A,"      ",A)') trim(params(i)%paramname),trim(params(i)%latex) 
        end do

        do i=1,size(derived_params)
            write(paramnames_unit,'(A,"      ",A)') trim(derived_params(i)%paramname),trim(derived_params(i)%latex) 
        end do


        close(paramnames_unit)
 
    end subroutine write_paramnames_file




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

    function posterior_file(settings,unnormalised,i) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        logical,intent(in),optional :: unnormalised
        integer,intent(in),optional :: i

        character(STR_LENGTH) :: file_name

        character(STR_LENGTH) :: cluster_num

        if(present(i)) then

            write(cluster_num,'(I5)') i
            file_name = trim(cluster_dir(settings)) // '/' // trim(settings%file_root) // '_'

            if(present(unnormalised)) then
                if(unnormalised) then
                    file_name = trim(file_name) // 'unnormalised_'
                else
                    file_name = trim(file_name) // 'unnormalised_temp_'
                endif
            end if

            file_name = trim(file_name) // trim(adjustl(cluster_num)) //'.txt'
        else 

            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root)

            if(present(unnormalised)) then
                if(unnormalised) then
                    file_name = trim(file_name) // '_unnormalised'
                else
                    file_name = trim(file_name) // '_unnormalised_temp' 
                endif
            end if

            file_name = trim(file_name) // '.txt'
        end if

    end function posterior_file

    function equals_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_eq.txt'

    end function equals_file

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

    function paramnames_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.paramnames'

    end function paramnames_file

end module
