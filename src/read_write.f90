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
    end subroutine


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
        write(write_resume_unit,'("=== Number of dead clusters ===")')
        write(write_resume_unit,fmt_int) RTI%ncluster_dead

        write(write_resume_unit,'("=== number of global weighted posterior points ===")')
        write(write_resume_unit,fmt_int) RTI%nposterior_global(1)
        write(write_resume_unit,'("=== number of global equally weighted posterior points ===")')
        write(write_resume_unit,fmt_int) RTI%nequals_global(1)

        ! Write out the grade information
        write(write_resume_unit,'("=== Number of grades ===")')
        write(write_resume_unit,fmt_int) size(settings%grade_dims)
        write(fmt_int,'("(",I0,A,")")') size(settings%grade_dims),INT_FMT   ! define the integer array format
        write(write_resume_unit,'("=== positions of grades ===")')
        write(write_resume_unit,fmt_int) settings%grade_dims
        write(write_resume_unit,'("=== number of repeats ===")')
        write(write_resume_unit,fmt_int) RTI%num_repeats
        write(write_resume_unit,'("=== Number of likelihood calls ===")')                    
        write(write_resume_unit,fmt_int) RTI%nlike


        ! write number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%ncluster,INT_FMT   ! define the integer array format

        write(write_resume_unit,'("=== Number of live points in each cluster ===")')
        write(write_resume_unit,fmt_int) RTI%nlive           
        write(write_resume_unit,'("=== Number of phantom points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nphantom         
        write(write_resume_unit,'("=== Number of weighted posterior points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nposterior
        write(write_resume_unit,'("=== Number of equally weighted posterior points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nequals
        write(write_resume_unit,'("=== Minimum loglikelihood positions ===")')                    
        write(write_resume_unit,fmt_int) RTI%i


        ! write number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%ncluster_dead,INT_FMT   ! define the integer array format

        write(write_resume_unit,'("=== Number of weighted posterior points in each dead cluster ===")')                    
        if(RTI%ncluster_dead>0) write(write_resume_unit,fmt_int) RTI%nposterior_dead
        write(write_resume_unit,'("=== Number of equally weighted posterior points in each dead cluster ===")')                    
        if(RTI%ncluster_dead>0) write(write_resume_unit,fmt_int) RTI%nequals_dead

        ! write evidences
        write(fmt_dbl,'("(",I0,A,")")') 1, DB_FMT              ! Initialise the double format
        write(write_resume_unit,'("=== global evidence -- log(<Z>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZ                
        write(write_resume_unit,'("=== global evidence^2 -- log(<Z^2>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZ2              
        write(write_resume_unit,'("=== posterior thin factor ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%thin_posterior


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

        write(write_resume_unit,'("=== maximum log weights -- log(w_p) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%maxlogweight

        if(RTI%ncluster_dead>0) write(fmt_dbl,'("(",I0,A,")")') RTI%ncluster_dead, DB_FMT   ! Initialise the double array format
        write(write_resume_unit,'("=== local dead evidence -- log(<Z_p>) ===")')
        if(RTI%ncluster_dead>0) write(write_resume_unit,fmt_dbl) RTI%logZp_dead
        write(write_resume_unit,'("=== local dead evidence^2 -- log(<Z_p^2>) ===")')
        if(RTI%ncluster_dead>0) write(write_resume_unit,fmt_dbl) RTI%logZp2_dead
        write(write_resume_unit,'("=== maximum dead log weights -- log(w_p) ===")')                    
        if(RTI%ncluster_dead>0) write(write_resume_unit,fmt_dbl) RTI%maxlogweight_dead


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
            write(write_resume_unit,'("--- live points ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_point=1,RTI%nlive(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%live(:,i_point,i_cluster)
            end do
        end do

        ! write phantom points
        write(write_resume_unit,'("=== phantom points ===")')                    
        do i_cluster=1,RTI%ncluster
            write(write_resume_unit,'("--- phantom points ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_point=1,RTI%nphantom(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%phantom(:,i_point,i_cluster)
            end do
        end do

        write(write_resume_unit,'("=== weighted posterior points ===")')                    
        write(fmt_dbl,'("(",I0,A,")")') settings%nposterior, DB_FMT   ! Initialise the double array format for weighted posterior points
        do i_cluster=1,RTI%ncluster
            write(write_resume_unit,'("--- weighted posterior points ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_point=1,RTI%nposterior(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%posterior(:,i_point,i_cluster)
            end do
        end do

        write(write_resume_unit,'("=== dead weighted posterior points ===")')                    
        do i_cluster=1,RTI%ncluster_dead
            write(write_resume_unit,'("--- dead weighted posterior points ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_point=1,RTI%nposterior_dead(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%posterior_dead(:,i_point,i_cluster)
            end do
        end do

        write(write_resume_unit,'("=== global weighted posterior points ===")')                    
        do i_point=1,RTI%nposterior_global(1)
            write(write_resume_unit,fmt_dbl) RTI%posterior_global(:,i_point,1)
        end do

        write(write_resume_unit,'("=== equally weighted posterior points ===")')                    
        write(fmt_dbl,'("(",I0,A,")")') settings%np, DB_FMT   ! Initialise the double array format for weighted posterior points
        do i_cluster=1,RTI%ncluster
            write(write_resume_unit,'("--- equally weighted posterior points ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_point=1,RTI%nequals(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%equals(:,i_point,i_cluster)
            end do
        end do

        write(write_resume_unit,'("=== dead equally weighted posterior points ===")')                    
        do i_cluster=1,RTI%ncluster_dead
            write(write_resume_unit,'("--- dead equally weighted posterior points ",I3,"/",I3," ---")') i_cluster, RTI%ncluster
            do i_point=1,RTI%nequals_dead(i_cluster)
                write(write_resume_unit,fmt_dbl) RTI%equals_dead(:,i_point,i_cluster)
            end do
        end do

        write(write_resume_unit,'("=== global equally weighted posterior points ===")')                    
        do i_point=1,RTI%nequals_global(1)
            write(write_resume_unit,fmt_dbl) RTI%equals_global(:,i_point,1)
        end do
        ! Close the writing file
        close(write_resume_unit)

    end subroutine write_resume_file


    subroutine read_resume_file(settings,RTI)
        use utils_module, only: DB_FMT,INT_FMT,fmt_len,read_resume_unit
        use run_time_module, only: run_time_info
        use settings_module, only: program_settings
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
        double precision, dimension(settings%nposterior) :: temp_posterior
        double precision, dimension(settings%np) :: temp_equal

        integer, allocatable,dimension(:) :: nlive
        integer, allocatable,dimension(:) :: nphantom
        integer, allocatable,dimension(:) :: nposterior
        integer, allocatable,dimension(:) :: nposterior_dead
        integer, allocatable,dimension(:) :: nposterior_global
        integer, allocatable,dimension(:) :: nequals
        integer, allocatable,dimension(:) :: nequals_dead
        integer, allocatable,dimension(:) :: nequals_global


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

        read(read_resume_unit,*)                         ! 
        read(read_resume_unit,fmt_int) RTI%ndead         ! number of dead points
        read(read_resume_unit,*)                         ! 
        read(read_resume_unit,fmt_int) RTI%ncluster      ! number of clusters
        read(read_resume_unit,*)                         ! 
        read(read_resume_unit,fmt_int) RTI%ncluster_dead ! number of dead clusters

        ! Allocate nlive and nphantom arrays based on these
        allocate(RTI%nlive(RTI%ncluster),RTI%nphantom(RTI%ncluster),RTI%nposterior(RTI%ncluster),RTI%nequals(RTI%ncluster),RTI%i(RTI%ncluster),nlive(RTI%ncluster),nphantom(RTI%ncluster),nposterior(RTI%ncluster),nequals(RTI%ncluster),RTI%nposterior_dead(RTI%ncluster_dead),RTI%nequals_dead(RTI%ncluster_dead),nposterior_dead(RTI%ncluster_dead),nequals_dead(RTI%ncluster_dead),RTI%nposterior_stack(RTI%ncluster),RTI%num_repeats(size(settings%grade_dims)),RTI%nequals_global(1),nequals_global(1),RTI%nposterior_global(1),nposterior_global(1))
        RTI%nphantom=0
        RTI%nlive=0
        RTI%nposterior=0
        RTI%nposterior_dead=0
        RTI%nposterior_global=0
        RTI%nequals=0
        RTI%nequals_dead=0
        RTI%nequals_global=0
        RTI%nposterior_stack=0

        read(read_resume_unit,*)                         ! 
        read(read_resume_unit,fmt_int) nposterior_global ! number of weighted posteriors
        read(read_resume_unit,*)                         ! 
        read(read_resume_unit,fmt_int) nequals_global    ! number of equally weighted posteriors

        ! read in out the grade information
        read(read_resume_unit,*)                    
        read(read_resume_unit,fmt_int) i_temp
        if(size(settings%grade_dims)/=i_temp) call halt_program('resume error: Grades do not match')
        allocate(RTI%nlike(i_temp))
        write(fmt_int,'("(",I0,A,")")') size(settings%grade_dims),INT_FMT   ! define the integer array format
        read(read_resume_unit,*)                    
        read(read_resume_unit,*)
        read(read_resume_unit,*)                    
        read(read_resume_unit,fmt_int) RTI%num_repeats
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) RTI%nlike    ! number of likelihood calls

        ! Read in number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%ncluster,INT_FMT  ! define the integer array format

        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) nlive        ! temporary number of live points
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) nphantom     ! temporary number of phantom points
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) nposterior   ! temporary number of weighted posteriors
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) nequals      ! temporary number of equally weighted posterors
        read(read_resume_unit,*)                    !
        read(read_resume_unit,fmt_int) RTI%i        ! minimum loglikelihood positions

        write(fmt_int,'("(",I0,A,")")') RTI%ncluster_dead,INT_FMT  ! define the integer array format
        read(read_resume_unit,*)                         ! 
        if(RTI%ncluster_dead>0) read(read_resume_unit,fmt_int) nposterior_dead   ! temporary number of weighted posteriors
        read(read_resume_unit,*)                         ! 
        if(RTI%ncluster_dead>0) read(read_resume_unit,fmt_int) nequals_dead      ! temporary number of equally weighted posterors


        ! Check to see if this is consistent with settings
        if(settings%nlive/=sum(nlive)) call halt_program('resume error: nlive does not match')


        ! Allocate the rest of the arrays
        allocate(                                                                     &
            RTI%live(settings%nTotal,settings%nlive,RTI%ncluster),                    &
            RTI%phantom(settings%nTotal,settings%nlive,RTI%ncluster),                 &
            RTI%posterior(settings%nposterior,settings%nlive,RTI%ncluster),           &
            RTI%posterior_dead(settings%nposterior,settings%nlive,RTI%ncluster_dead), &
            RTI%posterior_global(settings%nposterior,settings%nlive,1),               &
            RTI%equals(settings%np,settings%nlive,RTI%ncluster),                      &
            RTI%equals_dead(settings%np,settings%nlive,RTI%ncluster_dead),            &
            RTI%equals_global(settings%np,settings%nlive,1),                          &
            RTI%logLp(RTI%ncluster),                                                  &
            RTI%logXp(RTI%ncluster),                                                  &
            RTI%logZp(RTI%ncluster),                                                  &
            RTI%logZp_dead(RTI%ncluster_dead),                                        &
            RTI%logZXp(RTI%ncluster),                                                 &
            RTI%logZp2(RTI%ncluster),                                                 &
            RTI%logZp2_dead(RTI%ncluster_dead),                                       &
            RTI%logZpXp(RTI%ncluster),                                                &
            RTI%logXpXq(RTI%ncluster,RTI%ncluster),                                   &
            RTI%covmat(settings%nDims,settings%nDims,RTI%ncluster),                   &
            RTI%cholesky(settings%nDims,settings%nDims,RTI%ncluster),                 &
            RTI%maxlogweight(RTI%ncluster),                                           &
            RTI%maxlogweight_dead(RTI%ncluster_dead)                                  &
            )


        ! Read in evidences
        write(fmt_dbl,'("(",I0,A,")")') 1, DB_FMT              ! Initialise the double format
        read(read_resume_unit,*)                               !
        read(read_resume_unit,fmt_dbl) RTI%logZ                ! global evidence estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZ2               ! global evidence^2 estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%thin_posterior      ! what to thin the posterior by



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
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%maxlogweight        ! max log weights

        write(fmt_dbl,'("(",I0,A,")")') RTI%ncluster_dead, DB_FMT   ! Initialise the double array format
        read(read_resume_unit,*)                                    ! 
        if(RTI%ncluster_dead>0) read(read_resume_unit,fmt_dbl) RTI%logZp_dead               ! local dead evidence estimate
        read(read_resume_unit,*)                                    ! 
        if(RTI%ncluster_dead>0) read(read_resume_unit,fmt_dbl) RTI%logZp2_dead              ! local dead evidence^2 estimate
        read(read_resume_unit,*)                               ! 
        if(RTI%ncluster_dead>0) read(read_resume_unit,fmt_dbl) RTI%maxlogweight_dead   ! max dead log weights


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
            read(read_resume_unit,*)                               
            do i_point=1,nlive(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_point
                call add_point(temp_point,RTI%live,RTI%nlive,i_cluster)
            end do
        end do

        ! Read in phantom points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster
            read(read_resume_unit,*)                               
            do i_point=1,nphantom(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_point
                call add_point(temp_point,RTI%phantom,RTI%nphantom,i_cluster)
            end do
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nposterior, DB_FMT   ! Initialise the double array format for weighted posterior points

        ! Read in weighted posterior points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster
            read(read_resume_unit,*)                               
            do i_point=1,nposterior(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_posterior
                call add_point(temp_posterior,RTI%posterior,RTI%nposterior,i_cluster)
            end do
        end do

        ! read in dead weighted posterior points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster_dead
            read(read_resume_unit,*)                               
            do i_point=1,nposterior_dead(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_posterior
                call add_point(temp_posterior,RTI%posterior_dead,RTI%nposterior_dead,i_cluster)
            end do
        end do

        ! Read in global weighted posterior points
        read(read_resume_unit,*)                               
        do i_point=1,nposterior_global(1)
            read(read_resume_unit,fmt_dbl) temp_posterior
            call add_point(temp_posterior,RTI%posterior_global,RTI%nposterior_global,1)
        end do



        write(fmt_dbl,'("(",I0,A,")")') settings%np, DB_FMT   ! Initialise the double array format for equally weighted posterior points

        ! Read in equally weighted posterior points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster
            read(read_resume_unit,*)                               
            do i_point=1,nequals(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_equal
                call add_point(temp_equal,RTI%equals,RTI%nequals,i_cluster)
            end do
        end do

        ! read in dead equally weighted posterior points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%ncluster_dead
            read(read_resume_unit,*)                               
            do i_point=1,nequals_dead(i_cluster)
                read(read_resume_unit,fmt_dbl) temp_equal
                call add_point(temp_equal,RTI%equals_dead,RTI%nequals_dead,i_cluster)
            end do
        end do

        ! Read in global equally weighted posterior points
        read(read_resume_unit,*)                               
        do i_point=1,nequals_global(1)
            read(read_resume_unit,fmt_dbl) temp_equal
            call add_point(temp_equal,RTI%equals_global,RTI%nequals_global,1)
        end do

        ! Close the reading unit
        close(read_resume_unit)


        ! Allocate the posterior stack if we're calculating this
        allocate(RTI%posterior_stack(settings%nposterior,settings%nlive,RTI%ncluster))
        RTI%nposterior_stack = 0 ! Initialise number of posterior points at 0

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
        double precision :: weight

        integer, dimension(RTI%ncluster+RTI%ncluster_dead) :: ordering

        character(len=fmt_len) :: fmt_dbl

        ! Assign the writing format
        write(fmt_dbl,'("(",I0,A,")")') settings%np,DB_FMT 

        ! ============= Equally weighted posteriors ================
        if(settings%equals) then
            ! ------------- global equally weighted posteriors ----------------

            ! Open the equally weighted global posterior file
            open(write_equals_unit,file=trim(posterior_file(settings,.false.,.true.)))

            ! Print out the posteriors
            do i_post=1,RTI%nequals_global(1)
                write(write_equals_unit,fmt_dbl) 1d0,RTI%equals_global(settings%p_2l:,i_post,1)
            end do

            ! Close the equally weighted global posterior file
            close(write_equals_unit)


            ! Sort the indices into order
            ordering = sort_doubles([-RTI%logZp,-RTI%logZp_dead])

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
            do i_post=1,RTI%nposterior_global(1)
                weight = exp(RTI%posterior_global(settings%pos_w,i_post,1) + RTI%posterior_global(settings%pos_l,i_post,1) - RTI%maxlogweight_global) 
                if( weight>0d0  ) write(write_posterior_unit,fmt_dbl) weight,-2*RTI%posterior_global(settings%pos_l,i_post,1),RTI%posterior_global(settings%pos_p0:,i_post,1) 
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


    subroutine write_stats_file(settings,RTI,nlikesum)
        use utils_module, only: DB_FMT,fmt_len,write_stats_unit,logsubexp
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info,calculate_logZ_estimate
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: RTI
        integer,dimension(:),   intent(in) :: nlikesum    !> number of likelihood calls since last call

        double precision                           :: logZ       
        double precision                           :: varlogZ  
        double precision, dimension(RTI%ncluster) :: logZp      
        double precision, dimension(RTI%ncluster) :: varlogZp 
        double precision, dimension(RTI%ncluster_dead) :: logZp_dead      
        double precision, dimension(RTI%ncluster_dead) :: varlogZp_dead 

        integer :: p

        character(len=fmt_len) :: fmt_Z,fmt_nlike

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

        write(fmt_Z,'("(""log(Z_"",",A,","")  = "",", A, ","" +/- "",", A, """ (Still Active)"")")') 'I2',DB_FMT,DB_FMT
        do p=1,RTI%ncluster
            write(write_stats_unit,fmt_Z) p, logZp(p), sqrt(abs(varlogZp(p)))
        end do
        write(fmt_Z,'("(""log(Z_"",",A,","")  = "",", A, ","" +/- "",", A, ")")') 'I2',DB_FMT,DB_FMT
        do p=1,RTI%ncluster_dead
            write(write_stats_unit,fmt_Z) p+RTI%ncluster, logZp_dead(p), sqrt(abs(varlogZp_dead(p)))
        end do

        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Run-time information:")')
        write(write_stats_unit, '("---------------------")')
        write(write_stats_unit,'("")')
        write(write_stats_unit,'(" ncluster:   ", I8," /",I8           )') RTI%ncluster, RTI%ncluster+RTI%ncluster_dead
        write(write_stats_unit,'(" nposterior: ", I8                   )') RTI%nposterior_global(1)
        write(write_stats_unit,'(" nequals:    ", I8                   )') RTI%nequals_global(1)
        write(write_stats_unit,'(" ndead:      ", I8)') RTI%ndead
        write(write_stats_unit,'(" nlive:      ", I8)') settings%nlive

        write(fmt_nlike,'("("" nlike:      "",",I0,"I8)")') size(nlikesum)
        write(write_stats_unit,fmt_nlike) RTI%nlike

        write(fmt_nlike,'(  "("" <nlike>:    "","  ,I0,   "F8.2,""   (""",I0,"F8.2 "" per slice )"")")') size(nlikesum), size(nlikesum)
        write(write_stats_unit,fmt_nlike) dble(nlikesum)/dble(settings%update_files),dble(nlikesum)/dble(RTI%num_repeats*settings%update_files)


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

    function paramnames_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.paramnames'

    end function paramnames_file

end module
