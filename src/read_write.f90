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
        implicit none

        !> settings variable to get the base_dir and root_name out of
        type(program_settings), intent(in) :: settings

        integer :: i_cluster ! cluster iterator

        logical :: deleted

        deleted = delete_file( stats_file(settings) )           ! Delete stats file
        deleted = delete_file( phys_live_file(settings) )       ! Delete phys_live file
        deleted = delete_file( resume_file(settings,.false.) )  ! Delete temp resume file
        deleted = delete_file( resume_file(settings,.true.) )   ! Delete resume file


        ! Delete normalised .txt files
        deleted = delete_file( posterior_file(settings) )

        i_cluster = 1
        do while ( delete_file( posterior_file(settings,i_cluster) ) )
            i_cluster = i_cluster + 1
        end do

        ! Delete the unnormalised .txt files
        deleted = delete_file( posterior_file(settings,.true.) )

        i_cluster = 1
        do while ( delete_file( posterior_file(settings,.true.,i_cluster) ) )
            i_cluster = i_cluster + 1
        end do


    end subroutine delete_files


    subroutine write_resume_file(settings,RTI)
        use utils_module, only: DB_FMT,INT_FMT,fmt_len,write_resume_unit
        use run_time_module, only: run_time_info
        use settings_module, only: program_settings

        implicit none
        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(out) :: RTI

        integer :: i_cluster


        character(len=fmt_len) :: fmt_int
        character(len=fmt_len) :: fmt_dbl



        ! Open the .resume file
        open(write_resume_unit,file=trim(resume_file(settings,.true.)), action='write') 

        ! write integers 
        write(fmt_int,'("(",I0,A,")")') 1,INT_FMT   ! define the integer format

        ! number of dimensions
        write(write_resume_unit,'("=== Number of dimensions ===")')                    
        write(write_resume_unit,fmt_int) settings%nDims

        ! number of derived parameters
        write(write_resume_unit,'("=== Number of derived parameters ===")')                    
        write(write_resume_unit,fmt_int) settings%nDerived

        write(write_resume_unit,'("=== Number of dead points/iterations ===")')                    
        write(write_resume_unit,fmt_int) RTI%ndead
        write(write_resume_unit,'("=== Total number of likelihood calls ===")')                    
        write(write_resume_unit,fmt_int) RTI%nlike
        write(write_resume_unit,'("=== Number of clusters ===")')                    
        write(write_resume_unit,fmt_int) RTI%nclusters


        ! write number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%nclusters,INT_FMT   ! define the integer array format

        write(write_resume_unit,'("=== Number of live points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nlive               ! number of live points
        write(write_resume_unit,'("=== Number of phantom points in each cluster ===")')                    
        write(write_resume_unit,fmt_int) RTI%nphantom            ! number of phantom points

        ! write evidences
        write(fmt_dbl,'("(",I0,A,")")') 1, DB_FMT              ! Initialise the double format
        write(write_resume_unit,'("=== global evidence -- log(<Z>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZ                ! global evidence estimate
        write(write_resume_unit,'("=== global evidence^2 -- log(<Z^2>) ===")')                    
        write(write_resume_unit,fmt_dbl) RTI%logZ2               ! global evidence^2 estimate



        write(fmt_dbl,'("(",I0,A,")")') RTI%nclusters, DB_FMT   ! Initialise the double array format
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
        do i_cluster=1,RTI%nclusters
            write(write_resume_unit,fmt_dbl) RTI%logXpXq         ! local volume cross correlation
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nTotal, DB_FMT   ! Initialise the double array format for live points

        ! write live points
        write(write_resume_unit,'("=== live points ===")')                    
        do i_cluster=1,RTI%nclusters
            write(write_resume_unit,fmt_dbl) RTI%live(:,:,i_cluster)
        end do

        ! write phantom points
        write(write_resume_unit,'("=== phantom points ===")')                    
        do i_cluster=1,RTI%nclusters
            write(write_resume_unit,fmt_dbl) RTI%phantom(:,:,i_cluster)
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
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(out) :: RTI

        integer :: i_cluster
        integer :: i_temp


        character(len=fmt_len) :: fmt_int
        character(len=fmt_len) :: fmt_dbl


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
        read(read_resume_unit,fmt_int) RTI%nlike    ! number of likelihood calls
        read(read_resume_unit,*)                    ! 
        read(read_resume_unit,fmt_int) RTI%nclusters! number of clusters

        ! Allocate nlive and nphantom arrays based on these
        allocate(                           &
            RTI%nlive(RTI%nclusters),       &
            RTI%nphantom(RTI%nclusters),    &
            RTI%nposterior(RTI%nposterior)  &
            )

        ! Read in number of live and phantom points
        write(fmt_int,'("(",I0,A,")")') RTI%nclusters,INT_FMT  ! define the integer array format

        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_int) RTI%nlive               ! number of live points
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_int) RTI%nphantom            ! number of phantom points


        ! Check to see if this is consistent with settings
        if(settings%nlive/=sum(RTI%nlive)) call halt_program('resume error: nlive does not match')


        ! Allocate the rest of the arrays
        allocate(                                                        &
            RTI%live(settings%nTotal,settings%nlive,RTI%nclusters),      &
            RTI%phantom(settings%nTotal,settings%nlive,RTI%nclusters),   &
            RTI%logXp(RTI%nclusters),                                    &
            RTI%logZXp(RTI%nclusters),                                   &
            RTI%logZp(RTI%nclusters),                                    &
            RTI%logZp2(RTI%nclusters),                                   &
            RTI%logZpXp(RTI%nclusters),                                  &
            RTI%logXpXq(RTI%nclusters,RTI%nclusters),                    &
            )


        ! Read in evidences
        write(fmt_dbl,'("(",I0,A,")")') 1, DB_FMT              ! Initialise the double format
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZ                ! global evidence estimate
        read(read_resume_unit,*)                               ! 
        read(read_resume_unit,fmt_dbl) RTI%logZ2               ! global evidence^2 estimate



        write(fmt_dbl,'("(",I0,A,")")') RTI%nclusters, DB_FMT  ! Initialise the double array format
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
        do i_cluster=1,RTI%nclusters
            read(read_resume_unit,fmt_dbl) RTI%logXpXq         ! local volume cross correlation
        end do


        write(fmt_dbl,'("(",I0,A,")")') settings%nTotal, DB_FMT   ! Initialise the double array format for live points

        ! Read in live points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%nclusters
            read(read_resume_unit,fmt_dbl) RTI%live(:,:,i_cluster)
        end do

        ! Read in phantom points
        read(read_resume_unit,*)                               
        do i_cluster=1,RTI%nclusters
            read(read_resume_unit,fmt_dbl) RTI%phantom(:,:,i_cluster)
        end do

        ! Close the reading unit
        close(read_resume_unit)


        ! Allocate the posterior array if we're calculating this
        if(settings%calculate_posterior) allocate(RTI%posterior(settings%nTotal,settings%nlive,RTI%nclusters))
        RTI%nposterior = 0 ! Initialise number of posterior points at 0

        ! Allocate the cholesky and covmat arrays
        allocate(                                                      &
            RTI%cholesky(settings%nDims,settings%nDims,RTI%nclusters), &
            RTI%covmat(settings%nDims,settings%nDims,RTI%nclusters)    &
            )

    end subroutine read_resume_file



    subroutine write_posterior_file(settings,info,posterior_points,nposterior) 
        use utils_module, only: DB_FMT,fmt_len,write_txt_unit,write_untxt_unit,read_untxt_unit,logsigma
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: info
        double precision,intent(inout), dimension(settings%nposterior,settings%nstack,0:settings%ncluster*2) :: posterior_points
        integer,intent(inout),dimension(0:settings%ncluster*2) :: nposterior

        integer :: i_posterior
        integer :: i_cluster

        double precision :: logweight

        double precision, dimension(settings%nposterior) :: posterior_point


        integer :: io_stat

        character(len=fmt_len) :: fmt_dbl_nposterior
        character(len=fmt_len) :: fmt_dbl_nposterior_norm


        ! Initialise the formats
        write(fmt_dbl_nposterior,'("(",I0,A,")")') settings%nposterior, DB_FMT
        write(fmt_dbl_nposterior_norm,'("(",I0,A,")")') settings%nposterior-1, DB_FMT
        
        ! Open the .txt file for writing
        open(write_txt_unit,file=trim(posterior_file(settings)), action='write') 
        ! Open the old unormalised .txt file for reading
        open(read_untxt_unit,file=trim(posterior_file(settings,.true.)), action='read',iostat=io_stat) 
        ! Open the new unormalised .txt file for writing
        open(write_untxt_unit,file=trim(posterior_file(settings,.false.)), action='write') 



        ! Trim off points that are too low in weight
        do while(io_stat==0) 

            ! Read a point from the old unormalised text file
            read(read_untxt_unit,fmt_dbl_nposterior,iostat=io_stat) posterior_point

            ! If this is still above the cdf threshold
            if( posterior_point(settings%pos_Z) - info%logevidence > logsigma(settings%sigma_posterior) ) then

                ! Re-copy to the new unnormalised posterior file
                write(write_untxt_unit,fmt_dbl_nposterior) posterior_point

                ! And add it to the .txt file
                logweight = posterior_point(settings%pos_w) - info%logevidence 

                if(logweight < log(huge(1d0)) .and. logweight > log(tiny(1d0)) ) &
                    write(write_txt_unit,fmt_dbl_nposterior_norm) &
                    exp(logweight),-2*posterior_point(settings%pos_l),posterior_point(settings%pos_p0:settings%pos_d1)

            end if

        end do

        ! Now add the new posterior points
        do i_posterior=1,nposterior(0)


            ! Add to the new unnormalised posterior file
            write(write_untxt_unit,fmt_dbl_nposterior) posterior_points(:,i_posterior,0)

            logweight = posterior_points(settings%pos_w,i_posterior,0)-info%logevidence

            if(logweight < log(huge(1d0)) .and. logweight > log(tiny(1d0)) ) &
                write(write_txt_unit,fmt_dbl_nposterior_norm) &
                exp(logweight),-2*posterior_point(settings%pos_l),posterior_point(settings%pos_p0:settings%pos_d1)
        end do

        ! Close the files
        close(write_txt_unit)
        close(write_untxt_unit)
        close(read_untxt_unit)

        ! Move the temporary file to the new one
        call rename( trim(posterior_file(settings,.false.)), trim(posterior_file(settings,.true.)) )

        ! Re-set the posterior points
        nposterior(0) = 0


        if(settings%do_clustering) then

            do i_cluster=1,info%ncluster_A+info%ncluster_P

                if(nposterior(i_cluster) >= 1) then

                    ! Open the .txt file for writing
                    open(write_txt_unit,file=trim(posterior_file(settings,i=i_cluster)), action='write') 
                    ! Open the old unormalised .txt file for reading
                    open(read_untxt_unit,file=trim(posterior_file(settings,.true.,i_cluster)), action='read',iostat=io_stat) 
                    ! Open the new unormalised .txt file for writing
                    open(write_untxt_unit,file=trim(posterior_file(settings,.false.,i_cluster)), action='write') 

                    ! Trim off points that are too low in weight
                    do while(io_stat==0) 

                        ! Read a point from the old unormalised text file
                        read(read_untxt_unit,fmt_dbl_nposterior,iostat=io_stat) posterior_point

                        ! If this is still above the cdf threshold
                        if( posterior_point(settings%pos_Z) - info%logevidence > logsigma(settings%sigma_posterior) ) then

                            ! Re-copy to the new unnormalised posterior file
                            write(write_untxt_unit,fmt_dbl_nposterior) posterior_point

                            ! And add it to the .txt file
                            logweight = posterior_point(settings%pos_w) - info%logZ(i_cluster)

                            if(logweight < log(huge(1d0)) .and. logweight > log(tiny(1d0)) ) &
                                write(write_txt_unit,fmt_dbl_nposterior_norm) &
                                exp(logweight),-2*posterior_point(settings%pos_l),posterior_point(settings%pos_p0:settings%pos_d1)

                        end if

                    end do


                    ! Now add the new posterior points
                    do i_posterior=1,nposterior(i_cluster)

                        ! Add to the new unnormalised posterior file
                        write(write_untxt_unit,fmt_dbl_nposterior) posterior_points(:,i_posterior,i_cluster)

                        logweight = posterior_points(settings%pos_w,i_posterior,i_cluster)-info%logZ(i_cluster)

                        if(logweight < log(huge(1d0)) .and. logweight > log(tiny(1d0)) ) &
                            write(write_txt_unit,fmt_dbl_nposterior)   &
                            exp(logweight),-2*posterior_point(settings%pos_l),posterior_point(settings%pos_p0:settings%pos_d1)
                    end do

                    close(write_txt_unit)
                    close(write_untxt_unit)
                    close(read_untxt_unit)

                    ! Move the temporary file to the new one
                    call rename( trim(posterior_file(settings,.false.,i_cluster)), trim(posterior_file(settings,.true.,i_cluster)) )

                end if

            end do

            nposterior(1:) = 0

        end if



    end subroutine write_posterior_file

    subroutine write_phys_live_points(settings,info,live_points)
        use utils_module, only: DB_FMT,fmt_len,write_phys_unit,write_phys_cluster_unit
        use settings_module, only: program_settings 
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: info
        double precision, intent(in), dimension(settings%nTotal,settings%nlive,settings%ncluster) :: live_points

        integer i_err

        integer i_live
        integer i_cluster

        character(len=fmt_len) :: fmt_dbl_nphys


        ! Initialise the formats
        write(fmt_dbl_nphys,'("(",I0,A,")")') settings%nDims+settings%nDerived+1, DB_FMT


        ! Delete the old files
        open(write_phys_unit,file=trim(phys_live_file(settings)), action='write', iostat=i_err) 
        if(i_err.eq.0) close(write_phys_unit,status='delete')

        if(settings%do_clustering) then

            do i_cluster=1,settings%ncluster
                open(write_phys_cluster_unit,file=trim(phys_live_file(settings,i_cluster)), action='write', iostat=i_err) 
                if(i_err.eq.0) close(write_phys_cluster_unit,status='delete')
            end do
        end if

        ! Open a new file for appending to
        open(write_phys_unit,file=trim(phys_live_file(settings)), action='write', position="append")

        do i_cluster = 1,info%ncluster_A

            if(settings%do_clustering) open(write_phys_cluster_unit,file=trim(phys_live_file(settings,i_cluster)), action='write') 

            do i_live=1,info%n(i_cluster)
                if(settings%do_clustering) then
                    write(write_phys_unit,fmt_dbl_nphys) live_points(settings%p0:settings%d1,i_live,i_cluster), live_points(settings%l0,i_live,i_cluster)
                    write(write_phys_cluster_unit,fmt_dbl_nphys) live_points(settings%p0:settings%d1,i_live,i_cluster), live_points(settings%l0,i_live,i_cluster)
                else
                    write(write_phys_unit,fmt_dbl_nphys) live_points(settings%p0:settings%d1,i_live,i_cluster), live_points(settings%l0,i_live,i_cluster)
                end if
            end do

            if(settings%do_clustering) close(write_phys_cluster_unit)

        end do

        close(write_phys_unit)


    end subroutine write_phys_live_points


    subroutine write_stats_file(settings,info,ndead)
        use utils_module, only: DB_FMT,fmt_len,write_stats_unit,logzero,logsubexp
        use settings_module, only: program_settings
        use run_time_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: info
        integer,                intent(in) :: ndead

        double precision, dimension(info%ncluster_A + info%ncluster_P) :: mu
        double precision, dimension(info%ncluster_A + info%ncluster_P) :: sigma

        integer :: i

        character(len=fmt_len) :: fmt_1

        open(write_stats_unit,file=trim(stats_file(settings)), action='write') 


        write(write_stats_unit, '("Evidence estimates:")')
        write(write_stats_unit, '("===================")')
        write(write_stats_unit, '("  - The evidence Z is a log-normally distributed, with location and scale parameters mu and sigma.")')
        write(write_stats_unit, '("  - We denote this as log(Z) = mu +/- sigma.")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Global evidence:")')
        write(write_stats_unit, '("----------------")')
        write(write_stats_unit,'("")')
        if(info%logevidence>logzero .and. info%logevidence2>logzero) then
            mu(1)    = 2*info%logevidence - 0.5*info%logevidence2              
            sigma(1) = sqrt(info%logevidence2 - 2*info%logevidence)
            write(fmt_1,'("(""log(Z)       = "",", A, ","" +/- "",", A, ")")') DB_FMT,DB_FMT
            write(write_stats_unit,fmt_1) mu(1),sigma(1)
        else
            write(write_stats_unit, '(" Too early to produce a sensible evidence estimate ")')
        end if

        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Local evidences:")')
        write(write_stats_unit, '("----------------")')
        write(write_stats_unit, '(I2, " clusters found so far")') info%ncluster_A + info%ncluster_P
        write(write_stats_unit, '(" ", I2, " still active")') info%ncluster_A 
        write(write_stats_unit,'("")')

        mu    = 2*info%logZ(:info%ncluster_A+info%ncluster_P) - 0.5*info%logZ2(:info%ncluster_A+info%ncluster_P)              
        sigma = sqrt(info%logZ2(:info%ncluster_A+info%ncluster_P) - 2*info%logZ(:info%ncluster_A+info%ncluster_P))


        if(info%ncluster_A>=1) then
            do i=1,info%ncluster_A
                write(fmt_1,'("(""log(Z_"",",A,","")  = "",", A, ","" +/- "",", A, ")")') 'I2',DB_FMT,DB_FMT
                if(info%logZ(i)>logzero) then
                    write(write_stats_unit,fmt_1) i, mu(i),sigma(i)
                else
                    write(write_stats_unit,'("log(Z_",I2,")  = ?", "(still evaluating)")') i
                end if
            end do
        end if

        if(info%ncluster_P>=1) then
            do i=info%ncluster_A+1,info%ncluster_A+info%ncluster_P

                if(info%logZ(i)>logzero) then
                    write(write_stats_unit,fmt_1) i, mu(i),sigma(i)
                else
                    write(write_stats_unit,'("log(Z_",I2,")  = ?")') i
                end if
            end do
        end if

        write(write_stats_unit,'("")')
        write(write_stats_unit,'("")')
        write(write_stats_unit, '("Run-time information:")')
        write(write_stats_unit, '("---------------------")')
        write(write_stats_unit,'("")')
        write(write_stats_unit,'(" ndead:    ", I8)') ndead
        write(write_stats_unit,'(" nlive:    ", I8)') settings%nlive
        write(write_stats_unit,'(" active clusters: ", I8)') info%ncluster_A
        if(info%ncluster_A>1) then
            do i=1,info%ncluster_A
                write(write_stats_unit,'(" nlive (cluster", I3, "): ", I8)') i, info%n(i)
            end do
        end if
        if(info%ncluster_T>1) write(write_stats_unit, '( I4, " clusters found so far" )') info%ncluster_A + info%ncluster_P



        close(write_stats_unit)

    end subroutine write_stats_file





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

    function posterior_file(settings,reading,i) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings
        logical,intent(in),optional :: reading
        integer,intent(in),optional :: i

        character(STR_LENGTH) :: file_name

        character(STR_LENGTH) :: cluster_num

        if(present(i)) then

            write(cluster_num,'(I5)') i
            file_name = trim(cluster_dir(settings)) // '/' // trim(settings%file_root) // '_'

            if(present(reading)) then
                if(reading) then
                    file_name = trim(file_name) // 'unnormalised_'
                else
                    file_name = trim(file_name) // 'unnormalised_temp_'
                endif
            end if

            file_name = trim(file_name) // trim(adjustl(cluster_num)) //'.txt'
        else 

            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root)

            if(present(reading)) then
                if(reading) then
                    file_name = trim(file_name) // '_unnormalised'
                else
                    file_name = trim(file_name) // '_unnormalised_temp' 
                endif
            end if

            file_name = trim(file_name) // '.txt'
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
            file_name = trim(cluster_dir(settings)) // '/' // trim(settings%file_root) // '_phys_live_' // trim(adjustl(cluster_num)) //'.txt'
        else 
            file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '_phys_live.txt'
        end if

    end function phys_live_file


end module
