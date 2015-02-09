!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

    subroutine write_resume_file(settings,info,live_points,nphantom,phantom_points,ndead,total_likelihood_calls)
        use utils_module, only: DB_FMT,INT_FMT,fmt_len,write_resume_unit
        use evidence_module, only: run_time_info
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: info

        double precision,intent(in), dimension(settings%nTotal,settings%nlive,settings%ncluster) :: live_points

        integer,intent(in),dimension(settings%ncluster) :: nphantom
        double precision,intent(in), dimension(settings%nTotal,settings%nstack,settings%ncluster) :: phantom_points

        integer,intent(in) :: ndead
        integer,intent(in) :: total_likelihood_calls

        integer :: i_cluster

        character(len=fmt_len) :: fmt_dbl_ncluster_A
        character(len=fmt_len) :: fmt_dbl_ncluster_T
        character(len=fmt_len) :: fmt_dbl_nTotal
        character(len=fmt_len) :: fmt_dbl_2

        character(len=fmt_len) :: fmt_int_ncluster_A
        character(len=fmt_len) :: fmt_int_1


        ! Initialise the formats
        write(fmt_dbl_ncluster_A,'("(",I0,A,")")') info%ncluster_A, DB_FMT
        write(fmt_dbl_ncluster_T,'("(",I0,A,")")') info%ncluster_T, DB_FMT
        write(fmt_dbl_nTotal    ,'("(",I0,A,")")') settings%nTotal, DB_FMT
        write(fmt_dbl_2         ,'("(",I0,A,")")') 2,               DB_FMT

        write(fmt_int_ncluster_A,'("(",I0,A,")")') info%ncluster_A, INT_FMT
        write(fmt_int_1         ,'("(",I0,A,")")') 1,               INT_FMT
        
        
        ! Open the .resume file
        open(write_resume_unit,file=trim(resume_file(settings))//'_new', action='write') 


        ! Cluster information: 

        ! Number of Active clusters
        write(write_resume_unit,fmt_int_1) info%ncluster_A
        ! Number of Passive clusters
        write(write_resume_unit,fmt_int_1) info%ncluster_P
        ! Total amount of evidence information
        write(write_resume_unit,fmt_int_1) info%ncluster_T
        ! global log evidence and evidence^2
        write(write_resume_unit,fmt_dbl_2) info%logevidence, info%logevidence2
        ! number of live points in each cluster
        write(write_resume_unit,fmt_int_ncluster_A) info%n(:info%ncluster_A)
        ! Log likelihood contours
        write(write_resume_unit,fmt_dbl_ncluster_A) info%logL(:info%ncluster_A)
        ! Log volume
        write(write_resume_unit,fmt_dbl_ncluster_A) info%logX(:info%ncluster_A)
        ! Log evidence
        write(write_resume_unit,fmt_dbl_ncluster_A) info%logZ(:info%ncluster_A)
        ! Log evidence^2
        write(write_resume_unit,fmt_dbl_ncluster_A) info%logZ2(:info%ncluster_A)
        ! Correlations: log(X_iX_j)
        write(write_resume_unit,fmt_dbl_ncluster_A) info%logXX(:info%ncluster_A,:info%ncluster_A)
        ! Correlations: log(Z_iX_j)
        write(write_resume_unit,fmt_dbl_ncluster_T) info%logZX(:info%ncluster_T,:info%ncluster_A)


        ! Live points
        do i_cluster=1,info%ncluster_A
            write(write_resume_unit,fmt_dbl_nTotal) live_points(:,:info%n(i_cluster),i_cluster)
        end do

        ! number of phantom points
        write(write_resume_unit,fmt_int_ncluster_A) nphantom(:info%ncluster_A)
        ! Phantom points
        do i_cluster=1,info%ncluster_A
            write(write_resume_unit,fmt_dbl_nTotal) phantom_points(:,:nphantom(i_cluster),i_cluster)
        end do

        ! number of dead points
        write(write_resume_unit,fmt_int_1) ndead
        ! total likelihood calls
        write(write_resume_unit,fmt_int_1) total_likelihood_calls


        close(write_resume_unit)

        call rename(trim(resume_file(settings))//'_new',trim(resume_file(settings)))

    end subroutine write_resume_file


    subroutine read_resume_file(settings,info,live_points,nphantom,phantom_points,ndead,total_likelihood_calls)
        use utils_module, only: DB_FMT,INT_FMT,fmt_len,read_resume_unit
        use evidence_module, only: run_time_info
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(inout) :: info

        double precision,intent(out), dimension(settings%nTotal,settings%nlive,settings%ncluster) :: live_points

        integer,intent(out),dimension(settings%ncluster) :: nphantom
        double precision,intent(out), dimension(settings%nTotal,settings%nstack,settings%ncluster) :: phantom_points

        integer,intent(out) :: ndead
        integer,intent(out) :: total_likelihood_calls

        integer :: i_cluster

        character(len=fmt_len) :: fmt_dbl_ncluster_A
        character(len=fmt_len) :: fmt_dbl_ncluster_T
        character(len=fmt_len) :: fmt_dbl_nTotal
        character(len=fmt_len) :: fmt_dbl_2

        character(len=fmt_len) :: fmt_int_ncluster_A
        character(len=fmt_len) :: fmt_int_1


        ! Initialise the formats
        write(fmt_dbl_ncluster_A,'("(",I0,A,")")') info%ncluster_A, DB_FMT
        write(fmt_dbl_ncluster_T,'("(",I0,A,")")') info%ncluster_T, DB_FMT
        write(fmt_dbl_nTotal    ,'("(",I0,A,")")') settings%nTotal, DB_FMT
        write(fmt_dbl_2         ,'("(",I0,A,")")') 2,               DB_FMT

        write(fmt_int_ncluster_A,'("(",I0,A,")")') info%ncluster_A, INT_FMT
        write(fmt_int_1         ,'("(",I0,A,")")') 1,               INT_FMT
        
        
        ! Open the .resume file
        open(read_resume_unit,file=trim(resume_file(settings)), action='read') 


        ! Cluster information: 

        ! Number of Active clusters
        read(read_resume_unit,fmt_int_1) info%ncluster_A
        ! Number of Passive clusters
        read(read_resume_unit,fmt_int_1) info%ncluster_P
        ! Total amount of evidence information
        read(read_resume_unit,fmt_int_1) info%ncluster_T
        ! global log evidence and evidence^2
        read(read_resume_unit,fmt_dbl_2) info%logevidence, info%logevidence2
        ! number of live points in each cluster
        read(read_resume_unit,fmt_int_ncluster_A) info%n(:info%ncluster_A)
        ! Log likelihood contours
        read(read_resume_unit,fmt_dbl_ncluster_A) info%logL(:info%ncluster_A)
        ! Log volume
        read(read_resume_unit,fmt_dbl_ncluster_A) info%logX(:info%ncluster_A)
        ! Log evidence
        read(read_resume_unit,fmt_dbl_ncluster_A) info%logZ(:info%ncluster_A)
        ! Log evidence^2
        read(read_resume_unit,fmt_dbl_ncluster_A) info%logZ2(:info%ncluster_A)
        ! Correlations: log(X_iX_j)
        read(read_resume_unit,fmt_dbl_ncluster_A) info%logXX(:info%ncluster_A,:info%ncluster_A)
        ! Correlations: log(Z_iX_j)
        read(read_resume_unit,fmt_dbl_ncluster_T) info%logZX(:info%ncluster_T,:info%ncluster_A)


        ! Live points
        do i_cluster=1,info%ncluster_A
            read(read_resume_unit,fmt_dbl_nTotal) live_points(:,:info%n(i_cluster),i_cluster)
        end do

        ! number of phantom points
        read(read_resume_unit,fmt_int_ncluster_A) nphantom(:info%ncluster_A)
        ! Phantom points
        do i_cluster=1,info%ncluster_A
            read(read_resume_unit,fmt_dbl_nTotal) phantom_points(:,:nphantom(i_cluster),i_cluster)
        end do

        ! number of dead points
        read(read_resume_unit,fmt_int_1) ndead
        ! total likelihood calls
        read(read_resume_unit,fmt_int_1) total_likelihood_calls


        close(read_resume_unit)

    end subroutine read_resume_file



    subroutine write_posterior_file(settings,info,posterior_points,nposterior) 
        use utils_module, only: DB_FMT,fmt_len,write_txt_unit,write_untxt_unit,read_untxt_unit,logsigma
        use settings_module, only: program_settings
        use evidence_module, only: run_time_info 
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
        use evidence_module, only: run_time_info 
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
        use evidence_module, only: run_time_info 
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
    function resume_file(settings) result(file_name)
        use settings_module, only: program_settings
        use utils_module,    only: STR_LENGTH
        implicit none
        type(program_settings), intent(in) :: settings

        character(STR_LENGTH) :: file_name

        file_name = trim(settings%base_dir) // '/' // trim(settings%file_root) // '.resume'

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
