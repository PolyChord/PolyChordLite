!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

    subroutine write_resume_file(settings,info,live_points,nphantom,phantom_points,ndead,total_likelihood_calls,nposterior,posterior_points)
        use utils_module, only: DBL_FMT,write_resume_unit
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

        integer,intent(in),dimension(0:settings%ncluster) :: nposterior
        double precision,intent(in), dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior,0:settings%ncluster*2) :: posterior_points

        integer :: i_err

        integer :: i_cluster


        
        ! Open the .resume file, note the presence of iostat prevents program
        ! termination during this write
        open(write_resume_unit,file=trim(settings%file_root) // '.resume', action='write', iostat=i_err) 


        ! Cluster information: 

        ! Number of Active clusters
        write(write_resume_unit,'(I)') info%ncluster_A
        ! Number of Active + Passive clusters
        write(write_resume_unit,'(I)') info%ncluster_T
        ! global log evidence and evidence^2
        write(write_resume_unit,'(2E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logevidence, info%logevidence2
        ! number of live points in each cluster
        write(write_resume_unit,'(<info%ncluster_A>I)') info%n(:info%ncluster_A)
        ! Log likelihood contours
        write(write_resume_unit,'(<info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logL(:info%ncluster_A)
        ! Log volume
        write(write_resume_unit,'(<info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logX(:info%ncluster_A)
        ! Log evidence
        write(write_resume_unit,'(<info%ncluster_T>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logZ(:info%ncluster_A)
        ! Log evidence^2
        write(write_resume_unit,'(<info%ncluster_T>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logZ2(:info%ncluster_A)
        ! Correlations: log(X_iX_j)
        write(write_resume_unit,'(<info%ncluster_A*info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logXX(:info%ncluster_A,:info%ncluster_A)
        ! Correlations: log(Z_iX_j)
        write(write_resume_unit,'(<info%ncluster_T*info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logZX(:info%ncluster_T,:info%ncluster_A)


        ! Live points
        do i_cluster=1,info%ncluster_A
            write(write_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_points(:,:info%n(i_cluster),i_cluster)
        end do

        ! number of phantom points
        write(write_resume_unit,'(<info%ncluster_A>I)') nphantom(:info%ncluster_A)
        ! Phantom points
        do i_cluster=1,info%ncluster_A
            write(write_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') phantom_points(:,:nphantom(i_cluster),i_cluster)
        end do

        ! number of dead points
        write(write_resume_unit,'(I)') ndead
        ! total likelihood calls
        write(write_resume_unit,'(I)') total_likelihood_calls
        ! Number of saved posterior points
        write(write_resume_unit,'(<info%ncluster_A+1>I)') nposterior(0:info%ncluster_A) 
        ! posterior points
        do i_cluster=0,info%ncluster_A
            write(write_resume_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_points(:,:nposterior(i_cluster),i_cluster)
        end do

        close(write_resume_unit)

    end subroutine write_resume_file

    subroutine read_resume_file(settings,info,live_points,nphantom,phantom_points,ndead,total_likelihood_calls,nposterior,posterior_points)
        use utils_module, only: DBL_FMT,read_resume_unit
        use evidence_module, only: run_time_info
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(out) :: info

        double precision,intent(out), dimension(settings%nTotal,settings%nlive,settings%ncluster) :: live_points

        integer,intent(out),dimension(settings%ncluster) :: nphantom
        double precision,intent(out), dimension(settings%nTotal,settings%nstack,settings%ncluster) :: phantom_points

        integer,intent(out) :: ndead
        integer,intent(out) :: total_likelihood_calls

        integer,intent(out),dimension(0:settings%ncluster) :: nposterior
        double precision,intent(out), dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior,0:settings%ncluster*2) :: posterior_points

        integer :: i_err

        integer :: i_cluster


        
        ! Open the .resume file, note the presence of iostat prevents program
        ! termination during this read
        open(read_resume_unit,file=trim(settings%file_root) // '.resume', action='read', iostat=i_err) 


        ! Cluster information: 

        ! Number of Active clusters
        read(read_resume_unit,'(I)') info%ncluster_A
        ! Number of Active + Passive clusters
        read(read_resume_unit,'(I)') info%ncluster_T
        ! global log evidence and evidence^2
        read(read_resume_unit,'(2E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logevidence, info%logevidence2
        ! number of live points in each cluster
        read(read_resume_unit,'(<info%ncluster_A>I)') info%n(:info%ncluster_A)
        ! Log likelihood contours
        read(read_resume_unit,'(<info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logL(:info%ncluster_A)
        ! Log volume
        read(read_resume_unit,'(<info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logX(:info%ncluster_A)
        ! Log evidence
        read(read_resume_unit,'(<info%ncluster_T>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logZ(:info%ncluster_A)
        ! Log evidence^2
        read(read_resume_unit,'(<info%ncluster_T>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logZ2(:info%ncluster_A)
        ! Correlations: log(X_iX_j)
        read(read_resume_unit,'(<info%ncluster_A*info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logXX(:info%ncluster_A,:info%ncluster_A)
        ! Correlations: log(Z_iX_j)
        read(read_resume_unit,'(<info%ncluster_T*info%ncluster_A>E<DBL_FMT(1)>.<DBL_FMT(2)>)') info%logZX(:info%ncluster_T,:info%ncluster_A)


        do i_cluster=1,info%ncluster_A
            read(read_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_points(:,:info%n(i_cluster),i_cluster)
        end do

        ! number of phantom points
        read(read_resume_unit,'(<info%ncluster_A>I)') nphantom(:info%ncluster_A) 
        ! Phantom points
        do i_cluster=1,info%ncluster_A
            read(read_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') phantom_points(:,:nphantom(i_cluster),i_cluster)
        end do

        ! number of dead points
        read(read_resume_unit,'(I)') ndead
        ! total likelihood calls
        read(read_resume_unit,'(I)') total_likelihood_calls
        ! Number of saved posterior points
        read(read_resume_unit,'(<info%ncluster_A+1>I)') nposterior(0:info%ncluster_A)
        ! posterior points
        do i_cluster=0,info%ncluster_A
            read(read_resume_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_points(:,:nposterior(i_cluster),i_cluster)
        end do

        close(read_resume_unit)

    end subroutine read_resume_file



    subroutine write_posterior_file(settings,info,posterior_points,nposterior,active_clusters) 
        use utils_module, only: DBL_FMT,write_txt_unit,logzero,STR_LENGTH
        use settings_module, only: program_settings
        use evidence_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: info
        double precision,intent(in), dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior,0:settings%ncluster*2) :: posterior_points
        integer,intent(in),dimension(0:settings%ncluster) :: nposterior
        integer,intent(in)                :: active_clusters

        integer :: i_err

        integer :: i_posterior
        integer i_cluster
        character(STR_LENGTH) :: cluster_num
        
        ! Open the .txt file, note the presence of iostat prevents program
        ! termination during this write
        open(write_txt_unit,file=trim(settings%file_root) // '.txt' , action='write', iostat=i_err) 

        do i_posterior=1,nposterior(i_cluster)
            if(posterior_points(1,i_posterior,i_cluster)-info%logevidence < log(huge(1d0)) ) then
                write(write_txt_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)')   &
                    exp(posterior_points(1,i_posterior,i_cluster)-info%logevidence),posterior_points(2:,i_posterior,i_cluster)
            end if
        end do

        close(write_txt_unit)



        do i_cluster=1,settings%ncluster

            if(nposterior(i_cluster) >= 1) then
                write(cluster_num,'(I5)') i_cluster
                open(write_txt_unit,file=trim(settings%file_root) // "_cluster_" //trim(adjustl(cluster_num))// '.txt' , action='write', iostat=i_err) 

                do i_posterior=1,nposterior(i_cluster)
                    if(posterior_points(1,i_posterior,i_cluster)-info%logZ(i_cluster) < log(huge(1d0)) ) then
                        write(write_txt_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)')   &
                            exp(posterior_points(1,i_posterior,i_cluster)-info%logevidence),posterior_points(2:,i_posterior,i_cluster)
                    end if
                end do

                close(write_txt_unit)
            end if

        end do


    end subroutine write_posterior_file

    subroutine write_phys_live_points(settings,info,live_points)
        use utils_module, only: DBL_FMT,write_phys_unit,write_phys_cluster_unit,STR_LENGTH
        use settings_module, only: program_settings 
        use evidence_module, only: run_time_info 
        implicit none

        type(program_settings), intent(in) :: settings
        type(run_time_info),    intent(in) :: info
        double precision, intent(in), dimension(settings%nTotal,settings%nlive,settings%ncluster) :: live_points

        integer i_err

        integer i_live
        integer i_cluster
        character(STR_LENGTH) :: cluster_num

        ! Delete the old files
        open(write_phys_unit,file=trim(settings%file_root) // '_phys_live.txt' , action='write', iostat=i_err) 
        if(i_err.eq.0) close(write_phys_unit,status='delete')

        do i_cluster=1,settings%ncluster
            write(cluster_num,'(I5)') i_cluster
            open(write_phys_cluster_unit,file=trim(settings%file_root) // '_phys_live_' //trim(adjustl(cluster_num))//'.txt' , action='write', iostat=i_err) 
            if(i_err.eq.0) close(write_phys_cluster_unit,status='delete')
        end do

        ! Open a new file for appending to
        open(write_phys_unit,file=trim(settings%file_root) // '_phys_live.txt' , action='write', iostat=i_err, position="append")

        do i_cluster = 1,info%ncluster_A

            if(settings%do_clustering) then
                write(cluster_num,'(I5)') i_cluster
                open(write_phys_cluster_unit,file=trim(settings%file_root) // '_phys_live_' //trim(adjustl(cluster_num))//'.txt' , action='write', iostat=i_err) 
            end if

            do i_live=1,info%n(i_cluster)
                if(settings%do_clustering) then
                    write(write_phys_unit,'(<settings%nDims+settings%nDerived>E<DBL_FMT(1)>.<DBL_FMT(2)>, I5, E<DBL_FMT(1)>.<DBL_FMT(2)>)') &
                        live_points(settings%p0:settings%d1,i_live,i_cluster), i_cluster, live_points(settings%l0,i_live,i_cluster)
                    write(write_phys_cluster_unit,'(<settings%nDims+settings%nDerived>E<DBL_FMT(1)>.<DBL_FMT(2)>, I5, E<DBL_FMT(1)>.<DBL_FMT(2)>)') &
                        live_points(settings%p0:settings%d1,i_live,i_cluster), i_cluster, live_points(settings%l0,i_live,i_cluster)
                else
                    write(write_phys_unit,'(<settings%nDims+settings%nDerived+1>E<DBL_FMT(1)>.<DBL_FMT(2)>)') &
                        live_points(settings%p0:settings%d1,i_live,i_cluster), live_points(settings%l0,i_live,i_cluster)
                end if
            end do

            if(settings%do_clustering) close(write_phys_cluster_unit)

        end do

        close(write_phys_unit)


    end subroutine write_phys_live_points




end module read_write_module
