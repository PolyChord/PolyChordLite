!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

    subroutine write_resume_file(settings,M,live_data,evidence_vec,ndead)
        use utils_module, only: STR_LENGTH,DBL_FMT
        use model_module, only: model
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(model),            intent(in) :: M
        double precision, dimension(M%nTotal,settings%nlive) :: live_data
        double precision, dimension(6)             :: evidence_vec
        integer :: ndead


        
        open(100,file=trim(settings%file_root) // '.resume' ) 

        ! Live points
        write(100,'(<M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
        ! Evidence vector
        write(100,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
        ! number of dead points
        write(100,'(I)') ndead

        close(100)

    end subroutine write_resume_file




!    subroutine write_log_files(settings,M,live_data,evidence_vec,ndead)
!        use utils_module, only: STR_LENGTH,DBL_FMT
!        use model_module, only: model
!        use settings_module, only: program_settings
!
!        implicit none
!
!
!        type(program_settings), intent(in) :: settings
!        type(model),            intent(in) :: M
!        double precision, dimension(M%nTotal,settings%nlive) :: live_data
!        double precision, dimension(6)             :: evidence_vec
!        integer :: ndead
!        character(len=STR_LENGTH) :: temp_format
!
!
!        ! Write the live_points to file
!        open(100,file=trim(settings%file_root) // '_live.dat' ) 
!
!        temp_format ='(<M%nTotal>'//DBL_FMT//')'
!        write(100,fmt=temp_format) live_data
!
!        close(100)
!
!        ! Write the evidence vector to file
!        open(100,file=trim(settings%file_root) // '_stats.dat' ) 
!
!        write(100,'(A)') '# evidence vector'
!        write(100,'(6'//DBL_FMT//')') evidence_vec
!        write(100,'(A)') '# number of dead points'
!        write(100,'(I)') ndead
!
!        close(100)
!
!    end subroutine write_log_files
!
!
!
!    subroutine write_resume_files(settings,M,live_data,evidence_vec,ndead)
!        use utils_module, only: STR_LENGTH,DBL_FMT
!        use model_module, only: model
!        use settings_module, only: program_settings
!
!        implicit none
!
!
!        type(program_settings), intent(in) :: settings
!        type(model),            intent(in) :: M
!        double precision, dimension(M%nTotal,settings%nlive) :: live_data
!        double precision, dimension(6)             :: evidence_vec
!        integer :: ndead
!
!        character(len=STR_LENGTH) :: temp_format
!
!        ! Write the live_points to file
!        open(100,file=trim(settings%file_root) // '_live.dat' ) 
!
!        temp_format='(<M%nTotal>'//DBL_FMT//')'
!        write(100,fmt=temp_format) live_data
!
!        close(100)
!
!        ! Write the evidence vector to file
!        open(100,file=trim(settings%file_root) // '_stats.dat' ) 
!
!        write(100,'(A)') '# evidence vector'
!        write(100,'(8'//DBL_FMT//')') evidence_vec
!        write(100,'(A)') '# number of dead points'
!        write(100,'(I)') ndead
!
!        close(100)
!
!    end subroutine write_resume_files
!
!
!
!    !> Write all of the details of the algorithm settings and model for reproducibility purposes
!    subroutine write_settings_file(settings,M)
!        use utils_module, only: STR_LENGTH,DBL_FMT
!        use model_module, only: model
!        use settings_module, only: program_settings
!
!        implicit none
!
!
!        type(program_settings), intent(in) :: settings
!        type(model),            intent(in) :: M
!        integer :: ndead
!
!        character(len=STR_LENGTH) :: temp_format
!
!
!        ! Write the settings to file
!        open(100,file=trim(settings%file_root) // '_settings.dat' ) 
!
!        write(100,'(A)') '###### Algorithm Settings ###### '
!
!        write(100,'(A)') '# nlive '
!        write(100,'(A)') settings%file_root
!
!        write(100,'(A)') '# nlive '
!        write(100,'(I)') settings%nlive
!
!        write(100,'(A)') '# feedback '
!        write(100,'(I)') settings%feedback
!
!        write(100,'(A)')     '# precision_criterion '
!        write(100,'(6'//DBL_FMT//')') settings%precision_criterion
!
!        write(100,'(A)') '# max_ndead '
!        write(100,'(I)') settings%max_ndead
!
!        write(100,'(A)') '# save_dead '
!        write(100,'(L)') settings%save_dead
!
!        write(100,'(A)') '# num_chords '
!        write(100,'(I)') settings%num_chords
!
!        write(100,'(A)') '###### Model Settings ########## '
!
!        write(100,'(A)') '# nDims '
!        write(100,'(I)') M%nDims
!
!        write(100,'(A)') '# nDerived '
!        write(100,'(I)') M%nDerived
!
!        if(M%uniform_num>0) then
!            write(100,'(A)') '# uniform_num '
!            write(100,'(I)') M%uniform_num
!
!            temp_format='(<M%uniform_num>'//DBL_FMT//')'
!            write(100,'(A)') '# uniform_params '
!            write(100,fmt=temp_format) M%uniform_params
!        end if
!        if(M%log_uniform_num>0) then
!            write(100,'(A)') '# log_uniform_num '
!            write(100,'(I)') M%log_uniform_num
!
!            temp_format='(<M%log_uniform_num>'//DBL_FMT//')'
!            write(100,'(A)') '# log_uniform_params '
!            write(100,fmt=temp_format) M%log_uniform_params
!        end if
!        if(M%gaussian_num>0) then
!            write(100,'(A)') '# gaussian_num '
!            write(100,'(I)') M%gaussian_num
!
!            temp_format='(<M%gaussian_num>'//DBL_FMT//')'
!            write(100,'(A)') '# gaussian_params '
!            write(100,fmt=temp_format) M%gaussian_params
!        end if
!        if(M%sorted_uniform_num>0) then
!            write(100,'(A)') '# sorted_uniform_num '
!            write(100,'(I)') M%sorted_uniform_num
!
!            temp_format='(<M%sorted_uniform_num>'//DBL_FMT//')'
!            write(100,'(A)') '# sorted_uniform_params '
!            write(100,fmt=temp_format) M%sorted_uniform_params
!        end if
!
!        close(100)
!
!
!    end subroutine write_settings_file




end module read_write_module
