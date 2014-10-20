!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

    subroutine write_resume_file(settings,live_points,stack_size,phantom_points,evidence_vec,ndead,total_likelihood_calls,nposterior,posterior_array)
        use utils_module, only: DBL_FMT,write_resume_unit
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        integer,intent(in) :: stack_size
        double precision,intent(in), dimension(settings%nTotal,settings%nlive) :: live_points
        double precision,intent(in), dimension(settings%nTotal,settings%nstack) :: phantom_points
        integer :: nposterior
        integer :: total_likelihood_calls
        double precision, dimension(settings%nDims+settings%nDerived+2,settings%nmax_posterior) :: posterior_array
        double precision, dimension(:)             :: evidence_vec
        integer :: ndead

        integer :: i_err


        
        ! Open the .resume file, note the presence of iostat prevents program
        ! termination during this write
        open(write_resume_unit,file=trim(settings%file_root) // '.resume', action='write', iostat=i_err) 

        ! Live points
        write(write_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_points
        ! Stack size
        write(write_resume_unit,'(I)') stack_size
        ! Phantom points
        write(write_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') phantom_points(:,:stack_size)
        ! Evidence vector
        write(write_resume_unit,'(<size(evidence_vec)>E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
        ! number of dead points
        write(write_resume_unit,'(I)') ndead
        ! total likelihood calls
        write(write_resume_unit,'(I)') total_likelihood_calls
        ! Number of saved posterior points
        write(write_resume_unit,'(I)') nposterior
        ! posterior points
        write(write_resume_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)

        close(write_resume_unit)

    end subroutine write_resume_file

    subroutine write_posterior_file(settings,posterior_array,evidence,nposterior) 
        use utils_module, only: DBL_FMT,write_txt_unit,logzero
        use settings_module, only: program_settings

        implicit none


        integer :: nposterior
        type(program_settings), intent(in) :: settings
        double precision, dimension(settings%nDims+settings%nDerived+2,nposterior) :: posterior_array
        double precision                           :: evidence

        integer :: i_err

        integer :: i_posterior
        
        ! Open the .txt file, note the presence of iostat prevents program
        ! termination during this write
        open(write_txt_unit,file=trim(settings%file_root) // '.txt' , action='write', iostat=i_err) 

        do i_posterior=1,nposterior
            if(posterior_array(1,i_posterior)-evidence < log(huge(1d0)) ) then
                write(write_txt_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)')   &
                    exp(posterior_array(1,i_posterior)-evidence),posterior_array(2:,i_posterior)
            end if
        end do

        close(write_txt_unit)

    end subroutine write_posterior_file

    subroutine write_phys_live_points(settings,live_points)
        use utils_module, only: DBL_FMT,write_phys_unit
        use settings_module, only: program_settings
        implicit none

        type(program_settings), intent(in) :: settings
        double precision, intent(in), dimension(settings%nTotal,settings%nstack) :: live_points

        integer i_err

        integer i_live

        open(write_phys_unit,file=trim(settings%file_root) // '_phys_live.txt' , action='write', iostat=i_err) 

        do i_live=1,settings%nlive
            write(write_phys_unit,'(<settings%nDims+settings%nDerived+3>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_points(settings%p0:settings%d1,i_live), live_points(settings%cluster,i_live), live_points(settings%l0,i_live)
        end do

        close(write_phys_unit)

    end subroutine write_phys_live_points




end module read_write_module
