!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

    subroutine write_resume_file(settings,M,live_data,evidence_vec,ndead,mean_likelihood_calls,total_likelihood_calls,nposterior,posterior_array)
        use utils_module, only: DBL_FMT,write_resume_unit
        use model_module, only: model
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(model),            intent(in) :: M
        double precision, dimension(:,:) :: live_data
        integer :: nposterior
        double precision :: mean_likelihood_calls
        integer :: total_likelihood_calls
        double precision, dimension(settings%nDims+2,settings%nmax_posterior) :: posterior_array
        double precision, dimension(6)             :: evidence_vec
        integer :: ndead

        integer :: i_err


        
        ! Open the .resume file, note the presence of iostat prevents program
        ! termination during this write
        open(write_resume_unit,file=trim(settings%file_root) // '.resume', action='write', iostat=i_err) 

        ! Live points
        write(write_resume_unit,'(<settings%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
        ! Evidence vector
        write(write_resume_unit,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
        ! number of dead points
        write(write_resume_unit,'(I)') ndead
        ! mean likelihood calls
        write(write_resume_unit,'(E<DBL_FMT(1)>.<DBL_FMT(2)>)') mean_likelihood_calls
        ! total likelihood calls
        write(write_resume_unit,'(I)') total_likelihood_calls
        ! Number of saved posterior points
        write(write_resume_unit,'(I)') nposterior
        ! posterior points
        write(write_resume_unit,'(<M%nDims+M%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array(:,:nposterior)

        close(write_resume_unit)

    end subroutine write_resume_file

    subroutine write_posterior_file(settings,M,posterior_array,evidence,nposterior) 
        use utils_module, only: DBL_FMT,write_txt_unit,logzero
        use model_module, only: model
        use settings_module, only: program_settings

        implicit none


        integer :: nposterior
        type(program_settings), intent(in) :: settings
        type(model),            intent(in) :: M
        double precision, dimension(settings%nDims+settings%nDerived+2,nposterior) :: posterior_array
        double precision                           :: evidence

        integer :: i_err

        integer :: i_posterior

        double precision :: logminimum_weight

        logminimum_weight = log(settings%minimum_weight)
        
        ! Open the .txt file, note the presence of iostat prevents program
        ! termination during this write
        open(write_txt_unit,file=trim(settings%file_root) // '.txt' , action='write', iostat=i_err) 

        do i_posterior=1,nposterior

            if (posterior_array(1,i_posterior)-evidence > logminimum_weight)                 &
                write(write_txt_unit,'(<settings%nDims+settings%nDerived+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)')   &
                exp(posterior_array(1,i_posterior)-evidence),posterior_array(2:,i_posterior)
        end do

        close(write_txt_unit)

    end subroutine write_posterior_file




end module read_write_module
