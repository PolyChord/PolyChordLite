!> Module for reading and writing to files
module read_write_module
    implicit none

    contains

    subroutine write_resume_file(settings,M,first_index_live,live_data,evidence_vec,ndead,first_index_posterior,last_index_posterior,posterior_array)
        use utils_module, only: DBL_FMT,write_resume_unit
        use model_module, only: model
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(model),            intent(in) :: M
        integer :: first_index_live
        double precision, dimension(M%nTotal,settings%nlive) :: live_data
        integer :: first_index_posterior
        integer :: last_index_posterior
        double precision, dimension(M%nDims+2,settings%nmax_posterior) :: posterior_array
        double precision, dimension(6)             :: evidence_vec
        integer :: ndead

        integer :: i_err


        
        ! Open the .resume file, note the presence of iostat prevents program
        ! termination during this write
        open(write_resume_unit,file=trim(settings%file_root) // '.resume', action='write', iostat=i_err) 

        ! First index of the live points
        write(write_resume_unit,'(I)') first_index_live
        ! Live points
        write(write_resume_unit,'(<M%nTotal>E<DBL_FMT(1)>.<DBL_FMT(2)>)') live_data
        ! Evidence vector
        write(write_resume_unit,'(6E<DBL_FMT(1)>.<DBL_FMT(2)>)') evidence_vec
        ! number of dead points
        write(write_resume_unit,'(I)') ndead
        ! First index of the posterior_points
        write(write_resume_unit,'(I)') first_index_posterior
        write(write_resume_unit,'(I)') last_index_posterior
        ! posterior points
        write(write_resume_unit,'(<M%nDims+4>E<DBL_FMT(1)>.<DBL_FMT(2)>)') posterior_array

        close(write_resume_unit)

    end subroutine write_resume_file

    subroutine write_posterior_file(settings,M,posterior_array,evidence,first_index) 
        use utils_module, only: DBL_FMT,write_txt_unit,logzero
        use model_module, only: model
        use settings_module, only: program_settings

        implicit none


        type(program_settings), intent(in) :: settings
        type(model),            intent(in) :: M
        double precision, dimension(M%nDims+2,settings%nmax_posterior) :: posterior_array
        double precision                           :: evidence
        integer :: first_index

        integer :: next_index

        integer :: i_err

        
        ! Open the .txt file, note the presence of iostat prevents program
        ! termination during this write
        open(write_txt_unit,file=trim(settings%file_root) // '.txt' , action='write', iostat=i_err) 

        ! Initialise the index at the first index
        next_index = first_index
        do while(next_index/=-1 .or. posterior_array(3,next_index)>logzero) 
            ! Write this point to file
            write(write_txt_unit,'(<M%nDims+2>E<DBL_FMT(1)>.<DBL_FMT(2)>)')&
                exp(posterior_array(3,next_index)-evidence),posterior_array(4:,next_index)

            ! Advance to the next point in the linked list
            next_index= nint(posterior_array(2,next_index))
        end do

        close(write_txt_unit)

    end subroutine write_posterior_file




end module read_write_module
