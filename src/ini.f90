module ini_module


  implicit none

  character(len=2), parameter :: comment='#!'

contains

    subroutine read_priors(file_name,priors) 
        use priors_module, only: prior,prior_type_from_string,unknown_type
        use utils_module,  only: STR_LENGTH,params_unit
        use abort_module,  only: halt_program
        use params_module, only: param_type,add_parameter
        implicit none
        
        character(len=*)                                   :: file_name !> The name of the file
        type(prior), dimension(:), allocatable,intent(out) :: priors    !> The array of priors to be returned

        
        character(len=STR_LENGTH)                :: filename        ! The file name to read priors from
        character(len=STR_LENGTH)                :: line_buffer     ! Line buffer
        character(len=STR_LENGTH)                :: paramname       ! parameter name
        character(len=STR_LENGTH)                :: latex           ! latex name
        integer                                  :: speed           ! parameter speed
        character(len=STR_LENGTH)                :: prior_type_str  ! prior type string
        integer                                  :: prior_type      ! prior type integer
        integer                                  :: prior_block     ! prior block
        double precision,allocatable,dimension(:):: prior_params    ! prior parameters

        type(param_type),dimension(:),allocatable :: params         ! Parameter array


        integer :: i_comment=1 ! Comment index
        integer :: i_break     ! Params index
        integer :: io_stat=0   ! Error checker

        write(filename,'(A)') file_name
        open(unit=params_unit,file=trim(filename))

        ! Skip the leading comment lines
        do while(i_comment/=0.and.io_stat==0)
            read(params_unit,'(A)',iostat=io_stat) line_buffer
            i_comment = scan(line_buffer,comment)
        end do

        do while(io_stat==0) 

            !1) Parameter name
            read(line_buffer,*) paramname !read in string
            paramname=trim(paramname)     !trim string


            !2) Latex name
            call next_element(line_buffer,'|') ! advance
            read(line_buffer,*) latex          ! read in string
            latex=trim(latex)                  ! trim string


            !3) Parameter speed/grade
            call next_element(line_buffer,'|') ! advance
            read(line_buffer,*) speed          ! read in integer


            !4) Prior type
            call next_element(line_buffer,'|')                     ! advance
            read(line_buffer,*) prior_type_str                     ! read in string
            prior_type = prior_type_from_string(prior_type_str)    ! convert to integer

            ! Halt if we don't know this prior type
            if(prior_type==unknown_type) call halt_program('read_priors error: Unknown prior type for parameter '//trim(paramname)) 

            !5) Prior block
            call next_element(line_buffer,'|') ! advance
            read(line_buffer,*) prior_block    ! read in integer


            call next_element(line_buffer,'|')              ! advance
            call get_prior_params(prior_params,line_buffer) ! get the prior params

            ! Add this parameter to the array
            call add_parameter(params,paramname,latex,speed,prior_type,prior_block,prior_params) 




            ! Read the next line
            read(params_unit,'(A)',iostat=io_stat) line_buffer
        end do




        close(params_unit)
        stop


    end subroutine read_priors

    subroutine next_element(line_buffer,delimiter) 
        use utils_module,  only: STR_LENGTH
        implicit none
        character(len=STR_LENGTH),intent(inout)  :: line_buffer ! Line buffer
        character :: delimiter

        line_buffer = trim(line_buffer(scan(line_buffer,delimiter)+1:)) ! Find the next element
    end subroutine next_element

    subroutine get_prior_params(prior_params,line_buffer)
        use utils_module,  only: STR_LENGTH
        implicit none
        character(len=STR_LENGTH),intent(inout)               :: line_buffer ! Line buffer
        double precision,allocatable,dimension(:),intent(out) :: prior_params      ! prior parameters
        double precision,allocatable,dimension(:)             :: temp_params       ! prior parameters

        integer :: i

        ! deallocate it if it's already allocated
        if(allocated(prior_params)) deallocate(prior_params)

        ! Trim leading spaces
        line_buffer = adjustl(line_buffer)

        ! If it's an empty string, then allocate a zero-length array
        allocate(prior_params(0),temp_params(0)) 
        i=0

        do while(trim(line_buffer)/='') 

            deallocate(temp_params)
            allocate(temp_params(i))
            temp_params = prior_params

            deallocate(prior_params)
            allocate(prior_params(i+1))
            prior_params(:i) = temp_params

            read(line_buffer,*) prior_params(i+1)

            call next_element(line_buffer,' ')
            line_buffer = adjustl(line_buffer)
            i=i+1

        end do

        write(*,*) prior_params



    end subroutine get_prior_params


 
end module ini_module
