module ini_module
    use utils_module, only: dp


  implicit none

  character(len=2), parameter :: comment='#!'
  character(len=2), parameter :: equals='=:'
  character(len=1), parameter :: section_start='['
  character(len=1), parameter :: section_end=']'

contains


    subroutine initialise_program(settings,priors,params,derived_params)
        use priors_module, only: create_priors,prior
        use settings_module,   only: program_settings,initialise_settings
        use params_module, only: param_type
        use read_write_module, only: write_paramnames_file
        implicit none
        
        type(program_settings),intent(inout)            :: settings  !> Program settings
        type(prior),dimension(:),allocatable,intent(out):: priors

        type(param_type),dimension(:),allocatable,intent(in) :: params         ! Parameter array
        type(param_type),dimension(:),allocatable,intent(in) :: derived_params ! Derived parameter array

        settings%nDims = size(params)
        settings%nDerived = size(derived_params)

        if(settings%write_paramnames.and.(settings%equals.or.settings%posteriors)) call write_paramnames_file(settings,params,derived_params)

        call create_priors(priors,params,settings)

        ! Calculate all of the rest of the settings
        call initialise_settings(settings)   

    end subroutine initialise_program





    subroutine read_params(file_name,settings,params,derived_params)
        use settings_module,   only: program_settings
        use params_module, only: param_type
        implicit none
        
        character(len=*)                    :: file_name !> The name of the file
        type(program_settings),intent(out)  :: settings  !> Program settings

        type(param_type),dimension(:),allocatable,intent(out) :: params         ! Parameter array
        type(param_type),dimension(:),allocatable,intent(out) :: derived_params ! Derived parameter array


        settings%nlive              = get_integer(file_name,'nlive')
        settings%num_repeats        = get_integer(file_name,'num_repeats')
        settings%nprior             = get_integer(file_name,'nprior',-1)
        settings%do_clustering      = get_logical(file_name,'do_clustering',.false.)
        settings%feedback           = get_integer(file_name,'feedback',1)
        settings%precision_criterion= get_double(file_name,'precision_criterion',1d-3)
        settings%logzero            = get_double(file_name,'logzero',-1d30)
        settings%max_ndead          = get_integer(file_name,'max_ndead',-1)
        settings%boost_posterior    = get_double(file_name,'boost_posterior',0d0)
        settings%posteriors         = get_logical(file_name,'posteriors',.false.)
        settings%equals             = get_logical(file_name,'equals',.false.)
        settings%cluster_posteriors = get_logical(file_name,'cluster_posteriors',.false.)
        settings%write_resume       = get_logical(file_name,'write_resume',.false.)
        settings%write_paramnames   = get_logical(file_name,'write_paramnames',.false.)
        settings%read_resume        = get_logical(file_name,'read_resume',.false.)
        settings%write_stats        = get_logical(file_name,'write_stats',.true.)
        settings%write_live         = get_logical(file_name,'write_live',.false.)
        settings%write_dead         = get_logical(file_name,'write_dead',.true.)
        settings%write_prior        = get_logical(file_name,'write_prior',.false.)
        settings%compression_factor = get_double(file_name,'compression_factor',exp(-1d0))
        settings%base_dir           = get_string(file_name,'base_dir','chains')
        settings%file_root          = get_string(file_name,'file_root','test')
        settings%seed               = get_integer(file_name,'seed', -1)

        settings%generate_from_seed = get_logical(file_name,'generate_from_seed',.false.)
        settings%nprior_repeat      = get_integer(file_name,'nprior_repeat',-1)
        call get_doubles(file_name,'seed_point',settings%seed_point)

        call get_doubles(file_name,'grade_frac',settings%grade_frac)
        call get_integers(file_name,'nlives',settings%nlives)
        call get_doubles(file_name,'loglikes',settings%loglikes)

        call get_params(file_name,params,derived_params)  

        settings%nDims         = size(params)
        settings%nDerived      = size(derived_params)

    end subroutine read_params


    function default_params(nDims,root,latex)
        use priors_module, only: uniform_type
        use params_module, only: param_type,add_parameter
        implicit none
        integer, intent(in) :: nDims
        type(param_type),dimension(:),allocatable :: default_params

        integer :: i
        character(len=*) :: root
        character(len=*) :: latex
        character(len=500) :: i_string
        character(len=500) :: paramname
        character(len=500) :: latex_paramname

        allocate(default_params(0))

        do i=1,nDims
            write(i_string, '(I20)') i
            paramname = trim(adjustl(root))//trim(adjustl(i_string))
            latex_paramname = trim(adjustl(latex))//'_{'//trim(adjustl(i_string))//'}'
            call add_parameter(default_params,trim(adjustl(paramname)),trim(adjustl(latex_paramname)),1,uniform_type,1,[ 0d0 , 1d0 ])
        end do

    end function default_params



    function get_string(file_name,key_word,dflt,ith)
        use utils_module,  only: STR_LENGTH,params_unit
        use abort_module,  only: halt_program
        implicit none
        character(len=*),intent(in)  :: file_name !> The name of the file to search in
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        character(len=*),intent(in),optional  :: dflt  !> keyword to search for
        integer,intent(in), optional :: ith       !> Get the ith instance of this string

        character(len=STR_LENGTH) :: get_string  ! string following keyword

        character(len=STR_LENGTH) :: keyword    !> keyword to search for
        character(len=STR_LENGTH) :: filename   ! The fortran readable filename

        character(len=STR_LENGTH)                :: line_buffer     ! Line buffer

        integer :: io_stat   ! check to see if we've reached the end of the file

        integer :: i_equals ! placement of equals signs
        integer :: counter



        ! Convert filename to something fortran can read
        write(filename,'(A)') file_name
        open(unit=params_unit,file=trim(filename),iostat=io_stat)
        if(io_stat/=0) call halt_program('ini error: '//trim(file_name)//' does not exist')
        write(keyword,'(A)') key_word

        get_string = ''

        counter=1

        do 
            ! Read in the next line
            read(params_unit,'(A)',iostat=io_stat) line_buffer
            ! Exit if io_stat is less than zero
            if (io_stat < 0) exit

            ! Skip any comment lines
            if( scan(line_buffer,comment) /=0 ) cycle

            ! Search for equals signs
            i_equals = scan(line_buffer,equals)
            if(i_equals==0) cycle

            ! check to see if this matches our keyword
            if( trim(adjustl(line_buffer(:i_equals-1))) == trim(adjustl(keyword)) ) then


                if(present(ith)) then
                    
                    if(counter==ith) then
                        get_string = adjustl(line_buffer(i_equals+1:))
                        exit
                    else
                        counter=counter+1
                    end if
                else
                    get_string = adjustl(line_buffer(i_equals+1:))
                    exit
                end if

            end if

        end do

        if(trim(get_string)==''.and. present(dflt)) get_string=trim(dflt)

        ! close the file
        close(params_unit)


    end function get_string

    function get_double(file_name,key_word,dflt)
        use utils_module,  only: STR_LENGTH
        use abort_module,  only: halt_program
        character(len=*),intent(in)  :: file_name !> The name of the file to search in
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        real(dp),intent(in),optional :: dflt

        character(len=STR_LENGTH) :: string  ! string following keyword
        real(dp) :: get_double  ! double following keyword

        string = get_string(file_name,key_word)
        if(trim(string)/='') then
            read(string,*) get_double
        else if(present(dflt)) then
            get_double = dflt
        else
            call halt_program('ini error: no keyword '//trim(key_word))
        end if

    end function get_double

    subroutine get_doubles(file_name,key_word,doubles)
        use utils_module,  only: STR_LENGTH
        use array_module,  only: reallocate
        character(len=*),intent(in)  :: file_name !> The name of the file to search in
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        real(dp), intent(out), allocatable, dimension(:) :: doubles

        character(len=STR_LENGTH) :: string  ! string following keyword

        integer :: i

        ! Allocate a zero size array
        if(allocated(doubles)) deallocate(doubles)
        allocate(doubles(0))

        ! Get the string
        string = get_string(file_name,key_word)

        do while( trim(string)/='' )
            i=size(doubles)+1
            call reallocate(doubles,i)
            read(string,*) doubles(i)
            call next_element(string,' ')
        end do

    end subroutine get_doubles

    subroutine get_integers(file_name,key_word,integers)
        use utils_module,  only: STR_LENGTH
        use array_module,  only: reallocate
        character(len=*),intent(in)  :: file_name !> The name of the file to search in
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        integer, intent(out), allocatable, dimension(:) :: integers

        character(len=STR_LENGTH) :: string  ! string following keyword

        integer :: i

        ! Allocate a zero size array
        if(allocated(integers)) deallocate(integers)
        allocate(integers(0))

        ! Get the string
        string = get_string(file_name,key_word)

        do while( trim(string)/='' )
            i=size(integers)+1
            call reallocate(integers,i)
            read(string,*) integers(i)
            call next_element(string,' ')
        end do

    end subroutine get_integers

    function get_integer(file_name,key_word,dflt)
        use utils_module,  only: STR_LENGTH
        use abort_module,  only: halt_program
        character(len=*),intent(in)  :: file_name !> The name of the file to search in
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        integer,intent(in),optional :: dflt

        character(len=STR_LENGTH) :: string  ! string following keyword
        integer :: get_integer  ! integer following keyword

        if(present(dflt)) get_integer=dflt

        string = get_string(file_name,key_word)
        if(trim(string)/='') then
            read(string,*) get_integer
        else if(present(dflt)) then
            get_integer = dflt
        else
            call halt_program('ini error: no keyword '//trim(key_word))
        end if

    end function get_integer

    function get_logical(file_name,key_word,dflt)
        use utils_module,  only: STR_LENGTH
        use abort_module,  only: halt_program
        character(len=*),intent(in)  :: file_name !> The name of the file to search in
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        logical,intent(in),optional :: dflt

        character(len=STR_LENGTH) :: string  ! string following keyword
        logical :: get_logical  ! logical following keyword

        string = get_string(file_name,key_word)
        if(trim(string)/='') then
            read(string,*) get_logical
        else if(present(dflt)) then
            get_logical = dflt
        else
            call halt_program('ini error: no keyword '//trim(key_word))
        end if

    end function get_logical















    subroutine get_params(file_name,params,derived_params) 
        use priors_module, only: prior_type_from_string,unknown_type
        use utils_module,  only: STR_LENGTH
        use abort_module,  only: halt_program
        use params_module, only: param_type,add_parameter
        implicit none
        
        character(len=*),intent(in)               :: file_name !> The name of the file
        type(param_type),dimension(:),allocatable,intent(out) :: params         ! Parameter array
        type(param_type),dimension(:),allocatable,intent(out) :: derived_params ! Derived parameter array

        
        character(len=STR_LENGTH)                :: line_buffer     ! Line buffer
        character(len=STR_LENGTH)                :: paramname       ! parameter name
        character(len=STR_LENGTH)                :: latex           ! latex name
        integer                                  :: speed           ! parameter speed
        character(len=STR_LENGTH)                :: prior_type_str  ! prior type string
        integer                                  :: prior_type      ! prior type integer
        integer                                  :: prior_block     ! prior block
        real(dp),allocatable,dimension(:):: prior_params    ! prior parameters

        integer                                  :: sub_cluster     ! whether or not to do sub clustering on this parameter
        character(1),parameter                   :: sc='*'          ! indicator for sub clustering


        integer :: i_param

        allocate(params(0),derived_params(0))

        i_param = 1
        do 
            line_buffer = get_string(file_name,'P',ith=i_param)
            if(trim(line_buffer)=='') exit
            i_param=i_param+1

            !1) Parameter name
            read(line_buffer,*) paramname !read in string
            sub_cluster = index(paramname,sc)
            if(sub_cluster>0) then
                paramname=trim(paramname(:sub_cluster-1))
            else
                paramname=trim(paramname)     !trim string
            end if


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
            if(prior_type==unknown_type) call halt_program('get_priors error: Unknown prior type for parameter '//trim(paramname)) 

            !5) Prior block
            call next_element(line_buffer,'|') ! advance
            read(line_buffer,*) prior_block    ! read in integer


            call next_element(line_buffer,'|')              ! advance
            call get_prior_params(prior_params,line_buffer) ! get the prior params

            ! Add this parameter to the array
            if(sub_cluster==0) then
                call add_parameter(params,paramname,latex,speed,prior_type,prior_block,prior_params) 
            else
                call add_parameter(params,paramname,latex,speed,prior_type,prior_block,prior_params,.true.) 
            end if

        end do

        i_param = 1
        do 
            line_buffer = get_string(file_name,'D',ith=i_param)
            if(line_buffer=='') exit
            i_param=i_param+1

            !1) Parameter name
            read(line_buffer,*) paramname !read in string
            paramname=trim(paramname)     !trim string


            !2) Latex name
            call next_element(line_buffer,'|') ! advance
            read(line_buffer,*) latex          ! read in string
            latex=trim(latex)                  ! trim string

            deallocate(prior_params)
            allocate(prior_params(0))

            ! Add this parameter to the array
            call add_parameter(derived_params,paramname,latex,0,0,0,prior_params) 
        end do

    end subroutine get_params



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
        real(dp),allocatable,dimension(:),intent(out) :: prior_params      ! prior parameters
        real(dp),allocatable,dimension(:)             :: temp_params       ! prior parameters

        integer :: i

        ! deallocate it if it's already allocated
        if(allocated(prior_params)) deallocate(prior_params)
        if(allocated(temp_params)) deallocate(temp_params)

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

    end subroutine get_prior_params








    function split_string(string,separator)
        implicit none
        character(len=*),intent(in) :: string
        character(1),intent(in) :: separator
        character(len=len(string)), dimension(:), allocatable :: split_string

        character(len=len(string)), dimension(:), allocatable :: temp_split_string

        integer :: i,j,k


        allocate(split_string(0))
        temp_split_string = split_string
        i = 0
        k = 0
        do 
            j = index(string(i+1:),separator)
            if(trim(adjustl(string(i+1:i+j)))=='') exit

            k = k+1 
            deallocate(split_string)
            allocate(split_string(k))
            split_string(:k-1) = temp_split_string
            deallocate(temp_split_string)
            split_string(k) = trim(adjustl(string(i+1:i+j))) 
            temp_split_string = split_string

            i = i+j
        end do

    end function split_string



 
end module ini_module
