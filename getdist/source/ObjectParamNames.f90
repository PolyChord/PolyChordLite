    module ParamNames
    use MpiUtils
    use FileUtils
    use StringUtils
    use MiscUtils
    implicit none
    private
    integer, parameter :: ParamNames_maxlen = 128

    Type TParamNames
        integer :: nnames =0
        integer :: num_MCMC = 0
        integer :: num_derived = 0
        character(LEN=ParamNames_maxlen), dimension(:), allocatable ::  name
        character(LEN=ParamNames_maxlen), dimension(:), allocatable ::  label
        character(LEN=ParamNames_maxlen), dimension(:), allocatable ::  comment
        logical, dimension(:), allocatable ::  is_derived
    contains
    procedure :: AddNames => ParamNames_Add
    procedure :: AddFile => ParamNames_AddFile
    procedure :: Alloc => ParamNames_Alloc
    procedure :: AssignItem => ParamNames_AssignItem
    procedure :: AsString => ParamNames_AsString
    procedure :: Dealloc => ParamNames_Dealloc
    procedure :: HasReadIniForParam => ParamNames_HasReadIniForParam
    procedure :: Index => ParamNames_Index
    procedure :: Init => ParamNames_Init
    procedure :: LabelForName => ParamNames_label
    procedure :: MaxNameLen => ParamNames_MaxNameLen
    procedure :: NameAtIndex => ParamNames_name
    procedure :: NameOrNumber => ParamNames_NameOrNumber
    procedure :: ParseLine => ParamNames_ParseLine
    procedure :: ReadIndices => ParamNames_ReadIndices
    procedure :: ReadIniForParam => ParamNames_ReadIniForParam
    procedure :: SetLabels => ParamNames_SetLabels
    procedure :: SetUnnamed => ParamNames_SetUnnamed
    procedure :: WriteFile => ParamNames_WriteFile
    generic :: Add => AddNames, AddFile
    end Type TParamNames

    public TParamNames, ParamNames_maxlen
    contains

    function ParamNames_ParseLine(this,InLine,n) result(res)
    class(TParamNames) :: this
    character(LEN=*) :: InLine
    character(LEN=ParamNames_maxlen*3) :: name, label
    integer n
    logical res
    integer pos, alen,status

    alen = len_trim(InLIne)
    pos =1
    do while (pos < alen .and. IsWhiteSpace(InLIne(pos:pos)))
        pos = pos+1
    end do
    read(InLine(pos:), *, iostat=status) name
    res= (status==0)
    if (.not. res) return
    pos = pos + len_trim(name)
    do while (pos < alen .and. IsWhiteSpace(InLIne(pos:pos)))
        pos = pos+1
    end do
    label = trim(adjustl(InLine(pos:alen)))
    pos = scan(label,'#')
    if (pos/=0) then
        this%comment(n) = label(pos+1: len_trim(label))
        label = label(1:pos-1)
    else
        this%comment(n) = ''
    endif
    pos = scan(label,char(9))
    if (pos/=0) label = label(1:pos-1)
    name = trim(adjustl(name))
    alen = len_trim(name)
    if (name(alen:alen)=='*') then
        name(alen:alen)=' '
        this%is_derived(n) = .true.
    else
        this%is_derived(n) = .false.
    end if
    this%name(n) = trim(name)
    this%label(n) = trim(label)

    end function ParamNames_ParseLine

    subroutine ParamNames_Alloc(this,n)
    class(TParamNames) :: this
    integer,intent(in) :: n

    call this%Dealloc()
    allocate(this%name(n))
    allocate(this%label(n))
    allocate(this%comment(n))
    allocate(this%is_derived(n))
    this%nnames = n
    this%is_derived = .false.
    this%num_MCMC = 0
    this%num_derived = 0
    this%name = ''
    this%comment=''
    this%label=''

    end subroutine ParamNames_Alloc

    subroutine ParamNames_dealloc(this)
    class(TParamNames) :: this
    if (allocated(this%name)) &
        deallocate(this%name,this%label,this%comment,this%is_derived)

    end subroutine ParamNames_dealloc

    subroutine ParamNames_Init(this, filename)
    class(TParamNames) :: this
    character(Len=*), intent(in) :: filename
    integer n
    character(LEN=:), allocatable :: InLine
    Type(TTextFile) :: F

    call F%Open(filename)
    n = F%Lines()
    call this%Alloc(n)

    n=0
    do while (F%ReadLine(InLine))
        if (InLine=='') cycle
        n=n+1
        if (.not. this%ParseLine(InLine,n)) then
            call MpiStop(concat('ParamNames_Init: error parsing line: ',n))
        end if
    end do
    call F%Close()

    this%nnames = n
    this%num_derived = count(this%is_derived)
    this%num_MCMC = this%nnames - this%num_derived

    end subroutine ParamNames_Init

    subroutine ParamNames_AssignItem(this, Names2,n,i)
    class(TParamNames), target :: this, Names2
    integer n, i

    this%name(n) = Names2%name(i)
    this%label(n) = Names2%label(i)
    this%comment(n) = Names2%comment(i)
    this%is_derived(n) = Names2%is_derived(i)

    end subroutine ParamNames_AssignItem

    subroutine ParamNames_AddFile(this, filename, check_duplicates)
    class(TParamNames), target :: this
    character(Len=*), intent(in) :: filename
    logical, intent(in), optional :: check_duplicates
    Type(TParamNames) :: P

    call P%Init(filename)
    call this%Add(P,check_duplicates)

    end subroutine ParamNames_AddFile


    subroutine ParamNames_Add(this, Names, check_duplicates)
    class(TParamNames), target :: this, Names
    logical, intent(in), optional :: check_duplicates
    integer n,i, newold, derived
    class(TParamNames),pointer :: P, NamesOrig

    allocate(NamesOrig, source = this)

    n=0
    do i=1, Names%nnames
        if (NamesOrig%index(Names%name(i))==-1) then
            n=n+1
        else
            if (DefaultFalse(check_duplicates)) &
                & call MpiStop('ParamNames_Add: Duplicate name tag - '//trim(Names%name(i)))
        end if
    end do
    if (n==0) return

    call this%Alloc(NamesOrig%nnames + n)
    this%nnames = 0
    do derived=0,1
        P=> NamesOrig
        do newold=0,1
            do i=1, P%nnames
                if (this%index(P%name(i))==-1) then
                    if (derived==0 .and. .not. P%is_derived(i) .or.derived==1 .and. P%is_derived(i) ) then
                        this%nnames = this%nnames + 1
                        call this%AssignItem(P, this%nnames , i)
                    end if
                end if
            end do
            P=> Names
        enddo
    end do
    if (this%nnames/= NamesOrig%nnames + n) stop 'ParamNames_Add: duplicate parameters?'

    this%num_derived = count(this%is_derived)
    this%num_MCMC= this%nnames-this%num_derived

    call NamesOrig%Dealloc()
    deallocate(NamesOrig)

    end subroutine ParamNames_Add

    subroutine ParamNames_SetLabels(this,filename)
    class(TParamNames) :: this
    Type(TParamNames) :: LabNames
    character(Len=*), intent(in) :: filename
    integer i,ix

    call LabNames%init(filename)
    do i=1, LabNames%nnames
        ix = this%index(LabNames%name(i))
        if (ix/=-1) then
            this%label(ix) = LabNames%label(i)
        end if
    end do

    end subroutine ParamNames_SetLabels

    subroutine ParamNames_SetUnnamed(this, nparams, prefix)
    class(TParamNames) :: this
    integer, intent(in) :: nparams
    character(LEN=*), intent(in), optional :: prefix
    integer i

    call this%Alloc(nparams)
    this%nnames = nparams
    this%num_derived = 0
    this%num_MCMC = nparams
    do i=1, nparams
        this%name(i) = trim(PresentDefault('param',prefix))// IntToStr(i)
    end do

    end subroutine ParamNames_SetUnnamed

    function ParamNames_index(this,name) result(ix)
    class(TParamNames) :: this
    character(len=*), intent(in) :: name
    integer ix,i

    do i=1,this%nnames
        if (this%name(i) == name) then
            ix = i
            return
        end if
    end do
    ix = -1

    end function ParamNames_index


    function ParamNames_label(this,name) result(lab)
    class(TParamNames) :: this
    character(len=*), intent(in) :: name
    character(len = :), allocatable :: lab
    integer ix

    ix = this%index(name)
    if (ix>0) then
        lab = trim(this%label(ix))
    else
        lab = ''
    end if

    end function ParamNames_label

    function ParamNames_name(this,ix) result(name)
    class(TParamNames) :: this
    character(len=:), allocatable  :: name
    integer, intent(in) :: ix

    if (ix <= this%nnames) then
        name = trim(this%name(ix))
    else
        name = ''
    end if

    end function ParamNames_name


    subroutine ParamNames_ReadIndices(this,InLine, params, num, unknown_value)
    class(TParamNames) :: this
    character(LEN=*), intent(in) :: InLine
    integer, intent(out) :: params(*)
    integer, intent(in), optional :: unknown_value
    integer  :: num, status
    character(LEN=ParamNames_maxlen) part
    integer param,len,ix, pos, max_num, outparam, outvalue
    integer, parameter :: unknown_num = 1024
    character(LEN=:), allocatable :: skips

    skips=''
    if (num==0) return
    len = len_trim(InLine)
    pos = 1
    if (num==-1) then
        max_num = unknown_num
    else
        max_num = num
    end if
    outparam=0
    do param = 1, max_num
        do while (pos <= len)
            if (IsWhiteSpace(InLine(pos:pos))) then
                pos = pos+1
            else
                exit
            endif
        end do
        read(InLine(pos:), *, iostat=status) part
        if (status/=0) exit
        pos = pos + len_trim(part)
        ix = this%index(part)
        if (ix>0) then
            outvalue = ix
        else
            if (verify(trim(part),'0123456789') /= 0) then
                if (present(unknown_value)) then
                    skips = skips //' '//trim(part)
                    if (unknown_value/=-1) then
                        outvalue = unknown_value
                    else
                        cycle
                    end if
                else
                    call MpiStop( 'ParamNames: Unknown parameter name '//trim(part))
                end if
            else
                read(part,*) outvalue
            end if
        end if
        outparam = outparam +1
        if (max_num == unknown_num) num = outparam
        params(outparam) = outvalue
    end do
    if (status==0) return

    if (skips/='') write(*,'(a)') ' skipped unused params:'// skips
    if (max_num==unknown_num) return
    call MpiStop('ParamNames: Not enough this or numbers - '//trim(InLine))

    end subroutine ParamNames_ReadIndices

    function ParamNames_AsString(this, i, want_comment) result(line)
    class(TParamNames) :: this
    integer, intent(in) :: i
    logical ,intent(in), optional :: want_comment
    character(LEN=:), allocatable :: Line
    logical wantCom

    if (present(want_comment)) then
        wantCom = want_comment
    else
        wantCom = .false.
    end if

    if (i> this%nnames) call MpiStop('ParamNames_AsString: index out of range')
    Line = trim(this%name(i))
    if (this%is_derived(i)) Line = Line // '*'
    Line =  Line//char(9)//trim(this%label(i))
    if (wantCom .and. this%comment(i)/='') then
        Line = Line//char(9)//'#'//trim(this%comment(i))
    end if

    end function ParamNames_AsString

    subroutine ParamNames_WriteFile(this, fname, indices, add_derived)
    class(TParamNames) :: this
    character(LEN=*), intent(in) :: fname
    integer, intent(in), optional :: indices(:)
    logical, intent(in), optional :: add_derived
    integer i
    Type(TTextFile) :: F

    call F%CreateFile(fname)
    if (present(indices)) then
        do i=1, size(indices)
            call F%Write(this%AsString(indices(i)))
        end do
        if (present(add_derived)) then
            if (add_derived) then
                do i=1,this%num_derived
                    call F%Write(this%AsString(this%num_mcmc+i))
                end do
            end if
        end if
    else
        do i=1, this%nnames
            call F%Write(this%AsString(i))
        end do
    end if

    call F%Close()

    end subroutine ParamNames_WriteFile


    function ParamNames_NameOrNumber(this,ix, tag_derived) result(name)
    class(TParamNames) :: this
    character(len=:), allocatable  :: name
    logical, intent(in), optional :: tag_derived
    integer, intent(in) :: ix

    name = trim(this%name(ix))
    if (name == '') name = trim(IntToStr(ix))
    if (DefaultFalse(tag_derived) .and. this%is_derived(ix)) name = name //'*'

    end function ParamNames_NameOrNumber

    function ParamNames_MaxNameLen(this) result(alen)
    class(TParamNames) :: this
    integer alen, i

    alen = 0
    do i=1, this%nnames
        alen = max(alen, len_trim(this%NameOrNumber(i)))
    end do

    end function ParamNames_MaxNameLen


    function ParamNames_ReadIniForParam(this,Ini,Key, param) result(input)
    ! read Key[name] or Keyn where n is the parameter number
    use IniObjects
    class(TParamNames) :: this
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    integer, intent(in) :: param
    character(LEN=:), allocatable :: input

    input = ''
    if (this%nnames>0) then
        input = Ini%Read_String(trim(key)//'['//trim(this%name(param))//']')
    end if
    if (input=='') then
        input = Ini%Read_String(trim(Key)//trim(IntToStr(param)))
    end if

    end function ParamNames_ReadIniForParam

    function ParamNames_HasReadIniForParam(this,Ini,Key, param) result(B)
    ! read Key[name] or Keyn where n is the parameter number
    use IniObjects
    class(TParamNames) :: this
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    integer, intent(in) :: param
    logical B

    B = .false.
    if (this%nnames>0) then
        B = Ini%HasKey(trim(key)//'['//trim(this%name(param))//']')
    end if
    if (.not. B) then
        B = Ini%HasKey(trim(Key)//trim(IntToStr(param)))
    end if

    end function ParamNames_HasReadIniForParam



    end module ParamNames
