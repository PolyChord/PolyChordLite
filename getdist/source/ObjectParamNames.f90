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
    procedure :: Add => ParamNames_Add
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
    procedure :: WriteFile => ParamNames_WriteFile
    end Type TParamNames

    public TParamNames, ParamNames_maxlen
    contains

    function IsWhiteSpace(C)
    character, intent(in) :: C
    logical IsWhiteSpace

    IsWhiteSpace = (C==' ') .or. (C==char(9))

    end function IsWhiteSpace


    function ParamNames_ParseLine(Names,InLine,n) result(res)
    class(TParamNames) :: Names
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
        Names%comment(n) = label(pos+1: len_trim(label))
        label = label(1:pos-1)
    else
        Names%comment(n) = ''
    endif
    pos = scan(label,char(9))
    if (pos/=0) label = label(1:pos-1)
    name = trim(adjustl(name))
    alen = len_trim(name)
    if (name(alen:alen)=='*') then
        name(alen:alen)=' '
        Names%is_derived(n) = .true.
    else
        Names%is_derived(n) = .false.
    end if
    Names%name(n) = trim(name)
    Names%label(n) = trim(label)

    end function ParamNames_ParseLine

    subroutine ParamNames_Alloc(Names,n)
    class(TParamNames) :: Names
    integer,intent(in) :: n

    call Names%Dealloc()
    allocate(Names%name(n))
    allocate(Names%label(n))
    allocate(Names%comment(n))
    allocate(Names%is_derived(n))
    Names%nnames = n
    Names%is_derived = .false.
    Names%num_MCMC = 0
    Names%num_derived = 0
    Names%name = ''
    Names%comment=''
    Names%label=''

    end subroutine ParamNames_Alloc

    subroutine ParamNames_dealloc(Names)
    class(TParamNames) :: Names
    if (allocated(Names%name)) &
    deallocate(Names%name,Names%label,Names%comment,Names%is_derived)

    end subroutine ParamNames_dealloc

    subroutine ParamNames_Init(Names, filename)
    class(TParamNames) :: Names
    character(Len=*), intent(in) :: filename
    integer n
    character(LEN=:), allocatable :: InLine
    Type(TTextFile) :: F

    call F%Open(filename)
    n = F%Lines()
    call Names%Alloc(n)

    n=0
    do while (F%ReadLine(InLine))
        if (InLine=='') cycle
        n=n+1
        if (.not. Names%ParseLine(InLine,n)) then
            call MpiStop(concat('ParamNames_Init: error parsing line: ',n))
        end if
    end do
    call F%Close()

    Names%nnames = n
    Names%num_derived = count(Names%is_derived)
    Names%num_MCMC = Names%nnames - Names%num_derived

    end subroutine ParamNames_Init

    subroutine ParamNames_AssignItem(Names, Names2,n,i)
    class(TParamNames), target :: Names, Names2
    integer n, i

    Names%name(n) = Names2%name(i)
    Names%label(n) = Names2%label(i)
    Names%comment(n) = Names2%comment(i)
    Names%is_derived(n) = Names2%is_derived(i)

    end subroutine ParamNames_AssignItem

    subroutine ParamNames_Add(Names, Names2, check_duplicates)
    class(TParamNames), target :: Names, Names2
    logical, intent(in), optional :: check_duplicates
    integer n,i, newold, derived
    class(TParamNames),pointer :: P, NamesOrig

    allocate(NamesOrig, source = Names)

    n=0
    do i=1, names2%nnames
        if (NamesOrig%index(Names2%name(i))==-1) then
            n=n+1
        else
            if (PresentDefault(.false., check_duplicates)) &
            & call MpiStop('ParamNames_Add: Duplicate name tag - '//trim(Names2%name(i)))
        end if
    end do
    if (n==0) return

    call Names%Alloc(NamesOrig%nnames + n)
    Names%nnames = 0
    do derived=0,1
        P=> NamesOrig
        do newold=0,1
            do i=1, P%nnames
                if (Names%index(P%name(i))==-1) then
                    if (derived==0 .and. .not. P%is_derived(i) .or.derived==1 .and. P%is_derived(i) ) then
                        Names%nnames = Names%nnames + 1
                        call Names%AssignItem(P, Names%nnames , i)
                    end if
                end if
            end do
            P=> Names2
        enddo
    end do
    if (Names%nnames/= NamesOrig%nnames + n) stop 'ParamNames_Add: duplicate parameters?'

    Names%num_derived = count(Names%is_derived)
    Names%num_MCMC= Names%nnames-Names%num_derived

    call NamesOrig%Dealloc()
    deallocate(NamesOrig)

    end subroutine ParamNames_Add

    subroutine ParamNames_SetLabels(Names,filename)
    class(TParamNames) :: Names
    Type(TParamNames) :: LabNames
    character(Len=*), intent(in) :: filename
    integer i,ix

    call LabNames%init(filename)
    do i=1, LabNames%nnames
        ix = Names%index(LabNames%name(i))
        if (ix/=-1) then
            Names%label(ix) = LabNames%label(i)
        end if
    end do

    end subroutine ParamNames_SetLabels

    function ParamNames_index(Names,name) result(ix)
    class(TParamNames) :: Names
    character(len=*), intent(in) :: name
    integer ix,i

    do i=1,Names%nnames
        if (Names%name(i) == name) then
            ix = i
            return
        end if
    end do
    ix = -1

    end function ParamNames_index


    function ParamNames_label(Names,name) result(lab)
    class(TParamNames) :: Names
    character(len=*), intent(in) :: name
    character(len = :), allocatable :: lab
    integer ix

    ix = Names%index(name)
    if (ix>0) then
        lab = trim(Names%label(ix))
    else
        lab = ''
    end if

    end function ParamNames_label

    function ParamNames_name(Names,ix) result(name)
    class(TParamNames) :: Names
    character(len=:), allocatable  :: name
    integer, intent(in) :: ix

    if (ix <= Names%nnames) then
        name = trim(Names%name(ix))
    else
        name = ''
    end if

    end function ParamNames_name


    subroutine ParamNames_ReadIndices(Names,InLine, params, num, unknown_value)
    class(TParamNames) :: Names
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
        ix = Names%index(part)
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
    call MpiStop('ParamNames: Not enough names or numbers - '//trim(InLine))

    end subroutine ParamNames_ReadIndices

    function ParamNames_AsString(Names, i, want_comment) result(line)
    class(TParamNames) :: Names
    integer, intent(in) :: i
    logical ,intent(in), optional :: want_comment
    character(LEN=:), allocatable :: Line
    logical wantCom

    if (present(want_comment)) then
        wantCom = want_comment
    else
        wantCom = .false.
    end if

    if (i> Names%nnames) call MpiStop('ParamNames_AsString: index out of range')
    Line = trim(Names%name(i))
    if (Names%is_derived(i)) Line = Line // '*'
    Line =  Line//char(9)//trim(Names%label(i))
    if (wantCom .and. Names%comment(i)/='') then
        Line = Line//char(9)//'#'//trim(Names%comment(i))
    end if

    end function ParamNames_AsString

    subroutine ParamNames_WriteFile(Names, fname, indices, add_derived)
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: fname
    integer, intent(in), optional :: indices(:)
    logical, intent(in), optional :: add_derived
    integer i
    Type(TTextFile) :: F

    call F%CreateFile(fname)
    if (present(indices)) then
        do i=1, size(indices)
            call F%Write(Names%AsString(indices(i)))
        end do
        if (present(add_derived)) then
            if (add_derived) then
                do i=1,Names%num_derived
                    call F%Write(Names%AsString(Names%num_mcmc+i))
                end do
            end if
        end if
    else
        do i=1, Names%nnames
            call F%Write(Names%AsString(i))
        end do
    end if

    call F%Close()

    end subroutine ParamNames_WriteFile


    function ParamNames_NameOrNumber(Names,ix) result(name)
    class(TParamNames) :: Names
    character(len=:), allocatable  :: name
    integer, intent(in) :: ix

    name = trim(Names%name(ix))
    if (name == '') name = trim(IntToStr(ix))

    end function ParamNames_NameOrNumber

    function ParamNames_MaxNameLen(Names) result(alen)
    class(TParamNames) :: Names
    integer alen, i

    alen = 0
    do i=1, Names%nnames
        alen = max(alen, len_trim(Names%NameOrNumber(i)))
    end do

    end function ParamNames_MaxNameLen


    function ParamNames_ReadIniForParam(Names,Ini,Key, param) result(input)
    ! read Key[name] or Keyn where n is the parameter number
    use IniObjects
    class(TParamNames) :: Names
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    integer, intent(in) :: param
    character(LEN=:), allocatable :: input

    input = ''
    if (Names%nnames>0) then
        input = Ini%Read_String(trim(key)//'['//trim(Names%name(param))//']')
    end if
    if (input=='') then
        input = Ini%Read_String(trim(Key)//trim(IntToStr(param)))
    end if

    end function ParamNames_ReadIniForParam

    function ParamNames_HasReadIniForParam(Names,Ini,Key, param) result(B)
    ! read Key[name] or Keyn where n is the parameter number
    use IniObjects
    class(TParamNames) :: Names
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    integer, intent(in) :: param
    logical B

    B = .false.
    if (Names%nnames>0) then
        B = Ini%HasKey(trim(key)//'['//trim(Names%name(param))//']')
    end if
    if (.not. B) then
        B = Ini%HasKey(trim(Key)//trim(IntToStr(param)))
    end if

    end function ParamNames_HasReadIniForParam



    end module ParamNames
