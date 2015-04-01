    module settings
    use MiscUtils
    use FileUtils
    use StringUtils
    use MpiUtils
    use IniObjects
    use ParamNames
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit,error_unit
    implicit none

#ifdef SINGLE
    integer, parameter :: mcp= KIND(1.0)
#else
    integer, parameter :: mcp= KIND(1.d0)
#endif
#ifdef MPI
#ifdef SINGLE
    integer, parameter :: MPI_real_mcp = MPI_REAL
#else
    integer, parameter :: MPI_real_mcp = MPI_DOUBLE_PRECISION
#endif
#endif
    integer,parameter :: time_dp = KIND(1.d0)


    double precision, parameter :: pi=3.14159265358979323846264338328d0, &
        twopi=2*pi, fourpi=4*pi
    double precision, parameter :: root2 = 1.41421356237309504880168872421d0, sqrt2 = root2
    double precision, parameter :: log2 = 0.693147180559945309417232121458d0

    real, parameter :: pi_r = 3.141592653, twopi_r = 2*pi_r, fourpi_r = twopi_r*2

    real(mcp), parameter :: const_c = 2.99792458e8_mcp

    logical :: use_fast_slow = .false.

    character(LEN=*), parameter :: CosmoMC_Version = 'Feb2015'

    character(LEN=:), allocatable :: chisq_label

    Type, extends(TIniFile) :: TNameKeyIniFile
        !Allowing things like param[name] =
    contains
    procedure :: NamedKey => TNameKeyIniFile_NamedKey
    end type TNameKeyIniFile


    Type, extends(TNameKeyIniFile) :: TSettingIni
    contains
    procedure :: FailStop => TSettingIni_FailStop
    procedure :: ReadFilename => TSettingIni_ReadFilename
    procedure :: ReadRelativeFilename => TSettingIni_ReadRelativeFilename
    procedure :: ReplaceDirs => TSettingIni_ReplaceDirs
    procedure :: TagValuesForName => TSettingIni_TagValuesForName
    procedure :: SettingValuesForTagName => TSettingIni_SettingValuesForTagName
    end type

    real(mcp) :: AccuracyLevel = 1
    !Set to >1 to use theory calculation etc on higher accuracy settings.
    !Does not affect MCMC (except making it all slower)

    logical :: flush_write = .true.

    logical :: new_chains = .true.

    integer, parameter :: max_likelihood_functions = 50

    integer, parameter :: max_data_params = 200
    integer, parameter :: max_theory_params = 50
    integer, parameter :: max_num_params = max_theory_params + max_data_params

    !Set to false if using a slow likelihood function so no there's point is treating
    !'fast' parameters differently (in fact, doing so will make performance worse)

    integer, parameter :: sampling_metropolis = 1, sampling_slice = 2, sampling_fastslice =3, &
        sampling_slowgrid = 4,  sampling_multicanonical = 5,  sampling_wang_landau = 6, &
        sampling_fast_dragging = 7

    integer :: sampling_method = sampling_metropolis

    !For fast dragging method, baseline number of intermediate drag steps
    real(mcp) :: dragging_steps = 3._mcp

    !The rest are set up automatically
    logical  ::  generic_mcmc= .false.
    !set to true to not call CAMB, etc.
    !write GenericLikelihoodFunction in calclike.f90

    character(LEN=:), allocatable :: DataDir, LocalDir

    Type(TSettingIni), save :: CustomParams

    logical :: stop_on_error = .true. !whether to stop with error, or continue ignoring point

    integer :: num_theory_params, index_data, index_semislow=-1 !set later depending on datasets and theory parameterization

    integer, dimension(:), allocatable :: params_used
    integer num_params, num_params_used, num_data_params

    integer :: num_threads = 0
    integer :: instance = 0
    integer :: MPIchains = 1, MPIrank = 0
    real(time_dp) :: MPIRun_start_time

    logical :: checkpoint = .false.


    integer :: output_lines = 0
    Type(TTextFile), save :: ChainOutFile = TTextFile(RealFormat='(*(E16.7))')
    Type(TTextFile), save :: LogFile

    integer :: Feedback = 0

    real(mcp), parameter :: logZero = 1e30_mcp
    character (LEN =1024) :: FileChangeIni = '', FileChangeIniAll = ''
    character(LEN=:), allocatable :: baseroot, rootname

    integer, parameter :: stdout = output_unit

    contains

    subroutine InitializeGlobalSettingDefaults
    character, parameter :: backslash = char(92)

    DataDir='data/'
    LocalDir='./'
    chisq_label = backslash//'chi^2_{'//backslash//'rm %s}'

    end subroutine InitializeGlobalSettingDefaults

    subroutine DoStop(S, abort)
    character(LEN=*), intent(in), optional :: S
    logical, intent(in), optional :: abort
    logical wantbort
#ifdef MPI
    integer ierror
    real(time_dp) runTime
#endif

    call ChainOutFile%Close()

    if (present(abort)) then
        wantbort = abort
    else
        wantbort = .false.
    end if

    if (present(S) .and. (wantbort .or. MPIRank==0)) write (*,*) trim(S)
#ifdef MPI
    runTime = MPI_WTime() - MPIrun_Start_Time
    if (Feedback > 0 .and. MPIRank==0) &
        write (*,'("Total time:  ",I0,"  (",F10.5," hours  )")')nint(runTime),runTime/(60*60)
    ierror =0
    if (wantbort) then
        !Abort all in case other continuing chains want to communicate with us
        !in the case when max number of samples is reached
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
    else
        call mpi_finalize(ierror)
    end if
#endif

#ifdef DECONLY
    pause
#endif
    stop
    end subroutine DoStop

    subroutine TSettingIni_FailStop(this)
    class(TSettingIni) :: this

    call MpiStop()

    end subroutine TSettingIni_FailStop

    function TSettingIni_ReplaceDirs(this,inname, ADir) result(filename)
    class(TSettingIni) :: this
    character(LEN=*) :: inname
    character(LEN=:), allocatable :: filename
    character(LEN=*), optional, intent(in) :: ADir

    filename = inname
    call StringReplace('%DATASETDIR%',PresentDefault(DataDir,ADir),filename)
    call StringReplace('%LOCALDIR%',LocalDir,filename)
    if (allocated(this%Original_filename)) &
        & call StringReplace('%THISDIR%',File%ExtractPath(this%Original_filename),filename)

    end function TSettingIni_ReplaceDirs

    function TSettingIni_ReadRelativeFilename(this,key, ADir, NotFoundFail) result (filename)
    class(TSettingIni) :: this
    character(LEN=*), intent(in) :: Key
    character(LEN=*), optional, intent(in) :: ADir
    character(LEN=:), allocatable :: filename
    logical, optional :: NotFoundFail

    filename = this%ReadFileName(key,Adir,NotFoundFail,.true.)

    end function TSettingIni_ReadRelativeFilename

    function TSettingIni_ReadFilename(this,key, ADir, NotFoundFail, relative) result (filename)
    class(TSettingIni) :: this
    character(LEN=*), intent(in) :: Key
    character(LEN=*), optional, intent(in) :: ADir
    character(LEN=:), allocatable :: filename
    logical, optional :: NotFoundFail, relative
    integer i

    filename = this%Read_String(key, NotFoundFail)
    if (filename=='') return

    filename = this%ReplaceDirs(filename, ADir)

    do i=1, CustomParams%Count
        call StringReplace('%'//CustomParams%Name(i)//'%',&
            trim(this%ReplaceDirs(CustomParams%Value(i), ADir)) ,filename)
    end do
    if (DefaultFalse(relative) .and. .not. File%IsFullPath(filename)) then
        filename =  File%ExtractPath(this%Original_filename)//filename
    end if

    end function TSettingIni_ReadFilename

    subroutine TSettingIni_TagValuesForName(this, name, OutList, filename)
    !Reads all entries of the form "name[tag] = value", storing tag=value in OutList
    class(TSettingIni) :: this
    character(LEN=*), intent(in) :: name
    character(LEN=:), allocatable :: tag, value
    class(TNameValueList) :: OutList
    integer i, ix
    character(LEN=:), pointer :: KeyName
    logical, intent(in), optional :: filename

    call OutList%Clear()
    tag = trim(name) //'['

    do i=1, this%Count
        KeyName=>this%Name(i)
        if (StringStarts(KeyName, tag) .and. this%Value(i)/='') then
            ix = index(KeyName, ']')
            if (ix /= len(keyName))  call this%Error('Error reading tagged key', name)
            if (index(KeyName,',')/=0) cycle
            value =this%Read_String(KeyName)
            if (DefaultFalse(filename) .and. value/='') value = this%ReplaceDirs(value)
            call OutList%Add(KeyName(len(tag)+1:len(KeyName)-1), value)
            !!Note: Use Read_String so read values are stored
        end if
    end do

    end subroutine TSettingIni_TagValuesForName


    subroutine TSettingIni_SettingValuesForTagName(this, name, tag, OutList, filename)
    !Reads all entries of the form "name[tag,setting] = value", storing setting=value in OutList
    class(TSettingIni) :: this
    character(LEN=*), intent(in) :: name, tag
    character(LEN=:), allocatable :: value, stem
    class(TNameValueList) :: OutList
    integer i, ix
    character(LEN=:), pointer :: KeyName
    logical, intent(in), optional :: filename

    call OutList%Clear()
    stem = trim(name) //'['//trim(tag)//','

    do i=1, this%Count
        KeyName=>this%Name(i)
        if (StringStarts(KeyName, stem)) then
            ix = index(KeyName, ']')
            if (ix /= len(KeyName)) call this%Error('Error reading tagged key setting', name)
            value =this%Read_String(KeyName)
            if (DefaultFalse(filename) .and. value/='') value = this%ReplaceDirs(value)
            call OutList%Add(KeyName(len(stem)+1:len(KeyName)-1), value)
        end if
    end do

    end subroutine TSettingIni_SettingValuesForTagName


    function TNameKeyIniFile_NamedKey(this, Key, Name) result(NamedKey)
    class(TNameKeyIniFile) :: this
    character(LEN=*), intent(in) :: Key, Name
    character(LEN=:), allocatable :: NamedKey

    NamedKey = trim(key)//'['//trim(Name)//']'

    end function TNameKeyIniFile_NamedKey


    subroutine CheckParamChangeF(F)
    character(LEN=*), intent(in) ::  F
    logical bad, doexit
    Type(TSettingIni) :: Ini

    if (F /= '') then
        call Ini%Open(F, bad, .false.)
        if (bad) return
        doexit = (Ini%Read_Int('exit',0) == 1)
        FeedBack = Ini%Read_Int('feedback',Feedback)
        num_threads = Ini%Read_Int('num_threads',num_threads)
        call Ini%Close()
        if (F== FileChangeIni) call File%Delete(FileChangeini)
        if (doexit) call MpiStop('exit requested')
    end if

    end subroutine CheckParamChangeF

    subroutine CheckParamChange

    call CheckParamChangeF(FileChangeIni)
    if (FileChangeIni/=FileChangeIniAll) call CheckParamChangeF(FileChangeIniAll)

    end subroutine CheckParamChange

    subroutine Timer(Msg, start)
    character(LEN=*), intent(in), optional :: Msg
    real(time_dp), save :: timer_start
    real(time_dp), optional :: start
    real(time_dp) T

    if (present(start)) then
        T=start
    else
        T=timer_start
    end if

    if (present(Msg)) then
        write (*,*) trim(Msg)//': ', TimerTime() - T
    end if
    if (.not. present(start)) timer_start= TimerTime()

    end subroutine Timer

    subroutine DoAbort(S)
    character(LEN=*), intent(in), optional :: S
#ifdef MPI
    integer ierror
#endif
    if (present(S)) write (*,*) trim(S)
#ifdef MPI
    call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
#endif

#ifdef DECONLY
    pause
#endif
    stop
    end subroutine DoAbort


    end module settings
