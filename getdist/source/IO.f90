
    module IO
    !Module to wrap output functions
    !Replace file e.g. to store output in a pipline database (indexed by filename)
    !Valid handles must be /=0
    use settings
    use FileUtils
    use ParamNames
    implicit none

    contains

    subroutine IO_WriteProposeMatrix(pmat, prop_mat, comment)
    use MatrixUtils
    real(mcp) pmat(:,:)
    character(LEN=*), intent(in) :: prop_mat
    character(LEN=*), optional, intent(in) :: comment

    if (present(comment)) then
        call Matrix_write(prop_mat,pmat,.true.,commentline=comment)
    else
        call Matrix_write(prop_mat,pmat,.true.)
    endif

    end subroutine IO_WriteProposeMatrix

    subroutine IO_ReadProposeMatrix(NameMapping,pmat, prop_mat)
    class(TParamNames) :: NameMapping
    real(mcp) pmat(:,:)
    character(LEN=*), intent(in) :: prop_mat
    real(mcp), allocatable :: tmpMat(:,:)
    integer i,x,y
    integer status
    character(LEN=:), allocatable :: InLine, comment
    integer num, cov_params(max_num_params)
    Type(TTextFile) :: F

    call F%Open(prop_mat)
    comment=''
    if (F%ReadLineSkipEmptyAndComments(InLine, comment=comment) .and. comment/='') then
        !Have paramnames to identify
        num=-1
        call NameMapping%ReadIndices(comment, cov_params, num, unknown_value=0)
        allocate(tmpMat(num,num))
        pmat=0
        y=0
        do
            read(InLine,*,iostat=status) tmpMat(:,y+1)
            if (status/=0) exit
            y=y+1
            if (y==num) exit
            if (.not. F%ReadLineSkipEmptyAndComments(InLine)) exit
        end do
        if (y/=num) call mpiStop('ReadProposeMatrix: wrong number of rows/columns in .covmat')
        call F%Close()
        do y=1,num
            if (cov_params(y)/=0) then
                do x=1,num
                    if (cov_params(x)/=0) pmat(cov_params(x),cov_params(y)) = tmpMat(x,y)
                end do
            end if
        end do
        deallocate(tmpMat)
        return

    end if
    call F%Close()

    i=File%TxtNumberColumns(InLine)
    if (i==num_params) then
        call File%ReadTextMatrix(prop_mat,pmat, num_params, num_params)
    else if (i==num_theory_params) then
        allocate(tmpMat(num_theory_params,num_theory_params))
        call File%ReadTextMatrix(prop_mat,tmpmat, num_theory_params, num_theory_params)
        pmat=0
        pmat(1:num_theory_params,1:num_theory_params) = tmpMat
        deallocate(tmpMat)
    else
        call MpiStop('Propose matrix the wrong size: '//trim(prop_mat))
    end if

    end subroutine IO_ReadProposeMatrix


    subroutine IO_OutputChainRow(F, mult, like, values)
    class(TFileStream) :: F
    real(mcp) mult, like, values(:)

    call F%Write( [mult, like, values])

    if (flush_write) call F%Flush()

    end subroutine IO_OutputChainRow

    function IO_ReadChainRow(F, mult, like, values, params_used, chainOK, samples_chains) result(OK)
    !Returns OK=false if end of file or if not enough values on each line, otherwise OK = true
    !Returns chainOK = false if bad line or NaN, chainOK=false and OK=true for NaN (continue reading)
    logical OK
    class(TTextFile) :: F
    real(mcp), intent(out) :: mult, like, values(:)
    integer, intent(in) :: params_used(:)
    logical, optional, intent(out) :: ChainOK
    logical, optional, intent(in) :: samples_chains
    logical samples_are_chains
    character(LEN=:), allocatable :: InLine
    real(mcp) invals(size(params_used))
    integer status

    if (present(samples_chains)) then
        samples_are_chains=samples_chains
    else
        samples_are_chains = .true.
    endif

    if (present(ChainOK)) chainOK = .true.

    OK = F%ReadLine(InLine)
    if (.not. OK) return
    if (SCAN (InLine, 'N') /=0) then
        if (present(ChainOK)) chainOK = .false.
        return
    end if

    if (samples_are_chains) then
        read(InLine, *, iostat=status) mult, like, invals
    else
        mult=1
        like=1
        read(InLine, *, iostat=status) invals
    end if
    OK = status==0
    if (present(ChainOK)) chainOK = OK
    if (OK) values(params_used) =invals

    end function IO_ReadChainRow

    subroutine IO_ReadLastChainParams(name, mult, like, values, params_used)
    character(LEN=*), intent(in) :: name
    real(mcp), intent(out) :: mult, like, values(:)
    integer, intent(in) :: params_used(:)
    character(LEN=:), allocatable :: InLine

    InLine = File%LastLine(name)
    read(InLine, *) mult, like, values(params_used)

    end subroutine IO_ReadLastChainParams

    subroutine IO_OutputParamNames(Names, fname, indices, add_derived)
    class(TParamNames) :: Names
    character(len=*), intent(in) :: fname
    integer, intent(in), optional :: indices(:)
    logical, intent(in), optional :: add_derived

    call Names%WriteFile(trim(fname)//'.paramnames', indices, add_derived)

    end subroutine IO_OutputParamNames

    subroutine  IO_WriteBounds(Names, fname, limmin,limmax, limbot,limtop, indices)
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: fname
    real(mcp), intent(in) :: limmin(:), limmax(:)
    logical, intent(in) :: limbot(:), limtop(:)
    integer, intent(in) :: indices(:)
    integer i,ix
    Type(TTextFile) :: F
    character(LEN=17) :: lim1,lim2

    call F%CreateFile(fname)
    do i=1, size(indices)
        ix = indices(i)
        if (limbot(ix) .or. limtop(ix)) then
            if (limbot(ix)) then
                write(lim1, F%RealFormat) limmin(ix)
            else
                lim1='    N'
            end if
            if (limtop(ix)) then
                write(lim2, F%RealFormat) limmax(ix)
            else
                lim2='    N'
            end if
            write(F%unit,'(1A22,2A17)') Names%NameOrNumber(ix-2), lim1, lim2
        end if
    end do
    call F%Close()

    end subroutine IO_WriteBounds


    subroutine IO_ReadParamNames(Names, in_root, prior_ranges)
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: in_root
    character(LEN=ParamNames_maxlen) name
    character(LEN=:), allocatable :: infile
    real(mcp) :: prior_ranges(:,:), minmax(2)
    integer ix, status
    Type(TTextFile) :: F

    prior_ranges=0
    infile = trim(in_root) // '.paramnames'
    if (File%Exists(infile)) then
        call Names%Init(infile)
        infile = trim(in_root) // '.ranges'
        if (File%Exists(infile)) then
            call F%Open(infile)
            do
                read(F%unit, *, iostat=status) name, minmax
                if (status/=0) exit
                ix = Names%index(name)
                if (ix/=-1) prior_ranges(:,ix) = minmax
            end do
            call F%Close()
        end if
    end if

    end subroutine IO_ReadParamNames

    function IO_ReadChainRows(in_root, chain_ix,chain_num, ignorerows, nrows, &
    ncols,max_rows,coldata,samples_are_chains) result(OK)
    !OK = false if chain not found or not enough samples
    character(LEN=*), intent(in) :: in_root
    integer,intent(in) :: chain_ix, chain_num
    integer, intent(in) :: max_rows, ignorerows
    integer, intent(in) :: ncols
    real(KIND(1.d0)), intent(inout) :: coldata(ncols,0:max_rows) !(col_index, row_index)
    logical, intent(in) :: samples_are_chains
    integer, intent(inout) :: nrows
    logical OK, ChainOK
    real(mcp) invars(1:ncols)
    integer status
    integer row_start
    integer, allocatable :: indices(:)
    character(LEN=:), allocatable :: infile
    integer i
    Type(TTextFile) :: F

    row_start=nrows
    if (chain_num == 0) then
        infile = trim(in_root) // '.txt'
    else
        infile = trim(in_root) //'_'//IntToStr(chain_ix)// '.txt'
    end if

    write (*,*) 'reading ' // trim(infile)

    call F%Open(infile, status=status)
    if (status/=0) then
        write (*,'(" chain ",1I4," missing")') chain_ix
        OK = .false.
        return
    end if

    if (ignorerows >=1) then
        if (.not. F%SkipLines(ignorerows)) then
            OK = .false.
            return
        end if
    end if

    allocate(indices(ncols-2))
    indices=[ (I, I=3, ncols) ] ! [3:ncols]
    OK = .true.
    do
        if (.not. IO_ReadChainRow(F, invars(1), invars(2), &
        invars,indices,chainOK,samples_are_chains)) then
            if (.not. chainOK) then
                write (*,*) 'error reading line ', nrows -row_start + ignorerows ,' - skipping to next row'
                cycle
            endif
            return
        else
            if (.not. chainOK) then
                write (*,*) 'WARNING: skipping line with probable NaN'
                cycle
            end if
        end if

        coldata(1:ncols, nrows) = invars(1:ncols)
        nrows = nrows + 1
        if (nrows > max_rows) stop 'need to increase max_rows'
    end do

    end function IO_ReadChainRows


    subroutine IO_OutputMargeStats(Names, froot,num_vars,num_contours, contours,contours_str, &
    cont_lines, colix, mean, sddev, has_limits_bot, has_limits_top, labels)
    use ParamNames
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: froot
    integer, intent(in) :: num_vars, num_contours
    logical,intent(in) :: has_limits_bot(:,:),has_limits_top(:,:)
    real(mcp), intent(in) :: mean(*), sddev(*), contours(*), cont_lines(:,:,:)
    character(LEN=*), intent(in) :: contours_str
    integer,intent(in) :: colix(*)
    character(LEN=128) labels(*), tag, nameFormat
    integer i,j
    character(LEN=*), parameter :: txtFormat = '(1A15)'
    Type(TTextFile) F

    j = max(9,Names%MaxNameLen())
    nameFormat = concat('(1A',j+1,')')

    call F%CreateFile(trim(froot)//'.margestats')
    F%RealFormat = '(*(1E15.7))'
    call F%Write('Marginalized limits: ' // contours_str)
    call F%NewLine()
    call F%WriteLeftAligned(nameFormat,'parameter')
    call F%WriteInLine('  ')
    call F%WriteLeftAligned(txtFormat,'mean')
    call F%WriteLeftAligned(txtFormat,'sddev')
    do j=1, num_contours
        call F%WriteLeftAligned(txtFormat,concat('lower',j))
        call F%WriteLeftAligned(txtFormat,concat('upper',j))
        call F%WriteLeftAligned('(1A7)',concat('limit',j))
    end do
    call F%NewLine()

    do j=1, num_vars
        call F%WriteLeftAligned(nameFormat,Names%NameOrNumber(colix(j)-2))
        call F%WriteInLineItems(mean(j), sddev(j))
        do i=1, num_contours
            call F%WriteInLineItems(cont_lines(j,1,i),cont_lines(j,2,i)) !don't send array section to avoid ifort 14 bug
            if (has_limits_bot(i,colix(j)).and. has_limits_top(i,colix(j))) then
                tag='none'
            elseif (has_limits_bot(i,colix(j))) then
                tag= '>'
            elseif (has_limits_top(i,colix(j))) then
                tag = '<'
            else
                tag = 'two'
            end if
            call F%WriteInLine('  '//tag, '(1A7)')
        end do
        call F%WriteTrim('   '//labels(colix(j)))
    end do
    call F%Close()

    end subroutine IO_OutputMargeStats


    end module IO
