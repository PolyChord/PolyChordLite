    module Interpolation
    use FileUtils
    use MpiUtils
    implicit none

#ifdef SINGLE
    integer, parameter :: sp_acc = KIND(1.0)
#else
    integer, parameter :: sp_acc = KIND(1.d0)
#endif
    real(sp_acc), parameter :: SPLINE_DANGLE=1.e30_sp_acc
    integer, parameter :: GI = sp_acc

    Type :: TInterpolator
        private
        logical :: initialized =.false.
    contains
    procedure :: FirstUse => TInterpolator_FirstUse
    procedure :: Error => TInterpolator_Error
    end Type TInterpolator

    Type, extends(TInterpolator) :: TCubicSpline
        integer n
        real(sp_acc), dimension(:), allocatable :: X
        real(sp_acc), dimension(:), allocatable :: F,  ddF
    contains
    procedure, private :: TCubicSpline_Init
    procedure, private :: TCubicSpline_InitInt
    procedure, private :: TCubicSpline_Value
    procedure, private :: TCubicSpline_ArrayValue
    procedure, private :: TCubicSpline_IntRangeValue
    procedure :: InitFromFile => TCubicSpline_InitFromFile
    procedure :: InitInterp => TCubicSpline_InitInterp
    procedure :: Derivative => TCubicSpline_Derivative
    generic :: Value => TCubicSpline_Value
    generic :: Array => TCubicSpline_ArrayValue, TCubicSpline_IntRangeValue
    generic :: Init => TCubicSpline_Init, TCubicSpline_InitInt
    FINAL :: TCubicSpline_Free !not needed in standard, just for compiler bugs
    end Type

    Type, extends(TInterpolator) :: TInterpGrid2D
        !      ALGORITHM 760, COLLECTED ALGORITHMS FROM ACM.
        !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
        !      VOL. 22, NO. 3, September, 1996, P.  357--361.
        REAL(GI), private, allocatable :: wk(:,:,:)
        REAL(GI), allocatable :: x(:), y(:)
        REAL(GI), allocatable :: z(:,:)
        integer nx, ny
    contains
    procedure :: Init => TInterpGrid2D_Init
    procedure :: InitFromFile => TInterpGrid2D_InitFromFile
    procedure :: Value => TInterpGrid2D_Value !one point
    procedure :: Values => TInterpGrid2D_Values !array of points
    procedure :: Error => TInterpGrid2D_error
    procedure, private :: InitInterp => TInterpGrid2D_InitInterp
    FINAL :: TInterpGrid2D_Free
    end Type TInterpGrid2D


    public TCubicSpline, TInterpGrid2D, SPLINE_DANGLE

    contains

    subroutine TInterpolator_FirstUse(W)
    class(TInterpolator) W

    W%Initialized = .true.
    call W%error('TInterpolator not initialized')

    end subroutine TInterpolator_FirstUse

    subroutine TInterpolator_error(W,S)
    class(TInterpolator):: W
    character(LEN=*), intent(in) :: S

    call MpiStop('Interpolation error: '//trim(S))

    end subroutine TInterpolator_error


    subroutine TCubicSpline_Init(W, Xarr,  values, n, End1, End2 )
    class(TCubicSpline) :: W
    real(sp_acc), intent(in) :: Xarr(:), values(:)
    integer, intent(in), optional :: n
    real(sp_acc), intent(in), optional :: End1, End2

    if (present(n)) then
        W%n = n
    else
        W%n = size(Xarr)
    end if

    allocate(W%F(W%n))
    W%F = values(1:W%n)
    allocate(W%X(W%n))
    W%X = Xarr(1:W%n)
    call W%InitInterp(End1, End2)

    end subroutine TCubicSpline_Init

    subroutine TCubicSpline_InitInt(W, Xarr,  values, n, End1, End2 )
    class(TCubicSpline) :: W
    integer, intent(in) :: XArr(:)
    real(sp_acc), intent(in) :: values(:)
    real(sp_acc), allocatable :: XReal(:)
    integer, intent(in), optional :: n
    real(sp_acc), intent(in), optional :: End1, End2

    allocate(XReal(size(XArr)))
    XReal =  XArr
    call W%Init(XReal, values, n, End1, End2)

    end subroutine TCubicSpline_InitInt

    subroutine TCubicSpline_InitInterp(W,End1,End2)
    class(TCubicSpline):: W
    real(sp_acc), intent(in), optional :: End1, End2
    real(sp_acc) :: e1,e2

    if (present(End1)) then
        e1 = End1
    else
        e1 = SPLINE_DANGLE
    end if

    if (present(End2)) then
        e2 = End2
    else
        e2 = SPLINE_DANGLE
    end if

    allocate(W%ddF(W%n))
    call spline(W%X,W%F,W%n,e1,e2,W%ddF)
    W%Initialized = .true.

    end subroutine TCubicSpline_InitInterp

    subroutine TCubicSpline_InitFromFile(W, Filename, xcol, ycol)
    class(TCubicSpline):: W
    character(LEN=*), intent(in) :: Filename
    integer, intent(in), optional :: xcol, ycol
    integer :: ixcol=1,iycol=2
    integer :: nx,ny
    integer :: parse, status
    real(sp_acc), allocatable :: tmp(:)
    character(LEN=:), allocatable :: InLine
    Type(TTextFile) :: F
    if (present(xcol)) ixcol=xcol
    if (present(ycol)) iycol=ycol

    allocate(tmp(max(ixcol,iycol)))

    call F%Open(Filename)

    do parse=1,2
        ny=0
        nx=0
        do
            if (.not. F%ReadLineSkipEmptyAndComments(InLine)) exit

            read(InLine,*, iostat=status) tmp
            if (status/=0) call W%Error('Error reading line: '//trim(InLine))

            nx=nx+1
            if (parse==2) then
                W%X(nx) = tmp(ixcol)
                W%F(nx) = tmp(iycol)
            endif
        end do

        if (parse==2) exit

        if (nx<2) call W%Error('not enough values to interpolate')
        W%n = nx
        allocate(W%X(W%n))
        allocate(W%F(W%n))
        status=0
        call F%Rewind()
    end do

    call F%Close()

    call W%InitInterp()

    end subroutine TCubicSpline_InitFromFile


    subroutine TCubicSpline_Free(W)
    Type(TCubicSpline) :: W

    if (allocated(W%X)) then
        deallocate(W%X)
        deallocate(W%F)
        deallocate(W%ddF)
    end if
    W%Initialized = .false.

    end subroutine TCubicSpline_Free

    function TCubicSpline_Value(W, x, error )
    class(TCubicSpline) :: W
    real(sp_acc) :: TCubicSpline_Value
    real(sp_acc), intent(in) :: x
    integer, intent(inout), optional :: error !initialize to zero outside, changed if bad
    integer llo,lhi
    real(sp_acc) a0,b0,ho

    if (.not. W%Initialized) call W%FirstUse

    if (x< W%X(1) .or. x> W%X(W%n)) then
        if (present(error)) then
            error = -1
            TCubicSpline_Value=0
            return
        else
            write (*,*) 'TCubicSpline_Value: out of range ', x
            stop
        end if
    end if

    llo=1
    do while (W%X(llo+1) < x)  !could do binary search here if large
        llo = llo + 1
    end do

    lhi=llo+1
    ho=W%X(lhi)-W%X(llo)
    a0=(W%X(lhi)-x)/ho
    b0=(x-W%X(llo))/ho
    TCubicSpline_Value = a0*W%F(llo)+ b0*W%F(lhi)+((a0**3-a0)* W%ddF(llo) +(b0**3-b0)*W%ddF(lhi))*ho**2/6

    end function TCubicSpline_Value

    subroutine TCubicSpline_ArrayValue(W, x, y, error )
    !Get array of values y(x), assuming x is monotonically increasing
    class(TCubicSpline) :: W
    real(sp_acc), intent(in) :: x(1:)
    real(sp_acc), intent(out) :: y(1:)
    integer, intent(inout), optional :: error !initialize to zero outside, changed if bad
    integer llo,lhi, i
    real(sp_acc) a0,b0,ho

    if (.not. W%Initialized) call W%FirstUse

    if (x(1)< W%X(1) .or. x(size(x))> W%X(W%n)) then
        if (present(error)) then
            error = -1
            return
        else
            write (*,*) 'TCubicSpline_Value: out of range ', x
            stop
        end if
    end if

    llo=1
    do i=1, size(x)
        do while (W%X(llo+1) < x(i))  !could do binary search here if large
            llo = llo + 1
        end do

        lhi=llo+1
        ho=W%X(lhi)-W%X(llo)
        a0=(W%X(lhi)-x(i))/ho
        b0=(x(i)-W%X(llo))/ho
        y(i) = a0*W%F(llo)+ b0*W%F(lhi)+((a0**3-a0)* W%ddF(llo) +(b0**3-b0)*W%ddF(lhi))*ho**2/6
    end do

    end subroutine TCubicSpline_ArrayValue

    subroutine TCubicSpline_IntRangeValue(W, xmin, xmax, y, error )
    !Get array of values y(x), assuming x is monotonically increasing
    class(TCubicSpline) :: W
    integer, intent(in) :: xmin, xmax
    real(sp_acc), intent(out) :: y(xmin:)
    integer, intent(inout), optional :: error !initialize to zero outside, changed if bad
    integer llo,lhi, x
    real(sp_acc) a0,b0,ho

    if (.not. W%Initialized) call W%FirstUse

    if (xmin< W%X(1) .or. xmax> W%X(W%n)) then
        if (present(error)) then
            error = -1
            return
        else
            write (*,*) 'TCubicSpline_Value: out of range ', x
            stop
        end if
    end if

    llo=1
    do x=xmin, xmax
        do while (W%X(llo+1) < x)  !could do binary search here if large
            llo = llo + 1
        end do

        lhi=llo+1
        ho=W%X(lhi)-W%X(llo)
        a0=(W%X(lhi)-x)/ho
        b0=(x-W%X(llo))/ho
        y(x) = a0*W%F(llo)+ b0*W%F(lhi)+((a0**3-a0)* W%ddF(llo) +(b0**3-b0)*W%ddF(lhi))*ho**2/6
    end do

    end subroutine TCubicSpline_IntRangeValue

    ! Get derivative of spline
    function TCubicSpline_Derivative(W, x, error )
    class(TCubicSpline) :: W
    real(sp_acc) :: TCubicSpline_Derivative
    real(sp_acc), intent(in) :: x
    integer, intent(inout), optional :: error !initialize to zero outside, changed if bad
    integer llo,lhi
    real(sp_acc) a0,b0,ho,dely

    if (.not. W%Initialized) call W%FirstUse

    if (x< W%X(1) .or. x> W%X(W%n)) then
        if (present(error)) then
            error = -1
            TCubicSpline_Derivative=0
            return
        else
            write (*,*) 'TCubicSpline_Derivative: out of range ', x
            stop
        end if
    end if

    llo=1
    do while (W%X(llo+1) < x)  !could do binary search here if large
        llo = llo + 1
    end do

    lhi=llo+1
    ho=W%X(lhi)-W%X(llo)
    a0=(W%X(lhi)-x)/ho
    b0=(x-W%X(llo))/ho
    dely=W%F(lhi)-W%F(llo)
    TCubicSpline_Derivative = dely/ho-(3.d0*a0**2-1.0d0)*ho*W%ddF(llo)/6.0d0 + &
    (3.d0*b0**2-1.0d0)*ho*W%ddF(lhi)/6.0d0

    end function TCubicSpline_Derivative


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! calculates array of second derivatives used by cubic spline
    ! interpolation. y2 is array of second derivatives, yp1 and ypn are first
    ! derivatives at end points.

    !Thanks Martin Reinecke
    subroutine spline(x,y,n,d11,d1n,d2)
    integer, intent(in) :: n
    real(sp_acc), intent(in) :: x(n), y(n), d11, d1n
    real(sp_acc), intent(out) :: d2(n)
    integer i
    real(sp_acc) xp,qn,sig,un,xxdiv,u(n-1),d1l,d1r

    d1r= (y(2)-y(1))/(x(2)-x(1))
    if (d11==SPLINE_DANGLE) then
        d2(1)=0._sp_acc
        u(1)=0._sp_acc
    else
        d2(1)=-0.5_sp_acc
        u(1)=(3._sp_acc/(x(2)-x(1)))*(d1r-d11)
    endif

    do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1._sp_acc/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1._sp_acc/(sig*d2(i-1)+2._sp_acc)

        d2(i)=(sig-1._sp_acc)*xp

        u(i)=(6._sp_acc*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
    end do
    d1l=d1r

    if (d1n==SPLINE_DANGLE) then
        qn=0._sp_acc
        un=0._sp_acc
    else
        qn=0.5_sp_acc
        un=(3._sp_acc/(x(n)-x(n-1)))*(d1n-d1l)
    endif

    d2(n)=(un-qn*u(n-1))/(qn*d2(n-1)+1._sp_acc)
    do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
    end do
    end subroutine spline


    subroutine TInterpGrid2D_InitInterp(W)
    class(TInterpGrid2D):: W

    allocate(W%Wk(3,W%NX,W%NY))
    CALL rgpd3p(W%nx, W%ny, W%x, W%y, W%z, W%wk)
    W%Initialized = .true.

    end subroutine TInterpGrid2D_InitInterp

    !F2003 wrappers by AL, 2013
    subroutine TInterpGrid2D_Init(W, x, y, z)
    class(TInterpGrid2D):: W
    REAL(GI), INTENT(IN)      :: x(:)
    REAL(GI), INTENT(IN)      :: y(:)
    REAL(GI), INTENT(IN)      :: z(:,:)

    W%nx = size(x)
    W%ny = size(y)
    allocate(W%x(W%nx), source = x)
    allocate(W%y(W%ny), source = y)
    allocate(W%z(size(z,1),size(z,2)), source = z)

    call W%InitInterp()

    end subroutine TInterpGrid2D_Init


    subroutine TInterpGrid2D_error(W,S)
    class(TInterpGrid2D):: W
    character(LEN=*), intent(in) :: S

    write(*,*) 'InterGrid2D Error: '//trim(S)
    stop

    end subroutine TInterpGrid2D_error


    subroutine TInterpGrid2D_InitFromFile(W, Filename, xcol, ycol, zcol)
    class(TInterpGrid2D):: W
    character(LEN=*), intent(in) :: Filename
    integer, intent(in), optional :: xcol, ycol, zcol
    integer :: ixcol=1,iycol=2,izcol=3
    integer :: nx,ny
    integer :: parse, status
    real(sp_acc), allocatable :: tmp(:)
    real(sp_acc) lasty
    logical :: first
    character(LEN=:), allocatable :: InLine
    Type(TTextFile) :: F

    if (present(xcol)) ixcol=xcol
    if (present(ycol)) iycol=ycol
    if (present(zcol)) izcol=zcol

    allocate(tmp(max(ixcol,iycol,izcol)))

    call F%Open(FileName)

    do parse=1,2
        first = .true.
        ny=0
        nx=0
        do while(F%ReadLineSkipEmptyAndComments(InLine))

        read(InLine,*, iostat=status) tmp
        if (status/=0) call W%Error('Error reading line: '//trim(InLine))

        if (first .or. tmp(iycol)/=lasty) then
            lasty=tmp(iycol)
            ny=ny+1
            nx=0
            first = .false.
        end if

        nx=nx+1
        if (parse==2) then
            if (ny==1) then
                W%x(nx) = tmp(ixcol)
            else
                if (tmp(ixcol)/=W%x(nx)) call W%Error('Non-grid x values')
            end if
            if (nx==1) then
                W%y(ny) = tmp(iycol)
            else
                if (tmp(iycol)/=W%y(ny))call W%Error('Non-grid y values')
            end if
            W%z(nx, ny) = tmp(izcol)
        endif
        end do

        if (parse==2) exit

        if (nx<2 .or. ny<2) call W%Error('not enough values to interpolate')
        W%nx = nx
        W%ny = ny
        allocate(W%x(W%nx))
        allocate(W%y(W%ny))
        allocate(W%z(nx,ny))
        status=0
        call F%Rewind()
    end do

    call F%Close()

    call W%InitInterp()

    end subroutine TInterpGrid2D_InitFromFile

    subroutine TInterpGrid2D_Free(W)
    Type(TInterpGrid2D):: W

    if (allocated(W%Wk)) then
        deallocate(W%x)
        deallocate(W%y)
        deallocate(W%Wk)
        deallocate(W%z)
    end if
    W%Initialized = .false.

    end subroutine TInterpGrid2D_Free

    function TInterpGrid2D_Value(W,x,y,error) result (res)
    !Z matrix not stored internally to save mem, so must pass again
    class(TInterpGrid2D) :: W
    real(GI) res, z(1), xx(1),yy(1)
    real(GI), intent(in) :: x,y
    integer, intent(inout), optional :: error

    xx(1)=x
    yy(1)=y
    call W%Values(1,xx,yy,z,error)
    res = z(1)

    end function TInterpGrid2D_Value

    subroutine TInterpGrid2D_Values(W, nip, x,y,z, error)
    !Z matrix not stored internally to save mem, so must pass again
    class(TInterpGrid2D) :: W
    integer, intent(in) :: nip
    real(GI), intent(out):: z(*)
    real(GI), intent(in) :: x(*),y(*)
    integer, intent(inout), optional :: error
    integer md,ier

    md=2
    if (.not. W%Initialized)  call W%FirstUse

    call rgbi3p(W%Wk,md, W%nx, W%ny, W%x, W%y, W%z, nip, x, y, z, ier)
    if (present(error)) then
        error=ier
    elseif (ier/=0) then
        call W%Error('error interpolating value')
    end if

    end subroutine TInterpGrid2D_Values


    SUBROUTINE rgbi3p(Wk,md, nxd, nyd, xd, yd, zd, nip, xi, yi, zi, ier)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2003-06-11  Time: 10:11:03

    ! Rectangular-grid bivariate interpolation
    ! (a master subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine performs interpolation of a bivariate function, z(x,y), on a
    ! rectangular grid in the x-y plane.  It is based on the revised Akima method.

    ! In this subroutine, the interpolating function is a piecewise function
    ! composed of a set of bicubic (bivariate third-degree) polynomials, each
    ! applicable to a rectangle of the input grid in the x-y plane.
    ! Each polynomial is determined locally.

    ! This subroutine has the accuracy of a bicubic polynomial, i.e., it
    ! interpolates accurately when all data points lie on a surface of a
    ! bicubic polynomial.

    ! The grid lines can be unevenly spaced.

    ! The input arguments are
    !   MD  = mode of computation
    !       = 1 for new XD, YD, or ZD data (default)
    !       = 2 for old XD, YD, and ZD data,
    !   NXD = number of the input-grid data points in the x coordinate
    !         (must be 2 or greater),
    !   NYD = number of the input-grid data points in the y coordinate
    !         (must be 2 or greater),
    !   XD  = array of dimension NXD containing the x coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   YD  = array of dimension NYD containing the y coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   ZD  = two-dimensional array of dimension NXD*NYD
    !         containing the z(x,y) values at the input-grid data points,
    !   NIP = number of the output points at which interpolation
    !         of the z value is desired (must be 1 or greater),
    !   XI  = array of dimension NIP containing the x coordinates
    !         of the output points,
    !   YI  = array of dimension NIP containing the y coordinates
    !         of the output points.

    ! The output arguments are
    !   ZI  = array of dimension NIP where the interpolated z
    !         values at the output points are to be stored,
    !   IER = error flag
    !       = 0 for no errors
    !       = 1 for NXD = 1 or less
    !       = 2 for NYD = 1 or less
    !       = 3 for identical XD values or XD values out of sequence
    !       = 4 for identical YD values or YD values out of sequence
    !       = 5 for NIP = 0 or less.

    ! N.B. The workspace has been removed from the argument list.
    !   WK  = three dimensional array of dimension 3*NXD*NYD used internally
    !         as a work area.

    ! The very fisrt call to this subroutine and the call with a new XD, YD, and
    ! ZD array must be made with MD=1.  The call with MD=2 must be preceded by
    ! another call with the same XD, YD, and ZD arrays.  Between the call with
    ! MD=2 and its preceding call, the WK array must not be disturbed.

    ! The constant in the PARAMETER statement below is
    !   NIPIMX = maximum number of output points to be processed at a time.
    ! The constant value has been selected empirically.

    ! This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.


    ! Specification statements
    !     .. Parameters ..

    INTEGER, INTENT(IN)   :: md
    INTEGER, INTENT(IN)   :: nxd
    INTEGER, INTENT(IN)   :: nyd
    REAL(GI), INTENT(IN)      :: xd(nxd)
    REAL(GI), INTENT(IN)      :: yd(nyd)
    REAL(GI), INTENT(IN)  :: zd(nxd,nyd)
    INTEGER, INTENT(IN)   :: nip
    REAL(GI), INTENT(IN)  :: xi(nip)
    REAL(GI), INTENT(IN)  :: yi(nip)
    REAL(GI), INTENT(OUT)  :: zi(nip)
    INTEGER, INTENT(OUT)  :: ier
    REAL(GI), INTENT(INOUT)  :: wk(3,nxd,nyd)

    !     ..
    !     .. Local Scalars ..
    INTEGER, PARAMETER  :: nipimx=51

    INTEGER  :: iip, ix, iy, nipi
    !     ..
    !     .. Local Arrays ..
    INTEGER  :: inxi(nipimx), inyi(nipimx)

    !     ..
    !     .. External Subroutines ..
    ! EXTERNAL         rglctn, rgpd3p, rgplnl
    !     ..
    !     .. Intrinsic Functions ..
    ! INTRINSIC        MIN
    !     ..

    ! Preliminary processing
    ! Error check
    IF (nxd <= 1) GO TO 40
    IF (nyd <= 1) GO TO 50
    DO  ix = 2,nxd
        IF (xd(ix) <= xd(ix-1)) GO TO 60
    END DO
    DO  iy = 2,nyd
        IF (yd(iy) <= yd(iy-1)) GO TO 70
    END DO
    IF (nip <= 0) GO TO 80
    ier = 0

    ! Calculation
    ! Estimates partial derivatives at all input-grid data points (for MD=1).
    IF (md /= 2) THEN
        CALL rgpd3p(nxd, nyd, xd, yd, zd, wk)
    END IF

    ! DO-loop with respect to the output point
    ! Processes NIPIMX output points, at most, at a time.
    DO  iip = 1,nip,nipimx
        nipi = MIN(nip- (iip-1),nipimx)
        ! Locates the output points.
        CALL rglctn(nxd, nyd, xd, yd, nipi, xi(iip), yi(iip), inxi, inyi)

        ! Calculates the z values at the output points.
        CALL rgplnl(nxd, nyd, xd, yd, zd, wk, nipi, xi(iip), yi(iip), inxi, inyi, &
        zi(iip))
    END DO
    RETURN

    ! Error exit
40  WRITE (*,FMT=9000)
    ier = 1
    GO TO 90
50  WRITE (*,FMT=9010)
    ier = 2
    GO TO 90
60  WRITE (*,FMT=9020) ix,xd(ix)
    ier = 3
    GO TO 90
70  WRITE (*,FMT=9030) iy,yd(iy)
    ier = 4
    GO TO 90
80  WRITE (*,FMT=9040)
    ier = 5
90  WRITE (*,FMT=9050) nxd,nyd,nip
    RETURN

    ! Format statements for error messages
9000 FORMAT (/' *** RGBI3P Error 1: NXD = 1 or less')
9010 FORMAT (/' *** RGBI3P Error 2: NYD = 1 or less')
9020 FORMAT (/' *** RGBI3P Error 3: Identical XD values or',  &
    ' XD values out of sequence'/ '    IX =', i6, ',  XD(IX) =', e11.3)
9030 FORMAT (/' *** RGBI3P Error 4: Identical YD values or',  &
    ' YD values out of sequence',/,'    IY =',i6,',  YD(IY) =', e11.3)
9040 FORMAT (/' *** RGBI3P Error 5: NIP = 0 or less')
9050 FORMAT ('    NXD =', i5,',  NYD =', i5,',  NIP =', i5/)
    END SUBROUTINE rgbi3p



    SUBROUTINE rgsf3p(wk, md, nxd, nyd, xd, yd, zd, nxi, xi, nyi, yi, zi, ier)

    ! Rectangular-grid surface fitting
    ! (a master subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine performs surface fitting by interpolating values of a
    ! bivariate function, z(x,y), on a rectangular grid in the x-y plane.
    ! It is based on the revised Akima method.

    ! In this subroutine, the interpolating function is a piecewise function
    ! composed of a set of bicubic (bivariate third-degree) polynomials, each
    ! applicable to a rectangle of the input grid in the x-y plane.
    ! Each polynomial is determined locally.

    ! This subroutine has the accuracy of a bicubic polynomial, i.e., it fits the
    ! surface accurately when all data points lie on a surface of a bicubic
    ! polynomial.

    ! The grid lines of the input and output data can be unevenly spaced.

    ! The input arguments are
    !   MD  = mode of computation
    !       = 1 for new XD, YD, or ZD data (default)
    !       = 2 for old XD, YD, and ZD data,
    !   NXD = number of the input-grid data points in the x
    !         coordinate (must be 2 or greater),
    !   NYD = number of the input-grid data points in the y
    !         coordinate (must be 2 or greater),
    !   XD  = array of dimension NXD containing the x coordinates
    !         of the input-grid data points (must be in a
    !         monotonic increasing order),
    !   YD  = array of dimension NYD containing the y coordinates
    !         of the input-grid data points (must be in a
    !         monotonic increasing order),
    !   ZD  = two-dimensional array of dimension NXD*NYD
    !         containing the z(x,y) values at the input-grid data points,
    !   NXI = number of output grid points in the x coordinate
    !         (must be 1 or greater),
    !   XI  = array of dimension NXI containing the x coordinates
    !         of the output grid points,
    !   NYI = number of output grid points in the y coordinate
    !         (must be 1 or greater),
    !   YI  = array of dimension NYI containing the y coordinates
    !         of the output grid points.

    ! The output arguments are
    !   ZI  = two-dimensional array of dimension NXI*NYI, where the interpolated
    !         z values at the output grid points are to be stored,
    !   IER = error flag
    !       = 0 for no error
    !       = 1 for NXD = 1 or less
    !       = 2 for NYD = 1 or less
    !       = 3 for identical XD values or XD values out of sequence
    !       = 4 for identical YD values or YD values out of sequence
    !       = 5 for NXI = 0 or less
    !       = 6 for NYI = 0 or less.

    ! N.B. The workspace has been removed from the argument list.
    !   WK  = three-dimensional array of dimension 3*NXD*NYD used internally
    !         as a work area.

    ! The very first call to this subroutine and the call with a new XD, YD, or
    ! ZD array must be made with MD=1.  The call with MD=2 must be preceded by
    ! another call with the same XD, YD, and ZD arrays.  Between the call with
    ! MD=2 and its preceding call, the WK array must not be disturbed.

    ! The constant in the PARAMETER statement below is
    !   NIPIMX = maximum number of output points to be processed at a time.
    ! The constant value has been selected empirically.

    ! This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.


    ! Specification statements
    !     .. Parameters ..

    INTEGER, INTENT(IN)   :: md
    INTEGER, INTENT(IN)   :: nxd
    INTEGER, INTENT(IN)   :: nyd
    REAL(GI), INTENT(IN)      :: xd(nxd)
    REAL(GI), INTENT(IN)      :: yd(nyd)
    REAL(GI), INTENT(IN OUT)  :: zd(nxd,nyd)
    INTEGER, INTENT(IN)   :: nxi
    REAL(GI), INTENT(IN OUT)  :: xi(nxi)
    INTEGER, INTENT(IN)   :: nyi
    REAL(GI), INTENT(IN)      :: yi(nyi)
    REAL(GI), INTENT(IN OUT)  :: zi(nxi,nyi)
    INTEGER, INTENT(OUT)  :: ier
    REAL(GI), INTENT(INOUT)  :: wk(3,nxd,nyd)

    !     ..
    !     .. Local Scalars ..
    INTEGER, PARAMETER  :: nipimx=51

    INTEGER  :: ix, ixi, iy, iyi, nipi
    !     ..
    !     .. Local Arrays ..
    REAL(GI)     :: yii(nipimx)
    INTEGER  :: inxi(nipimx), inyi(nipimx)

    !     ..
    !     .. External Subroutines ..
    ! EXTERNAL         rglctn,rgpd3p,rgplnl
    !     ..
    !     .. Intrinsic Functions ..
    ! INTRINSIC        MIN
    !     ..

    ! Preliminary processing
    ! Error check
    IF (nxd <= 1) GO TO 60
    IF (nyd <= 1) GO TO 70
    DO  ix = 2,nxd
        IF (xd(ix) <= xd(ix-1)) GO TO 80
    END DO
    DO  iy = 2,nyd
        IF (yd(iy) <= yd(iy-1)) GO TO 90
    END DO
    IF (nxi <= 0) GO TO 100
    IF (nyi <= 0) GO TO 110
    ier = 0

    ! Calculation
    ! Estimates partial derivatives at all input-grid data points
    ! (for MD=1).
    IF (md /= 2) THEN
        CALL rgpd3p(nxd, nyd, xd, yd, zd, wk)
    END IF

    ! Outermost DO-loop with respect to the y coordinate of the output grid points
    DO  iyi = 1,nyi
        DO  ixi = 1,nipimx
            yii(ixi) = yi(iyi)
        END DO

        ! Second DO-loop with respect to the x coordinate of the output grid points
        ! Processes NIPIMX output-grid points, at most, at a time.
        DO  ixi = 1,nxi,nipimx
            nipi = MIN(nxi- (ixi-1), nipimx)
            ! Locates the output-grid points.
            CALL rglctn(nxd, nyd, xd, yd, nipi, xi(ixi), yii, inxi, inyi)

            ! Calculates the z values at the output-grid points.
            CALL rgplnl(nxd, nyd, xd, yd, zd, wk, nipi, xi(ixi), yii, inxi, inyi,  &
            zi(ixi,iyi))
        END DO
    END DO
    RETURN

    ! Error exit
60  WRITE (*,FMT=9000)
    ier = 1
    GO TO 120
70  WRITE (*,FMT=9010)
    ier = 2
    GO TO 120
80  WRITE (*,FMT=9020) ix,xd(ix)
    ier = 3
    GO TO 120
90  WRITE (*,FMT=9030) iy,yd(iy)
    ier = 4
    GO TO 120
100 WRITE (*,FMT=9040)
    ier = 5
    GO TO 120
110 WRITE (*,FMT=9050)
    ier = 6
120 WRITE (*,FMT=9060) nxd,nyd,nxi,nyi
    RETURN

    ! Format statements for error messages
9000 FORMAT (/' *** RGSF3P Error 1: NXD = 1 or less')
9010 FORMAT (/' *** RGSF3P Error 2: NYD = 1 or less')
9020 FORMAT (/' *** RGSF3P Error 3: Identical XD values or',  &
    ' XD values out of sequence',/,'    IX =',i6,',  XD(IX) =', e11.3)
9030 FORMAT (/' *** RGSF3P Error 4: Identical YD values or',  &
    ' YD values out of sequence',/,'    IY =',i6,',  YD(IY) =', e11.3)
9040 FORMAT (/' *** RGSF3P Error 5: NXI = 0 or less')
9050 FORMAT (/' *** RGSF3P Error 6: NYI = 0 or less')
9060 FORMAT ('    NXD =', i5, ',  NYD =', i5, ',  NXI =', i5,',  NYI =', i5 /)
    END SUBROUTINE rgsf3p



    !     ..
    ! Statement Function definitions
    ! z2f(xx1,xx2,zz0,zz1) = (zz1-zz0)*xx2/xx1 + zz0
    ! z3f(xx1,xx2,xx3,zz0,zz1,zz2) = ((zz2-zz0)* (xx3-xx1)/xx2 -  &
    !    (zz1-zz0)* (xx3-xx2)/xx1)* (xx3/ (xx2-xx1)) + zz0

    FUNCTION z2f(xx1, xx2, zz0, zz1) RESULT(fn_val)

    REAL(GI), INTENT(IN)  :: xx1, xx2, zz0, zz1
    REAL(GI)              :: fn_val

    fn_val = (zz1-zz0)*xx2/xx1 + zz0
    RETURN
    END FUNCTION z2f



    FUNCTION z3f(xx1, xx2, xx3, zz0, zz1, zz2) RESULT(fn_val)

    REAL(GI), INTENT(IN)  :: xx1, xx2, xx3, zz0, zz1, zz2
    REAL(GI)              :: fn_val

    fn_val = ((zz2-zz0)*(xx3-xx1)/xx2 - (zz1-zz0)*(xx3-xx2)/xx1) *  &
    (xx3/(xx2-xx1)) + zz0
    RETURN
    END FUNCTION z3f



    SUBROUTINE rgpd3p(nxd, nyd, xd, yd, zd, pdd)

    ! Partial derivatives of a bivariate function on a rectangular grid
    ! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine estimates three partial derivatives, zx, zy, and
    ! zxy, of a bivariate function, z(x,y), on a rectangular grid in
    ! the x-y plane.  It is based on the revised Akima method that has
    ! the accuracy of a bicubic polynomial.

    ! The input arguments are
    !   NXD = number of the input-grid data points in the x
    !         coordinate (must be 2 or greater),
    !   NYD = number of the input-grid data points in the y
    !         coordinate (must be 2 or greater),
    !   XD  = array of dimension NXD containing the x coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   YD  = array of dimension NYD containing the y coordinates of the
    !         input-grid data points (must be in a monotonic increasing order),
    !   ZD  = two-dimensional array of dimension NXD*NYD
    !         containing the z(x,y) values at the input-grid data points.

    ! The output argument is
    !   PDD = three-dimensional array of dimension 3*NXD*NYD,
    !         where the estimated zx, zy, and zxy values at the
    !         input-grid data points are to be stored.


    ! Specification statements
    !     .. Scalar Arguments ..

    INTEGER, INTENT(IN)  :: nxd
    INTEGER, INTENT(IN)  :: nyd
    REAL(GI), INTENT(IN)     :: xd(nxd)
    REAL(GI), INTENT(IN)     :: yd(nyd)
    REAL(GI), INTENT(IN)     :: zd(nxd,nyd)
    REAL(GI), INTENT(OUT)    :: pdd(3,nxd,nyd)

    !     ..
    !     .. Local Scalars ..
    REAL(GI) :: b00, b00x, b00y, b01, b10, b11, cx1, cx2, cx3, cy1, cy2,  &
    cy3, disf, dnm, dz00, dz01, dz02, dz03, dz10, dz11, dz12,  &
    dz13, dz20, dz21, dz22, dz23, dz30, dz31, dz32, dz33,  &
    dzx10, dzx20, dzx30, dzxy11, dzxy12, dzxy13, dzxy21,  &
    dzxy22, dzxy23, dzxy31, dzxy32, dzxy33, dzy01, dzy02,  &
    dzy03, epsln, pezx, pezxy, pezy, smpef, smpei, smwtf,  &
    smwti, sx, sxx, sxxy, sxxyy, sxy, sxyy, sxyz, sxz, sy, syy,  &
    syz, sz, volf, wt, x0, x1, x2, x3, y0, y1, y2,  &
    y3, z00, z01, z02, z03, z10, z11, z12, z13, z20, z21, z22,  &
    z23, z30, z31, z32, z33, zxdi, zxydi, zydi
    INTEGER :: ipex, ipey, ix0, ix1, ix2, ix3, iy0, iy1, iy2, iy3, nx0, ny0
    !     ..
    !     .. Local Arrays ..
    REAL(GI)    :: b00xa(4), b00ya(4), b01a(4), b10a(4), cxa(3,4), cya(3,4),   &
    sxa(4), sxxa(4), sya(4), syya(4), xa(3,4), ya(3,4),   &
    z0ia(3,4), zi0a(3,4)

    INTEGER, parameter :: idlt(3,4) = RESHAPE([-3, -2, -1, -2, -1, 1,-1, 1,2, 1 ,2, 3 ], [ 3, 4 ] )
    !AL Jun 14: Fixed F90 translation

    ! Calculation
    ! Initial setting of some local variables
    nx0 = MAX(4,nxd)
    ny0 = MAX(4,nyd)

    ! Double DO-loop with respect to the input grid points
    DO  iy0 = 1,nyd
        DO  ix0 = 1,nxd
            x0 = xd(ix0)
            y0 = yd(iy0)
            z00 = zd(ix0,iy0)

            ! Part 1.  Estimation of ZXDI
            ! Initial setting
            smpef = 0.0
            smwtf = 0.0
            smpei = 0.0
            smwti = 0.0
            ! DO-loop with respect to the primary estimate
            DO  ipex = 1,4
                ! Selects necessary grid points in the x direction.
                ix1 = ix0 + idlt(1,ipex)
                ix2 = ix0 + idlt(2,ipex)
                ix3 = ix0 + idlt(3,ipex)
                IF ((ix1 < 1) .OR. (ix2 < 1) .OR. (ix3 < 1) .OR.  &
                (ix1 > nx0) .OR. (ix2 > nx0) .OR. (ix3 > nx0)) CYCLE
                ! Selects and/or supplements the x and z values.
                x1 = xd(ix1) - x0
                z10 = zd(ix1,iy0)
                IF (nxd >= 4) THEN
                    x2 = xd(ix2) - x0
                    x3 = xd(ix3) - x0
                    z20 = zd(ix2,iy0)
                    z30 = zd(ix3,iy0)
                ELSE IF (nxd == 3) THEN
                    x2 = xd(ix2) - x0
                    z20 = zd(ix2,iy0)
                    x3 = 2*xd(3) - xd(2) - x0
                    z30 = z3f(x1,x2,x3,z00,z10,z20)
                ELSE IF (nxd == 2) THEN
                    x2 = 2*xd(2) - xd(1) - x0
                    z20 = z2f(x1,x2,z00,z10)
                    x3 = 2*xd(1) - xd(2) - x0
                    z30 = z2f(x1,x3,z00,z10)
                END IF
                dzx10 = (z10-z00)/x1
                dzx20 = (z20-z00)/x2
                dzx30 = (z30-z00)/x3
                ! Calculates the primary estimate of partial derivative zx as
                ! the coefficient of the bicubic polynomial.
                cx1 = x2*x3/ ((x1-x2)* (x1-x3))
                cx2 = x3*x1/ ((x2-x3)* (x2-x1))
                cx3 = x1*x2/ ((x3-x1)* (x3-x2))
                pezx = cx1*dzx10 + cx2*dzx20 + cx3*dzx30
                ! Calculates the volatility factor and distance factor in the x
                ! direction for the primary estimate of zx.
                sx = x1 + x2 + x3
                sz = z00 + z10 + z20 + z30
                sxx = x1*x1 + x2*x2 + x3*x3
                sxz = x1*z10 + x2*z20 + x3*z30
                dnm = 4.0*sxx - sx*sx
                b00 = (sxx*sz-sx*sxz)/dnm
                b10 = (4.0*sxz-sx*sz)/dnm
                dz00 = z00 - b00
                dz10 = z10 - (b00+b10*x1)
                dz20 = z20 - (b00+b10*x2)
                dz30 = z30 - (b00+b10*x3)
                volf = dz00**2 + dz10**2 + dz20**2 + dz30**2
                disf = sxx

                ! Calculates the EPSLN value, which is used to decide whether or
                ! not the volatility factor is essentially zero.
                epsln = (z00**2+z10**2+z20**2+z30**2)*1.0E-12
                ! Accumulates the weighted primary estimates of zx and their weights.
                IF (volf > epsln) THEN
                    ! - For a finite weight.
                    wt = 1.0/ (volf*disf)
                    smpef = smpef + wt*pezx
                    smwtf = smwtf + wt
                ELSE
                    ! - For an infinite weight.
                    smpei = smpei + pezx
                    smwti = smwti + 1.0
                END IF

                ! Saves the necessary values for estimating zxy
                xa(1,ipex) = x1
                xa(2,ipex) = x2
                xa(3,ipex) = x3
                zi0a(1,ipex) = z10
                zi0a(2,ipex) = z20
                zi0a(3,ipex) = z30
                cxa(1,ipex) = cx1
                cxa(2,ipex) = cx2
                cxa(3,ipex) = cx3
                sxa(ipex) = sx
                sxxa(ipex) = sxx
                b00xa(ipex) = b00
                b10a(ipex) = b10
            END DO

            ! Calculates the final estimate of zx.
            IF (smwti < 0.5) THEN
                ! - When no infinite weights exist.
                zxdi = smpef/smwtf
            ELSE
                ! - When infinite weights exist.
                zxdi = smpei/smwti
            END IF
            ! End of Part 1.

            ! Part 2.  Estimation of ZYDI
            ! Initial setting
            smpef = 0.0
            smwtf = 0.0
            smpei = 0.0
            smwti = 0.0
            ! DO-loop with respect to the primary estimate
            DO  ipey = 1,4
                ! Selects necessary grid points in the y direction.
                iy1 = iy0 + idlt(1,ipey)
                iy2 = iy0 + idlt(2,ipey)
                iy3 = iy0 + idlt(3,ipey)
                IF ((iy1 < 1) .OR. (iy2 < 1) .OR. (iy3 < 1) .OR.  &
                (iy1 > ny0) .OR. (iy2 > ny0) .OR. (iy3 > ny0)) CYCLE
                ! Selects and/or supplements the y and z values.
                y1 = yd(iy1) - y0
                z01 = zd(ix0,iy1)
                IF (nyd >= 4) THEN
                    y2 = yd(iy2) - y0
                    y3 = yd(iy3) - y0
                    z02 = zd(ix0,iy2)
                    z03 = zd(ix0,iy3)
                ELSE IF (nyd == 3) THEN
                    y2 = yd(iy2) - y0
                    z02 = zd(ix0,iy2)
                    y3 = 2*yd(3) - yd(2) - y0
                    z03 = z3f(y1,y2,y3,z00,z01,z02)
                ELSE IF (nyd == 2) THEN
                    y2 = 2*yd(2) - yd(1) - y0
                    z02 = z2f(y1,y2,z00,z01)
                    y3 = 2*yd(1) - yd(2) - y0
                    z03 = z2f(y1,y3,z00,z01)
                END IF
                dzy01 = (z01-z00)/y1
                dzy02 = (z02-z00)/y2
                dzy03 = (z03-z00)/y3
                ! Calculates the primary estimate of partial derivative zy as
                ! the coefficient of the bicubic polynomial.
                cy1 = y2*y3/ ((y1-y2)* (y1-y3))
                cy2 = y3*y1/ ((y2-y3)* (y2-y1))
                cy3 = y1*y2/ ((y3-y1)* (y3-y2))
                pezy = cy1*dzy01 + cy2*dzy02 + cy3*dzy03
                ! Calculates the volatility factor and distance factor in the y
                ! direction for the primary estimate of zy.
                sy = y1 + y2 + y3
                sz = z00 + z01 + z02 + z03
                syy = y1*y1 + y2*y2 + y3*y3
                syz = y1*z01 + y2*z02 + y3*z03
                dnm = 4.0*syy - sy*sy
                b00 = (syy*sz-sy*syz)/dnm
                b01 = (4.0*syz-sy*sz)/dnm
                dz00 = z00 - b00
                dz01 = z01 - (b00+b01*y1)
                dz02 = z02 - (b00+b01*y2)
                dz03 = z03 - (b00+b01*y3)
                volf = dz00**2 + dz01**2 + dz02**2 + dz03**2
                disf = syy

                ! Calculates the EPSLN value, which is used to decide whether or
                ! not the volatility factor is essentially zero.
                epsln = (z00**2+z01**2+z02**2+z03**2)*1.0E-12
                ! Accumulates the weighted primary estimates of zy and their weights.
                IF (volf > epsln) THEN
                    ! - For a finite weight.
                    wt = 1.0/ (volf*disf)
                    smpef = smpef + wt*pezy
                    smwtf = smwtf + wt
                ELSE
                    ! - For an infinite weight.
                    smpei = smpei + pezy
                    smwti = smwti + 1.0
                END IF
                ! Saves the necessary values for estimating zxy
                ya(1,ipey) = y1
                ya(2,ipey) = y2
                ya(3,ipey) = y3
                z0ia(1,ipey) = z01
                z0ia(2,ipey) = z02
                z0ia(3,ipey) = z03
                cya(1,ipey) = cy1
                cya(2,ipey) = cy2
                cya(3,ipey) = cy3
                sya(ipey) = sy
                syya(ipey) = syy
                b00ya(ipey) = b00
                b01a(ipey) = b01
            END DO

            ! Calculates the final estimate of zy.
            IF (smwti < 0.5) THEN
                ! - When no infinite weights exist.
                zydi = smpef/smwtf
            ELSE
                ! - When infinite weights exist.
                zydi = smpei/smwti
            END IF
            ! End of Part 2.

            ! Part 3.  Estimation of ZXYDI
            ! Initial setting
            smpef = 0.0
            smwtf = 0.0
            smpei = 0.0
            smwti = 0.0
            ! Outer DO-loops with respect to the primary estimates in the x direction
            DO  ipex = 1,4
                ix1 = ix0 + idlt(1,ipex)
                ix2 = ix0 + idlt(2,ipex)
                ix3 = ix0 + idlt(3,ipex)
                IF ((ix1 < 1) .OR. (ix2 < 1) .OR. (ix3 < 1) .OR.  &
                (ix1 > nx0) .OR. (ix2 > nx0) .OR. (ix3 > nx0)) CYCLE
                ! Retrieves the necessary values for estimating zxy in the x direction.
                x1 = xa(1,ipex)
                x2 = xa(2,ipex)
                x3 = xa(3,ipex)
                z10 = zi0a(1,ipex)
                z20 = zi0a(2,ipex)
                z30 = zi0a(3,ipex)
                cx1 = cxa(1,ipex)
                cx2 = cxa(2,ipex)
                cx3 = cxa(3,ipex)
                sx = sxa(ipex)
                sxx = sxxa(ipex)
                b00x = b00xa(ipex)
                b10 = b10a(ipex)

                ! Inner DO-loops with respect to the primary estimates in the y direction
                DO  ipey = 1,4
                    iy1 = iy0 + idlt(1,ipey)
                    iy2 = iy0 + idlt(2,ipey)
                    iy3 = iy0 + idlt(3,ipey)
                    IF ((iy1 < 1) .OR. (iy2 < 1) .OR. (iy3 < 1) .OR. (iy1 > ny0) .OR.  &
                    (iy2 > ny0) .OR. (iy3 > ny0)) CYCLE
                    ! Retrieves the necessary values for estimating zxy in the y direction.
                    y1 = ya(1,ipey)
                    y2 = ya(2,ipey)
                    y3 = ya(3,ipey)
                    z01 = z0ia(1,ipey)
                    z02 = z0ia(2,ipey)
                    z03 = z0ia(3,ipey)
                    cy1 = cya(1,ipey)
                    cy2 = cya(2,ipey)
                    cy3 = cya(3,ipey)
                    sy = sya(ipey)
                    syy = syya(ipey)
                    b00y = b00ya(ipey)
                    b01 = b01a(ipey)
                    ! Selects and/or supplements the z values.
                    IF (nyd >= 4) THEN
                        z11 = zd(ix1,iy1)
                        z12 = zd(ix1,iy2)
                        z13 = zd(ix1,iy3)
                        IF (nxd >= 4) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z23 = zd(ix2,iy3)
                            z31 = zd(ix3,iy1)
                            z32 = zd(ix3,iy2)
                            z33 = zd(ix3,iy3)
                        ELSE IF (nxd == 3) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z23 = zd(ix2,iy3)
                            z31 = z3f(x1,x2,x3,z01,z11,z21)
                            z32 = z3f(x1,x2,x3,z02,z12,z22)
                            z33 = z3f(x1,x2,x3,z03,z13,z23)
                        ELSE IF (nxd == 2) THEN
                            z21 = z2f(x1,x2,z01,z11)
                            z22 = z2f(x1,x2,z02,z12)
                            z23 = z2f(x1,x2,z03,z13)
                            z31 = z2f(x1,x3,z01,z11)
                            z32 = z2f(x1,x3,z02,z12)
                            z33 = z2f(x1,x3,z03,z13)
                        END IF
                    ELSE IF (nyd == 3) THEN
                        z11 = zd(ix1,iy1)
                        z12 = zd(ix1,iy2)
                        z13 = z3f(y1,y2,y3,z10,z11,z12)
                        IF (nxd >= 4) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z31 = zd(ix3,iy1)
                            z32 = zd(ix3,iy2)
                        ELSE IF (nxd == 3) THEN
                            z21 = zd(ix2,iy1)
                            z22 = zd(ix2,iy2)
                            z31 = z3f(x1,x2,x3,z01,z11,z21)
                            z32 = z3f(x1,x2,x3,z02,z12,z22)
                        ELSE IF (nxd == 2) THEN
                            z21 = z2f(x1,x2,z01,z11)
                            z22 = z2f(x1,x2,z02,z12)
                            z31 = z2f(x1,x3,z01,z11)
                            z32 = z2f(x1,x3,z02,z12)
                        END IF
                        z23 = z3f(y1,y2,y3,z20,z21,z22)
                        z33 = z3f(y1,y2,y3,z30,z31,z32)
                    ELSE IF (nyd == 2) THEN
                        z11 = zd(ix1,iy1)
                        z12 = z2f(y1,y2,z10,z11)
                        z13 = z2f(y1,y3,z10,z11)
                        IF (nxd >= 4) THEN
                            z21 = zd(ix2,iy1)
                            z31 = zd(ix3,iy1)
                        ELSE IF (nxd == 3) THEN
                            z21 = zd(ix2,iy1)
                            z31 = z3f(x1,x2,x3,z01,z11,z21)
                        ELSE IF (nxd == 2) THEN
                            z21 = z2f(x1,x2,z01,z11)
                            z31 = z2f(x1,x3,z01,z11)
                        END IF
                        z22 = z2f(y1,y2,z20,z21)
                        z23 = z2f(y1,y3,z20,z21)
                        z32 = z2f(y1,y2,z30,z31)
                        z33 = z2f(y1,y3,z30,z31)
                    END IF
                    ! Calculates the primary estimate of partial derivative zxy as
                    ! the coefficient of the bicubic polynomial.
                    dzxy11 = (z11-z10-z01+z00)/ (x1*y1)
                    dzxy12 = (z12-z10-z02+z00)/ (x1*y2)
                    dzxy13 = (z13-z10-z03+z00)/ (x1*y3)
                    dzxy21 = (z21-z20-z01+z00)/ (x2*y1)
                    dzxy22 = (z22-z20-z02+z00)/ (x2*y2)
                    dzxy23 = (z23-z20-z03+z00)/ (x2*y3)
                    dzxy31 = (z31-z30-z01+z00)/ (x3*y1)
                    dzxy32 = (z32-z30-z02+z00)/ (x3*y2)
                    dzxy33 = (z33-z30-z03+z00)/ (x3*y3)
                    pezxy = cx1* (cy1*dzxy11+cy2*dzxy12+cy3*dzxy13) +  &
                    cx2* (cy1*dzxy21+cy2*dzxy22+cy3*dzxy23) +  &
                    cx3* (cy1*dzxy31+cy2*dzxy32+cy3*dzxy33)
                    ! Calculates the volatility factor and distance factor in the x
                    ! and y directions for the primary estimate of zxy.
                    b00 = (b00x+b00y)/2.0
                    sxy = sx*sy
                    sxxy = sxx*sy
                    sxyy = sx*syy
                    sxxyy = sxx*syy
                    sxyz = x1* (y1*z11+y2*z12+y3*z13) + x2* (y1*z21+y2*z22+y3*z23) +  &
                    x3* (y1*z31+y2*z32+y3*z33)
                    b11 = (sxyz-b00*sxy-b10*sxxy-b01*sxyy)/sxxyy
                    dz00 = z00 - b00
                    dz01 = z01 - (b00+b01*y1)
                    dz02 = z02 - (b00+b01*y2)
                    dz03 = z03 - (b00+b01*y3)
                    dz10 = z10 - (b00+b10*x1)
                    dz11 = z11 - (b00+b01*y1+x1* (b10+b11*y1))
                    dz12 = z12 - (b00+b01*y2+x1* (b10+b11*y2))
                    dz13 = z13 - (b00+b01*y3+x1* (b10+b11*y3))
                    dz20 = z20 - (b00+b10*x2)
                    dz21 = z21 - (b00+b01*y1+x2* (b10+b11*y1))
                    dz22 = z22 - (b00+b01*y2+x2* (b10+b11*y2))
                    dz23 = z23 - (b00+b01*y3+x2* (b10+b11*y3))
                    dz30 = z30 - (b00+b10*x3)
                    dz31 = z31 - (b00+b01*y1+x3* (b10+b11*y1))
                    dz32 = z32 - (b00+b01*y2+x3* (b10+b11*y2))
                    dz33 = z33 - (b00+b01*y3+x3* (b10+b11*y3))
                    volf = dz00**2 + dz01**2 + dz02**2 + dz03**2 +  &
                    dz10**2 + dz11**2 + dz12**2 + dz13**2 +  &
                    dz20**2 + dz21**2 + dz22**2 + dz23**2 +  &
                    dz30**2 + dz31**2 + dz32**2 + dz33**2
                    disf = sxx*syy
                    ! Calculates EPSLN.
                    epsln = (z00**2 + z01**2 + z02**2 + z03**2 + z10**2 +   &
                    z11**2 + z12**2 + z13**2 + z20**2 + z21**2 + z22**2 +   &
                    z23**2 + z30**2 + z31**2 + z32**2 + z33**2)* 1.0E-12
                    ! Accumulates the weighted primary estimates of zxy and their weights.
                    IF (volf > epsln) THEN
                        ! - For a finite weight.
                        wt = 1.0/ (volf*disf)
                        smpef = smpef + wt*pezxy
                        smwtf = smwtf + wt
                    ELSE
                        ! - For an infinite weight.
                        smpei = smpei + pezxy
                        smwti = smwti + 1.0
                    END IF
                END DO
            END DO

            ! Calculates the final estimate of zxy.
            IF (smwti < 0.5) THEN
                ! - When no infinite weights exist.
                zxydi = smpef/smwtf
            ELSE
                ! - When infinite weights exist.
                zxydi = smpei/smwti
            END IF
            ! End of Part 3

            pdd(1,ix0,iy0) = zxdi
            pdd(2,ix0,iy0) = zydi
            pdd(3,ix0,iy0) = zxydi
        END DO
    END DO
    RETURN
    END SUBROUTINE rgpd3p



    SUBROUTINE rglctn(nxd, nyd, xd, yd, nip, xi, yi, inxi, inyi)

    ! Location of the desired points in a rectangular grid
    ! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine locates the desired points in a rectangular grid
    ! in the x-y plane.

    ! The grid lines can be unevenly spaced.

    ! The input arguments are
    !   NXD  = number of the input-grid data points in the x
    !          coordinate (must be 2 or greater),
    !   NYD  = number of the input-grid data points in the y
    !          coordinate (must be 2 or greater),
    !   XD   = array of dimension NXD containing the x coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   YD   = array of dimension NYD containing the y coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   NIP  = number of the output points to be located (must be 1 or greater),
    !   XI   = array of dimension NIP containing the x coordinates
    !          of the output points to be located,
    !   YI   = array of dimension NIP containing the y coordinates
    !          of the output points to be located.

    ! The output arguments are
    !   INXI = integer array of dimension NIP where the interval
    !          numbers of the XI array elements are to be stored,
    !   INYI = integer array of dimension NIP where the interval
    !          numbers of the YI array elements are to be stored.
    ! The interval numbers are between 0 and NXD and between 0 and NYD,
    ! respectively.


    ! Specification statements
    !     .. Scalar Arguments ..

    INTEGER, INTENT(IN)   :: nxd
    INTEGER, INTENT(IN)   :: nyd
    REAL(GI), INTENT(IN)      :: xd(nxd)
    REAL(GI), INTENT(IN)      :: yd(nyd)
    INTEGER, INTENT(IN)   :: nip
    REAL(GI), INTENT(IN)      :: xi(nip)
    REAL(GI), INTENT(IN)      :: yi(nip)
    INTEGER, INTENT(OUT)  :: inxi(nip)
    INTEGER, INTENT(OUT)  :: inyi(nip)

    !     ..
    !     .. Local Scalars ..
    REAL(GI)     :: xii, yii
    INTEGER  :: iip, imd, imn, imx, ixd, iyd, nintx, ninty

    !     ..
    ! DO-loop with respect to IIP, which is the point number of the output point
    DO  iip = 1,nip
        xii = xi(iip)
        yii = yi(iip)
        ! Checks if the x coordinate of the IIPth output point, XII, is
        ! in a new interval.  (NINTX is the new-interval flag.)
        IF (iip == 1) THEN
            nintx = 1
        ELSE
            nintx = 0
            IF (ixd == 0) THEN
                IF (xii > xd(1)) nintx = 1
            ELSE IF (ixd < nxd) THEN
                IF ((xii < xd(ixd)) .OR. (xii > xd(ixd+1))) nintx = 1
            ELSE
                IF (xii < xd(nxd)) nintx = 1
            END IF
        END IF

        ! Locates the output point by binary search if XII is in a new interval.
        ! Determines IXD for which XII lies between XD(IXD) and XD(IXD+1).
        IF (nintx == 1) THEN
            IF (xii <= xd(1)) THEN
                ixd = 0
            ELSE IF (xii < xd(nxd)) THEN
                imn = 1
                imx = nxd
                imd = (imn+imx)/2
10              IF (xii >= xd(imd)) THEN
                    imn = imd
                ELSE
                    imx = imd
                END IF
                imd = (imn+imx)/2
                IF (imd > imn) GO TO 10
                ixd = imd
            ELSE
                ixd = nxd
            END IF
        END IF
        inxi(iip) = ixd

        ! Checks if the y coordinate of the IIPth output point, YII, is
        ! in a new interval.  (NINTY is the new-interval flag.)
        IF (iip == 1) THEN
            ninty = 1
        ELSE
            ninty = 0
            IF (iyd == 0) THEN
                IF (yii > yd(1)) ninty = 1
            ELSE IF (iyd < nyd) THEN
                IF ((yii < yd(iyd)) .OR. (yii > yd(iyd+1))) ninty = 1
            ELSE
                IF (yii < yd(nyd)) ninty = 1
            END IF
        END IF

        ! Locates the output point by binary search if YII is in a new interval.
        ! Determines IYD for which YII lies between YD(IYD) and YD(IYD+1).
        IF (ninty == 1) THEN
            IF (yii <= yd(1)) THEN
                iyd = 0
            ELSE IF (yii < yd(nyd)) THEN
                imn = 1
                imx = nyd
                imd = (imn+imx)/2
20              IF (yii >= yd(imd)) THEN
                    imn = imd
                ELSE
                    imx = imd
                END IF
                imd = (imn+imx)/2
                IF (imd > imn) GO TO 20
                iyd = imd
            ELSE
                iyd = nyd
            END IF
        END IF
        inyi(iip) = iyd
    END DO
    RETURN
    END SUBROUTINE rglctn



    SUBROUTINE rgplnl(nxd, nyd, xd, yd, zd, pdd, nip, xi, yi, inxi, inyi, zi)

    ! Polynomials for rectangular-grid bivariate interpolation and surface fitting
    ! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

    ! Hiroshi Akima
    ! U.S. Department of Commerce, NTIA/ITS
    ! Version of 1995/08

    ! This subroutine determines a polynomial in x and y for a rectangle of the
    ! input grid in the x-y plane and calculates the z value for the desired
    ! points by evaluating the polynomial for rectangular-grid bivariate
    ! interpolation and surface fitting.

    ! The input arguments are
    !   NXD  = number of the input-grid data points in the x
    !          coordinate (must be 2 or greater),
    !   NYD  = number of the input-grid data points in the y
    !          coordinate (must be 2 or greater),
    !   XD   = array of dimension NXD containing the x coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   YD   = array of dimension NYD containing the y coordinates of the
    !          input-grid data points (must be in a monotonic increasing order),
    !   ZD   = two-dimensional array of dimension NXD*NYD
    !          containing the z(x,y) values at the input-grid data points,
    !   PDD  = three-dimensional array of dimension 3*NXD*NYD
    !          containing the estimated zx, zy, and zxy values
    !          at the input-grid data points,
    !   NIP  = number of the output points at which interpolation
    !          is to be performed,
    !   XI   = array of dimension NIP containing the x coordinates
    !          of the output points,
    !   YI   = array of dimension NIP containing the y coordinates
    !          of the output points,
    !   INXI = integer array of dimension NIP containing the
    !          interval numbers of the input grid intervals in the
    !          x direction where the x coordinates of the output points lie,
    !   INYI = integer array of dimension NIP containing the
    !          interval numbers of the input grid intervals in the
    !          y direction where the y coordinates of the output points lie.

    ! The output argument is
    !   ZI   = array of dimension NIP, where the interpolated z
    !          values at the output points are to be stored.


    ! Specification statements
    !     .. Scalar Arguments ..

    INTEGER, INTENT(IN)  :: nxd
    INTEGER, INTENT(IN)  :: nyd
    REAL(GI), INTENT(IN)     :: xd(nxd)
    REAL(GI), INTENT(IN)     :: yd(nyd)
    REAL(GI), INTENT(IN)     :: zd(nxd,nyd)
    REAL(GI), INTENT(IN)     :: pdd(3,nxd,nyd)
    INTEGER, INTENT(IN)  :: nip
    REAL(GI), INTENT(IN)     :: xi(nip)
    REAL(GI), INTENT(IN)     :: yi(nip)
    INTEGER, INTENT(IN)  :: inxi(nip)
    INTEGER, INTENT(IN)  :: inyi(nip)
    REAL(GI), INTENT(OUT)    :: zi(nip)

    !     ..
    !     .. Local Scalars ..
    REAL(GI) :: a, b, c, W, dx, dxsq, dy, dysq, p00, p01, p02, p03, p10, p11,  &
    p12, p13, p20, p21, p22, p23, p30, p31, p32, p33, q0, q1, q2,  &
    q3, u, v, x0, xii, y0, yii, z00, z01, z0dx, z0dy, z10, z11,  &
    z1dx, z1dy, zdxdy, zii, zx00, zx01, zx0dy, zx10, zx11,  &
    zx1dy, zxy00, zxy01, zxy10, zxy11, zy00, zy01, zy0dx, zy10, zy11, zy1dx
    INTEGER :: iip, ixd0, ixd1, ixdi, ixdipv, iyd0, iyd1, iydi, iydipv
    !     ..
    !     .. Intrinsic Functions ..
    ! INTRINSIC        MAX
    !     ..

    ! Calculation
    ! Outermost DO-loop with respect to the output point
    DO  iip = 1,nip
        xii = xi(iip)
        yii = yi(iip)
        IF (iip == 1) THEN
            ixdipv = -1
            iydipv = -1
        ELSE
            ixdipv = ixdi
            iydipv = iydi
        END IF
        ixdi = inxi(iip)
        iydi = inyi(iip)

        ! Retrieves the z and partial derivative values at the origin of
        ! the coordinate for the rectangle.
        IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
            ixd0 = MAX(1,ixdi)
            iyd0 = MAX(1,iydi)
            x0 = xd(ixd0)
            y0 = yd(iyd0)
            z00 = zd(ixd0,iyd0)
            zx00 = pdd(1,ixd0,iyd0)
            zy00 = pdd(2,ixd0,iyd0)
            zxy00 = pdd(3,ixd0,iyd0)
        END IF

        ! Case 1.  When the rectangle is inside the data area in both the
        ! x and y directions.
        IF ((ixdi > 0 .AND. ixdi < nxd) .AND. (iydi > 0 .AND. iydi < nyd)) THEN
            ! Retrieves the z and partial derivative values at the other three
            ! vertices of the rectangle.
            IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
                ixd1 = ixd0 + 1
                dx = xd(ixd1) - x0
                dxsq = dx*dx
                iyd1 = iyd0 + 1
                dy = yd(iyd1) - y0
                dysq = dy*dy
                z10 = zd(ixd1,iyd0)
                z01 = zd(ixd0,iyd1)
                z11 = zd(ixd1,iyd1)
                zx10 = pdd(1,ixd1,iyd0)
                zx01 = pdd(1,ixd0,iyd1)
                zx11 = pdd(1,ixd1,iyd1)
                zy10 = pdd(2,ixd1,iyd0)
                zy01 = pdd(2,ixd0,iyd1)
                zy11 = pdd(2,ixd1,iyd1)
                zxy10 = pdd(3,ixd1,iyd0)
                zxy01 = pdd(3,ixd0,iyd1)
                zxy11 = pdd(3,ixd1,iyd1)
                ! Calculates the polynomial coefficients.
                z0dx = (z10-z00)/dx
                z1dx = (z11-z01)/dx
                z0dy = (z01-z00)/dy
                z1dy = (z11-z10)/dy
                zx0dy = (zx01-zx00)/dy
                zx1dy = (zx11-zx10)/dy
                zy0dx = (zy10-zy00)/dx
                zy1dx = (zy11-zy01)/dx
                zdxdy = (z1dy-z0dy)/dx
                a = zdxdy - zx0dy - zy0dx + zxy00
                b = zx1dy - zx0dy - zxy10 + zxy00
                c = zy1dx - zy0dx - zxy01 + zxy00
                W = zxy11 - zxy10 - zxy01 + zxy00
                p00 = z00
                p01 = zy00
                p02 = (2.0* (z0dy-zy00)+z0dy-zy01)/dy
                p03 = (-2.0*z0dy+zy01+zy00)/dysq
                p10 = zx00
                p11 = zxy00
                p12 = (2.0* (zx0dy-zxy00)+zx0dy-zxy01)/dy
                p13 = (-2.0*zx0dy+zxy01+zxy00)/dysq
                p20 = (2.0* (z0dx-zx00)+z0dx-zx10)/dx
                p21 = (2.0* (zy0dx-zxy00)+zy0dx-zxy10)/dx
                p22 = (3.0* (3.0*a-b-c)+W)/ (dx*dy)
                p23 = (-6.0*a+2.0*b+3.0*c-W)/ (dx*dysq)
                p30 = (-2.0*z0dx+zx10+zx00)/dxsq
                p31 = (-2.0*zy0dx+zxy10+zxy00)/dxsq
                p32 = (-6.0*a+3.0*b+2.0*c-W)/ (dxsq*dy)
                p33 = (2.0* (2.0*a-b-c)+W)/ (dxsq*dysq)
            END IF

            ! Evaluates the polynomial.
            u = xii - x0
            v = yii - y0
            q0 = p00 + v* (p01+v* (p02+v*p03))
            q1 = p10 + v* (p11+v* (p12+v*p13))
            q2 = p20 + v* (p21+v* (p22+v*p23))
            q3 = p30 + v* (p31+v* (p32+v*p33))
            zii = q0 + u* (q1+u* (q2+u*q3))
            ! End of Case 1

            ! Case 2.  When the rectangle is inside the data area in the x
            ! direction but outside in the y direction.
        ELSE IF ((ixdi > 0.AND.ixdi < nxd) .AND.  &
        (iydi <= 0.OR.iydi >= nyd)) THEN
            ! Retrieves the z and partial derivative values at the other
            ! vertex of the semi-infinite rectangle.
            IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
                ixd1 = ixd0 + 1
                dx = xd(ixd1) - x0
                dxsq = dx*dx
                z10 = zd(ixd1,iyd0)
                zx10 = pdd(1,ixd1,iyd0)
                zy10 = pdd(2,ixd1,iyd0)
                zxy10 = pdd(3,ixd1,iyd0)
                ! Calculates the polynomial coefficients.
                z0dx = (z10-z00)/dx
                zy0dx = (zy10-zy00)/dx
                p00 = z00
                p01 = zy00
                p10 = zx00
                p11 = zxy00
                p20 = (2.0* (z0dx-zx00)+z0dx-zx10)/dx
                p21 = (2.0* (zy0dx-zxy00)+zy0dx-zxy10)/dx
                p30 = (-2.0*z0dx+zx10+zx00)/dxsq
                p31 = (-2.0*zy0dx+zxy10+zxy00)/dxsq
            END IF
            ! Evaluates the polynomial.
            u = xii - x0
            v = yii - y0
            q0 = p00 + v*p01
            q1 = p10 + v*p11
            q2 = p20 + v*p21
            q3 = p30 + v*p31
            zii = q0 + u* (q1+u* (q2+u*q3))
            ! End of Case 2

            ! Case 3.  When the rectangle is outside the data area in the x
            ! direction but inside in the y direction.
        ELSE IF ((ixdi <= 0.OR.ixdi >= nxd) .AND.  &
        (iydi > 0 .AND. iydi < nyd)) THEN
            ! Retrieves the z and partial derivative values at the other
            ! vertex of the semi-infinite rectangle.
            IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
                iyd1 = iyd0 + 1
                dy = yd(iyd1) - y0
                dysq = dy*dy
                z01 = zd(ixd0,iyd1)
                zx01 = pdd(1,ixd0,iyd1)
                zy01 = pdd(2,ixd0,iyd1)
                zxy01 = pdd(3,ixd0,iyd1)
                ! Calculates the polynomial coefficients.
                z0dy = (z01-z00)/dy
                zx0dy = (zx01-zx00)/dy
                p00 = z00
                p01 = zy00
                p02 = (2.0*(z0dy-zy00)+z0dy-zy01)/dy
                p03 = (-2.0*z0dy+zy01+zy00)/dysq
                p10 = zx00
                p11 = zxy00
                p12 = (2.0*(zx0dy-zxy00) + zx0dy - zxy01)/dy
                p13 = (-2.0*zx0dy + zxy01 + zxy00)/dysq
            END IF

            ! Evaluates the polynomial.
            u = xii - x0
            v = yii - y0
            q0 = p00 + v* (p01 + v*(p02+v*p03))
            q1 = p10 + v* (p11 + v*(p12+v*p13))
            zii = q0 + u*q1
            ! End of Case 3

            ! Case 4.  When the rectangle is outside the data area in both the
            ! x and y direction.
        ELSE IF ((ixdi <= 0 .OR. ixdi >= nxd) .AND.  &
        (iydi <= 0 .OR. iydi >= nyd)) THEN
            ! Calculates the polynomial coefficients.
            IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
                p00 = z00
                p01 = zy00
                p10 = zx00
                p11 = zxy00
            END IF
            ! Evaluates the polynomial.
            u = xii - x0
            v = yii - y0
            q0 = p00 + v*p01
            q1 = p10 + v*p11
            zii = q0 + u*q1
        END IF
        ! End of Case 4
        zi(iip) = zii
    END DO

    RETURN
    END SUBROUTINE rgplnl

    end module Interpolation
