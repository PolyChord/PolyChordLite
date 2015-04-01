
    module Samples
    use ObjectLists
    use settings
    use Interpolation
    implicit none

    integer, parameter :: sample_prec = mcp

    Type, extends(TRealArrayList):: TSampleList
    contains
    procedure :: ConfidVal => TSampleList_ConfidVal
    end Type TSampleList

    Type TKernel1D
        real(mcp), allocatable :: Win(:), xWin(:)
        real(mcp) :: h
        integer :: winw
        real(mcp), allocatable :: a0(:), boundary_K(:), boundary_xK(:)
    contains
    procedure ::  Init => Kernel1D_Init
    end Type TKernel1D
    
    !Spline interpolated density function
    Type TDensity1D
        integer n
        real(mcp) :: spacing
        real(mcp), dimension(:), allocatable :: X
        real(mcp), dimension(:), allocatable :: P, ddP
    contains
    procedure :: Prob => Density1d_prob
    procedure :: Init => Density1D_Init
    procedure :: Free => Density1d_Free
    procedure :: InitSpline => Density1D_initSpline
    procedure :: Limits => Density1D_Limits
    end Type TDensity1D


    contains

    function GelmanRubinEvalues(cov, meanscov, evals, num) result(OK)
    use MatrixUtils
    integer, intent(in) :: num
    real(sample_prec) :: cov(num,num), meanscov(num,num), evals(num)
    integer jj, error
    real(sample_prec) rot(num,num), rotmeans(num,num)
    real(sample_prec) :: sc
    logical OK

    rot = cov
    rotmeans = meanscov
    do jj=1,num
        sc = sqrt(cov(jj,jj))
        rot(jj,:) = rot(jj,:) / sc
        rot(:,jj) = rot(:,jj) / sc
        rotmeans(jj,:) = rotmeans(jj,:) /sc
        rotmeans(:,jj) = rotmeans(:,jj) /sc
    end do

    call Matrix_CholeskyRootInverse(rot, error=error)
    OK = error==0
    if (OK) then
        rotmeans =  matmul(matmul(rot, rotmeans), transpose(rot))
        call Matrix_Diagonalize(rotmeans, evals, num)
    end if

    end function GelmanRubinEvalues


    subroutine TSampleList_ConfidVal(L, ix, limfrac, ix1, ix2, Lower, Upper)
    !Taking the ix'th entry in each array to be a sample, value for which
    !limfrac of the items between ix1 and ix2 (inc) are above or below
    !e.g. if limfrac = 0.05 get two tail 90% confidence limits
    Class(TSampleList) :: L
    Type(TSampleList) :: SortItems
    integer, intent(IN) :: ix
    real(sample_prec), intent(IN) :: limfrac
    real(sample_prec), intent(OUT), optional :: Lower, Upper
    integer, intent(IN), optional :: ix1,ix2
    integer b,t,samps
    real(sample_prec) pos, d

    b=1
    t=L%Count
    if (present(ix1)) b = ix1
    if (present(ix2)) t = ix2
    samps = t - b + 1
    call SortItems%AssignPointers(L, b, t)
    call SortItems%SortArr(ix)

    if (present(Lower)) then
        pos = (samps-1)*limfrac + 1
        b = max(int(pos),1)
        Lower = SortItems%Value(b, ix)
        if (b < samps .and. pos>b) then
            d = pos - b
            Lower = Lower*(1 - d) + d * SortItems%Value(b+1,ix)
        end if
    end if
    if (present(Upper)) then
        pos = (samps-1)*(1.-limfrac) + 1
        b = max(int(pos),1)
        Upper = SortItems%Value(b,ix)
        if (b < samps .and. pos>b) then
            d = pos - b
            Upper = Upper*(1 - d) + d * SortItems%Value(b+1,ix)
        end if
    end if
    call SortItems%Clear(itemsOnly=.true.)

    end subroutine TSampleList_ConfidVal

    function Density1D_prob(D, x)
    Class(TDensity1D) :: D
    real(sample_prec) :: Density1D_prob
    real(sample_prec), intent(in) :: x
    integer llo,lhi
    real(sample_prec) a0,b0

    if (x>D%X(D%n) - D%spacing/1d6) then
        if (x > D%X(D%n) + D%spacing/1d6) then
            write (*,*) 'Density: x too big ', x
            stop
        end if
        Density1D_prob=D%P(D%n)
        return
    end if
    if (x< D%X(1)- D%spacing/1d6) then
        write (*,*) 'Density: out of range ', x, D%X(1), D%X(D%n)
        stop
    else
    end if
    llo = 1 + max(0,int((x-D%X(1))/D%spacing))
    lhi=llo+1
    a0=(D%X(lhi)-x)/D%spacing
    b0=(x-D%X(llo))/D%spacing
    Density1D_prob = a0*D%P(llo)+ b0*D%P(lhi)+((a0**3-a0)* D%ddP(llo) +(b0**3-b0)*D%ddP(lhi))*D%spacing**2/6

    end function Density1D_prob

    subroutine Density1D_Init(D,n, spacing)
    Class(TDensity1D) :: D
    integer, intent(in) :: n
    real(sample_prec), intent(in) :: spacing

    call D%Free()
    D%n=n
    allocate(D%X(n))
    allocate(D%P(n))
    allocate(D%ddP(n))
    D%spacing = spacing

    end subroutine Density1D_Init

    subroutine Density1D_Free(D)
    Class(TDensity1D) :: D

    if (allocated(D%X)) then
        deallocate(D%X,D%P,D%ddP)
    end if
    end subroutine Density1D_Free

    subroutine Density1D_initSpline(D)
    Class(TDensity1D) :: D

    call spline(D%x,D%P,D%n,SPLINE_DANGLE,SPLINE_DANGLE,D%ddP)
    end subroutine Density1D_initSpline

    subroutine Density1D_Limits(D, p, mn, mx, lim_bot,lim_top)
    class(TDensity1D) :: D
    integer i, bign
    real(mcp), intent(in) :: p
    logical,intent(out) :: lim_bot,lim_top
    real(mcp), intent(out) :: mn,mx
    integer, parameter :: factor = 100
    real(mcp) norm, try, try_sum, try_t,try_b,try_last
    real(mcp), allocatable :: grid(:)

    bign=(D%n-1)*factor + 1
    allocate(grid(bign))
    do i=1, bign
        grid(i) = D%Prob(D%X(1) + (i-1)*D%spacing/factor)
    end do
    norm = sum(grid)
    norm = norm - 0.5*D%P(D%n) - 0.5*D%P(1)

    try_t = maxval(grid)
    try_b = 0
    try_last = -1
    do
        try = (try_b + try_t)/2
        try_sum = sum(grid,mask = grid >=  try)
        if (try_sum < p*norm) then
            try_t = (try_b+try_t)/2
        else
            try_b = (try_b+try_t)/2
        end if
        if (abs(try_sum/try_last - 1) < 0.0001) exit
        try_last = try_sum
    end do
    try = (try_b+try_t)/2
    lim_bot = grid(1) >= try
    if (lim_bot) then
        mn = D%P(1)
    else
        do i=1, bign
            if (grid(i) > try) then
                mn = D%X(1) + (i-1)*D%spacing/factor
                exit
            end if
        end do
    end if
    lim_top=grid(bign) >= try
    if (lim_top) then
        mx = D%P(D%n)
    else
        do i=bign,1,-1
            if (grid(i) > try) then
                mx = D%X(1) + (i-1)*D%spacing/factor
                exit
            end if
        end do
    end if

    end subroutine Density1D_Limits


    subroutine Kernel1D_Init(this,winw, h, has_boundary)    
    class(TKernel1D) :: this
    integer, intent(in) :: winw
    real(mcp), intent(in) :: h
    logical, intent(in) :: has_boundary
    real(mcp), allocatable :: a1(:), a2(:)
    integer i
    
    this%winw=winw
    this%h = h
    allocate(this%Win(-winw:winw))
    do i=-winw,winw
        this%Win(i)=exp(-(i/h)**2/2)
    end do
    this%Win=this%Win/sum(this%win)
    if (has_boundary) then
        allocate(this%xWin(-winw:winw))
        allocate(this%a0(0:winw))
        allocate(a1(0:winw), a2(0:winw))
        allocate(this%boundary_K(0:winw))
        allocate(this%boundary_xK(0:winw))
        a2(0)=0
        do i=-winw,0
            this%xWin(i) = this%Win(i)*i
            a2(0) = a2(0) + this%Win(i)*i**2
        end do
        this%a0(0) = sum(this%Win(-winw:0))
        a1(0) = sum(this%xWin(-winw:0)) 
        do i=1,winw
            this%xWin(i) = this%Win(i)*i
            this%a0(i) = this%a0(i-1)+ this%Win(i)
            a1(i) = a1(i-1) + this%xWin(i)
            a2(i) = a2(i-1) + this%xWin(i)*i
        end do
        this%boundary_K = a2/ (this%a0*a2-a1**2)
        this%boundary_xK = a1/ (this%a0*a2-a1**2)
    end if
    
    end subroutine Kernel1D_Init 

    end module Samples
