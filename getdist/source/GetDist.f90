    ! Program to process .txt chain files produced by CosmoMC (cosmologist.info/cosmomc)
    ! Calculates statistics, thins, produces 1D and 2D marginalized plots and scatter plots
    ! and writes out useful information.
    ! Option to adjust priors and map parameters - write code below

    !Output files are (depending on options)

    ! * file_root.margestats - the mean, limits and stddev from the 1D marginalized distributions
    ! * file_root.likestats  - the best fit, and limits from the N-D confidence region
    ! * file_root.m          - MatLab .m file to produce 1D marginalized plots
    ! * file_root_2D.m       - MatLab .m file to produce 2D marginalized plots
    ! * file_root_3D.m       - MatLab .m file to produce 2D sample plots colored by a third parameter
    ! * file_root_tri.m      - MatLab .m file to produce triangle plots (1D on diagonal, 2D off-diagonal)
    ! * file_root_thin.txt   - thinned merged versions of input files
    ! * file_root.covmat     - covariance matrix for parameters
    ! * file_root.corr       - covariance matrix for parameters normalized to have unit variance
    ! * file_root.PCA        - human-readable file giving details of PCA, constrained parameters, etc
    ! * file_root.converge   - Various statistics to help assess chain convergence/sampling error

    ! You can easily edit the .m files produced to produce custom layouts of plots, etc.
    ! Data for plots are exported to the plot_data_dir folder, other files to out_dir

    !March 03: fixed bug computing the limits in the .margestats file
    !May 03: Added support for triangle plots
    !Dec 03: Support for non-chain samples, auto-correlation convergence, auto_label,
    !         split tests to quanify sampling error on upper and lower quantiles
    !Jan 05: MatLab 7 options, various minor changes
    !Jul 05: Added limits info to .margestats output
    !Apr 06: Fixed some version confusions
    !Aug 06: speeded 2D plotting
    !Oct 06: Added plot_data_dir and out_dir
    !Nov 07: uses Matrix_utils, fixed .ps file output filenames,
    !        added markerxxx for vertical lines in 1D matlabl plots
    !May 08: Option to process WMAP5-formatted chains (thanks to Mike Nolta)
    !        Added font_scale to scale default font sizes, num_contours input parameter (allows more than 2)
    !Sept 09: use of .paramnames and parameter_names file for names and labels; referencing by name as alternative to number
    !         allowed map_params with 1-column input format
    !         plotparams parameter ordering is preserved (e.g. to change order of 1D plots)
    !Oct 09:  fixed bug in credible intervals with prior cutoffs (Jan Hammann)
    !May 10:  added finish_run_command to run system command when getdist finished (e.g. "matlab < %ROOTNAME%.m")
    !          added line_labels to write out roots of lines being plotted (matlab)
    !          added no_triangle_axis_labels
    !Jan 11:  fix to confidence interval calcaultion for obscure cases with very empty tails (thanks Andrew Fowlie)
    !Jan 12:  increased precision of output files
    !Oct 12:  plot_output for pdf,eps,ps output support; subplot_size_inch
    !         fix for more than 100 parameters
    !Nov-Dec 12: removed sm support; use specific parameter in triangle plots; new R-1 definition, etc...
    !.. Mar 13: numerous changes..
    !   Apr 13: output .py files, converge_test_limit, etc.
    !Oct 13: add other mean likelihoods to .likestats
    module MCSamples
    use settings
    use MatrixUtils
    use ParamNames
    use Samples
    use ArrayUtils
    use RandUtils
    use FileUtils
    use ObjectLists
    implicit none

    Type(TParamNames), save :: NameMapping

    integer, parameter :: gp = KIND(1.d0)
    character(LEN=*), parameter :: float_format = '(*(E16.7))'


    ! #define coldata(a, b) coldata(b, a)
    !uncomment to switch order for compilation so confid_val is faster - suggested by Reijo Keskitalo March07

    real(gp), dimension(:,:), allocatable :: coldata  !(col_index, row_index)
    integer, dimension(:), allocatable :: thin_ix
    integer thin_rows

    integer, parameter :: max_rows = 500000
    integer, parameter :: max_cols = 500
    integer, parameter :: max_chains = 100
    integer, parameter :: max_split_tests = 4
    integer, parameter :: max_contours = 5
    integer :: corr_length_thin = 0
    integer :: corr_length_steps = 15

    Type(TDensity1D) :: Density1D

    integer chain_indices(max_chains), num_chains_used
    integer bestfit_ix
    real(mcp) meanlike, maxlike

    integer nrows, ncols, num_bins, num_bins_2D
    real(mcp) numsamp, max_mult, mean_mult
    integer ND_cont1, ND_cont2
    real(mcp), allocatable :: ND_limit_top(:,:), ND_limit_bot(:,:)
    integer covmat_dimension
    integer colix(max_cols), num_vars  !Parameters with non-blank labels
    real(mcp) prior_ranges(2,max_cols)
    real(mcp) mean(max_cols), sddev(max_cols)
    real(mcp), dimension(:,:), allocatable :: corrmatrix
    character(LEN=:), allocatable :: plot_data_dir, out_dir, rootdirname, in_root
    character(LEN=128) labels(max_cols)
    character(LEN=128) pname(max_cols)
    logical has_limits(max_cols), has_limits_top(max_cols),has_limits_bot(max_cols)
    !has limits for plotting, to avoid underweighting at edges
    logical marge_limits_top(max_contours, max_cols),marge_limits_bot(max_contours, max_cols)
    !has limits for 1 vs 2 tail outputs in margestats
    logical has_markers(max_cols)
    real(mcp) markers(max_cols)
    integer num_contours
    character(LEN=3) :: plot_ext = 'm'
    real(mcp) matlab_version
    character(LEN=3) :: plot_output = 'pdf'
    real :: subplot_size_inch = 3.0
    real :: subplot_size_inch2, subplot_size_inch3

    logical :: matlab_latex = .false.
    real(mcp) :: font_scale = 1.
    real(mcp) contours(max_contours)
    real(mcp) :: max_frac_twotail(max_contours) !how small the end bin must be relative to max to use two tail
    integer :: indep_thin = 0
    integer chain_numbers(max_chains)
    logical isused(max_cols)
    logical force_twotail
    logical  plot_meanlikes
    logical :: mean_loglikes = .false.
    logical shade_meanlikes, make_single_samples
    integer single_thin,cust2DPlots(max_cols**2)
    integer ix_min(max_cols),ix_max(max_cols)
    real(mcp) limmin(max_cols),limmax(max_cols)
    real(mcp) center(max_cols), param_min(max_cols), param_max(max_cols), range_min(max_cols), range_max(max_cols)
    logical BW,do_shading
    Type(TStringList), save :: ComparePlots
    logical :: prob_label = .false.
    logical :: plots_only, no_plots
    real(mcp) :: smooth_scale_1D=-1.d0, smooth_scale_2D = 1.d0
    real(mcp) :: credible_interval_threshold = 0.05d0
    integer :: boundary_correction_method = 1  !0 old basic norm correction, 1 linear boundary kernel

    contains

    subroutine AdjustPriors
    !Can adjust the multiplicity of each sample in coldata(1, rownum) for new priors
    !Be careful as this code is parameterisation dependent
    !  integer i, ix
    !   real(mcp) ombh2, chisq

    stop 'You need to write the AdjustPriors subroutine in GetDist.f90 first!'

    !   write (*,*) 'Adjusting priors'
    !   ix = NameMapping%Index('omegabh2')
    !   do i=0, nrows-1
    !!E.g. ombh2 prior
    !      ombh2 = coldata(ix+2,i)
    !      chisq = (ombh2 - 0.0213)**2/0.001 **2
    !      coldata(1,i) = coldata(1,i)*exp(-chisq/2)
    !      coldata(2,i) = coldata(2,i) + chisq/2
    !
    !   end do

    end subroutine AdjustPriors

    subroutine MapParameters(invars)
    real(gp) invars(1:ncols)
    ! map parameters in invars: eg. invars(3)=invars(3)*invars(4)

    !    invars(2+13)=invars(17+2)*exp(-invars(2+4))
    stop 'Need to write MapParameters routine first'

    end subroutine MapParameters

    subroutine CoolChain(cool)
    real(mcp), intent(in) :: cool
    integer i
    real(mcp) maxL, newL

    write (*,*) 'Cooling chains by ', cool
    MaxL = minval(coldata(2,0:nrows-1))
    do i=0, nrows-1
        newL = coldata(2,i)*cool
        coldata(1,i) = coldata(1,i)*exp(-(newL - coldata(2,i)) - MaxL*(1-cool) )
        coldata(2,i) = newL
    end do

    end subroutine CoolChain

    subroutine DeleteZeros
    integer i,ii

    ii=0
    do i=0, nrows-1
        if (coldata(1,i)/=0) then
            coldata(:,ii) = coldata(:,i)
            ii=ii+1
        end if
    end do
    if (ii==0) stop 'Prior has removed all models!'
    if (ii /= nrows) write (*,*) 'Prior removed ',nrows-ii,'models'
    nrows = ii

    end subroutine DeleteZeros

    subroutine SortColData(bycol)
    !Sort coldata in order of likelihood
    integer, intent(in) :: bycol
    integer i
    Type(TSampleList) :: S

    call S%SetCapacity(nrows)
    do i = 0, nrows -1
        call S%Add(coldata(1:ncols,i))
    end do

    call S%SortArr(bycol)
    do i = 0, nrows-1
        coldata(1:ncols,i) =  S%Item(i+1)
    end do
    call S%Clear()
    end  subroutine SortColData


    subroutine MakeSingleSamples(single_thin)
    !Make file of weight-1 samples by choosing samples with probability given by their weight
    integer  i, single_thin
    real(mcp) maxmult
    Type(TTextFile) :: F
    Feedback = 0
    call initRandom()

    call F%CreateFile(plot_data_dir//trim(rootname)//'_single.txt')
    maxmult = maxval(coldata(1,0:nrows-1))
    do i= 0, nrows -1
        if (ranmar() <= coldata(1,i)/maxmult/single_thin) &
            write (F%unit,float_format) 1.0, coldata(2,i), coldata(colix(1:num_vars),i)
    end do
    call F%Close()

    end subroutine MakeSingleSamples

    subroutine WriteThinData(fname,cool)
    character(LEN=*), intent(in) :: fname
    real(mcp),intent(in) :: cool
    integer i
    real(mcp) MaxL, NewL
    Type(TTextFile) :: F

    if (cool /= 1) write (*,*) 'Cooled thinned output with temp: ', cool

    MaxL = minval(coldata(2,0:nrows-1))

    call F%CreateFile(fname)
    do i=0, thin_rows -1
        if (cool/=1) then
            newL = coldata(2,thin_ix(i))*cool
            write (F%unit,float_format) exp(-(newL - coldata(2,thin_ix(i))) - MaxL*(1-cool) ), newL, coldata(3:ncols,thin_ix(i))
        else
            write (F%unit,float_format) 1., coldata(2:ncols,thin_ix(i))
        end if
    end do
    write (*,*) 'Wrote ',thin_rows, 'thinned samples'
    call F%Close()

    end subroutine WriteThinData


    subroutine ThinData(fac, ix1,ix2)
    !Make thinned samples

    integer, intent(in) :: fac
    integer, intent(in), optional :: ix1,ix2
    integer i, tot, nout, nend, mult

    if (allocated(thin_ix)) deallocate(thin_ix)
    allocate(thin_ix(0:nint(numsamp)/fac))

    tot= 0
    nout=0
    i = 0
    if (present(ix1)) i=ix1
    nend = nrows
    if (present(ix2)) nend = ix2+1

    mult = coldata(1,i)
    do while (i< nend)
        if (abs(nint(coldata(1,i)) - coldata(1,i)) > 1e-4) &
            stop 'non-integer weights in ThinData'

        if (mult + tot < fac) then
            tot = tot + mult
            i=i+1
            if (i< nend) mult = nint(coldata(1,i))
        else
            thin_ix(nout) = i
            nout= nout+1
            if (mult == fac - tot) then
                i=i+1
                if (i< nend) mult = nint(coldata(1,i))
            else
                mult = mult - (fac -tot)
            end if
            tot = 0
        end if
    end do

    thin_rows = nout

    end subroutine ThinData


    subroutine GetCovMatrix
    use IO
    integer i, j, nused
    real(mcp) mean(max_cols)
    real(mcp) scale
    real(mcp), dimension(:,:), allocatable :: covmatrix
    integer used_ix(max_cols)
    character(LEN=:), allocatable :: outline

    allocate(corrmatrix(ncols-2,ncols-2))

    corrmatrix = 0
    nused = 0
    do i=1, ncols-2
        if (isused(i+2)) then
            if (i <= NameMapping%num_MCMC) then
                nused = nused + 1
                used_ix(nused)=i
            end if
            mean(i) = sum(coldata(1,0:nrows-1)*coldata(i+2,0:nrows-1))/numsamp
        end if
    end do

    do i = 1, ncols-2
        if (isused(i+2)) then
            do j = i, ncols-2
                if (isused(j+2)) then
                    corrmatrix(i,j) = sum(coldata(1,0:nrows-1)*((coldata(2+i,0:nrows-1)-mean(i))* &
                        (coldata(2+j,0:nrows-1)-mean(j))))/numsamp
                    corrmatrix(j,i) = corrmatrix(i,j)
                end if
            end do
        end if
    end do

    if (NameMapping%nnames/=0) then
        outline=''
        do i=1, nused
            outline = outline//' '//NameMapping%NameAtIndex(used_ix(i))
        end do
        allocate(covmatrix(nused,nused))
        covmatrix = corrmatrix(used_ix(1:nused),used_ix(1:nused))
        call IO_WriteProposeMatrix(covmatrix,trim(rootdirname) //'.covmat', outline)
        !               call Matrix_write(trim(rootdirname) //'.covmat',covmatrix,.true.,commentline=outline)
        deallocate(covmatrix)
    else
        if (covmat_dimension /= 0 .and. .not. plots_only) then
            if (covmat_dimension > ncols -2) stop 'covmat_dimension larger than number of parameters'
            write (*,*) 'Writing covariance matrix for ',covmat_dimension,' parameters'
            allocate(covmatrix(covmat_dimension,covmat_dimension))
            covmatrix=corrmatrix(1:covmat_dimension,1:covmat_dimension)
            call IO_WriteProposeMatrix(covmatrix,trim(rootdirname) //'.covmat')
            !               call Matrix_write(trim(rootdirname) //'.covmat',covmatrix,.true.)
            deallocate(covmatrix)
        end if
    end if

    do i=1, ncols-2
        if (corrmatrix(i,i) > 0) then
            scale = sqrt(corrmatrix(i,i))
            corrmatrix(i,:) = corrmatrix(i,:)/scale
            corrmatrix(:,i) = corrmatrix(:,i)/scale
        end if
    end do
    if (.not. plots_only) call Matrix_write(trim(rootdirname) //'.corr',corrmatrix, .true.)

    end subroutine GetCovMatrix


    function MostCorrelated2D(i1,i2, direc)
    integer, intent(in) :: i1,i2, direc
    integer MostCorrelated2D
    !Find which parameter is most correllated with the degeneracy in ix1, ix2
    integer pars(2)
    real(mcp) mat2D(2,2), evals(2), u(2,2)
    real(mcp) corrs(2,ncols-2)

    if (direc /= 0 .and. direc /= -1) stop 'Invalid 3D color parameter'

    pars(1)=i1
    pars(2)=i2
    mat2D = corrmatrix(pars,pars)
    u = mat2D

    call Matrix_Diagonalize(u, evals, 2)
    corrs = matmul(transpose(u), corrmatrix(pars,:))
    corrs(:,pars) = 0

    MostCorrelated2D = MaxIndex(real(abs(corrs(2+direc,colix(1:num_vars)-2))), num_vars)
    MostCorrelated2D = colix(MostCorrelated2D) -2

    end function MostCorrelated2D

    subroutine GetFractionIndices(fraction_indices,n)
    integer, intent(in) :: n
    integer num, fraction_indices(*), i
    real(mcp) tot, aim

    tot = 0
    aim=numsamp/n
    num = 1
    fraction_indices(1) = 0
    fraction_indices(n+1) = nrows
    do i=0, nrows-1
        tot = tot + coldata(1,i)
        if (tot > aim) then
            num=num+1
            fraction_indices(num) = i
            if (num==n) exit
            aim=aim+numsamp/n
        end if
    end do

    end subroutine GetFractionIndices

    function ConfidVal(ix,limfrac,upper,ix1,ix2)
    integer, intent(IN) :: ix
    real(mcp), intent(IN) :: limfrac
    logical, intent(IN) :: upper
    integer, intent(IN), optional :: ix1,ix2
    real(mcp) ConfidVal, samps
    real(mcp) try, lasttry, try_t, try_b
    integer l,t
    !find upper and lower bounds

    l=0
    t=nrows-1
    if (present(ix1)) l=ix1
    if (present(ix2)) t = ix2
    try_b = minval(coldata(ix,l:t))
    try_t = maxval(coldata(ix,l:t))

    samps = sum(coldata(1,l:t))

    lasttry = -1
    if (upper) then
        do
            try = sum(coldata(1,l:t),mask = coldata(ix,l:t) > (try_b + try_t)/2)
            if (try > samps*limfrac) then
                try_b = (try_b+try_t)/2
            else
                try_t = (try_b+try_t)/2
            end if
            if (try == lasttry .and. (abs(try_b + try_t)< 1e-10_mcp &
                .or. abs((try_b - try_t)/(try_b + try_t)) < 0.05_mcp )) exit
            lasttry = try
        end do
    else
        do
            try = sum(coldata(1,l:t),mask = coldata(ix,l:t) < (try_b + try_t)/2)
            if (try > samps*limfrac) then
                try_t = (try_b+try_t)/2
            else
                try_b = (try_b+try_t)/2
            end if
            if (try == lasttry .and. (abs(try_b + try_t)< 1e-10_mcp &
                .or. abs((try_b - try_t)/(try_b + try_t)) < 0.05_mcp )) exit
            lasttry = try
        end do
    end if

    ConfidVal = try_t
    end function ConfidVal

    subroutine WriteS(S)
    character(LEN=*), intent(in) :: S

    write (*,*) trim(S)

    end subroutine WriteS

    subroutine PCA(pars,n,param_map,normparam_num)
    !Perform principle component analysis
    !In other words, get eigenvectors and eigenvalues for normalized variables
    !with optional (log) mapping
    integer, intent(in) :: n, pars(n), normparam_num
    character(LEN=*), intent(IN) ::  param_map
    integer i, j, normparam, locarr(1:1)
    character(LEN =60) fmt, tmpS

    real(mcp) PCmean(n), sd(n), newmean(n), newsd(n)
    real(mcp) evals(n)
    real(mcp) corrmatrix(n,n), u(n,n)
    real(mcp), dimension(:,:), allocatable :: PCdata
    character(LEN=100) PClabs(n), div, expo
    logical doexp
    Type(TTextFile) :: F

    write (*,*) 'Doing PCA for ',n,' parameters'
    call F%CreateFile(trim(rootdirname) //'.PCA')
    write (F%unit,*) 'PCA for parameters: '

    if (normparam_num /=0) then
        normparam = IndexOf(normparam_num,pars,n)
        if (normparam ==0) stop 'Invalid PCA normalization parameter'
    else
        normparam = 0
    end if

    allocate(PCdata(n,0:nrows-1))

    PCdata(:,:) = coldata(2+pars,0:nrows-1)
    fmt = trim(numcat('',n))//'f8.4)'

    doexp =.false.
    do i = 1, n
        if (param_map(i:i)=='L') then
            doexp = .true.
            PCdata(i,:) = log(PCdata(i,:))
            PClabs(i) = 'ln(' //trim(labels(pars(i)+2))//')'
        elseif (param_map(i:i)=='M') then
            doexp = .true.
            PCdata(i,:) = log(-PCdata(i,:))
            PClabs(i) = 'ln(-' //trim(labels(pars(i)+2))//')'
        else
            PClabs(i) = trim(labels(pars(i)+2))
        end if
        write (F%unit,*) pars(i),':',trim(PClabs(i))
        PCmean(i) = sum(coldata(1,0:nrows-1)*PCdata(i,:))/numsamp
        PCdata(i,:) = PCdata(i,:) - PCmean(i)
        sd(i) = sqrt(sum(coldata(1,0:nrows-1)*PCdata(i,:)**2)/numsamp)
        if (sd(i)/=0) PCdata(i,:) = PCdata(i,:)/sd(i)
        corrmatrix(i,i) = 1
    end do

    call F%NewLine()
    call F%Write('Correlation matrix for reduced parameters')
    do i = 1, n
        do j = i, n
            corrmatrix(i,j) = sum(coldata(1,0:nrows-1)*PCdata(i,:)*PCdata(j,:))/numsamp
            corrmatrix(j,i) = corrmatrix(i,j)
        end do
        write (F%unit,'(1I4,'': '','//trim(fmt)) pars(i), corrmatrix(i,:)
    end do

    u = corrmatrix
    call Matrix_Diagonalize(u, evals, n)

    call F%NewLine()
    call F%Write('e-values of correlation matrix')
    do i = 1, n
        write (F%unit,'(''PC'',1I2,'': '','//trim(fmt)) i, evals(i)
    end do

    call F%NewLine()
    call F%Write('e-vectors')

    do i = 1, n
        write (F%unit,'(1I3,'': '','//trim(fmt)) pars(i), u(i,:)
    end do

    if (normparam /= 0) then
        !Set so parameter normparam has exponent 1
        do i=1, n
            u(:,i) = u(:,i)/u(normparam,i)*sd(normparam)
        end do
    else
        !Normalize so main component has exponent 1
        do i=1, n
            locarr(1:1) = maxloc(abs(u(:,i)))
            u(:,i) = u(:,i)/u(locarr(1),i)*sd(locarr(1))
        end do
    end if

    do i = 0, nrows -1
        PCdata(:,i) = matmul(transpose(u), PCdata(:,i))
        if (doexp) PCdata(:,i) = exp(PCdata(:,i))
    end do

    call F%NewLine()
    call F%Write('Principle components')

    do i = 1, n
        tmpS = trim(numcat('PC',i))//' (e-value: '//trim(RealToStr(evals(i)))//')'
        !ifc gives recursive IO error if you use a write within a write,
        !even if separate or internal files
        write (F%unit,*) trim(tmpS)
        do j=1,n
            if (param_map(j:j)=='L' .or. param_map(j:j)=='M' ) then
                expo = RealToStr(1/sd(j)*u(j,i) )
                if (param_map(j:j)=='M') then
                    div = RealToStr( -exp(PCmean(j)))
                else
                    div = RealToStr( exp(PCmean(j)))
                end if
                write(F%unit,*) '['//trim(RealToStr(u(j,i)))//']   ('//trim(labels(pars(j)+2)) &
                    //'/'//trim(div)//')^{'//trim(expo)//'}'
            else
                expo = RealToStr(sd(j)/u(j,i))
                if (doexp) then
                    write(F%unit,*) '['//trim(RealToStr(u(j,i)))//']    exp(('// &
                        trim(labels(pars(j)+2))//'-'//trim(RealToStr(PCmean(j)))// ')/'// trim(expo)//')'
                else
                    write(F%unit,*) '['//trim(RealToStr(u(j,i)))//']   ('// &
                        trim(labels(pars(j)+2))//'-'//trim(RealToStr(PCmean(j)))// ')/'// trim(expo)
                end if
            end if
        end do
        newmean(i) = sum(coldata(1,0:nrows-1)*PCdata(i,:))/numsamp
        newsd(i) = sqrt(sum(coldata(1,0:nrows-1)*(PCdata(i,:)-newmean(i))**2)/numsamp)
        write (F%unit,*) '          = '//trim(RealToStr(newmean(i)))// ' +- '//trim(RealToStr(newsd(i)))
        write (F%unit,'('' ND limits: '',4f9.3)') minval(PCdata(i,0:ND_cont1)),maxval(PCdata(i,0:ND_cont1)), &
            minval(PCdata(i,0:ND_cont2)),maxval(PCdata(i,0:ND_cont2))
        write (F%unit,*) ''
    end do

    !Find out how correlated these components are with other parameters
    write (F%unit,*) 'Correlations of principle components'

    write (F%unit,trim(numcat('(''    '',',n))//'I8)') (I, I=1, n)

    do i=1, n
        PCdata(i,:) = (PCdata(i,:) - newmean(i)) /newsd(i)
    end do

    do j=1,n
        write (F%unit,'(''PC'',1I2)', advance = 'no') j
        do i=1,n
            write (F%unit,'(1f8.3)', advance ='no') sum(coldata(1,0:nrows-1)*PCdata(i,:)* &
                PCdata(j,:))/numsamp
        end do
        write (F%unit,*) ''
    end do

    do j=1, num_vars
        write (F%unit,'(1I4)', advance = 'no') colix(j)-2
        do i=1,n
            write (F%unit,'(1f8.3)', advance ='no') sum(coldata(1,0:nrows-1)*PCdata(i,:)* &
                (coldata(colix(j),0:nrows-1)-mean(j))/sddev(j))/numsamp
        end do
        write (F%unit,*) '  ('//trim(labels(colix(j)))//')'
    end do

    call F%Close()
    deallocate(PCdata)

    end subroutine PCA

    subroutine GetUsedCols
    integer j

    do j = 3, ncols
        isused(j) = any(coldata(j,0:nrows-1)/=coldata(j,0))
    end do
    end subroutine GetUsedCols

    subroutine DoConvergeTests(limfrac)
    real(mcp), intent(in) ::limfrac !the limit to use for split tests and Raftery-Lewis
    real(mcp) chain_means(max_chains,max_cols),chain_samp(max_chains)
    real(mcp) between_chain_var(max_cols), in_chain_var(max_cols)
    integer frac(max_split_tests+1), split_n, chain_start(max_chains)
    integer i,j,k,jj,kk,ix, maxoff
    real(mcp) split_tests(max_split_tests)
    real(mcp) mean(max_cols), fullmean(max_cols),fullvar(max_cols)
    real(mcp), parameter :: cutfrac = 0._mcp
    real(mcp) usedsamps,evals(max_cols),R, maxsamp
    real(mcp), dimension(:,:), allocatable :: cov, meanscov
    integer usedvars(max_cols), num, thin_fac(max_chains), markov_thin(max_chains)
    real(gp) u, g2, fitted, focus, alpha, beta, probsum, tmp1
    real(mcp) confid
    integer i1,i2,i3, nburn(max_chains), hardest, endb,hardestend
    integer tran(2,2,2), tran2(2,2), off
    integer, dimension(:), allocatable :: binchain
    real(mcp), parameter :: epsilon = 0.001
    double precision, dimension(:,:), allocatable :: corrs
    character(LEN=10) :: typestr
    integer autocorr_thin
    logical :: invertible
    Type(TTextFile) :: F
    ! Get statistics for individual chains, and do split tests on the samples

    call F%CreateFile(trim(rootdirname) //'.converge')

    if (num_chains_used > 1) write (*,*) 'Number of chains used =  ',num_chains_used

    chain_indices(num_chains_used+1) = nrows
    do i=1, num_chains_used
        chain_start(i) =  chain_indices(i) + nint((chain_indices(i+1)-chain_indices(i))*cutfrac)
        chain_samp(i) = sum(coldata(1,chain_start(i):chain_indices(i+1)-1))
    end do
    usedsamps = sum(chain_samp(1:num_chains_used))
    maxsamp = maxval(chain_samp(1:num_chains_used))
    num = 0
    do j = 3, ncols
        if (isused(j)) then
            mean(j)=0
            fullmean(j) =  sum(coldata(1,0:nrows-1)*coldata(j,0:nrows-1))/numsamp

            if (num_chains_used> 1) then
                num=num+1
                usedvars(num)=j
                do i=1, num_chains_used
                    mean(j) = mean(j) + sum(coldata(1,chain_start(i):chain_indices(i+1)-1)*&
                        coldata(j,chain_start(i):chain_indices(i+1)-1))
                end do
                mean(j) = mean(j)/usedsamps
            end if
        end if
    end do

    do j = 3, ncols
        if (isused(j)) &
            fullvar(j)=  sum(coldata(1,0:nrows-1)*(coldata(j,0:nrows-1)-fullmean(j))**2)/numsamp
    end do


    if (num_chains_used > 1) then
        write (F%unit,*) ''
        write(F%unit,*)  'Variance test convergence stats using remaining chains'
        write (F%unit,*) 'param var(chain mean)/mean(chain var)'
        write (F%unit,*) ''

        do j = 3, ncols
            between_chain_var(j) = 0
            in_chain_var(j) = 0

            if (isused(j)) then
                if (num_chains_used > 1) then
                    !Get stats for individual chains - the variance of the means over the mean of the variances
                    do i=1, num_chains_used
                        chain_means(i,j) =   sum(coldata(1,chain_start(i):chain_indices(i+1)-1)* &
                            coldata(j,chain_start(i):chain_indices(i+1)-1))/chain_samp(i)

                        between_chain_var(j) = between_chain_var(j) + &
                            ! chain_samp(i)/maxsamp* & !Weight for different length chains
                            (chain_means(i,j) - mean(j))**2

                        in_chain_var(j) = in_chain_var(j) +  & !chain_samp(i)/maxsamp *&
                            sum(coldata(1,chain_start(i):chain_indices(i+1)-1)* &
                            (coldata(j,chain_start(i):chain_indices(i+1)-1)-chain_means(i,j))**2)
                    end do
                    between_chain_var(j) = between_chain_var(j)/(num_chains_used-1) !(usedsamps/maxsamp -1)
                    in_chain_var(j) = in_chain_var(j)/usedsamps
                    write (F%unit,'(1I3,f13.5,"  '//trim(labels(j))//'")') j-2, &
                        between_chain_var(j) /in_chain_var(j)
                end if
            end if
        end do
    end if

    if (num_chains_used > 1 .and. covmat_dimension>0) then
        !Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
        !c.f. Brooks and Gelman 1997

        do while (usedvars(num) > covmat_dimension+2)
            num= num - 1
        end do

        allocate(meanscov(num,num))
        allocate(cov(num,num))

        do jj=1,num
            j=usedvars(jj)
            do kk=jj, num
                k = usedvars(kk)
                meanscov(jj,kk) =0
                cov(jj,kk)= 0
                do i= 1, num_chains_used
                    cov(jj,kk) = cov(jj,kk) + &
                        !sqrt(chain_samp(jj)*chain_samp(kk))/maxsamp * &
                        sum(coldata(1,chain_start(i):chain_indices(i+1)-1)* &
                        (coldata(j,chain_start(i):chain_indices(i+1)-1)-chain_means(i,j))* &
                        (coldata(k,chain_start(i):chain_indices(i+1)-1)-chain_means(i,k)))

                    meanscov(jj,kk) = meanscov(jj,kk)+ & !sqrt(chain_samp(jj)*chain_samp(kk))/maxsamp * &
                        (chain_means(i,j)-mean(j))*(chain_means(i,k)-mean(k))
                end do
                meanscov(kk,jj) = meanscov(jj,kk)
                cov(kk,jj) = cov(jj,kk)
            end do
        end do
        meanscov = meanscov/(num_chains_used-1) !(usedsamps/maxsamp -1)
        cov = cov / usedsamps

        invertible = GelmanRubinEvalues(cov, meanscov, evals, num)
        if (invertible) then
            write (F%unit,*) ''
            write (F%unit,'(a)') 'var(mean)/mean(var) for eigenvalues of covariance of means of orthonormalized parameters'
            R = 0
            do jj=1,num
                write (F%unit,'(1I3,f13.5)') jj,evals(jj)
                R = max(R,evals(jj))
            end do
            !R is essentially the Gelman and Rubin statistic
            write (*,'(" var(mean)/mean(var), remaining chains, worst e-value: R-1 = ",f13.5)') R
            deallocate(cov,meanscov)
        else
            write(*,*) 'WARNING: Gelman-Rubin covariance not invertible'
        end if
    end if


    !Do tests for robustness under using splits of the samples
    !Return the rms ([change in upper/lower quantile]/[standard deviation])
    !when data split into 2, 3,.. sets
    write (F%unit,*) ''
    write (F%unit,*)  'Split tests: rms_n([delta(upper/lower quantile)]/sd) n={2,3,4}:'
    write(F%unit,*) 'i.e. mean sample splitting change in the quantiles in units of the st. dev.'
    write (F%unit,*) ''
    do j = 3, ncols
        if (isused(j)) then
            do endb =0,1
                !Just test all now
                !               if (endb==0 .and. has_limits_top(j))cycle
                !               if (endb==1 .and. has_limits_bot(j))cycle
                do split_n = 2,max_split_tests
                    call GetFractionIndices(frac,split_n)
                    split_tests(split_n) = 0
                    confid =  ConfidVal(j,(1-limfrac)/2,endb==0,0,nrows-1)
                    do i=1,split_n
                        split_tests(split_n) = split_tests(split_n) + &
                            (ConfidVal(j,(1-limfrac)/2,endb==0,frac(i),frac(i+1)-1)-confid)**2
                    end do !i
                    split_tests(split_n) = sqrt(split_tests(split_n)/split_n/fullvar(j))
                end do
                if (endb==0) then
                    typestr = 'upper'
                else
                    typestr = 'lower'
                end if
                write (F%unit,'(1I3,'//trim(IntToStr(max_split_tests-1)) // 'f9.4,"  ' &
                    //trim(labels(j))//' '//trim(typestr)//'")') j-2, split_tests(2:max_split_tests)

            end do !endb
        end if
    end do

    ! Now do Raftery and Lewis method
    ! See http://www.stat.washington.edu/tech.reports/raftery-lewis2.ps
    ! Raw non-importance sampled chains only
    if (all(abs(coldata(1,0:nrows-1) - nint(max(0.6_gp,coldata(1,0:nrows-1))))<1e-4)) then
        nburn = 0
        hardest=-1
        hardestend=0
        do ix=1, num_chains_used
            thin_fac(ix) = nint(maxval(coldata(1,chain_indices(ix):chain_indices(ix+1)-1)))

            do j = 3, covmat_dimension+2
                if (isused(j) .and. (force_twotail .or. .not. has_limits(j))) then
                    do endb =0,1
                        !Get binary chain depending on whether above or below confidence value
                        u = ConfidVal(j,(1-limfrac)/2,endb==0,&
                            chain_indices(ix),chain_indices(ix+1)-1)
                        do !thin_fac
                            call ThinData(thin_fac(ix),chain_indices(ix),chain_indices(ix+1)-1)
                            if (thin_rows < 2) exit
                            allocate(binchain(0:thin_rows-1))
                            where (coldata(j,thin_ix(0:thin_rows-1)) >= u)
                                binchain = 1
                            elsewhere
                                binchain = 2
                            endwhere
                            tran = 0
                            !Estimate transitions probabilities for 2nd order process
                            do i = 2, thin_rows-1
                                tran(binchain(i-2),binchain(i-1),binchain(i)) = &
                                    tran(binchain(i-2),binchain(i-1),binchain(i)) +1
                            end do
                            deallocate(binchain)

                            !Test whether 2nd order is better than Markov using BIC statistic
                            g2 = 0
                            do i1=1,2
                                do i2=1,2
                                    do i3=1,2
                                        if (tran(i1,i2,i3)/=0) then
                                            fitted = dble(tran(i1,i2,1) + tran(i1,i2,2)) * &
                                                (tran(1,i2,i3) + tran(2,i2,i3))  / dble( tran(1,i2,1) + &
                                                tran(1,i2,2) + tran(2,i2,1) + tran(2,i2,2) )
                                            focus = dble( tran(i1,i2,i3) )
                                            g2 = g2 + log( focus / fitted ) * focus
                                        end if
                                    end do !i1
                                end do !i2
                            end do !i3
                            g2 = g2 * 2
                            if (g2 - log( dble(thin_rows-2) ) * 2 < 0) exit
                            thin_fac(ix) = thin_fac(ix) + 1
                        end do !thin_fac

                        !Get Markov transition probabilities for binary processes
                        if (sum(tran(:,1,2))==0 .or. sum(tran(:,2,1))==0) then
                            thin_fac(ix) = 0
                            goto 203
                        end if
                        alpha = sum(tran(:,1,2))/dble(sum(tran(:,1,1))+sum(tran(:,1,2)))
                        beta =  sum(tran(:,2,1))/dble(sum(tran(:,2,1))+sum(tran(:,2,2)))
                        probsum = alpha + beta
                        tmp1 = log(probsum * epsilon / max(alpha,beta))/ log( dabs(1.0d0 - probsum) )
                        if (int( tmp1 + 1 ) * thin_fac(ix) > nburn(ix)) then
                            nburn(ix) = int( tmp1 + 1 ) * thin_fac(ix)
                            hardest = j
                            hardestend = endb
                        end if
                    end do
                end if
            end do !j

            markov_thin(ix) = thin_fac(ix)

            !Get thin factor to have independent samples rather than Markov
            hardest = max(hardest,1)
            u = ConfidVal(hardest,(1-limfrac)/2,hardestend==0)
            thin_fac(ix) = thin_fac(ix) + 1
            do !thin_fac
                call ThinData(thin_fac(ix),chain_indices(ix),chain_indices(ix+1)-1)
                if (thin_rows < 2) exit
                allocate(binchain(0:thin_rows-1))
                where (coldata(hardest,thin_ix(0:thin_rows-1)) > u)
                    binchain = 1
                elsewhere
                    binchain = 2
                endwhere
                tran2 = 0
                !Estimate transitions probabilities for 2nd order process
                do i = 1, thin_rows-1
                    tran2(binchain(i-1),binchain(i)) = &
                        tran2(binchain(i-1),binchain(i)) +1
                end do
                deallocate(binchain)

                !Test whether independence is better than Markov using BIC statistic
                g2 = 0
                do i1=1,2
                    do i2=1,2
                        if (tran2(i1,i2)/=0) then
                            fitted = dble( (tran2(i1,1) + tran2(i1,2)) * (tran2(1,i2) + &
                                tran2(2,i2)) ) / dble(thin_rows -1)
                            focus = dble( tran2(i1,i2) )
                            if (fitted <= 0 .or. focus <= 0) then
                                write (*,*) 'Raftery and Lewis estimator had problems'
                                return
                            end if
                            g2 = g2 + dlog( focus / fitted ) * focus
                        end if
                    end do !i1
                end do !i2
                g2 = g2 * 2

                if (g2 - log( dble(thin_rows-1) ) < 0) exit
                thin_fac(ix) = thin_fac(ix) + 1
            end do !thin_fac

203         if (thin_rows < 2) thin_fac(ix) = 0
        end do !chains

        write (F%unit,*) ''
        write (F%unit,*) 'Raftery&Lewis statistics'
        write (F%unit,*) ''
        write (F%unit,*) 'chain  markov_thin  indep_thin    nburn'
        do ix = 1, num_chains_used
            if (thin_fac(ix)==0) then
                write (F%unit,'(1I4,"      Not enough samples")') chain_numbers(ix)
            else
                write (F%unit,'(1I4,3I12)') chain_numbers(ix), markov_thin(ix), &
                    thin_fac(ix), nburn(ix)
            end if
        end do

        if (any(thin_fac(1:num_chains_used)==0)) then
            write (*,*) 'RL: Not enough samples to estimate convergence stats'
        else
            call writeS('RL: Thin for Markov:         '//&
                Trim(IntToStr(maxval(markov_thin(1:num_chains_used)))))
            indep_thin = maxval(thin_fac(1:num_chains_used))
            call writeS('RL: Thin for indep samples:  '// &
                trim(IntToStr(indep_thin)))
            call WriteS('RL: Estimated burn in steps: '//&
                trim(IntToStr(maxval(nburn(1:num_chains_used))))//' ('//&
                trim(IntToStr(nint(maxval(nburn(1:num_chains_used))/mean_mult)))//' rows)')
        end if

        !!Get correlation lengths
        write (F%unit,*) ''
        write (F%unit,*) 'Parameter auto-correlations as function of step separation'
        write (F%unit,*) ''
        if (corr_length_thin/=0) then
            autocorr_thin = corr_length_thin
        else
            if (indep_thin ==0) then
                autocorr_thin = 20
            elseif (indep_thin <= 30) then
                autocorr_thin =  5
            else
                autocorr_thin = 5* (indep_thin/30)
            end if
        end if

        call ThinData(autocorr_thin,0,nrows-1)
        maxoff = min(corr_length_steps,thin_rows/(autocorr_thin*num_chains_used))
        allocate(Corrs(ncols,maxoff))
        corrs = 0
        do off =1,maxoff
            do i=off, thin_rows-1
                do j = 3, ncols
                    if (isused(j)) &
                        corrs(j,off) = corrs(j,off) + (coldata(j,thin_ix(i))-fullmean(j))* &
                        (coldata(j,thin_ix(i-off)) - fullmean(j))
                end do
            end do
            do j = 3, ncols
                if (isused(j)) corrs(j,off) = corrs(j,off)/(thin_rows-off)/fullvar(j)
            end do
        end do

        if (maxoff>0) then
            write (F%unit,'("   ",'//trim(IntToStr(maxoff)) // 'I8)')  &
                (/(I, I=autocorr_thin, maxoff*autocorr_thin,autocorr_thin)/)
            do j = 3, ncols
                if (isused(j)) then
                    write (F%unit,'(1I3,'//trim(IntToStr(maxoff)) // 'f8.3,"  ' &
                        //trim(labels(j))//'")') j-2, corrs(j,:)
                end if
            end do
        end if

        deallocate(Corrs)
    end if

    call F%Close()

    end subroutine DoConvergeTests

    subroutine Get1DDensity(j)
    integer j,i,ix2
    real(mcp), allocatable :: binlikes(:), bincounts(:), binsraw(:), prior_mask(:)
    real(mcp), allocatable :: finebins(:),finebinlikes(:)
    integer ix
    integer imin, imax, winw, end_edge, fine_edge
    real(mcp) width,fine_width, edge_fac, maxbin
    logical :: has_prior
    character(LEN=:), allocatable :: fname, filename
    integer, parameter :: fine_fac = 10
    real(mcp) :: smooth_1D, opt_width
    Type(TTextFile) :: F
    Type(TKernel1D) :: Kernel
    real(mcp) :: normed, corrected, xP, binnorm

    ix = colix(j)

    param_min(j) = minval(coldata(ix,0:nrows-1))
    param_max(j) = maxval(coldata(ix,0:nrows-1))
    !Want sensible range, but also need to check that wide enough for 2D plots where projection may look much more peaked
    range_min(j) = min(ND_limit_bot(2,j), ConfidVal(ix,0.001_mcp,.false.))
    range_max(j) = max(ND_limit_top(2,j), ConfidVal(ix,0.001_mcp,.true.))

    width = (range_max(j)-range_min(j))/(num_bins+1)
    if (width==0) return

    if (smooth_scale_1D<=0._mcp ) then
        !Automatically set smoothing scale from rule of thumb for Gaussian, e.g. see
        !http://en.wikipedia.org/wiki/Kernel_density_estimation
        !1/5 power is insensitive so just use v crude estimate of effective number
        opt_width = 1.06/max(1.d0,numsamp/max_mult)**0.2d0 *sddev(j)
        smooth_1D = opt_width/width*abs(smooth_scale_1D)
        if (smooth_1d<0.5) write(*,*) 'Warning: num_bins not large enough for optimal density'
        smooth_1D=max(1.d0, smooth_1d)
    elseif (smooth_scale_1D<1.0_mcp) then
        smooth_1D=smooth_scale_1D*sddev(j)/width
        if (smooth_1d< 1) write(*,*) 'Warning: num_bins not large enough to well sample smoothed density: '&
            //trim(NameMapping%NameOrNumber(ix-2))
    else
        smooth_1d = smooth_scale_1D
    end if
    end_edge = nint(smooth_1D*2)

    if (has_limits_bot(ix)) then
        if ( range_min(j)-limmin(ix) > width*end_edge .and. param_min(j)-limmin(ix)>width*smooth_1d) then
            !long way from limit
            has_limits_bot(ix) = .false.
        else
            range_min(j) = limmin(ix)
        end if
    end if

    if (has_limits_top(ix)) then
        if ( limmax(ix) - range_max(j) > width*end_edge .and. limmax(ix)-param_max(j)>width*smooth_1d) then
            has_limits_top(ix) = .false.
        else
            range_max(j) = limmax(ix)
        end if
    end if
    has_limits(ix)= has_limits_top(ix) .or. has_limits_bot(ix)

    if (has_limits_top(ix)) then
        center(j) = range_max(j)
    else
        center(j) = range_min(j)
    end if

    ix_min(j) = nint((range_min(j) - center(j))/width)
    ix_max(j) = nint((range_max(j) - center(j))/width)

    if (.not. has_limits_bot(ix)) ix_min(j) = ix_min(j)-end_edge
    if (.not. has_limits_top(ix)) ix_max(j) = ix_max(j)+end_edge

    allocate(binsraw(ix_min(j):ix_max(j)))
    binsraw = 0

    winw = nint(2.5*fine_fac*smooth_1D)
    fine_edge = winw+fine_fac*end_edge
    fine_width = width/fine_fac

    imin = nint((param_min(j) - center(j))/fine_width)
    imax = nint((param_max(j) - center(j))/fine_width)
    imin = min(ix_min(j)*fine_fac, imin)
    imax = max(ix_max(j)*fine_fac, imax)

    allocate(finebins(imin-fine_edge:imax+fine_edge))
    finebins=0
    if (plot_meanlikes) allocate(finebinlikes(imin-fine_edge:imax+fine_edge))
    if (plot_meanlikes) finebinlikes=0


    do i = 0, nrows-1
        ix2=nint((coldata(ix,i)-center(j))/width)
        if (ix2<=ix_max(j) .and. ix2>= ix_min(j)) binsraw(ix2) = binsraw(ix2) + coldata(1,i)
        ix2=nint((coldata(ix,i)-center(j))/fine_width)
        finebins(ix2) = finebins(ix2) + coldata(1,i)
        if (plot_meanlikes) then
            if (mean_loglikes) then
                finebinlikes(ix2) = finebinlikes(ix2) + coldata(1,i)*coldata(2,i)
            else
                finebinlikes(ix2) = finebinlikes(ix2) + coldata(1,i)*exp(meanlike - coldata(2,i))
            end if
        end if
    end do

    if (ix_min(j) /= ix_max(j)) then
        !account for underweighting near edges
        if (.not. has_limits_bot(ix) .and. binsraw(ix_min(j)+end_edge-1)==0  .and. &
            binsraw(ix_min(j)+end_edge) >  maxval(binsraw)/15) then
        call EdgeWarning(ix-2)
        end if
        if (.not. has_limits_top(ix) .and. binsraw(ix_max(j)-end_edge+1)==0  .and. &
            binsraw(ix_max(j)-end_edge) >  maxval(binsraw)/15) then
        call EdgeWarning(ix-2)
        end if
    end if
    deallocate(binsraw)

    has_prior = has_limits_bot(ix) .or. has_limits_top(ix)
    call Kernel%Init(winw, fine_fac*smooth_1D, has_prior)

    if (has_prior .and. boundary_correction_method==0) then
        allocate(prior_mask(imin-fine_edge:imax+fine_edge))
        prior_mask =1
        if (has_limits_bot(ix)) then
            prior_mask(ix_min(j)*fine_fac) = 0.5
            prior_mask(imin-fine_edge:ix_min(j)*fine_fac-1) = 0
        end if
        if (has_limits_top(ix)) then
            prior_mask(ix_max(j)*fine_fac) = 0.5
            prior_mask(ix_max(j)*fine_fac+1:imax+fine_edge) = 0
        end if
    end if

    !High resolution density (sampled many times per smoothing scale)
    if (has_limits_bot(ix)) imin=ix_min(j)*fine_fac
    if (has_limits_top(ix)) imax=ix_max(j)*fine_fac
    call Density1D%Init(imax-imin+1,fine_width)
    do i = imin,imax
        Density1D%P(i-imin+1) = sum(Kernel%Win* finebins(i-winw:i+winw))
        Density1D%X(i-imin+1) = center(j) + i*fine_width
        if (boundary_correction_method==0 .and. has_prior .and. Density1D%P(i-imin+1)>0) then
            !correct for normalization of window where it is cut by prior boundaries
            edge_fac=1/sum(Kernel%win*prior_mask(i-winw:i+winw))
            Density1D%P(i-imin+1)=Density1D%P(i-imin+1)*edge_fac
        end if
    end do
    if (boundary_correction_method==1) then
        if (has_limits_bot(ix)) then
            do i=0, winw
                normed = Density1D%P(i+1)/Kernel%a0(i)
                xP =  sum(Kernel%xWin* finebins(imin+i-winw:imin+i+winw))
                corrected = Density1D%P(i+1) * Kernel%boundary_K(i) + xP* Kernel%boundary_xK(i)
                Density1D%P(i+1) = normed * exp(corrected/normed -1)
            end do
        end if
        if (has_limits_top(ix)) then
            do i = imax-winw+1, imax
                normed = Density1D%P(i-imin+1)/Kernel%a0(imax-i)
                xP =  sum(Kernel%xWin* finebins(i-winw:i+winw))
                corrected = Density1D%P(i-imin+1) * Kernel%boundary_K(imax-i) - xP* Kernel%boundary_xK(imax-i)
                Density1D%P(i-imin+1) = normed * exp(corrected/normed -1)
            end do
        end if
    end if

    maxbin = maxval(Density1D%P)
    if (maxbin==0) then
        write (*,*) 'no samples in bin, param: '//trim(NameMapping%NameOrNumber(colix(j)-2))
        stop
    end if
    Density1D%P=Density1D%P/maxbin
    call Density1D%InitSpline()

    if (.not. no_plots) then
        !Output values for plots
        allocate(binCounts(ix_min(j):ix_max(j)))
        bincounts = density1D%P( ix_min(j)*fine_fac - imin+1:ix_max(j)*fine_fac - imin + 1: fine_fac)
        if (plot_meanlikes ) then
            allocate(binLikes(ix_min(j):ix_max(j)))
            binlikes = 0
            if (mean_loglikes) binlikes=logZero
            do ix2=ix_min(j), ix_max(j)
                binnorm = sum(Kernel%win* finebins(ix2*fine_fac-winw:ix2*fine_fac+winw))
                if (binnorm>0) then
                    binlikes(ix2)=  sum(Kernel%win* finebinlikes(ix2*fine_fac-winw:ix2*fine_fac+winw))/binnorm
                end if
            end do
            if (mean_loglikes) then
                maxbin = minval(binlikes)
                where (binlikes - maxbin < 30)
                    binlikes = exp(-(binlikes- maxbin))
                elsewhere
                    binlikes = 0
                end where
            endif
        end if

        fname = trim(dat_file_name(rootname,j))
        filename = plot_data_dir// fname
        call F%CreateFile(Filename//'.dat')
        do i = ix_min(j), ix_max(j)
            write (F%unit,float_format) center(j) + i*width, bincounts(i)
        end do
        if (ix_min(j) == ix_max(j)) write (49,'(1E16.7,'' 0'')') center(j) + ix_min*width
        call F%Close()

        if (plot_meanlikes) then
            call F%CreateFile(filename//'.likes')
            maxbin = maxval(binlikes(ix_min(j):ix_max(j)))
            do i = ix_min(j), ix_max(j)
                write (F%unit,float_format) center(j) + i*width, binlikes(i)/maxbin
            end do
            call F%Close()
        end if
    end if

    end subroutine Get1DDensity

    subroutine Get2DPlotData(j,j2)
    integer, intent(in) :: j,j2
    integer i,ix1,ix2
    real(mcp) norm, maxbin
    real(mcp) try_b, try_t,try_sum, try_last
    character(LEN=:), allocatable :: plotfile, filename
    character(LEN=256) :: numstr
    real(mcp), dimension(:,:), allocatable :: finebins, finebinlikes, Win
    integer :: fine_fac_base = 5, fine_fac
    real(mcp) widthj,widthj2, contour_levels(max_contours)
    integer imin,imax,jmin,jmax
    real(mcp) :: corr, edge_fac
    real(mcp), allocatable :: prior_mask(:,:)
    real(mcp), dimension(:,:), allocatable :: bins2D, bin2Dlikes
    real(mcp) widthx, widthy
    integer col1, col2
    integer ixmax,iymax,ixmin,iymin
    integer winw, nbin2D
    logical has_prior
    real(mcp) smooth_scale
    Type(TTextFile) :: F

    has_prior=has_limits(colix(j)) .or. has_limits(colix(j2))

    corr = corrmatrix(colix(j)-2,colix(j2)-2)
    if (abs(corr)<0.3_mcp) corr=0._mcp !keep things simple unless obvious degeneracy
    corr=max(-0.95_mcp,corr)
    corr=min(0.95_mcp,corr)
    nbin2D = min(4*num_bins_2D,nint(num_bins_2D/(1-abs(corr)))) !for tight degeneracies increase bin density
    widthx =  (range_max(j)-range_min(j))/(nbin2D+1)
    widthy =  (range_max(j2)-range_min(j2))/(nbin2D+1)
    smooth_scale = (smooth_scale_2D*nbin2D)/num_bins_2D
    fine_fac = max(2,nint(fine_fac_base/smooth_scale))

    ixmin = nint((range_min(j) - center(j))/widthx)
    ixmax = nint((range_max(j) - center(j))/widthx)

    iymin = nint((range_min(j2) - center(j2))/widthy)
    iymax = nint((range_max(j2) - center(j2))/widthy)

    if ( .not. has_limits_bot(colix(j))) ixmin = ixmin-1
    if ( .not. has_limits_bot(colix(j2))) iymin = iymin-1
    if ( .not. has_limits_top(colix(j))) ixmax = ixmax+1
    if ( .not. has_limits_top(colix(j2))) iymax = iymax+1

    allocate(bins2D(ixmin:ixmax,iymin:iymax))
    allocate(bin2Dlikes(ixmin:ixmax,iymin:iymax))
    bins2D = 0
    bin2Dlikes = 0

    winw = nint(fine_fac*smooth_scale)
    imin = (ixmin-3)*winw+1
    imax = (ixmax+3)*winw-1
    jmin = (iymin-3)*winw+1
    jmax = (iymax+3)*winw-1
    allocate(finebins(imin:imax,jmin:jmax))
    if (shade_meanlikes) allocate(finebinlikes(imin:imax,jmin:jmax))
    finebins = 0
    if (shade_meanlikes) finebinlikes=0

    widthj = widthx/fine_fac
    widthj2 = widthy/fine_fac
    col1 = colix(j)
    col2 = colix(j2)
    do i = 0, nrows-1
        ix1=nint((coldata(col1,i)-center(j))/widthj)
        ix2=nint((coldata(col2,i)-center(j2))/widthj2)
        if (ix1>=imin .and. ix1<=imax .and. ix2>=jmin .and. ix2 <=jmax) then
            finebins(ix1,ix2) = finebins(ix1,ix2) + coldata(1,i)
            if (shade_meanlikes) finebinlikes(ix1,ix2) = finebinlikes(ix1,ix2) + coldata(1,i)*exp(meanlike - coldata(2,i))
        end if
    end do

    winw = nint(2*fine_fac*smooth_scale)
    allocate(Win(-winw:winw,-winw:winw))

    do ix1=-winw,winw
        do ix2=-winw,winw
            !                Win(ix1,ix2) = exp(-(ix1**2+ix2**2)/real(fine_fac**2,mcp)/2)
            Win(ix1,ix2) = exp(-(ix1**2+ix2**2 - 2*corr*ix1*ix2)/(2*fine_fac**2*smooth_scale**2*(1-corr**2)))
        end do
    end do

    if (has_prior) then
        norm=sum(win)
        allocate(prior_mask(imin:imax,jmin:jmax))
        prior_mask =1
        if (has_limits_bot(colix(j))) then
            prior_mask(ixmin*fine_fac,:) = prior_mask(ixmin*fine_fac,:)/2
            prior_mask(imin:ixmin*fine_fac-1,:) = 0
        end if
        if (has_limits_top(colix(j))) then
            prior_mask(ixmax*fine_fac,:) = prior_mask(ixmax*fine_fac,:)/2
            prior_mask(ixmax*fine_fac+1:imax,:) = 0
        end if
        if (has_limits_bot(colix(j2))) then
            prior_mask(:,iymin*fine_fac) = prior_mask(:,iymin*fine_fac)/2
            prior_mask(:,jmin:iymin*fine_fac-1) = 0
        end if
        if (has_limits_top(colix(j2))) then
            prior_mask(:,iymax*fine_fac) = prior_mask(:,iymax*fine_fac)/2
            prior_mask(:,iymax*fine_fac+1:jmax) = 0
        end if
    end if

    do ix1=ixmin, ixmax
        do ix2=iymin,iymax
            bins2D(ix1,ix2) = sum(win* finebins(ix1*fine_fac-winw:ix1*fine_fac+winw, ix2*fine_fac-winw:ix2*fine_fac+winw))
            if (shade_meanlikes) bin2Dlikes(ix1,ix2)= &
                sum(win* finebinlikes(ix1*fine_fac-winw:ix1*fine_fac+winw,ix2*fine_fac-winw:ix2*fine_fac+winw ))

            if (has_prior) then
                !correct for normalization of window where it is cut by prior boundaries
                edge_fac=norm/sum(win*prior_mask(ix1*fine_fac-winw:ix1*fine_fac+winw, ix2*fine_fac-winw:ix2*fine_fac+winw))
                bins2D(ix1,ix2) = bins2D(ix1,ix2)*edge_fac
                if (shade_meanlikes) bin2Dlikes(ix1,ix2)=bin2Dlikes(ix1,ix2)*edge_fac
            end if
        end do
    end do

    deallocate(Win,finebins)
    if (has_prior)  deallocate(prior_mask)

    if (shade_meanlikes) then
        deallocate(finebinlikes)
        do ix1=ixmin,ixmax
            do ix2 =iymin,iymax
                if (bins2D(ix1,ix2) >0) bin2Dlikes(ix1,ix2) = bin2Dlikes(ix1,ix2)/bins2D(ix1,ix2)
            end do
        end do
    end if

    bins2D=bins2D/maxval(bins2D)
    ! Get contour containing contours(:) of the probability
    norm = sum(bins2D)

    do ix1 = 1, num_contours
        try_t = maxval(bins2D)
        try_b = 0
        try_last = -1
        do
            try_sum = sum(bins2D,mask = bins2D < (try_b + try_t)/2)
            if (try_sum > (1-contours(ix1))*norm) then
                try_t = (try_b+try_t)/2
            else
                try_b = (try_b+try_t)/2
            end if
            if (try_sum == try_last) exit
            try_last = try_sum
        end do
        contour_levels(ix1) = (try_b+try_t)/2
    end do

    where (bins2D < 1e-30)
        bins2D=0
    end where

    plotfile = dat_file_2D(rootname, j, j2)
    filename = plot_data_dir//trim(plotfile)
    call F%CreateFile(filename)
    do ix1 = ixmin, ixmax
        write (F%unit,float_format) bins2D(ix1,iymin:iymax)
    end do

    call F%CreateFile( filename //'_y')
    do i = ixmin, ixmax
        write (F%unit,float_format) center(j) + i*widthx
    end do

    call F%CreateFile( filename//'_x')
    do i = iymin, iymax
        write (F%unit,float_format) center(j2) + i*widthy
    end do

    call F%CreateFile(filename//'_cont')
    write(numstr,*) contour_levels(1:num_contours)
    if (num_contours==1) numstr = trim(numstr)//' '//trim(numstr)
    write (F%unit,*) trim(numstr)

    if (shade_meanlikes) then
        call F%CreateFile(filename //'_likes')
        maxbin = maxval(bin2Dlikes(ixmin:ixmax,iymin:iymax))
        do ix1 = ixmin, ixmax
            write (F%unit,float_format)  bin2Dlikes(ix1,iymin:iymax)/maxbin
        end do
    end if
    call F%Close()

    end subroutine Get2DPlotData

    function PlotContMATLAB(aunit,aroot,j,j2, DoShade)
    character(LEN=*), intent(in) :: aroot
    integer, intent(in) :: j, j2, aunit
    logical, intent(in) :: DoShade
    character(LEN=:), allocatable :: plotfile
    logical PlotContMATLAB

    PlotContMATLAB= .false.
    plotfile = dat_file_2D(aroot, j, j2)
    if (.not. File%Exists(plot_data_dir//trim(plotfile))) return
    PlotContMATLAB= .true.
    write (aunit,'(a)') "pts=load(fullfile(plotdir,'" // trim(plotfile) //"'));"

    write (aunit,'(a)') "tmp = load(fullfile(plotdir,'" // trim(plotfile) //'_x''));'
    write (aunit,*) 'x1 = tmp(:,1);'

    write (aunit,'(a)') "tmp = load(fullfile(plotdir,'" // trim(plotfile) //'_y''));'
    write (aunit,*) 'x2 = tmp(:,1);'
    if (DoShade) then
        if (shade_meanlikes) then
            write (aunit,'(a)') "ptsL=load(fullfile(plotdir,'" // trim(plotfile) //'_likes''));'
            write (aunit,'(a)') 'contourf(x1,x2,ptsL,64);'
        else
            write (aunit,*) 'contourf(x1,x2,pts,64);'
        end if

        write (aunit,*) 'set(gca,''climmode'',''manual''); shading flat; hold on;'
    end if

    if (num_contours /= 0) then
        write (aunit,'(a)') "cnt = load(fullfile(plotdir,'" // trim(plotfile) //'_cont''));'
    end if

    end function PlotContMATLAB

    function matlabLabel(i)
    integer, intent(in) :: i
    character(LEN=:), allocatable :: matlabLabel

    if (matlab_latex) then
        matlabLabel = '$$'//trim(labels(i))//'$$'
    else
        matlabLabel= labels(i)
    end if
    end function matlabLabel

    subroutine Write2DPlotMATLAB(aunit,j,j2, DoLabelx,DoLabely, hide_ticklabels)
    integer, intent(in) :: aunit,j,j2
    logical, intent(in) :: DoLabelx,DoLabely
    logical, intent(in), optional :: hide_ticklabels
    character(LEN=120) fmt
    integer i
    logical hide_ticks

    if (present(hide_ticklabels)) then
        hide_ticks = hide_ticklabels
    else
        hide_ticks = .false.
    endif

    if (PlotContMATLAB(aunit,rootname,j,j2,do_shading) .and. num_contours /= 0) then
        write (aunit,'(a)') '[C h] = contour(x1,x2,pts,cnt,lineM{1});'
        write (aunit,'(a)') 'set(h,''LineWidth'',lw1);'
        write (aunit,*) 'hold on; axis manual; '
        do i = 1, ComparePlots%Count
            fmt = trim(numcat('lineM{',i+1)) // '}'
            if (PlotContMATLAB(aunit,ComparePlots%Item(i),j,j2,.false.)) &
                write (aunit,'(a)') '[C h] = contour(x1,x2,pts,cnt,'//trim(fmt)//');'
            write (aunit,'(a)') 'set(h,''LineWidth'',lw2);'
        end do
    end if

    write (aunit,'(a)') 'hold off; set(gca,''Layer'',''top'',''FontSize'',axes_fontsize);'
    fmt = ''',''FontSize'',lab_fontsize);'
    if (DoLabelx) then
        write (aunit,'(a)') 'xlabel('''//trim(matlabLabel(colix(j2)))//trim(fmt)
    else
        if (hide_ticks) write(aunit,*) 'set(gca,''xticklabel'',[]);'
    end if
    if (DoLabely) then
        write (aunit,'(a)') 'ylabel('''//trim(matlabLabel(colix(j)))//trim(fmt)
    else
        if (hide_ticks) write(aunit,*) 'set(gca,''yticklabel'',[]);'
    end if
    end subroutine  Write2DPlotMATLAB


    subroutine WritePlotFileInit(unit,sm, subplot_size)
    integer, intent(in) :: unit
    logical, intent(in) :: sm
    real, intent(in) :: subplot_size
    integer sz, i

    if (plot_ext=='py') then
        write(unit,'(a)') 'import GetDistPlots, os'
        write(unit,'(a)') 'g=GetDistPlots.GetDistPlotter(plot_data='''// plot_data_dir//''')'
        write(unit,'(a)') 'g.settings.setWithSubplotSize('//trim(RealToStr(subplot_size))//')'
        write(unit,'(a)') 'outdir='''//out_dir//''''
        write(unit,'(a)', advance='NO') 'roots=['''//trim(rootname)//''''
        do i = 1, ComparePlots%Count
            write(unit,'(a)', advance='NO') ','''// ComparePlots%Item(i)//''''
        end do
        write(unit,'(a)') ']'
    else
        write(unit,*) 'plot_size_inch = '//trim(RealToStr(subplot_size))//';'
        write(unit,'(a)') 'plotdir='''//plot_data_dir//''';'
        write(unit,'(a)') 'outdir='''//out_dir//''';'
        sz = 12
        if (sm .and. subplot_size < 4) sz=9
        sz = nint(sz*font_scale)
        write(unit,*) trim(concat('lab_fontsize = ',sz,'; axes_fontsize = ',sz,';'))
        if (matlab_latex) write(unit,*) 'set(0,''DefaultTextInterpreter'',''Latex'');'
        write(unit,*) 'clf'
        if (BW) then
            if (plot_meanlikes) then
                write(unit,*) 'lineM = {''-k'',''-r'',''-b'',''-m'',''-g'',''-c'',''-y''};';
                write(unit,*) 'lineL = {'':k'','':r'','':b'','':m'','':g'','':c'','':y''};';
                write(unit,*) 'lw1=3;lw2=1;' !Line Widths
            else
                write(unit,*) 'lineM = {''-k'',''--r'',''-.b'','':m'',''--g'',''-.c'',''-y''};';
                write(unit,*) 'lw1=1;lw2=1;' !Line Widths
            end if
        else
            write(unit,*) 'lineM = {''-k'',''-r'',''-b'',''-m'',''-g'',''-c'',''-y'',''--k'',''--r'',''--b''};';
            write(unit,*) 'lineL = {'':k'','':r'','':b'','':m'','':g'','':c'','':y'',''-.k'',''-.r'',''-.b''};';
            write(unit,*) 'lw1=1;lw2=1;' !Line Widths
        end if
        write(unit,*) 'colstr=''krbmgcykrb'';'
    end if

    end subroutine WritePlotFileInit

    subroutine WritePlotFileExport(unit, tag, plot_col, plot_row)
    integer, intent(in) :: unit, plot_col, plot_row
    character(LEN=*), intent(in) :: tag
    character(LEN=10) command
    character(LEN=:), allocatable:: outname

    outname=''''//trim(rootname)//trim(tag)//'.'//trim(plot_output)//''''
    if (plot_ext=='m') then
        write(unit,*)  'set(gcf, ''PaperUnits'',''inches'');'
        write(unit,*) 'x=',plot_col,'*plot_size_inch; y=',plot_row,'*plot_size_inch;'
        write(unit,*) 'set(gcf, ''PaperPosition'',[0 0 x y]); set(gcf, ''PaperSize'',[x y]);'
        if (plot_output == 'ps')  command='-dpsc2'
        if (plot_output == 'pdf')  command='-dpdf'
        if (plot_output == 'eps') command='-depsc2'
        write (unit,'(a)') 'print('''//trim(command)//''',fullfile(outdir,'// outname //'));'
    elseif (plot_ext=='py') then
        write (unit,'(a)') 'g.export(os.path.join(outdir,'// outname //'))'
    end if

    end subroutine WritePlotFileExport

    function quoted_param_name(j) result(res)
    integer, intent(in) :: j
    character(LEN=:), allocatable :: res

    res=''''//trim(NameMapping%NameOrNumber(j))//''''

    end function quoted_param_name

    function quoted_param_name_used(j) result(res)
    character(LEN=:), allocatable :: res
    integer, intent(in) :: j

    res=quoted_param_name(colix(j)-2)

    end function quoted_param_name_used

    function python_param_array(params,num) result(res)
    integer i,j
    integer, intent(in) :: params(:), num
    character(LEN=:), allocatable :: res

    res='['
    do i=1, num
        j= params(i)
        if (i>1) res=res//','
        res=res // quoted_param_name_used(j)
    end do
    res= res //']'

    end function python_param_array


    subroutine Write1DplotMatLab(aunit,j)
    integer, intent(in) :: aunit, j
    character(LEN=:), allocatable :: fname
    character(LEN=60) fmt
    integer ix1

    fname = dat_file_name(rootname, j)
    write (aunit,'(a)') "pts=load(fullfile(plotdir,'" // trim(fname)//  '.dat''));'
    write (aunit,*) 'plot(pts(:,1),pts(:,2),lineM{1},''LineWidth'',lw1);'
    write (aunit,*) 'axis([-Inf,Inf,0,1.1]);axis manual;'
    write (aunit,*) 'set(gca,''FontSize'',axes_fontsize); hold on;'

    if (plot_meanlikes) then
        write (aunit,'(a)') "pts=load(fullfile(plotdir,'" // trim(fname)//  '.likes''));'
        write (aunit,*) 'plot(pts(:,1),pts(:,2),lineL{1},''LineWidth'',lw1);'
    end if

    if (has_markers(colix(j))) then
        write (aunit,'("line([",1e15.6," ",1e15.6,"],[0 2],''Color'',''k'');")') markers(colix(j)),markers(colix(j))
    end if

    do ix1 = 1, ComparePlots%Count
        fname = dat_file_name(ComparePlots%Item(ix1),j)
        if (File%Exists(plot_data_dir// trim(fname)//'.dat')) then
            write (aunit,'(a)') "pts = load(fullfile(plotdir,'" // trim(fname)// '.dat''));'
            fmt = trim(numcat('lineM{',ix1+1)) // '}'
            write (aunit,'(a)') 'plot(pts(:,1),pts(:,2),'//trim(fmt)//',''LineWidth'',lw2);'

            if (plot_meanlikes) then
                write (aunit,'(a)') "pts=load(fullfile(plotdir,'" //trim(fname)//'.likes''));'
                fmt = trim(numcat('lineL{',ix1+1)) // '}'
                write (aunit,'(a)') 'plot(pts(:,1),pts(:,2),'//trim(fmt)//',''LineWidth'',lw2);'
            end if
        end if
    end do

    end subroutine Write1DplotMatLab

    subroutine WriteMatlabLineLabels(F)
    class(TTextFile) :: F
    integer ix1
    call F%WriteFormat('set(gcf,''Units'',''Normal''); frac=1/%u;',ComparePlots%Count+1)
    call F%Write('ax=axes(''Units'',''Normal'',''Position'',[0.1,0.95,0.85,0.05],''Visible'',''off'');')
    call F%Write('text(0,0,'''//trim(rootname)//''',''color'', colstr(1),''Interpreter'',''none'');')
    do ix1 = 1, ComparePlots%Count
        call F%WriteFormat('text(frac*%u,0,''' // ComparePlots%Item(ix1) // &
            ''',''Interpreter'',''none'',''color'',colstr(%u));', ix1, ix1+1)
    end do

    end  subroutine WriteMatlabLineLabels

    subroutine EdgeWarning(param)
    integer, intent(in) :: param

    if (NameMapping%nnames==0 .or. NameMapping%name(param)=='') then
        call WriteS('Warning: sharp edge in parameter '//trim(intToStr(param))// &
            ' - check limits'//trim(intToStr(param)))
    else
        call WriteS('Warning: sharp edge in parameter '//trim(NameMapping%name(param))// &
            ' - check limits['//trim(NameMapping%name(param))//'] or limits'//trim(intToStr(param)))
    end if

    end subroutine EdgeWarning

    function dat_file_name(rootname,j)
    integer, intent(in) :: j
    character(LEN=*) :: rootname
    character(LEN= :), allocatable :: dat_file_name
    dat_file_name =  trim(rootname)//'_p_'//trim(NameMapping%NameOrNumber(colix(j)-2))
    end function dat_file_name

    function dat_file_2D(rootname,j, j2)
    integer, intent(in) :: j, j2
    character(LEN=*) :: rootname
    character(LEN=:), allocatable :: dat_file_2D
    dat_file_2D =  trim(rootname)//'_2D_'//trim(NameMapping%NameOrNumber(colix(j)-2)) &
        //'_'//trim(NameMapping%NameOrNumber(colix(j2)-2))
    end function dat_file_2D

    subroutine CheckMatlabAxes(afile)
    integer, intent(in) :: afile

    write (afile,*) 'ls =get(gca,''XTick'');sz=size(ls,2);'
    write (afile,*) 'if(sz>2 && abs(ls(sz)-ls(1))<0.01)'
    write (afile,*) ' set(gca,''XTick'',ls(:,1:round(sz/2):sz));'
    write (afile,*) 'elseif (sz>4)'
    write (afile,*) ' set(gca,''XTick'',ls(:,1:2:sz));'
    write (afile,*) 'end;'

    end subroutine CheckMatlabAxes

    ! normal inverse translate from
    !http://home.online.no/~pjacklam/notes/invnorm
    ! a routine written by john herrero
    function dinvnorm(p)
    real(mcp) dinvnorm, p,p_low,p_high
    real(mcp) a1,a2,a3,a4,a5,a6
    real(mcp) b1,b2,b3,b4,b5
    real(mcp) c1,c2,c3,c4,c5,c6
    real(mcp) d1,d2,d3,d4
    real(mcp) q,r
    a1=-39.6968302866538
    a2=220.946098424521
    a3=-275.928510446969
    a4=138.357751867269
    a5=-30.6647980661472
    a6=2.50662827745924
    b1=-54.4760987982241
    b2=161.585836858041
    b3=-155.698979859887
    b4=66.8013118877197
    b5=-13.2806815528857
    c1=-0.00778489400243029
    c2=-0.322396458041136
    c3=-2.40075827716184
    c4=-2.54973253934373
    c5=4.37466414146497
    c6=2.93816398269878
    d1=0.00778469570904146
    d2=0.32246712907004
    d3=2.445134137143
    d4=3.75440866190742
    p_low=0.02425
    p_high=1-p_low
    if(p < p_low) then
        q=dsqrt(-2*dlog(p))
        dinvnorm=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
    else if((p.le.p_high)) then
        q=p-0.5
        r=q*q
        dinvnorm=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
    else
        q=dsqrt(-2*dlog(1-p))
        dinvnorm=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ ((((d1*q+d2)*q+d3)*q+d4)*q+1)
    end if
    end function dinvnorm

    subroutine GetChainLikeSummary(unit)
    integer unit

    bestfit_ix = 0 !Since we have sorted the lines

    maxlike = coldata(2,bestfit_ix)
    write (unit,*) 'Best fit sample -log(Like) = ', maxlike

    if (coldata(2,nrows-1) - maxlike < 30) then
        meanlike = log(sum(exp((coldata(2,0:nrows-1) -maxlike))*coldata(1,0:nrows-1)) &
            / numsamp) + maxlike
        write (unit,*) 'Ln(mean 1/like) = ', meanlike
    end if

    meanlike = sum(coldata(2,0:nrows-1)*coldata(1,0:nrows-1)) / numsamp
    write (unit,*) 'mean(-Ln(like)) = ', meanlike

    meanlike = -log(sum(exp(-(coldata(2,0:nrows-1) -maxlike))*coldata(1,0:nrows-1)) / numsamp) + maxlike
    write (unit,*) '-Ln(mean like)  = ', meanlike

    end subroutine GetChainLikeSummary


    end module MCSamples

    program GetDist
    use MCSamples
    use IO
    use settings
    implicit none
    Type(TSettingIni) :: Ini
    character(LEN=:), allocatable ::  InputFile, parameter_names_file, parameter_names_labels
    character(LEN=:) ,allocatable :: InLine
    integer ix, ix1,i,ix2,ix3
    real(mcp) ignorerows

    logical bad, adjust_priors

    character(LEN=:), allocatable :: filename, infile, out_root
    character(LEN=120) matlab_col, InS1,InS2,fmt
    integer plot_row, plot_col, chain_num, first_chain,chain_ix, plot_num
    integer, parameter :: ReadAllChainsNum = 1024

    integer x,y,j,j2, i2, num_2D_plots, num_cust2D_plots
    integer num_3D_plots
    character(LEN=:), allocatable :: contours_str
    character(LEN=80) plot_3D(20), bin_limits
    integer thin_factor, outliers
    integer plot_2D_param, j2min
    real(mcp) try_b, try_t
    real(mcp) LowerUpperLimits(max_cols,2,max_contours), limfrac

    integer chain_exclude(max_chains), num_exclude
    logical map_params
    logical :: triangle_plot = .false.
    logical :: no_triangle_axis_labels = .false.
    real(mcp) counts
    real(mcp) cool
    integer PCA_num, PCA_normparam
    integer PCA_params(max_cols)
    integer plotparams_num, plotparams(max_cols)
    character(LEN=max_cols) PCA_func
    logical Done2D(max_cols,max_cols)
    logical line_labels

    real(mcp) thin_cool
    logical no_tests, auto_label, no_names

    logical :: single_column_chain_files, samples_are_chains
    !for single_colum_chain_files
    integer :: first_haschain, stat, ip, itmp, idx, nrows2(max_cols), tmp_params(max_cols)
    real(mcp) :: rtmp
    character(LEN=:), allocatable :: finish_run_command
    Type(TTextFile) :: ChainFile
    integer triangle_num, triangle_params(max_cols), max_scatter_points
    real(mcp) tail_limit_bot,tail_limit_top, tail_confid_bot, tail_confid_top
    logical make_scatter_samples
    real(mcp) :: converge_test_limit
    Type (TTextFile) FileMatlab, File2D, File3d, FileTri, LikeFile
    integer c

    NameMapping%nnames = 0

    InputFile = GetParam(1)
    if (InputFile == '') stop 'No parameter input file'

    call Ini%Open(InputFile,  bad, .false.)

    if (bad) stop 'Error opening parameter file'

    parameter_names_file = Ini%Read_String('parameter_names')
    if (parameter_names_file/='') call NameMapping%Init(parameter_names_file)
    parameter_names_labels = Ini%Read_String('parameter_names_labels')
    matlab_latex = Ini%read_logical('matlab_latex',.false.)

    if (Ini%HasKey('nparams')) then
        ncols = Ini%Read_Int('nparams') + 2
        if (Ini%HasKey('columnnum')) stop 'specify only one of nparams or columnnum'
    else
        ncols = Ini%Read_Int('columnnum',0)
    end if

    if (NameMapping%nnames/=0 .and. ncols==0) then
        ncols = NameMapping%nnames+2
    end if

    in_root = GetParam(2)
    if (in_root=='') in_root = Ini%Read_String('file_root', notFoundFail=.true.)
    rootname = File%ExtractName(in_root)
    chain_num = Ini%Read_Int('chain_num',-1) !-1 means keep reading until one not found
    if (chain_num==-1) chain_num=ReadAllChainsNum

    single_column_chain_files = Ini%Read_Logical( 'single_column_chain_files',.false.)
    prior_ranges = 0

    if ( single_column_chain_files ) then
        pname(1) = 'weight'
        pname(2) = 'lnlike'
        do ix=3, ncols
            pname(ix) = Ini%Read_String(concat('pname',ix-2))
            if (pname(ix)=='' .and. NameMapping%nnames/=0) then
                pname(ix) = NameMapping%name(ix-2)
            end if
        end do
    else
        if (parameter_names_file=='') then
            call IO_ReadParamNames(NameMapping,in_root,prior_ranges)
            if (ncols==0 .and. NameMapping%nnames/=0) ncols = NameMapping%nnames+2
        end if
        if (parameter_names_labels/='') call NameMapping%SetLabels(parameter_names_labels)

        if (ncols==0) then
            if (chain_num == 0) then
                infile = in_root // '.txt'
            else
                infile = in_root //'_1.txt'
            end if

            ncols = File%TxtFileColumns(infile)
            write (*,*) 'Reading ',ncols, 'columns'
        end if
    end if

    allocate(coldata(ncols,0:max_rows))

    num_bins = Ini%Read_Int('num_bins')
    num_bins_2D = Ini%Read_Int('num_bins_2D', num_bins)
    smooth_scale_1D = Ini%read_Double('smooth_scale_1D',smooth_scale_1D)
    smooth_scale_2D = Ini%read_Double('smooth_scale_2D',smooth_scale_2D) !smoothing scale in terms of bin scale
    if (smooth_scale_1D>0 .and. smooth_scale_1D>1) write(*,*) 'WARNING: smooth_scale_1D>1 is oversmoothed'
    if (smooth_scale_1D>0 .and. smooth_scale_1D>1.9) stop 'smooth_scale_1D>1 is now in stdev units'
    credible_interval_threshold =Ini%Read_Double('credible_interval_threshold',credible_interval_threshold)
    boundary_correction_method = Ini%Read_Int('boundary_correction_method',boundary_correction_method)

    ignorerows = Ini%Read_Double('ignore_rows',0.d0)

    adjust_priors = Ini%Read_Logical('adjust_priors',.false.)

    plot_ext=Ini%Read_String_Default('plot_ext','py')
    if (plot_ext=='m') matlab_version = Ini%Read_Real('matlab_version',8.)

    plot_output = Ini%Read_String_Default('plot_output',plot_output)
    call Ini%Read('subplot_size_inch', subplot_size_inch)
    subplot_size_inch2 = Ini%Read_Real('subplot_size_inch2', subplot_size_inch)
    subplot_size_inch3 = Ini%Read_Real('subplot_size_inch3', subplot_size_inch)

    font_scale  = Ini%Read_Real('font_scale',1.)
    finish_run_command = Ini%Read_String('finish_run_command')

    labels(1) = 'mult'
    labels(2) = 'likelihood'
    auto_label = Ini%Read_Logical('auto_label',.false.)
    no_names = NameMapping%nnames==0
    if (auto_label .and. no_names) then
        call NameMapping%Alloc(ncols-2)
    end if
    do ix=3, ncols
        if (auto_label) then
            labels(ix) = IntToStr(ix-2)
            if (no_names) then
                NameMapping%name(ix-2) = trim(labels(ix))
                NameMapping%label(ix-2) = trim(labels(ix))
            end if
        else
            if (NameMapping%nnames/=0 .and. .not. &
                NameMapping%HasReadIniForParam(Ini, 'lab',ix-2)) then
            labels(ix) = trim(NameMapping%label(ix-2))
            else
                labels(ix) = NameMapping%ReadIniForParam(Ini,'lab',ix-2)
            end if
        end if
    end do

    prob_label = Ini%Read_Logical('prob_label',.false.)

    samples_are_chains = Ini%Read_Logical('samples_are_chains',.true.)

    no_plots = Ini%Read_Logical('no_plots',.false.)
    plots_only = Ini%Read_Logical('plots_only',.false.)
    no_tests = plots_only .or. Ini%Read_Logical('no_tests',.false.)
    line_labels = Ini%Read_Logical('line_labels',.false.)

    thin_factor = Ini%Read_Int('thin_factor',0)
    thin_cool = Ini%read_Real('thin_cool',1.)

    first_chain = Ini%Read_Int('first_chain',1)

    make_single_samples = Ini%Read_logical('make_single_samples',.false.)
    single_thin = Ini%Read_Int('single_thin',1)
    cool = Ini%Read_Real('cool',1.)

    do ix = 1, Ini%Read_Int('compare_num',0)
        call ComparePlots%Add(File%ExtractName(Ini%Read_String(numcat('compare',ix))))
    end do

    has_limits_top = .false.
    has_limits_bot = .false.

    bin_limits = Ini%Read_String('all_limits')
    markers=0
    has_markers=.false.
    do ix=3, ncols
        if (prior_ranges(1,ix-2)/=prior_ranges(2,ix-2)) then
            !Whether we actually use these limits later depends on if there are samples near them
            limmin(ix) = prior_ranges(1,ix-2)
            limmax(ix) = prior_ranges(2,ix-2)
            has_limits_top(ix) = .true.
            has_limits_bot(ix) = .true.
        end if
        if (bin_limits /= '') then
            InLine = bin_limits
        else
            InLine = NameMapping%ReadIniForParam(Ini,'limits',ix-2)
        end if
        if (InLine /= '') then
            read (InLine,*) InS1, InS2
            if (trim(adjustl(InS1)) /= 'N') then
                has_limits_bot(ix) = .true.
                read(InS1,*) limmin(ix)
            end if
            if (trim(adjustl(InS2)) /= 'N') then
                has_limits_top(ix) = .true.
                read(InS2,*) limmax(ix)
            end if
        end if

        InLine = NameMapping%ReadIniForParam(Ini,'marker',ix-2)
        if (InLine /= '') then
            has_markers(ix) = .true.
            read(InLine,*) markers(ix)
        end if
    end do

    if (Ini%HasKey('plotparams_num')) stop 'plotparams_num deprectated; just use plot_params'
    InLine = Ini%Read_String('plot_params')
    if (InLine/='') then
        plotparams_num=-1
        call NameMapping%ReadIndices(InLine, plotparams, plotparams_num, unknown_value=-1)
    else
        plotparams_num = 0
    end if

    InLine = Ini%Read_String('plot_2D_param')
    if (InLine=='') then
        plot_2D_param = 0
    else
        call NameMapping%ReadIndices(InLine, tmp_params, 1)
        plot_2D_param = tmp_params(1)
        if (plot_2D_param/=0 .and. plotparams_num/=0 .and. &
            count(plotparams(1:plotparams_num)==plot_2D_param)==0) &
            stop 'plot_2D_param not in plotparams'
    end if

    if (plot_2D_param /= 0) then
        plot_2D_param = plot_2D_param + 2
        num_cust2D_plots = 0
    else
        !Use custom array of specific plots
        num_cust2D_plots = Ini%Read_Int('plot_2D_num',0)
        do ix = 1, num_cust2D_plots
            InLine = Ini%Read_String(numcat('plot',ix))
            call NameMapping%ReadIndices(InLine, tmp_params, 2)
            if (plotparams_num/=0 .and. (count(plotparams(1:plotparams_num)==tmp_params(1))==0 .or. &
                count(plotparams(1:plotparams_num)==tmp_params(2))==0)) then
            write(*,*) trim(numcat('plot',ix)) //': parameter not in plotparams'
            stop
            end if
            cust2DPLots(ix) = tmp_params(1)+2 + (tmp_params(2)+2)*1000
        end do
    end if

    triangle_plot = Ini%Read_Logical('triangle_plot',.false.)
    if (triangle_plot) then
        no_triangle_axis_labels = Ini%read_Logical('no_triangle_axis_labels',.false.)
        InLine = Ini%Read_String('triangle_params')
        triangle_num=-1
        if (InLine/='') then
            call NameMapping%ReadIndices(InLine, triangle_params, triangle_num, unknown_value=-1)
        end if
    end if


    InS1 =Ini%Read_String('exclude_chain')
    num_exclude = 0
    chain_exclude = 0
    read (InS1, *, end =20) chain_exclude

20  do while (chain_exclude(num_exclude+1)/=0)
        num_exclude = num_exclude + 1
    end do


    map_params = Ini%Read_Logical('map_params',.false.)
    if (map_params)  &
        write (*,*) 'WARNING: Mapping params - .covmat file is new params.'

    shade_meanlikes = Ini%Read_Logical('shade_meanlikes',.false.)

    matlab_col = Ini%Read_String('matlab_colscheme')

    out_dir = Ini%Read_String('out_dir')

    out_root = Ini%Read_String('out_root')
    if (out_root /= '') then
        rootname = out_root
        write (*,*) 'producing files with with root '//trim(out_root)
    end if

    plot_data_dir = Ini%Read_String( 'plot_data_dir' )
    if ( plot_data_dir == '' ) then
        plot_data_dir = 'plot_data/'
    end if

    plot_data_dir = File%CheckTrailingSlash(plot_data_dir)

    !    if (.not. no_plots .and. .not. DirectoryExists(plot_data_dir)) &
    !    stop 'plot_data directory does not exist'

    if (out_dir /= '') then
        out_dir = File%CheckTrailingSlash(out_dir)
        write (*,*) 'producing files in directory '//out_dir
    end if

    rootdirname = concat(out_dir,rootname)

    num_contours = Ini%Read_Int('num_contours',2)
    contours_str = ''
    do i=1, num_contours
        contours(i) = Ini%Read_Double(numcat('contour',i))
        if (i>1) contours_str = contours_str // '; '
        contours_str = contours_str // Ini%Read_String(numcat('contour',i))
        max_frac_twotail(i) = Ini%Read_Double(numcat('max_frac_twotail',i), exp(-dinvnorm((1-contours(i))/2)**2/2))
    end do
    if (.not. no_tests) then
        converge_test_limit = Ini%Read_Double('converge_test_limit',contours(num_contours))
        corr_length_thin = Ini%Read_Int('corr_length_thin',corr_length_thin)
        corr_length_steps = Ini%Read_Int('corr_length_steps',corr_length_steps)
    end if
    force_twotail = Ini%Read_Logical('force_twotail',.false.)
    if (force_twotail) write (*,*) 'Computing two tail limits'

    if (Ini%Read_String('cov_matrix_dimension')=='') then
        if (NameMapping%nnames/=0) covmat_dimension = NameMapping%num_MCMC
    else
        covmat_dimension = Ini%Read_Int('cov_matrix_dimension',0)
        if (covmat_dimension == -1) covmat_dimension = ncols-2
    end if

    plot_meanlikes = Ini%Read_Logical('plot_meanlikes',.false.)

    if (Ini%HasKey('do_minimal_1d_intervals')) &
        stop 'do_minimal_1d_intervals no longer used; set credible_interval_threshold instead'

    PCA_num = Ini%Read_Int('PCA_num',0)
    if (PCA_num /= 0) then
        if (PCA_num <2) stop 'Can only do PCA for 2 or more parameters'
        InLine = Ini%Read_String('PCA_params')
        PCA_func =Ini%Read_String('PCA_func')
        !Characters representing functional mapping
        if (PCA_func == '') PCA_func(1:PCA_num) = 'N' !no mapping
        if (InLine == 'ALL' .or. InLine == 'all') then
            PCA_params(1:PCA_num) = (/ (i, i=1,PCA_num)/)
        else
            call NameMapping%ReadIndices(InLine, PCA_params, PCA_num)
        end if
        InLIne = Ini%Read_String('PCA_normparam')
        if (InLine=='') then
            PCA_NormParam = 0
        else
            call NameMapping%ReadIndices(InLine, tmp_params, 1)
            PCA_NormParam = tmp_params(1)
        end if
    end if

    num_3D_plots = Ini%Read_Int('num_3D_plots',0)
    do ix =1, num_3D_plots
        plot_3D(ix) = Ini%Read_String(numcat('3D_plot',ix))
    end do
    make_scatter_samples = Ini%Read_Logical('make_scatter_samples',.false.)
    max_scatter_points = Ini%Read_int('max_scatter_points',2000)

    BW = Ini%Read_Logical('B&W',.false.)
    do_shading = Ini%Read_Logical('do_shading',.true.)
    call Ini%Close()

    !Read in the chains
    nrows = 0
    num_chains_used =0

    do chain_ix = first_chain, first_chain-1 + max(1,chain_num)
        if (any(chain_exclude(1:num_exclude)==chain_ix)) cycle

        num_chains_used = num_chains_used + 1
        if (num_chains_used > max_chains) stop 'Increase max_chains in GetDist'
        chain_indices(num_chains_used) = nrows
        chain_numbers(num_chains_used) = chain_ix

        if ( single_column_chain_files ) then
            !Use used for WMAP 5-year chains suppled on LAMBDA; code from Mike Nolta
            !Standard CosmoMC case below

            first_haschain=0
            do ip = 1,ncols
                infile = concat(File%CheckTrailingSlash(concat(in_root,chain_ix)), pname(ip))
                if (.not. File%Exists(infile)) then
                    write (*,'(a)') 'skipping missing ' // trim(infile)
                    coldata(ip,:) = 0
                    nrows2(ip) = -1
                else
                    write (*,'(a)') 'reading ' // trim(infile)
                    call ChainFile%Open(infile)
                    if (first_haschain==0) first_haschain=ip
                    nrows2(ip) = 0
                    idx = 0 !Jo -1
                    do while (ChainFile%ReadLine(InLine))
                        idx = idx + 1

                        if ( ignorerows >= 1 .and. idx <= nint(ignorerows) ) then
                            print *, 'ignoring row'
                            cycle
                        end if

                        if (SCAN (InLine, 'N') /=0) then
                            write (*,*) 'WARNING: skipping line with probable NaN'
                            cycle
                        end if

                        read(InLine,*,iostat=stat) itmp, rtmp
                        if ( stat /= 0 ) then
                            write (*,*) 'error reading line ', nrows2(ip) + int(ignorerows), ' - skipping rest of file'
                            exit
                        end if
                        if ( idx /= itmp ) then
                            print *, "*** out of sync", idx, itmp
                            stop
                        end if

                        coldata(ip,nrows2(ip)) = rtmp
                        nrows2(ip) = nrows2(ip) + 1
                        if (nrows2(ip) > max_rows) stop 'need to increase max_rows'
                    end do
                    call ChainFile%Close()
                end if
            end do
            if (first_haschain==0) stop 'no chain parameter files read!'
            do ip = 2,ncols
                if ( nrows2(ip)/=-1 .and. nrows2(ip) /= nrows2(first_haschain) ) then
                    print *, '*** nrows mismatch:'
                    print *, nrows2(1:ncols)
                    stop
                end if
            end do
            nrows = nrows2(first_haschain)
            print *, 'all columns match, nrows = ', nrows

        else  !Not single column chain files (usual cosmomc format)
            !This increments nrows by number read in
            if (.not. IO_ReadChainRows(in_root, chain_ix, chain_num, int(ignorerows),nrows,ncols,max_rows, &
                coldata,samples_are_chains)) then
            num_chains_used = num_chains_used - 1
            if (chain_num==ReadAllChainsNum) then
                chain_num = chain_ix-1
                exit
            end if
            cycle
            endif
        end if

        if (map_params) then
            do ip =chain_indices(num_chains_used),nrows-1
                call MapParameters(coldata(1:ncols, ip))
            end do
        end if

        if (ignorerows<1 .and. ignorerows/=0) then
            i = chain_indices(num_chains_used)
            j = nint((nrows-i-1)*ignorerows)
            do ix = i,nrows-j-1
                coldata(:,ix) = coldata(:,ix+j)
            end do
            nrows = nrows - j
        end if
    end do

    if (nrows == 0) stop 'No un-ignored rows! (check number of chains/burn in)'


    if (cool /= 1) call CoolChain(cool)
    !Adjust weights if requested
    if (adjust_priors) then
        call AdjustPriors
    end if

    call GetUsedCols !See which parameters are fixed

    mean_mult = sum(coldata(1,0:nrows-1))/nrows

    max_mult = (mean_mult*nrows)/min(nrows/2,500)
    outliers = count(coldata(1,0:nrows-1) > max_mult)
    if (outliers /=0) write (*,*) 'outlier fraction ', real(outliers)/nrows

    max_mult = maxval(coldata(1,0:nrows-1))
    numsamp = sum(coldata(1,0:nrows-1))

    if (.not. no_tests) call DoConvergeTests(converge_test_limit)
    if (adjust_priors) call DeleteZeros

    write (*,*) 'mean input multiplicity = ',mean_mult

    !Output thinned data if requested
    !Must do this with unsorted output
    if (thin_factor /= 0) then
        call ThinData(thin_factor)
        call WriteThinData(trim(rootdirname)//'_thin.txt',thin_cool)
    end if

    !Produce file of weight-1 samples if requested

    if ((num_3D_plots/=0 .and. .not. make_single_samples .or. make_scatter_samples) .and. .not. no_plots) then
        make_single_samples = .true.
        single_thin = max(1,nint(numsamp/max_mult)/max_scatter_points)
    end if

    !Only use variables whose labels are not empty (and in list of plotparams if plotparams_num /= 0)
    num_vars = 0

    if (plotparams_num/=0) then
        do j=1, plotparams_num
            ix = plotparams(j)+2
            if (ix <=ncols .and. labels(ix) /= '' .and. isused(ix)) then
                num_vars = num_vars + 1
                colix(num_vars) = ix
            end if
        end do
    else
        do ix = 3,ncols
            if (labels(ix) /= '' .and. isused(ix)) then
                if (plotparams_num == 0 .or. any(plotparams(1:plotparams_num)==ix-2)) then
                    num_vars = num_vars + 1
                    colix(num_vars) = ix
                end if
            end if
        end do
    end if
    do j = 1, num_vars
        mean(j) = sum(coldata(1,0:nrows-1)*coldata( colix(j),0:nrows-1))/numsamp
        sddev(j)  = sqrt(sum(coldata(1,0:nrows-1)*(coldata(colix(j),0:nrows-1) -mean(j))**2)/numsamp)
    end do

    if (make_single_samples) call MakeSingleSamples(single_thin)

    call IO_WriteBounds(NameMapping, plot_data_dir//trim(rootname)//'.bounds', &
        limmin,limmax,has_limits_bot,has_limits_top, colix(1:num_vars))

    !Sort data in order of likelihood of points

    call SortColData(2)

    numsamp = sum(coldata(1,0:nrows-1))

    !Get ND confidence region (index into sorted coldata)
    counts = 0
    ND_cont1=-1; ND_cont2=-1
    do j=0, nrows -1
        counts = counts + coldata(1,j)
        if (counts > numsamp*contours(1) .and. ND_cont1==-1) then
            ND_cont1 = j
            if (ND_cont2 /= -1) exit
        end if
        if (counts > numsamp*contours(2) .and. ND_cont2==-1) then
            ND_cont2 = j
            if (ND_cont1 /= -1) exit
        end if
    end do
    allocate(ND_limit_top(2,num_vars), ND_limit_bot(2,num_vars))
    do j=1,num_vars
        ND_limit_bot(1,j) = minval(coldata(colix(j),0:ND_cont1))
        ND_limit_bot(2,j) = minval(coldata(colix(j),0:ND_cont2))
        ND_limit_top(1,j) = maxval(coldata(colix(j),0:ND_cont1))
        ND_limit_top(2,j) = maxval(coldata(colix(j),0:ND_cont2))
    end do

    triangle_plot = triangle_plot .and. (num_vars > 1)
    if (triangle_plot) then
        if(triangle_num==-1) then
            triangle_num=num_vars
            triangle_params(1:triangle_num) = [ (c, c= 1, num_vars) ]
        else
            ix=triangle_num
            do j=ix,1,-1
                ix2=indexOf(triangle_params(j)+2,colix,num_vars)
                if (ix2==0) then
                    triangle_num=triangle_num-1
                    triangle_params(j:triangle_num) = triangle_params(j+1:triangle_num+1)
                else
                    triangle_params(j)= ix2
                end if
            end do
            triangle_plot = triangle_num > 1
        end if
    end if

    write (*,*) 'using ',nrows,' rows, processing ',num_vars,' parameters'
    if (indep_thin/=0) then
        write (*,*) 'Approx indep samples: ', nint(numsamp/indep_thin)
    else
        write (*,*) 'effective number of samples (assuming indep): ', nint(numsamp/max_mult)
    end if
    !Get covariance matrix and correlation matrix

    call GetCovMatrix
    if (PCA_num>0 .and. .not. plots_only) call PCA(PCA_params,PCA_num,PCA_func, PCA_NormParam)

    !Find best fit, and mean likelihood
    call GetChainLikeSummary(stdout)

    if (.not. no_plots) then
        !Output files for 1D plots
        plot_col =  nint(sqrt(num_vars/1.4))
        plot_row = (num_vars +plot_col-1)/plot_col
        call FileMatlab%CreateFile(trim(rootdirname) // '.'//trim(plot_ext))
        !MatLab file for 1D plots
        call WritePlotFileInit(FileMatlab%unit,num_vars>3, subplot_size_inch)
    end if

    LowerUpperLimits = 0

    !Do 1D bins
    do j = 1,num_vars
        ix = colix(j)

        call Get1DDensity(j)

        !Get limits, one or two tail depending on whether posterior goes to zero at the limits or not
        do ix1 = 1, num_contours
            marge_limits_bot(ix1,ix) =  has_limits_bot(ix) .and. .not. force_twotail &
                .and. Density1D%P(1) > max_frac_twotail(ix1)
            marge_limits_top(ix1,ix) =  has_limits_top(ix) .and. .not. force_twotail &
                .and. Density1D%P(Density1D%n) > max_frac_twotail(ix1)
            if (.not. marge_limits_bot(ix1,ix) .or. .not. marge_limits_top(ix1,ix)) then
                !give limit
                call Density1D%Limits(contours(ix1), tail_limit_bot, tail_limit_top, &
                    marge_limits_bot(ix1,ix), marge_limits_top(ix1,ix))
                limfrac = 1-contours(ix1)
                if (marge_limits_bot(ix1,ix)) then !fix to end of prior range
                    tail_limit_bot = range_min(j)
                elseif (marge_limits_top(ix1,ix)) then !1 tail limit
                    tail_limit_bot = ConfidVal(ix,limfrac,.false.)
                else !2 tail limit
                    tail_confid_bot = ConfidVal(ix,limfrac/2,.false.)
                end if
                if (marge_limits_top(ix1,ix)) then
                    tail_limit_top = range_max(j)
                elseif (marge_limits_bot(ix1,ix)) then
                    tail_limit_top = ConfidVal(ix,limfrac,.true.)
                else
                    tail_confid_top = ConfidVal(ix,limfrac/2,.true.)
                end if
                if (.not. marge_limits_bot(ix1,ix) .and. .not. marge_limits_top(ix1,ix)) then
                    !Two tail, check if limits are at very different density
                    if (abs(Density1D%Prob(tail_confid_top) - Density1D%Prob(tail_confid_bot)) < credible_interval_threshold) then
                        tail_limit_top=tail_confid_top
                        tail_limit_bot=tail_confid_bot
                    end if
                end if
                LowerUpperLimits(j,2,ix1) = tail_limit_top
                LowerUpperLimits(j,1,ix1) = tail_limit_bot
            else !no limit
                LowerUpperLimits(j,1,ix1) = range_min(j)
                LowerUpperLimits(j,2,ix1) = range_max(j)
            end if
        end do

        if (.not. no_plots .and. plot_ext=='m') then
            call FileMatlab%WriteFormat('subplot(%u,%u,%u);',plot_row,plot_col,j)
            call Write1DplotMatLab(FileMatlab%unit,j);
            if (prob_label)  write (FileMatlab%unit,*) 'ylabel(''Probability'')'
            write(FileMatlab%unit,*)  'xlabel('''//   trim(matlabLabel(ix))//''',''FontSize'',lab_fontsize);'
            write (FileMatlab%unit,*) 'set(gca,''ytick'',[]);hold off;'
        end if !no plots
    end do

    if (.not. no_plots) then
        if (line_labels .and. plot_ext=='m') call WriteMatlabLineLabels(FileMatlab)
        if (plot_ext=='py') write(FileMatlab%unit,'(a)') 'g.plots_1d(roots)'
        call WritePlotFileExport(FileMatlab%unit, '', plot_col, plot_row)
        call FileMatlab%close()

        if (triangle_plot) then
            call FileTri%CreateFile(trim(rootdirname) // '_tri.'//trim(plot_ext))
            call WritePlotFileInit(FileTri%unit,num_vars>4, subplot_size_inch)
            if (plot_ext=='py') then
                write(FileTri%unit,'(a)') 'g.triangle_plot(roots, '//trim(python_param_array(triangle_params,triangle_num))//')'
            elseif (plot_ext=='m') then
                do i=1, triangle_num
                    j=triangle_params(i)
                    call FileTri%WriteFormat('subplot(%u,%u,%u);',triangle_num,triangle_num,(i-1)*triangle_num+i)
                    call Write1DplotMatLab(FileTri%unit,j)
                    if (triangle_num > 2 .and. subplot_size_inch<3.5) call CheckMatlabAxes(FileTri%unit)
                    if (prob_label)  write (FileTri%unit,*) 'ylabel(''Probability'')'
                    if (j==num_vars) then
                        write(FileTri%unit,*)  'xlabel('''//   trim(matlabLabel(colix(j)))//''',''FontSize'',lab_fontsize);'
                    else if (no_triangle_axis_labels) then
                        write(FileTri%unit,*) 'set(gca,''xticklabel'',[]);'
                    end if
                    write (FileTri%unit,*) 'set(gca,''ytick'',[]);hold off;'
                end do
            end if
        end if
    end if


    !do 2D bins and produce matlab file
    if (plot_2D_param == 0 .and. num_cust2D_plots==0 .and. .not. no_plots) then
        !In this case output the most correlated variable combinations
        write (*,*) 'doing 2D plots for most correlated variables'
        try_t = 1d5

        x=0;y=0;

        num_cust2D_plots = 12
        do j = 1, num_cust2D_plots
            try_b = -1d5
            do ix1 =1 ,num_vars
                do ix2 = ix1+1,num_vars
                    if (abs(corrmatrix(colix(ix1)-2,colix(ix2)-2)) < try_t .and. &
                        abs(corrmatrix(colix(ix1)-2,colix(ix2)-2)) > try_b) then
                    try_b = abs(corrmatrix(colix(ix1)-2,colix(ix2)-2))
                    x = ix1; y = ix2;
                    end if
                end do
            end do
            if (try_b == -1d5) then
                num_cust2D_plots = j-1
                exit
            end if
            try_t = try_b

            cust2Dplots(j) = colix(x) + colix(y)*1000
        end do
    end if


    if (num_cust2D_plots == 0) then
        num_2D_plots = 0

        do j= 1, num_vars
            if (ix_min(j) /= ix_max(j)) then
                do j2 = j+1, num_vars
                    if (ix_min(j2) /= ix_max(j2)) then
                        if (plot_2D_param==0 .or. plot_2D_param == colix(j) .or. plot_2D_param == colix(j2)) &
                            num_2D_plots = num_2D_plots + 1
                    end if
                end do
            end if
        end do
    else
        num_2D_plots = num_cust2D_plots
    end if

    done2D= .false.
    if (num_2D_plots > 0 .and. .not. no_plots) then
        write (*,*) 'Producing ',num_2D_plots,' 2D plots'
        filename = trim(rootdirname)//'_2D.'//trim(plot_ext)
        call File2d%CreateFile(filename)
        call WritePlotFileInit(File2d%unit,num_2D_plots >=7,subplot_size_inch2)
        if (plot_ext=='py') write(File2d%unit,'(a)') 'pairs=[]'

        plot_col = nint(sqrt(num_2D_plots/1.4))
        plot_row = (num_2D_plots +plot_col-1)/plot_col
        plot_num = 0

        do j= 1, num_vars
            if (ix_min(j) /= ix_max(j)) then
                if (plot_2D_param/=0 .or. num_cust2D_plots /= 0) then
                    if (colix(j) == plot_2D_param) cycle
                    j2min = 1
                else
                    j2min = j+1
                end if

                do j2 = j2min, num_vars
                    if (ix_min(j2) /= ix_max(j2)) then
                        if (plot_2D_param/=0 .and. colix(j2) /= plot_2D_param) cycle
                        if (num_cust2D_plots /= 0 .and.  &
                            count(cust2Dplots(1:num_cust2D_plots) == colix(j)*1000 + colix(j2))==0) cycle

                        plot_num = plot_num + 1
                        done2D(j,j2) = .true.
                        if (.not. plots_only) call Get2DPlotData(j,j2)
                        if (plot_ext=='m') then
                            call File2d%WriteFormat('subplot(%u,%u,%u);',plot_row,plot_col,plot_num)
                            call Write2DPlotMATLAB(File2d%unit,j,j2,.true.,.true.)
                            if (plot_row*plot_col > 4 .and. subplot_size_inch<3.5) call CheckMatlabAxes(File2d%unit)
                        elseif (plot_ext=='py') then
                            call File2d%WriteFormat('pairs.append([%s,%s])',quoted_param_name_used(j),quoted_param_name_used(j2))
                        end if
                    end if
                end do
            end if
        end do
        if (plot_ext=='m') then
            if (line_labels) call WriteMatlabLineLabels(File2d)
            if (matlab_col/='') call File2D%Write(matlab_col)
        elseif (plot_ext=='py') then
            call File2D%Write('g.plots_2d(roots,param_pairs=pairs)')
        end if
        call WritePlotFileExport(File2d%unit, '_2D', plot_col, plot_row)
        call File2D%Close()
    end if

    if (triangle_plot .and. .not. no_plots) then
        !Add the off-diagonal 2D plots
        do i = 1, triangle_num
            do i2 = i+1, triangle_num
                j= triangle_params(i)
                j2=triangle_params(i2)
                if (.not. Done2D(j2,j) .and. .not. plots_only) call Get2DPlotData(j2,j)
                if (plot_ext=='m') then
                    call FileTri%WriteFormat('subplot(%u,%u,%u);',triangle_num,triangle_num,(i2-1)*triangle_num + i)
                    call Write2DPlotMATLAB(FileTri%unit,j2,j,i2==triangle_num, i==1,no_triangle_axis_Labels)
                    if (triangle_num > 2 .and. subplot_size_inch<3.5) call CheckMatlabAxes(FileTri%unit)
                end if
            end do
        end do
        if (plot_ext=='m') then
            if (no_triangle_axis_labels) then
                write(FileTri%unit,*) 'h = get(gcf,''Children'');'
                write(FileTri%unit,*) 'for i=1:length(h)'
                write(FileTri%unit,*) 'p=get(h(i),''position'');'
                write(FileTri%unit,*) 'sc=1.2;w=max(p(3)*sc,p(4)*sc);'
                write(FileTri%unit,*) 'p(1)=p(1)-(w-p(3))/2; p(2)=p(2)-(w-p(4))/2;p(3)=w;p(4)=w;'
                write(FileTri%unit,*) 'set(h(i),''position'',p);'
                write(FileTri%unit,*) 'end;'
            end if
            if (matlab_col/='') write (FileTri%unit,*) trim(matlab_col)
        end if
        call WritePlotFileExport(FileTri%unit, '_tri', triangle_num+1, triangle_num+1)
        call FileTri%Close()
    end if

    !Do 3D plots (i.e. 2D scatter plots with coloured points)

    if (num_3D_plots /=0 .and. .not. no_plots) then
        write (*,*) 'producing ',num_3D_plots, '2D colored scatter plots'
        filename = trim(rootdirname)//'_3D.'//trim(plot_ext)
        call File3D%CreateFile(filename)
        call WritePlotFileInit(File3D%unit,num_3D_plots>1,subplot_size_inch3)

        if (plot_ext=='m') then
            write (File3D%unit,*) 'clf;colormap(''jet'');'
            if (mod(num_3D_plots,2)==0 .and. num_3D_plots < 11) then
                plot_col = num_3D_plots/2
            else
                plot_col =  nint(sqrt(1.*num_3D_plots))
            end if
            plot_row = (num_3D_plots +plot_col-1)/plot_col
            write(File3D%unit,'(a)') 'pts=load(fullfile(plotdir,''' // trim(rootname)//'_single.txt''));'
        elseif (plot_ext=='py') then
            write(File3D%unit,'(a)') 'sets=[]'
        end if
        do j=1, num_3D_plots
            call NameMapping%ReadIndices(plot_3D(j), tmp_params, 3)
            !x, y, color
            ix1 = indexOf(tmp_params(1)+2,colix,num_vars)+2
            ix2 = indexOf(tmp_params(2)+2,colix,num_vars)+2
            ix3 = indexOf(tmp_params(3)+2,colix,num_vars)+2
            if (any([ix1,ix2,ix3]<=0)) then
                write (*,*) 'Warning: 3D plot parameters not in used parameters'
                continue
            end if
            !            if (ix3<1) ix3 = MostCorrelated2D(ix1,ix2,ix3)
            if (plot_ext=='m') then
                call File3D%WriteFormat('subplot(%u,%u,%u);', plot_row,plot_col,j)
                write (File3D%unit,*) '%Do params ',tmp_params(1:3)
                call File3D%WriteFormat('scatter(pts(:,%u),pts(:,%u),3,pts(:,%u));', ix1,ix2,ix3)
                fmt = ''',''FontSize'',lab_fontsize);'
                write (File3D%unit,*) 'xlabel('''//trim(matlabLabel(tmp_params(1)+2))//trim(fmt)
                write (File3D%unit,*) 'ylabel('''//trim(matlabLabel(tmp_params(2)+2))//trim(fmt)
                write (File3D%unit,*) 'set(gca,''FontSize'',axes_fontsize); ax = gca;'
                write (File3D%unit,*) 'hbar = colorbar(''horiz'');axes(hbar);'

                write (File3D%unit,*) 'xlabel('''//trim(matlabLabel(tmp_params(3)+2))//trim(fmt)
                write (File3D%unit,*) 'set(gca,''FontSize'',axes_fontsize);'
                if (num_3D_plots > 2 .and. matlab_version < 7) then
                    write (File3D%unit,*) ' p = get(ax,''Position'');'
                    write (File3D%unit,*) 'set(ax,''Position'',[p(1) (p(2)+p(4)/8) p(3) p(4)]);'
                elseif (matlab_version==7) then
                    !workaround for colorbar/label overlap bug
                    write (File3D%unit,*) 'fix_colorbar(hbar,ax); axes(ax);'
                end if
            elseif (plot_ext=='py') then
                call File3D%WriteFormat('sets.append([%s,%s,%s])',  quoted_param_name(tmp_params(1)), &
                    & quoted_param_name(tmp_params(2)), quoted_param_name(tmp_params(3)))
            end if
        end do
        if (plot_ext=='py') then
            write(File3D%unit,'(a)') 'g.plots_3d(roots,sets)'
        end if
        call WritePlotFileExport(File3D%unit, '_3D', plot_col, plot_row)
        call File3D%close()
    end if

    !write out stats
    !Marginalized
    if (.not. plots_only) &
        call IO_OutputMargeStats(NameMapping, rootdirname, num_vars,num_contours,contours_str, &
        LowerUpperLimits, colix, mean, sddev, marge_limits_bot, marge_limits_top, labels)

    call NameMapping%WriteFile(plot_data_dir//trim(rootname)//'.paramnames', colix(1:num_vars)-2)

    !Limits from global likelihood
    if (.not. plots_only) then
        call LikeFile%CreateFile(trim(rootdirname)//'.likestats')
        call GetChainLikeSummary(LikeFile%unit)
        call LikeFile%NewLine()
        write(LikeFile%unit,'(a)') 'param  bestfit        lower1         upper1         lower2         upper2'
        do j=1, num_vars
            write(LikeFile%unit,'(1I5,5E15.7,"   '//trim(labels(colix(j)))//'")') colix(j)-2, coldata(colix(j),bestfit_ix),&
                ND_limit_bot(1,j), ND_limit_top(1,j), ND_limit_bot(2,j), ND_limit_top(2,j)
        end do
        call LikeFile%Close()
    end if

    !Comment this out if your compiler doesn't support "system"
    if (finish_run_command /='') then
        call StringReplace('%ROOTNAME%',rootname,finish_run_command)
        call StringReplace('%PLOTDIR%',plot_data_dir,finish_run_command)
        call StringReplace('%PLOTROOT%',plot_data_dir//rootname,finish_run_command)
        call system(finish_run_command)
    end if

    end program GetDist
