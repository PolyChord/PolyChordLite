module loglikelihood_module

        integer :: nDims
        integer :: n_knots
        integer :: nStats

        double precision, allocatable, dimension(:)   :: theta_saved
        double precision, allocatable, dimension(:,:) :: spline_data

        double precision, allocatable, dimension(:)   :: x0
        double precision, allocatable, dimension(:)   :: y0
        double precision, allocatable, dimension(:)   :: sigmax
        double precision, allocatable, dimension(:)   :: sigmay

        double precision :: x_min_int, x_max_int
        double precision :: x_min, x_max

    contains

    function loglikelihood(theta,phi)
        use utils_module, only: logzero, logincexp,logTwoPi
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision                            :: loglikelihood ! loglikelihood value to output

        double precision :: loglikelihood_temp
        integer :: i_stats
        integer :: i_int
        double precision :: x,y


        ! Read in the spline array if necessary
        if(any(theta_saved /= theta ) ) then

            ! Calculate the spline array
            spline_data(:,1) = theta(1:n_knots)
            spline_data(:,2) = theta(n_knots+1:nDims)

            ! Save the current theta
            theta_saved = theta

        end if


        ! Calculate the likelihood
        loglikelihood_temp = logzero
        loglikelihood      = 0d0

        ! Iterate over all points in the file
        do i_stats = 1,nStats


            ! If there's no error in the x variable then this is a lot simpler
            if(sigmax(i_stats) <=0d0) then

                ! calculate the likelihood
                x = x0(i_stats)
                y = linear_interpolate(x,spline_data,0)
                loglikelihood_temp = - log( sigmay(i_stats) ) - logTwoPi/2d0 - ((y-y0(i_stats))/sigmay(i_stats))**2/2d0 

            else
                ! Integrate the spline
                loglikelihood_temp = log_exp_int(x0(i_stats),y0(i_stats),sigmax(i_stats),sigmay(i_stats),x_min,x_max)

                ! Normalise via the normalising constants, and the width of the
                ! integral
                loglikelihood_temp = loglikelihood_temp   &
                    - log( sigmay(i_stats) ) - log( sigmax(i_stats) ) - logTwoPi &
                    - log(x_max-x_min)

            end if


            ! Add the likelihood from this point to the total likelihood (note
            ! that this is in fact multiplication)
            loglikelihood = loglikelihood + loglikelihood_temp

        end do

    end function loglikelihood





    ! This function computes the logarithm of:
    !
    ! /xmax                                                
    ! |     exp( - (x-x0)^2/2sx^2 - (y0-f(x))^2/2sy^2 ) dx 
    ! /xmin                                                
    !
    ! where f(x) is a linear spline defined by spline_arr.
    !
    ! This is effectively a piecewise linear integral:
    ! 
    ! /xmax                                                
    ! |     exp( - (x-x0)^2/2sx^2 - (y0-(m x + c))^2/2sy^2 ) dx 
    ! /xmin                                                
    !
    ! where m and c change depending on which interval one is in.
    !
    ! Completing the square on the integrand yields:
    !     1   /           \2     1        1   
    ! - ----- | x - e s^2 |   - --- f  + --- e
    !   2 s^2 \           /      2        2   
    !
    ! where 
    !
    ! s = (sx^-2 + m^2 sy^-2)^-(1/2)
    ! e = x0/sx^2 + (y-c) m/ sy^2
    ! f = x0^2 /sx^2 + (y-c)^2 / sy^2
    !
    ! The integral can then be computed with error functions, since
    !
    ! /x2
    ! |      exp( (x-x0)^2/2s^2 ) dx = 
    ! /x1
    !
    ! sqrt(pi/2) s ( erf( (x2-x0)/sqrt2 s) - erf( (x1-x0)/sqrt2 s) )
    !
    function log_exp_int(x0,y0,sx,sy,xmin,xmax)
        use utils_module, only: logzero,logincexp
        implicit none

        double precision, intent(in) :: x0,y0,sx,sy,xmin,xmax
        double precision :: log_exp_int

        integer :: i,n
        double precision :: m,c,x1,x2,y1,y2

        double precision :: s,e,f

        double precision,parameter :: logsqrtpiby2 = log(sqrt(atan(1d0)*2d0))

        log_exp_int = logzero

        n = size(spline_data,1)

        do i=1,n-1

            ! Here are some useful local variables
            x1 = spline_data(i,1)
            y1 = spline_data(i,2)
            x2 = spline_data(i+1,1)
            y2 = spline_data(i+1,2)

            m = (y2-y1)/(x2-x1)
            c =  y1 - m*x1

            ! Check to see our integration range
            if(x2<xmin) cycle
            if(x1<xmin) x1 = xmin

            if(x1>xmax) cycle
            if(x2>xmax) x2 = xmax

            ! Define the reduced variables
            s  = (1/sx**2 + m**2/sy**2)**(-0.5)
            e  = x0/sx**2 + (y0-c)*m/sy**2
            f  = x0**2 /sx**2 + (y0-c)**2/sy**2


            call logincexp(log_exp_int,&
                logsqrtpiby2 + log(s) +&
                logderf(&
                (x1-e*s**2)/sqrt(2d0)/s ,&
                (x2-e*s**2)/sqrt(2d0)/s  &
                )&
                - f/2 + e**2*s**2/2 &
                )

        end do

    end function log_exp_int

    function logderf(a,b)
        use utils_module, only: logzero
        implicit none
        double precision,intent(in) :: a,b
        double precision logderf

        double precision :: erfa,erfb

        erfb = erf(b)
        erfa = erf(a)

        if(erfb<=erfa) then
            logderf = logzero
        else
            logderf = log(erfb-erfa)
        end if

    end function logderf

    function erfapprox(x)
        implicit none
        double precision, intent(in) :: x
        double precision erfapprox

        double precision, parameter :: pi = atan(1d0)*4d0
        double precision, parameter :: a = 8*(pi-3)/(3*pi*(4-pi))

        erfapprox = sign(1d0,x)*sqrt(1-exp(-x**2*(4/pi + a*x**2)/(1+a*x**2)))
    end function erfapprox


    subroutine setup_loglikelihood(settings,mpi_communicator)
        use settings_module,   only: program_settings
        use abort_module,      only: halt_program
        implicit none
        type(program_settings), intent(in) :: settings
        integer,intent(in) :: mpi_communicator

        integer, parameter :: stats_unit = 1000

        integer :: ioerror
        integer :: i_stats


        ! Save the current theta for use next call
        nDims = settings%nDims
        n_knots = nDims/2


        ! Allocate any arrays that need to be allocated
        if(allocated(theta_saved)) deallocate(theta_saved)
        if(allocated(spline_data)) deallocate(spline_data)
        if(allocated(x0)) deallocate(x0)
        if(allocated(y0)) deallocate(y0)
        if(allocated(sigmax)) deallocate(sigmax)
        if(allocated(sigmay)) deallocate(sigmay)

        allocate(theta_saved(nDims),spline_data(n_knots,2))


        ! Find out the number of points
        open(stats_unit,file='data/data.dat')

        nStats = 0
        do 
            read(stats_unit,*,iostat=ioerror)
            if( ioerror==0 ) then
                nStats = nStats+1
            else
                exit
            end if
        end do

        close(stats_unit)

        if (nStats == 0) call halt_program('Error reading data/data.dat')



        ! Allocate the data arrays
        allocate(x0(nStats),y0(nStats),sigmax(nStats),sigmay(nStats))

        ! Read in the data array
        open(stats_unit,file='data/data.dat')
        do i_stats=1,nStats
            read(stats_unit,*) x0(i_stats), y0(i_stats), sigmax(i_stats), sigmay(i_stats) 
        end do
        close(stats_unit)

        ! Read in the mins and maxs array
        open(stats_unit,file='data/data_min_max.dat')
        read(stats_unit,*) x_min
        read(stats_unit,*) x_max
        close(stats_unit)



    end subroutine setup_loglikelihood




    function linear_interpolate(x,spline_array,deriv) result(y)

        implicit none 

        double precision            :: y                  ! dependent variable
        double precision,intent(in) :: x                  ! independent variable
        double precision,intent(in) :: spline_array(:,:)  ! x coords of interpolation
                                                          ! y coords of interpolation
        integer,         intent(in) :: deriv              ! whether 0 or 1 derivatives

        integer i
        integer :: num_nodes
        double precision, dimension(size(spline_array,1)) :: x_node
        double precision, dimension(size(spline_array,1)) :: y_node

        !------------------------------------------------------------------------------
        !                                    
        !                     y_node(2)^                                           +
        !                            <-+->                         ^       +
        !       y_node(1)^      +      v   +                     <-+->
        !              <-+->                   +   ^       +       v y_node(num_nodes)
        !         +      v                       <-+->                 
        !  +                                       v y_node(3)  
        !                                                
        !   
        !----------------+-------------+-----------+---------------+-------------------
        !            x_node(1)    x_node(2)      x_node(3)      x_node(num_nodes)
        !
        ! This function constructs a linear interpolation: 
        ! http://en.wikipedia.org/wiki/Linear_interpolation
        !
        ! The interpolation points are
        ! {  ( x_node(i) , y_node(i) )  : i=1..num_nodes }
        ! 
        ! Outside of the range, the function is continued with the previous
        ! gradient
        !

        num_nodes = size(spline_array,1)
        x_node = spline_array(:,1)
        y_node = spline_array(:,2)


        ! We need a positive number of nodes
        if (num_nodes < 1) then
            write(*,*) 'linear_interpolate: cannot have less than 1 node'
            stop
        else if (num_nodes==1) then 
            ! If we have only a single node, then this is a 'flat function'
            ! defined only by the y_node
            select case(deriv)
            case(0) 
                y=y_node(1)
            case(1)
                y=0d0
            end select
            return
        endif

        ! We should also check that the x values are ordered properly
        do i=2,num_nodes
            if( x_node(i-1)>x_node(i) ) then
                write(*,*) 'linear_interpolate: nodes are not ordered correctly'
                write(*,*) 'x_nodes:', x_node
                stop
            endif
        end do


        ! We proceed from low x to high x.

        ! First consider the case when we're below the bottom x_node
        if ( x <= x_node(1) ) then
            ! we want a 'line continuation' on this side
            select case(deriv)
            case(0)
                y = y_node(1) + (y_node(2)-y_node(1)) * (x-x_node(1)) / (x_node(2)-x_node(1))
            case(1)
                y = (y_node(2)-y_node(1)) / (x_node(2)-x_node(1))
            end select
            return

        endif

        ! Now consider the cases when we're in between
        do i=2,num_nodes
            if ( x <= x_node(i) ) then
                select case(deriv)
                case(0)
                    y = y_node(i-1) + (y_node(i)-y_node(i-1)) * (x-x_node(i-1)) / (x_node(i)-x_node(i-1))
                case(1)
                    y = (y_node(i)-y_node(i-1)) / (x_node(i)-x_node(i-1))
                end select
                return
            endif
        enddo



        ! Finally consider the case when we're above the top x_node
        if ( x >= x_node(num_nodes) ) then
            ! we want a 'line continuation' on this side
            select case(deriv)
            case(0)
                y = y_node(num_nodes-1)                          &
                    + (y_node(num_nodes)-y_node(num_nodes-1))    &
                    * (x-x_node(num_nodes-1))                    &
                    / (x_node(num_nodes)-x_node(num_nodes-1))
            case(1)
                y =   (y_node(num_nodes)-y_node(num_nodes-1))    &
                    / (x_node(num_nodes)-x_node(num_nodes-1))
            end select
            return

        endif


    end function linear_interpolate



end module loglikelihood_module
