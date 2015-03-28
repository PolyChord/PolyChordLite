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
        integer,parameter :: n_int = 100
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
                ! Calculate the minimum and maximum x's so that we can integrate
                x_min_int = x0(i_stats)-5d0*sigmax(i_stats)
                x_max_int = x0(i_stats)+5d0*sigmax(i_stats)

                ! loglikelihood_temp denotes the likelihood associated with the point
                loglikelihood_temp = logzero

                ! do a simple integral by summing on equal intervals 
                do i_int = 1,n_int
                    x = x_min_int + (i_int-1d0)/(n_int-1d0) * (x_max_int-x_min_int)
                    y = linear_interpolate(x,spline_data,0)
                    call logincexp(loglikelihood_temp, - ((y-y0(i_stats))/sigmay(i_stats))**2/2d0 - ((x-x0(i_stats))/sigmax(i_stats))**2/2d0 )
                end do

                ! Normalise via the normalising constants, and the width of the
                ! integral
                ! note that the factor of (x_max_int-x_min_int) that should be here cancels out
                ! with the prior
                loglikelihood_temp = loglikelihood_temp   &
                    - log( sigmay(i_stats) ) - log( sigmax(i_stats) ) - logTwoPi &
                    - log( n_int-1d0 )

            end if


            ! Add the likelihood from this point to the total likelihood (note
            ! that this is in fact multiplication)
            loglikelihood = loglikelihood + loglikelihood_temp

        end do

    end function loglikelihood



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
