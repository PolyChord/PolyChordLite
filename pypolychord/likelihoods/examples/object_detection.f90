module loglikelihood_module

        integer :: nx,ny
        double precision :: xmin,xmax,ymin,ymax
        double precision :: sigma

        double precision, dimension(:,:), allocatable :: dat,xarr,yarr



    contains

    function loglikelihood(theta,phi)
        use utils_module, only: logzero, logincexp,logTwoPi
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision                            :: loglikelihood ! loglikelihood value to output

        loglikelihood = -sum(&
        (dat-signal(theta))**2/2/sigma**2&
        )
        
        loglikelihood = loglikelihood &
            - log(sigma**2*atan(1d0)*8)*nx*ny/2d0

    end function loglikelihood

    function signal(theta)
        implicit none
        double precision, intent(in),  dimension(:) :: theta         !> Input parameters
        double precision,dimension(nx,ny) :: signal
        integer :: i,Nobj
        double precision :: A,x,y,R

        signal = 0
        Nobj = size(theta)/4

        do i=1,Nobj
            A = theta(4*(i-1)+1)
            x = theta(4*(i-1)+2)
            y = theta(4*(i-1)+3)
            R = theta(4*(i-1)+4)
            signal = signal + A*exp( -((x-xarr)**2 + (y-yarr)**2)/2/R/R)
        end do

    end function signal




    subroutine setup_loglikelihood(settings)
        use settings_module,   only: program_settings
        use abort_module,      only: halt_program
        implicit none
        type(program_settings), intent(in) :: settings

        integer, parameter :: stats_unit = 1000

        integer :: i


        ! Allocate any arrays that need to be allocated

        ! Find out the number of points
        open(stats_unit,file='data/obj_info.dat')

        read(stats_unit,*) nx
        read(stats_unit,*) xmin
        read(stats_unit,*) xmax
        read(stats_unit,*) ny
        read(stats_unit,*) ymin
        read(stats_unit,*) ymax
        read(stats_unit,*) sigma

        close(stats_unit)

        allocate(dat(nx,ny))

        open(stats_unit,file='data/obj.dat')
        do i=1,ny
            read(stats_unit,*) dat(:,i)
        end do
        close(stats_unit)

        allocate(xarr(nx,ny),yarr(nx,ny))

        xarr = spread( &
            [(xmin + (xmax-xmin)/(nx-1)*(i-1),i=1,nx)],&
            1,ny )
        yarr = spread( &
            [(ymax + (ymin-ymax)/(ny-1)*(i-1),i=1,ny)],&
            2,nx )

    end subroutine setup_loglikelihood






end module loglikelihood_module
