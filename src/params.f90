module params_module
    use utils_module, only: STR_LENGTH
    implicit none

    type param_type
        character(len=STR_LENGTH) :: paramname   ! Name of parameter
        character(len=STR_LENGTH) :: latex       ! LaTeX name
        integer                   :: speed       ! grade of parameter
        integer                   :: prior_type  ! type of prior
        integer                   :: prior_block ! prior block
        logical                   :: sub_cluster ! sub_clustering?

        ! Parameters in the prior
        double precision, dimension(:), allocatable :: prior_params
    end type param_type

    contains

    subroutine assign_parameter(param,paramname,latex,speed,prior_type,prior_block,prior_params,sub_cluster)
        implicit none
        type(param_type),intent(out):: param       ! Parameter to be returned
        character(len=*),intent(in) :: paramname   ! Name of parameter
        character(len=*),intent(in) :: latex       ! LaTeX name
        integer         ,intent(in) :: speed       ! grade of parameter
        integer         ,intent(in) :: prior_type  ! type of prior
        integer         ,intent(in) :: prior_block ! prior block
        logical,optional,intent(in) :: sub_cluster ! sub clustering?

        ! Parameters in the prior
        double precision, dimension(:),intent(in) :: prior_params


        write(param%paramname,'(A)') paramname
        write(param%latex,'(A)') latex
        param%speed = speed
        param%prior_type = prior_type
        param%prior_block = prior_block
        allocate( param%prior_params(size(prior_params)) )
        param%prior_params = prior_params
        if(present(sub_cluster)) then
            param%sub_cluster = sub_cluster
        else 
            param%sub_cluster = .false.
        end if

    end subroutine assign_parameter

    subroutine add_parameter(params,paramname,latex,speed,prior_type,prior_block,prior_params,sub_cluster)
        implicit none
        ! Parameter array
        type(param_type),dimension(:),allocatable,intent(inout) :: params

        character(len=*),intent(in) :: paramname   ! Name of parameter
        character(len=*),intent(in) :: latex       ! LaTeX name
        integer         ,intent(in),optional :: speed       ! grade of parameter
        integer         ,intent(in),optional :: prior_type  ! type of prior
        integer         ,intent(in),optional :: prior_block ! prior block
        logical         ,intent(in),optional :: sub_cluster ! sub cluster on this?

        ! Parameters in the prior
        double precision, dimension(:),intent(in),optional :: prior_params 
        double precision, dimension(0) :: blank_params 

        ! expand parameter array
        type(param_type), dimension(:),allocatable :: temp_params

        integer :: num_params 


        num_params = size(params)

        allocate(temp_params(num_params))
        temp_params = params

        deallocate(params)
        allocate(params(num_params+1))

        params(1:num_params) = temp_params


        if(present(speed) .and. present(prior_type) .and. present(prior_block) .and. present(prior_params) ) then
            if(present(sub_cluster)) then
                call assign_parameter(params(num_params+1),paramname,latex,speed,prior_type,prior_block,prior_params,.true.) 
            else
                call assign_parameter(params(num_params+1),paramname,latex,speed,prior_type,prior_block,prior_params) 
            end if
        else
            call assign_parameter(params(num_params+1),paramname,latex,1,0,0,blank_params) 
        end if

    end subroutine add_parameter


end module params_module
