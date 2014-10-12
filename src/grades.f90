module grades_module
    use mpi
    implicit none

    type :: parameter_grades
        integer, allocatable, dimension(:) :: grade_index
        integer :: min_grade
        integer :: max_grade

        integer, allocatable, dimension(:) :: num_repeats
        integer, allocatable, dimension(:) :: grade_nDims



    end type parameter_grades


    contains

    subroutine allocate_grades(grades, grade_information)
        use utils_module, only: stdout_unit
        implicit none
        type(parameter_grades),intent(out) :: grades
        integer,intent(in),optional,dimension(:) :: grade_information
        integer :: errorcode,mpierror

        integer :: i


        if(present(grade_information)) then
            ! Minimum grade
            grades%min_grade=minval(grade_information)
            ! Maximum grade
            grades%max_grade=maxval(grade_information)

            allocate(&
                grades%grade_index(grades%min_grade:grades%max_grade),&
                grades%num_repeats(grades%min_grade:grades%max_grade),& 
                grades%grade_nDims(grades%min_grade:grades%max_grade)& 
                )

            do i=2,size(grade_information)
                if( grade_information(i)<grade_information(i-1) ) then
                    write(stdout_unit,'(" Parameters must be ordered in terms of grade, lowest to highest")')
                    call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
                end if
            end do
            do i=grades%min_grade,grades%max_grade
                grades%grade_index(i:i)=minloc(grade_information,mask=grade_information==i)
            end do
            do i=grades%min_grade,grades%max_grade
                grades%grade_nDims(i)=count(grade_information==i)
            end do


        else
            ! A null array for the grade information
            grades%min_grade=1
            grades%max_grade=1

        end if

        ! Set some default chain lengths to start off.
        grades%num_repeats=0
        where(grades%grade_nDims/=0) grades%num_repeats=1

    end subroutine allocate_grades

    function calc_num_babies(grades) result(num_babies)
        implicit none
        type(parameter_grades),intent(in) :: grades

        double precision :: num_babies

        num_babies = sum(grades%num_repeats*grades%grade_nDims)

    end function






end module grades_module
