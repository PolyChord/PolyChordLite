module grades_module
    implicit none

    type :: parameter_grades
        integer, allocatable, dimension(:) :: grade_information
        integer :: min_grade
        integer :: max_grade
        integer :: num_grades

        integer, allocatable, dimension(:) :: chain_lengths



    end type parameter_grades


    contains

    subroutine allocate_grades(grades, grade_information)
        implicit none
        type(parameter_grades),intent(out) :: grades
        integer,intent(in),optional,dimension(:) :: grade_information

        integer :: i


        if(present(grade_information)) then

            ! Allocate the grade information
            allocate(grades%grade_information(size(grade_information)))
            ! Assign the grade information
            grades%grade_information=grade_information

            ! Minimum grade
            grades%min_grade=minval(grade_information)
            ! Maximum grade
            grades%max_grade=maxval(grade_information)
            ! Number of grades
            grades%num_grades = grades%max_grade-grades%min_grade+1
        else
            ! A null array for the grade information
            allocate(grades%grade_information(0))
            grades%min_grade=1
            grades%max_grade=1
            grades%num_grades=1


        end if
        
        ! Allocate an array of chain lengths (note the bounding does not
        ! start from 1 necessarily
        allocate(grades%chain_lengths(grades%min_grade:grades%max_grade))

        ! Set some default chain lengths to start off.
        ! These are chosen as the number of parameters within each grade
        do i=grades%min_grade,grades%max_grade
            grades%chain_lengths(i) = count(grade_information==i)
        end do

    end subroutine allocate_grades

    function calc_chain_length(grades) result(chain_length)
        implicit none
        type(parameter_grades),intent(in) :: grades
        integer :: i

        double precision :: chain_length

        chain_length=1
        do i=grades%max_grade,grades%min_grade,-1
            if(grades%chain_lengths(i)/=0) then
                chain_length = chain_length * grades%chain_lengths(i) +1
            end if
        end do
        chain_length=chain_length-1

    end function


    function calc_graded_choleskys(covmat,nDims,grades) result(choleskys)
        use utils_module, only: calc_cholesky
        implicit none
        integer,intent(in) :: nDims
        double precision,intent(in), dimension(nDims,nDims) :: covmat
        type(parameter_grades),intent(in) :: grades

        double precision, dimension(nDims,nDims,grades%min_grade:grades%max_grade) :: choleskys

        integer i_grade
        integer i_dims,j_dims
        integer num
        integer n_fix
        integer n_var

        integer :: errcode

        integer, dimension(nDims) :: fix_indices
        integer, dimension(nDims) :: var_indices

        double precision, dimension(nDims,nDims) :: A,B,C,Binv,subcovmat

        choleskys=0d0

        do i_grade=grades%min_grade,grades%max_grade

            num=count(grades%grade_information>=i_grade)

            ! Calculate the indices for each grade
            if(num==nDims) then

                choleskys(:,:,i_grade) = calc_cholesky(covmat,nDims)

            else if(num>0) then
                n_fix=0
                n_var=0
                fix_indices=0
                var_indices=0
                do i_dims=1,nDims
                    if(grades%grade_information(i_dims)<i_grade) then
                        n_fix=n_fix+1
                        fix_indices(n_fix) = i_dims
                    else
                        n_var=n_var+1
                        var_indices(n_var) = i_dims
                    end if
                end do

                A=0d0
                B=0d0
                Binv=0d0
                C=0d0
                A(:n_var,:n_var)=covmat(var_indices(:n_var),var_indices(:n_var))
                B(:n_fix,:n_fix)=covmat(fix_indices(:n_fix),fix_indices(:n_fix))
                C(:n_fix,:n_var)=covmat(fix_indices(:n_fix),var_indices(:n_var))

                Binv=B
                call dpotrf ('U',n_fix,Binv,nDims,errcode)
                call dpotri ('U',n_fix,Binv,nDims,errcode)
                do i_dims=1,n_fix
                    do j_dims=i_dims+1,n_fix
                        Binv(j_dims,i_dims)= Binv(i_dims,j_dims)
                    end do
                end do

                subcovmat(:n_var,:n_var) = A(:n_var,:n_var) - &
                    matmul(transpose(C(:n_fix,:n_var)),matmul(Binv(:n_fix,:n_fix),C(:n_fix,:n_var)))

                ! form cholesky factor
                call dpotrf ('U',n_var,subcovmat,nDims,errcode)

                ! Hand it back to the cholesky array
                choleskys(var_indices(:n_var),var_indices(:n_var),i_grade) = subcovmat
                
                


            end if

        end do



    end function calc_graded_choleskys




end module grades_module
