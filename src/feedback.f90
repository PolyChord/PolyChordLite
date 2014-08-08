!> Module containing all of the statements written to the command line during
!! compilation
!! 
!! There are several levels of feedback, indicated by a variable stored in
!! settings%feedback  :
!!
!! | integer  | Feedback level  
!! |----------|--------------------
!! |    0     | Only output the final evidence and error
!! |    1     | Indicate when entering various major sections of the program

module feedback_module
    implicit none

    contains

    !> Called before running the program
    subroutine write_opening_statement(M,settings)
        use model_module,    only: model,logzero
        use settings_module, only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings  ! The program settings 
        type(model),            intent(in) :: M         ! The model details

        double precision, dimension(M%nTotal) :: temp
        double precision, dimension(M%nTotal,settings%nlive) :: temp2


        if(settings%feedback >=1) then
            write(*,'("Nested Sampling Algorithm")')
            write(*,'("  author: Will Handley")')
            write(*,'("   email: wh260@cam.ac.uk")')
            write(*,*)
        end if

        if(settings%feedback >=0) then
            write(*,'("nlive      :",I8)')   settings%nlive
            write(*,'("nDims      :",I8)')   M%nDims
            write(*,'("nDerived   :",I8)')   M%nDerived
            temp    = M%loglikelihood(temp(M%p0:M%p1),settings%feedback) ! Write out the likelihood
            temp    = settings%sampler(temp2,logzero,M,settings%feedback) ! Write out the sampler
        end if
       



    end subroutine write_opening_statement



    !> Called before generating the live points
    subroutine write_started_generating(feedback)
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 


        if (feedback>=1) then
            write(*,*) 'generating live points' 
        end if

    end subroutine write_started_generating

    !> Called during generation of the live points
    !!
    !! Prints a neat progress bar to show you how far its got
    subroutine write_generating_live_points(feedback,i_live,nlive)
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 

        !> The current number of live points generate
        integer, intent(in) :: i_live
        !> The number of live points to be generated
        integer, intent(in) :: nlive

        if (feedback>=1) then
            call progress(dble(i_live)/dble(nlive))
        end if

    end subroutine write_generating_live_points

    !> Called at the end of generating the live points
    subroutine write_finished_generating(feedback)
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 


        if (feedback>=1) then
            write(*,*) 'all live points generated' 
        end if

    end subroutine write_finished_generating


    !> Subroutine to produce a progress bar of (optional) length brz, with frac
    !! filled
    !!
    !! Note that you shouldn't put any write statements in between using this
    !! subroutine.
    subroutine progress(frac,brsz)
        implicit none

        !> fraction completed
        double precision,intent(in) :: frac
        !> The size of the progress bar
        integer,intent(in),optional :: brsz

        
        integer :: percent                    ! the percentage completed
        character(:), allocatable ::bar       ! the bounds on the progress bar
        integer :: bar_size                   ! bar size
        integer :: i                          ! loop variable


        if(present(brsz))then
            ! If a specific bar size is requested set bar_size to brsz
            bar_size=brsz
        else
            ! Otherwise default to 100
            bar_size=100
        endif


        ! Convert the fraction to an integer percentage
        percent = 100 * frac

        ! Create the progress bar structure:
        ! e.g.  bar = "???% |                     |"
        ! The question marks indicate where the percentage completed will be
        ! written to.
        !
        ! Start:
        bar="???% |"
        do i= 1, bar_size
            ! Put in bar_size number of spaces
            bar = bar//' '
        end do
        ! Finish
        bar = bar//'|'

        ! Write the percentage to the bar variable, so it now reads
        ! e.g.  bar = " 70% |                     |"
        write(unit=bar(1:3),fmt="(i3)") percent

        ! Produce a number of *'s to indicade the progress, so it now reads
        ! e.g.  bar = " 70% |***************      |"
        do i=1, int(frac*bar_size)
            bar(6+i:6+i)="*"
        enddo

        ! print the progress bar to stdout (unit 6)
        write(unit=6,fmt="(a1,a<bar_size+7>)",advance="no") char(13), bar

        if (percent/=100) then
            ! If we're not at the end, then we should flush the command line so
            ! that we re-start from the beginning next time
            flush(unit=6)
        else
            ! Otherwise, we should just write a blank line 
            write(unit=6,fmt=*)
        endif
        return
    end subroutine progress

    !> Called before starting sampling
    subroutine write_started_sampling(feedback)
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 


        if (feedback>=1) then
            write(*,*) 'started sampling' 
        end if

    end subroutine write_started_sampling


    !> Nicely formatted final output statement
    subroutine write_final_results(M,evidence_vec,ndead,feedback)
        use model_module,    only: model, prior_log_volume
        implicit none
        !> The model details
        type(model),            intent(in) :: M 
        !> The degree of feedback required
        integer, intent(in) :: feedback 
        !> the evidence information
        double precision, intent(in), dimension(2) :: evidence_vec
        !> the number of dead points
        integer,intent(in) :: ndead

        if (feedback>=0) then
            write(*,'(A42)')                                        ' ________________________________________ '
            write(*,'(A42)')                                        '|                                        |'
            write(*,'("| ndead  = ", I12, "                  |"  )') ndead
            write(*,'("| Z      = ", E12.5, " +/- ", E12.5,  " |")') evidence_vec(1:2)
            write(*,'("| log(Z) = ", F12.5, " +/- ", F12.5,  " |")') log(evidence_vec(1)), evidence_vec(2)/evidence_vec(1) 
            write(*,'("| check  = ", F12.5, " +/- ", F12.5,  " |")') evidence_vec(1:2) * exp(prior_log_volume(M))
            write(*,'(A42)')                                        '|________________________________________|'
        endif

    end subroutine write_final_results

end module feedback_module
