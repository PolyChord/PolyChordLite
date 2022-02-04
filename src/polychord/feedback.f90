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
    use utils_module, only: dp
    implicit none

    contains

    !> Called before running the program
    subroutine write_opening_statement(settings)
        use utils_module,    only: stdout_unit,title_fb,normal_fb
        use read_write_module, only: resume_file
        use settings_module, only: program_settings
        implicit none
        type(program_settings), intent(in) :: settings  ! The program settings 


        if(settings%feedback >=title_fb) then
            write(stdout_unit,'("")')
            write(stdout_unit,'("PolyChord: Next Generation Nested Sampling")')
            write(stdout_unit,'("copyright: Will Handley, Mike Hobson & Anthony Lasenby")')
            write(stdout_unit,'("  version: 1.20.0")')
            write(stdout_unit,'("  release: 1st June 2021")')
            write(stdout_unit,'("    email: wh260@mrao.cam.ac.uk")')
            write(stdout_unit,'("")')
        end if

        if(settings%feedback >=normal_fb) then
            write(stdout_unit,'("Run Settings"   )')
            write(stdout_unit,'("nlive    :",I8)')   settings%nlive
            write(stdout_unit,'("nDims    :",I8)')   settings%nDims
            write(stdout_unit,'("nDerived :",I8)')   settings%nDerived
            if(settings%do_clustering) write(stdout_unit,'("Doing Clustering")')
            if(settings%synchronous) write(stdout_unit,'("Synchronous parallelisation")')
            if(settings%equals) write(stdout_unit,'("Generating equally weighted posteriors")')
            if(settings%posteriors) write(stdout_unit,'("Generating weighted posteriors")')
            if((settings%equals.or.settings%posteriors).and.settings%cluster_posteriors.and.settings%do_clustering) write(stdout_unit,'("Clustering on posteriors")')
            if(settings%write_resume) write(stdout_unit,'("Writing a resume file to ",A)') trim(resume_file(settings,.false.))
            if(allocated(settings%sub_clustering_dimensions)) then
                if(size(settings%sub_clustering_dimensions)==1) then
                    write(stdout_unit,'("Sub clustering on ",I4," dimension")') size(settings%sub_clustering_dimensions)
                else 
                    write(stdout_unit,'("Sub clustering on ",I4," dimensions")') size(settings%sub_clustering_dimensions)
                end if
                write(*,*) settings%sub_clustering_dimensions
            end if

            write(stdout_unit,'("")')
        end if




    end subroutine write_opening_statement


    subroutine write_resuming(feedback)
        use utils_module,    only: stdout_unit,normal_fb
        implicit none
        integer, intent(in) :: feedback

        if(feedback>=normal_fb) then
            write(stdout_unit,'("Resuming from previous run")')
        end if

    end subroutine write_resuming


    !> Called before generating the live points
    subroutine write_started_generating(feedback)
        use utils_module,    only: stdout_unit,normal_fb
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 


        if (feedback>=normal_fb) then
            write(stdout_unit,'("generating live points")')
            write(stdout_unit,'("")')
        end if

    end subroutine write_started_generating

    !> Called during generation of the live points
    !!
    !! Prints a neat progress bar to show you how far its got
    subroutine write_generating_live_points(feedback,i_live,nlive)
        use utils_module,    only: fancy_fb
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 

        !> The current number of live points generate
        integer, intent(in) :: i_live
        !> The number of live points to be generated
        integer, intent(in) :: nlive

        if (feedback>=fancy_fb) then
            call progress(dble(i_live)/dble(nlive),100)
        end if

    end subroutine write_generating_live_points

    !> Called at the end of generating the live points
    subroutine write_finished_generating(feedback)
        use utils_module,    only: stdout_unit,normal_fb
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 


        if (feedback>=normal_fb) then
            write(stdout_unit,'("")')
            write(stdout_unit,'("all live points generated")')
            write(stdout_unit,'("")')
        end if

    end subroutine write_finished_generating


    !> Subroutine to produce a progress bar of (optional) length brz, with frac
    !! filled
    !!
    !! Note that you shouldn't put any write statements in between using this
    !! subroutine.
    subroutine progress(frac,bar_size)
        use utils_module,    only: stdout_unit
        implicit none

        !> fraction completed
        real(dp),intent(in) :: frac
        !> The size of the progress bar
        integer,intent(in) :: bar_size


        integer :: percent                    ! the percentage completed
        character(bar_size+7) :: bar          ! the bounds on the progress bar
        integer :: i                          ! loop variable

        character(len=40) :: bar_fmt  ! format string for outputting


        write(bar_fmt,'("(a1,a",i0,")")') bar_size+7

        ! Convert the fraction to an integer percentage
        percent = nint(100 * frac)

        ! Create the progress bar structure:
        ! e.g.  bar = "???% |                     |"
        ! The question marks indicate where the percentage completed will be
        ! written to.
        !
        ! Start:
        bar="???% |"
        do i= 1, bar_size-1
            ! Put in bar_size number of spaces
            bar(6+i:6+i) = ' '
        end do
        ! Finish
        bar(bar_size+7:bar_size+7) = '|'

        ! Write the percentage to the bar variable, so it now reads
        ! e.g.  bar = " 70% |                     |"
        write(unit=bar(1:3),fmt="(i3)") percent

        ! Produce a number of *'s to indicade the progress, so it now reads
        ! e.g.  bar = " 70% |***************      |"
        do i=1, int(frac*bar_size)
            bar(6+i:6+i)="*"
        enddo

        ! print the progress bar to stdout (unit 6)
        write(stdout_unit,fmt=bar_fmt,advance="no") char(13), bar
        flush(stdout_unit)

        return
    end subroutine progress


    subroutine write_num_repeats(num_repeats,feedback)
        use utils_module,    only: stdout_unit,fmt_len,INT_FMT,normal_fb
        implicit none
        integer, dimension(:) :: num_repeats
        integer, intent(in) :: feedback 

        character(len=fmt_len) fmt_int

        if(feedback>=normal_fb) then
            write(fmt_int,'("( ""number of repeats: "",",I0,A,")")') size(num_repeats),INT_FMT
            write(stdout_unit,fmt_int) num_repeats
        end if

    end subroutine write_num_repeats


    !> Called before starting sampling
    subroutine write_started_sampling(feedback)
        use utils_module,    only: stdout_unit,normal_fb
        implicit none
        !> The degree of feedback required
        integer, intent(in) :: feedback 


        if (feedback>=normal_fb) then
            write(stdout_unit,'("started sampling")')
            write(stdout_unit,'("")')
        end if

    end subroutine write_started_sampling

    !> Intermediate results
    subroutine write_intermediate_results(settings,RTI,nlikesum)
        use run_time_module, only: run_time_info,calculate_logZ_estimate
        use settings_module, only: program_settings
        use utils_module,    only: stdout_unit,fmt_len,normal_fb,sort_doubles
        implicit none
        type(program_settings), intent(in)    :: settings    !> program settings
        type(run_time_info),    intent(in)    :: RTI         !> run time info
        integer,dimension(:),   intent(in)    :: nlikesum    !> number of likelihood calls since last call

        integer :: p

        real(dp)                          :: logZ       
        real(dp)                          :: varlogZ  
        real(dp), dimension(RTI%ncluster) :: logZp      
        real(dp), dimension(RTI%ncluster) :: varlogZp 
        real(dp), dimension(RTI%ncluster_dead) :: logZp_dead      
        real(dp), dimension(RTI%ncluster_dead) :: varlogZp_dead 

        integer, dimension(RTI%ncluster+RTI%ncluster_dead) :: ordering

        character(len=fmt_len) fmt_head,fmt_live,fmt_phantom,fmt_posterior,fmt_equals,fmt_tail,fmt_nlike,fmt_Z

        integer :: int_width, i_char

        int_width = ceiling(log10(dble(maxval([RTI%nlive(:RTI%ncluster),RTI%nphantom(:RTI%ncluster),RTI%nposterior(:RTI%ncluster),RTI%nequals(:RTI%ncluster)]))))

        ! Get the number of active clusters
        write(fmt_head,'("(",I0,"(""_""))")') (int_width+2)*RTI%ncluster + 11
        write(fmt_live,'("(""lives      |"",",I0,"(I",I0,","" |""))")') RTI%ncluster,int_width
        write(fmt_phantom,'("(""phantoms   |"",",I0,"(I",I0,","" |""))")') RTI%ncluster,int_width
        write(fmt_posterior,'("(""posteriors |"",",I0,"(I",I0,","" |""))")') RTI%ncluster,int_width
        write(fmt_equals,'("(""equals     |"",",I0,"(I",I0,","" |""))")') RTI%ncluster,int_width
        write(fmt_tail,'("(",I0,"(""‾""))")') (int_width+2)*RTI%ncluster + 11


        if (settings%feedback>=normal_fb) then
            write(stdout_unit,fmt_head)
            write(stdout_unit,fmt_live)  RTI%nlive
            write(stdout_unit,fmt_phantom)  RTI%nphantom
            write(stdout_unit,fmt_posterior)  RTI%nposterior
            write(stdout_unit,fmt_equals)  RTI%nequals
            write(stdout_unit,fmt_tail)
            write(stdout_unit,'("ncluster   =",  I8," /",I8            )') RTI%ncluster, RTI%ncluster+RTI%ncluster_dead
            write(stdout_unit,'("ndead      =",  I20                   )') RTI%ndead
            write(stdout_unit,'("nposterior =",  I20                   )') RTI%nposterior_global
            write(stdout_unit,'("nequals    =",  I20                   )') RTI%nequals_global

            write(fmt_nlike,'("(""nlike      ="",",I0,"I20)")') size(nlikesum)
            write(stdout_unit,fmt_nlike) RTI%nlike

            write(fmt_nlike,'(  "(""<nlike>    ="","  ,I0,   "F15.2,""   (""",I0,"F15.2 "" per slice )"")")') size(nlikesum), size(nlikesum)
            write(stdout_unit,fmt_nlike) dble(nlikesum)/dble(settings%nlive),dble(nlikesum)/dble(RTI%num_repeats*settings%nlive)


            call calculate_logZ_estimate(RTI,logZ,varlogZ,logZp,varlogZp,logZp_dead,varlogZp_dead)            

            if(abs(logZ) < 1e9) then
                write(stdout_unit,'("log(Z)     = ", F15.2, " +/- ", F5.2)') logZ,sqrt(abs(varlogZ))
            else
                write(stdout_unit,'("log(Z)     = ?")')
            end if

            ordering = sort_doubles([-RTI%logZp,-RTI%logZp_dead])

            do p=1,RTI%ncluster+RTI%ncluster_dead

                i_char = ceiling(log10(float(p+1)))
                if(ordering(p)<=RTI%ncluster) then

                    if(abs(logZp(ordering(p)))<1e9) then
                        write(fmt_Z,'("(""log(Z_"",I",I1,","")"",A",I1,",""= "",F15.2,"" +/- "",F5.2"" (still evaluating)"")")') i_char,4-i_char
                        write(stdout_unit,fmt_Z) p, '', logZp(ordering(p)),sqrt(abs(varlogZp(ordering(p))))
                    else
                        write(fmt_Z,'("(""log(Z_"",I",I1,","")"",A",I1,",""= ? (still evaluating)"")")') i_char,4-i_char
                        write(stdout_unit,fmt_Z) p, ''
                    end if

                else
                    if(abs(logZp_dead(ordering(p)-RTI%ncluster))<1e9) then
                        write(fmt_Z,'("(""log(Z_"",I",I1,","")"",A",I1,",""= "",F15.2,"" +/- "",F5.2)")') i_char,4-i_char
                        write(stdout_unit,fmt_Z) p, '', logZp_dead(ordering(p)-RTI%ncluster),sqrt(abs(varlogZp_dead(ordering(p)-RTI%ncluster)))
                    else
                        write(fmt_Z,'("(""log(Z_"",I",I1,","")"",A",I1,",""= ?"")")') i_char,4-i_char
                        write(stdout_unit,fmt_Z) p, ''
                    end if
                end if

            end do

            write(stdout_unit,'("")')
            write(stdout_unit,'("")')
            write(stdout_unit,'("")')
        end if

    end subroutine write_intermediate_results



    !> Nicely formatted final output statement
    subroutine write_final_results(output_info,feedback)
        use utils_module,    only: stdout_unit,title_fb
        implicit none
        !> Output of the program.
        !! # log(evidence)
        !! # error(log(evidence))
        !! # ndead
        !! # number of likelihood calls
        real(dp), dimension(4) :: output_info
        !> The degree of feedback required
        integer, intent(in) :: feedback 

        if (feedback>=title_fb) then
            write(stdout_unit,'(A54)')                                        ' ____________________________________________________ '
            write(stdout_unit,'(A54)')                                        '|                                                    |'
            write(stdout_unit,'("| ndead  = ", I12, "                              |"  )') nint(output_info(3))
            write(stdout_unit,'("| log(Z) = ", F18.5, " +/- ", F18.5,  " |")') output_info(1),sqrt(abs(output_info(2)))
            write(stdout_unit,'(A54)')                                        '|____________________________________________________|'
        endif
    end subroutine write_final_results

end module feedback_module
