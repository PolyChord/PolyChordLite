!> This is the main driving routine of the nested sampling algorithm
program PolyChord

    ! ~~~~~~~ Loaded Modules ~~~~~~~
    use interfaces_module,        only: run_polychord
    use loglikelihood_module,     only: loglikelihood, setup_loglikelihood
    use utils_module,             only: STR_LENGTH
    use abort_module,             only: halt_program

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none
    character(len=STR_LENGTH)                 :: input_file     ! input file

    if(iargc()==1) then
        call getarg(1,input_file) 
    else
        call halt_program('PolyChord should be called with at most one argument, the input file')
    end if

    call run_polychord(loglikelihood, setup_loglikelihood, input_file)

end program PolyChord
