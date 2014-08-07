!> This module encodes the type 'program_settings' which contains all of the
!! details required to perform a nested sampling run.
module settings_module
    use model_module, only: model
    implicit none

    !> Type to contain all of the parameters involved in a nested sampling run
    Type :: program_settings

        !> The number of live points
        integer :: nlive =512   

        !> The degree of feedback to provide
        integer :: feedback = 1

        !> The degree of precision in the final answer
        double precision :: precision_criterion = 1d-5

        !> The maximum number of dead points/samples
        !!
        !! Set equal to -1 for no maximum number
        integer :: max_ndead = -1

        !> Whether or not to save the dead points
        !! 
        !! It may not be worth saving them when doing extremely high dimensional problems
        logical :: save_dead = .true.

        !> Pointer to the sampling procedure.
        !!
        !! e.g: MultiNest, Galilean Sampling, Hamiltonian sampling ...
        !! 
        !! If you wish to write a new sampling procedure, the best way is to
        !! create a new module 'my_sampling_procedure_module' and program it in
        !! there with the interface specified in samp. Once this is done you can
        !! point to in it the main program by writing
        !!
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !! settings%sampler => my_sampling_procedure
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !!
        !! where settings is of type program_settings: 
        !!
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !! type(program settings) :: settings
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !!
        !! and should be passed to the NestedSampling algorithm 
        !!
        procedure(samp), pass(settings),       pointer :: sampler


        !> Pointer to the evidence calculator
        !!
        !! e.g: KeetonEvidence, SkillingEvidence
        !!
        !! @todo write a skilling evidence algorithm
        !!
        !! If you wish to write a evidence calculator, the best way is to
        !! write a new subroutine 'my_evidence_calculator' with the interface
        !! specified in ev. Once this is done you can point to in it the main
        !! program by writing
        !!
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !! settings%evidence_calculator => my_evidence_calculator
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !!
        !! where settings is of type program_settings: 
        !!
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !! type(program settings) :: settings
        !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !!
        !! and should be passed to the NestedSampling algorithm 
        !!
        procedure(ev),   pass(settings), pointer :: evidence_calculator

    end type program_settings

    interface
        !> Interface to the sampling procedure
        !!
        !! The sampling procedure takes the details of the model (M) and the current
        !! set of live points (live_data) in order to generate a new_point
        !! uniformly sampled from within the loglikelihood contour specifed by
        !! loglikelihood_bound
        function samp(settings, live_data, loglikelihood_bound, M,feedback) result(new_point)

            import :: model
            import :: program_settings
            implicit none

            ! ------- Inputs -------
            !> program settings (mostly useful to pass on the number of live points)
            class(program_settings), intent(in) :: settings

            !> The current set of live points. 2D array:
            !!
            !! First index ranges over ( hypercube coords, physical coords, derived params, loglikelihood),
            !!
            !! Second index ranges over all the live points
            double precision, intent(in), dimension(:,:) :: live_data

            !> The current loglikelihood bound
            double precision, intent(in) :: loglikelihood_bound

            !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
            type(model),            intent(in) :: M

            !> Optional argument to cause the sampler to print out relevent information
            integer, intent(in), optional :: feedback

            ! ------- Outputs -------
            !> The newly generated point
            double precision,    dimension(M%nTotal)   :: new_point
        end function samp
    end interface

    interface
        !> Interface to an evidence calculator
        !!
        !! The evidence calculator recieves the relevent information, namely 
        !! * the loglikelihood of the newly created point ( new_loglikelihood )
        !! * the loglikelihood of the dying point  ( old_loglikelihood )
        !! * number of iterations/dead points ( ndead )
        !!
        !! It ouputs 
        !! * a length 2 vector ( evidence_vec ) with the [evidence, evidence error] in the value of the function
        !! * whether more samples are needed in the logical variable more_samples_needed
        !!
        function ev(settings,new_loglikelihood,old_loglikelihood,ndead,more_samples_needed) result (evidence_vec)

            import :: program_settings

            implicit none

            ! ------- Inputs ------- 
            !> program settings (mostly useful to pass on the number of live points)
            class(program_settings), intent(in) :: settings
            
            !> loglikelihood of the newest created point
            double precision,       intent(in) :: new_loglikelihood

            !> loglikelihood of the most recently dead point
            double precision,       intent(in) :: old_loglikelihood

            !> number of dead points/ number of iterations
            integer,                intent(in) :: ndead

            ! ------- Outputs ------- 
            !> Whether we have obtained enough samples for an accurate evidence
            logical,intent(out) :: more_samples_needed

            ! vector containing [evidence, evidence error]
            double precision, dimension(2) :: evidence_vec

        end function ev
    end interface

end module settings_module 
