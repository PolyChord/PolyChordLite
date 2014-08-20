!> This module encodes the type 'program_settings' which contains all of the
!! details required to perform a nested sampling run.
module settings_module
    use model_module,   only: model
    use utils_module,   only: STR_LENGTH
    implicit none

    !> Type to contain all of the parameters involved in a nested sampling run
    Type :: program_settings

        !> The file root for outputs
        character(STR_LENGTH) :: file_root='chains/test'

        !> The number of live points
        integer :: nlive =500

        !> The degree of feedback to provide
        integer :: feedback = 1

        !> The degree of precision in the final answer
        double precision :: precision_criterion = 1d-5

        !> The maximum number of dead points/samples
        !!
        !! Set equal to -1 for no maximum number
        integer :: max_ndead = -1

        !> The maximum number of posterior points
        !!
        !! This is for memory allocation purposes, it won't necessarily have
        !! this many points if they're not 'good enough'
        integer :: nmax_posterior = 100000

        !> The minimum weight of the posterior points
        double precision :: minimum_weight = 1d-15

        !> Whether or not to save the dead points
        !! 
        !! It may not be worth saving them when doing extremely high dimensional problems
        logical :: save_dead = .false.

        !> The number of chords to draw
        integer :: num_chords = 6

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
        !! set of live points (live_data) in order to generate a baby_point
        !! uniformly sampled from within the loglikelihood contour specifed by
        !! loglikelihood bound contained at the M%l1 index of seed_point
        function samp(settings, seed_point, min_max_array, M,feedback) result(baby_point)

            import :: model
            import :: program_settings
            implicit none

            ! ------- Inputs -------
            !> program settings (mostly useful to pass on the number of live points)
            class(program_settings), intent(in) :: settings

            !> The details of the model (e.g. number of dimensions,loglikelihood,etc)
            type(model),            intent(in) :: M

            !> The seed point
            double precision, intent(in), dimension(M%nTotal)   :: seed_point

            !> The minimum and maximum values from each of the live points
            double precision, intent(in),    dimension(:,:)   :: min_max_array

            !> Optional argument to cause the sampler to print out relevent information
            integer, intent(in), optional :: feedback

            ! ------- Outputs -------
            !> The newly generated point
            double precision,    dimension(M%nTotal)     :: baby_point


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
        subroutine ev(settings,new_loglikelihood,old_loglikelihood,ndead,more_samples_needed,evidence_vec)

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

            ! vector containing [evidence, evidence error]
            double precision, dimension(6) :: evidence_vec

            ! ------- Outputs ------- 
            !> Whether we have obtained enough samples for an accurate evidence
            logical,intent(out) :: more_samples_needed


        end subroutine ev
    end interface

end module settings_module 
