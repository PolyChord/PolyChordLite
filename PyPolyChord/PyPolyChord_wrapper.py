from PyPolyChord.output import PolyChordOutput
import numpy

try:
    from PyPolyChord import _PyPolyChord
except ImportError as e:
    print(ImportError)
    if str(e) == 'No module named _PyPolyChord':
        print('')
        print('   Could not load Python Extension _PyPolyChord.so')
        print('')
        print('   You have to build it first with:')
        print('')
        print('   $   make PyPolyChord')                  
        print('')
        print('   In the base PolyChord directory.')
        print('')

    elif str(e) == 'libchord.so: cannot open shared object file: No such file or directory':
        print('')
        print('   Could not load PolyChord library "libchord.so"')
        print('')
        print('   You have to build it first,')
        print('   and point the LD_LIBRARY_PATH environment variable to it:')
        print('')
        print('   /-- BASH: --------------------------------------------\\')
        print('   |                                                     |')
        print('   |$   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib |')
        print('   |                                                     |')
        print('   +-- CSH: ---------------------------------------------+')
        print('   |                                                     |')
        print('   |$   setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH:$PWD/lib |')
        print('   |                                                     |')
        print('   \\-----------------------------------------------------/')
        print('')
    else:
        raise e

    import sys
    sys.exit(1)

def mpi_notification():
    print('/.===============================================================================.\\')
    print('||                                                                               ||')
    print('||    IMPORTANT: PyPolyChord settings                                            ||')
    print('||    -------------------------------                                            ||')
    print('||                                                                               ||')
    print('||    If you get MPI errors, try setting LD_PRELOAD to your mpi installation     ||')
    print('||                                                                               ||')
    print('||   /-- BASH: --------------------------------------------------------\\         ||')
    print('||   |                                                                 |         ||')
    print('||   |$   export LD_PRELOAD=/usr/lib/openmpi/lib/libmpi.so:$LD_PRELOAD |         ||')
    print('||   |                                                                 |         ||')
    print('||   +-- CSH: ---------------------------------------------------------+         ||')
    print('||   |                                                                 |         ||')
    print('||   |$   setenv LD_PRELOAD /usr/lib/openmpi/lib/libmpi.so:$LD_PRELOAD |         ||')
    print('||   |                                                                 |         ||')
    print('||   \\--------------------------------------------------------------- /          ||')
    print('||                                                                               ||')
    print('||   Where /usr/lib/openmpi/lib/libmpi.so should be replaced with the            ||')
    print('||   appropriate loctaion of libmpi.so on your system.                           ||')
    print('||                                                                               ||')
    print('\\.===============================================================================./  ')
    print('')

def default_prior(cube):
    theta = cube
    return theta 


def run_polychord(loglikelihood, nDims, nDerived, settings, prior=default_prior):
    """
    Runs PolyChord.

    For full details of nested sampling and PolyChord, please refer to:

    * PolyChord paper: http://arxiv.org/abs/1506.00171
    * Nested Sampling paper: http://projecteuclid.org/euclid.ba/1340370944

    To run in mpi mode, just run your script with mpirun as usual.
    Make sure that PyPolyChord is compiled with MPI:
    $ make veryclean
    $ make PyPolyChord MPI=1

    If MPI fails with some kind of library error, please see the text in
    mpi_notification() above

    Users are also required to cite the PolyChord papers:
    arXiv:1502.01856
    arXiv:1506.00171
    in their publications.


    Parameters
    ----------

    loglikelihood: function
        This function computes the log-likelihood of the model and derived
        parameters (phi) from the physical coordinates (theta).

        Parameters
        ----------
        theta: float list
            physical coordinate. A list of of length nDims.

        Returns
        -------
        (logL, phi): (float, float list)
            return is a 2-tuple of the log-likelihood (logL) and the derived
            parameters (phi). phi is a list of length nDerived.

    nDims: int
        Dimensionality of the model, i.e. the number of physical parameters.

    nDerived: int
        The number of derived parameters (can be 0).

    settings: settings.Settings
        Object containing polychord settings

    Optional Arguments
    ------------------

    prior: function
        This function converts from the unit hypercube to the physical
        parameters.
        (Default: identity function => prior(cube) == cube )

        Parameters
        ----------
        cube: float list
            coordinates in the unit hypercube. A list of length nDims.

        Returns
        -------
        theta: float list
            physical coordinates. A list of length nDims.

    Returns
    -------
    None. (in Python)

    All output is currently produced in the form of text files in <base_dir>
    directory. If you would like to contribute to PyPolyChord and improve this,
    please get in touch:

    Will Handley: wh260@cam.ac.uk

    In general the contents of <base_dir> is a set of getdist compatible files.

    <root> = <base_dir>/<file_root>

    <root>.txt                                               (posteriors = True)
        Weighted posteriors. Compatible with getdist. Each line contains:
            -2*loglikelihood, weight, <physical parameters>, <derived parameters>
        Note that here the weights are not normalised, but instead are chosen
        so that the maximum weight is 1.0.

    <root>_equal_weights.txt                                     (equals = True)
        Weighted posteriors. Compatible with getdist. Each line contains:
            -2*loglikelihood, 1.0, <physical parameters>, <derived parameters>

    <root>_dead.txt                                          (write_dead = True)
        Dead points. Each line contains:
            loglikelihood, <physical parameters>, <derived parameters>

    <root>_phys_live.txt                                     (write_live = True)
        Live points. Each line contains:
            <physical parameters>, <derived parameters>, loglikelihood
        Mostly useful for run-time monitoring.

    <root>.resume
        Resume file. Human readable.

    <root>.stats
        Final output evidence statistics

    """

    # Test likelihood
    print("testing prior...")
    x = list(numpy.random.rand(nDims))
    theta = prior(x)
    if not isinstance(theta,list):
        raise TypeError("Return from prior must be a list")
    elif not all(isinstance(item,float) for item in theta):
        raise TypeError("Return from prior must be a list of floats")
    ret = loglikelihood(theta)
    if not isinstance(ret,tuple):
        raise TypeError("Return from likelihood must be a tuple of (<likelihood value>, <derived parameters>)")
    elif not len(ret) == 2:
        raise TypeError("Return from likelihood must be a tuple of (<likelihood value>, <derived parameters>)")
    elif not isinstance(ret[0],float):
        raise TypeError("Likelihood value (first of tuple returned from loglikelihood) must be a float")
    elif not isinstance(ret[1],list):
        raise TypeError("Derived parameters (second of tuple returned from loglikelihood) must be a list")
    elif not all(isinstance(item,float) for item in ret[1]):
        raise TypeError("Derived parameters (second of tuple returned from loglikelihood) must be a list of floats")

    # Run polychord from module library
    _PyPolyChord.run(loglikelihood,
                     prior,
                     nDims,
                     nDerived,
                     settings.nlive,
                     settings.num_repeats,
                     settings.do_clustering,
                     settings.feedback,
                     settings.precision_criterion,
                     settings.max_ndead,
                     settings.boost_posterior,
                     settings.posteriors,
                     settings.equals,
                     settings.cluster_posteriors,
                     settings.write_resume,
                     settings.write_paramnames,
                     settings.read_resume,
                     settings.write_stats,
                     settings.write_live,
                     settings.write_dead,
                     settings.update_files,
                     settings.base_dir,
                     settings.file_root,
                     settings.grade_frac,
                     settings.grade_dims)

    return PolyChordOutput(settings.base_dir,settings.file_root)


def run_nested_sampling(loglikelihood, nDims, nDerived, **kwargs):
    """
    Runs PolyChord (legacy interface).

    For full details of nested sampling and PolyChord, please refer to:

    * PolyChord paper: http://arxiv.org/abs/1506.00171
    * Nested Sampling paper: http://projecteuclid.org/euclid.ba/1340370944

    To run in mpi mode, just run your script with mpirun as usual.
    Make sure that PyPolyChord is compiled with MPI:
    $ make veryclean
    $ make PyPolyChord MPI=1

    If MPI fails with some kind of library error, please see the text in
    mpi_notification() above

    Users are also required to cite the PolyChord papers:
    arXiv:1502.01856
    arXiv:1506.00171
    in their publications.


    Parameters
    ----------

    loglikelihood: function
        This function computes the log-likelihood of the model and derived
        parameters (phi) from the physical coordinates (theta).

        Parameters
        ----------
        theta: float list
            physical coordinate. A list of of length nDims.

        Returns
        -------
        (logL, phi): (float, float list)
            return is a 2-tuple of the log-likelihood (logL) and the derived
            parameters (phi). phi is a list of length nDerived.

    nDims: int
        Dimensionality of the model, i.e. the number of physical parameters.

    nDerived: int
        The number of derived parameters (can be 0).


    Keyword arguments
    -----------------

    prior: function
        This function converts from the unit hypercube to the physical
        parameters.
        (Default: identity function => prior(cube) == cube )

        Parameters
        ----------
        cube: float list
            coordinates in the unit hypercube. A list of length nDims.

        Returns
        -------
        theta: float list
            physical coordinates. A list of length nDims.


    nlive: int
        (Default: nDims*25)
        The number of live points.
        Increasing nlive increases the accuracy of posteriors and evidences,
        and proportionally increases runtime ~ O(nlive).

    num_repeats : int
        (Default: nDims*5)
        The number of slice slice-sampling steps to generate a new point.
        Increasing num_repeats increases the reliability of the algorithm.
        Typically
        * for reliable evidences need num_repeats ~ O(5*nDims).
        * for reliable posteriors need num_repeats ~ O(nDims)

    do_clustering : boolean
        (Default: True)
        Whether or not to use clustering at run time.

    feedback : {0,1,2,3}
        (Default: 1)
        How much command line feedback to give

    precision_criterion : float
        (Default: 0.001)
        Termination criterion. Nested sampling terminates when the evidence
        contained in the live points is precision_criterion fraction of the
        total evidence.

    max_ndead : int
        (Default: -1)
        Alternative termination criterion. Stop after max_ndead iterations.
        Set negative to ignore (default).

    boost_posterior : float
        (Default: 0.0)
        Increase the number of posterior samples produced.  This can be set
        arbitrarily high, but you won't be able to boost by more than
        num_repeats
        Warning: in high dimensions PolyChord produces _a lot_ of posterior
        samples. You probably don't need to change this

    posteriors : boolean
        (Default: True)
        Produce (weighted) posterior samples. Stored in <root>.txt.

    equals : boolean
        (Default: True)
        Produce (equally weighted) posterior samples. Stored in
        <root>_equal_weights.txt

    cluster_posteriors : boolean
        (Default: True)
        Produce posterior files for each cluster?
        Does nothing if do_clustering=False.

    write_resume : boolean
        (Default: True)
        Create a resume file.

    read_resume : boolean
        (Default: True)
        Read from resume file.

    write_stats : boolean
        (Default: True)
        Write an evidence statistics file.

    write_live : boolean
        (Default: True)
        Write a live points file.

    write_dead : boolean
        (Default: True)
        Write a dead points file.

    update_files : int
        (Default: nlive)
        How often to update the files in <base_dir>.

    base_dir : string
        (Default: 'chains')
        Where to store output files.

    file_root : string
        (Default: 'test')
        Root name of the files produced.

    grade_frac : List[float]
        (Default: 1)
        The amount of time to spend in each speed.

    grade_dims : List[int]
        (Default: 1)
        The number of parameters within each speed.

    Returns
    -------
    None. (in Python)

    All output is currently produced in the form of text files in <base_dir>
    directory. If you would like to contribute to PyPolyChord and improve this,
    please get in touch:

    Will Handley: wh260@cam.ac.uk

    In general the contents of <base_dir> is a set of getdist compatible files.

    <root> = <base_dir>/<file_root>

    <root>.txt                                               (posteriors = True)
        Weighted posteriors. Compatible with getdist. Each line contains:
            -2*loglikelihood, weight, <physical parameters>, <derived parameters>
        Note that here the weights are not normalised, but instead are chosen
        so that the maximum weight is 1.0.

    <root>_equal_weights.txt                                     (equals = True)
        Weighted posteriors. Compatible with getdist. Each line contains:
            -2*loglikelihood, 1.0, <physical parameters>, <derived parameters>

    <root>_dead.txt                                          (write_dead = True)
        Dead points. Each line contains:
            loglikelihood, <physical parameters>, <derived parameters>

    <root>_phys_live.txt                                     (write_live = True)
        Live points. Each line contains:
            <physical parameters>, <derived parameters>, loglikelihood
        Mostly useful for run-time monitoring.

    <root>.resume
        Resume file. Human readable.

    <root>.stats
        Final output evidence statistics

    """

    prior = kwargs.pop('prior', default_prior)
    nlive = kwargs.pop('nlive', nDims*25)
    num_repeats = kwargs.pop('num_repeats', nDims*5)
    do_clustering = kwargs.pop('do_clustering', True)
    feedback = kwargs.pop('feedback', 1)
    precision_criterion = kwargs.pop('precision_criterion', 0.001)
    max_ndead = kwargs.pop('max_ndead', -1)
    boost_posterior = kwargs.pop('boost_posterior', 0.0)
    posteriors = kwargs.pop('posteriors', True)
    equals = kwargs.pop('equals', True)
    cluster_posteriors = kwargs.pop('cluster_posteriors', True)
    write_resume = kwargs.pop('write_resume', True)
    write_paramnames = kwargs.pop('write_paramnames', False)
    read_resume = kwargs.pop('read_resume', True)
    write_stats = kwargs.pop('write_stats', True)
    write_live = kwargs.pop('write_live', True)
    write_dead = kwargs.pop('write_dead', True)
    update_files = kwargs.pop('update_files', nlive)
    base_dir = kwargs.pop('base_dir', 'chains')
    file_root = kwargs.pop('file_root', 'test')
    grade_dims = kwargs.pop('grade_dims', [nDims])
    grade_frac = kwargs.pop('grade_frac', [1.0]*len(grade_dims))

    if kwargs:
        raise TypeError('Unexpected **kwargs in Contours constructor: %r' % kwargs)

    if len(grade_frac) != len(grade_dims):
        raise ValueError('grade_dims and grade_frac must be the same length')

    if sum(grade_dims) != nDims:
        raise ValueError('grade_dims must sum to the total dimensionality: sum(' + str(grade_dims) + ') /= %i' % nDims)


    # Run polychord from module library
    _PyPolyChord.run(loglikelihood, prior, nDims, nDerived, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, base_dir, file_root, grade_frac, grade_dims)

    return None

