from .output import PolyChordOutput
import numpy
import pkg_resources
import sys
import os

ld_library_path = pkg_resources.resource_string(__name__,'.ld_library_path.sh').decode('utf-8') 
ld_preload = pkg_resources.resource_string(__name__,'.ld_preload.sh').decode('utf-8') 

try:
    import _PyPolyChord
except ImportError as e:
    print(ImportError)
    if str(e) == 'libchord.so: cannot open shared object file: No such file or directory':
        print('PolyChord: You need to configure your LD_LIBRARY_PATH')
        print('           to point to your PolyChord installation.')
        print('           e.g. add this to your .bashrc: ')
        print('')
        print(ld_library_path)
        sys.exit(1)
    else:
        raise e

try:
    sys_ld_preload = os.environ['LD_PRELOAD']
except KeyError:
    sys_ld_preload = ''

if ld_preload and 'libmpi.so' not in sys_ld_preload:
    print('PolyChord: You need to configure your LD_PRELOAD')
    print('           to point to your mpi installation.')
    print('           e.g. add this to your .bashrc:')
    print('')
    print(ld_preload)
    sys.exit(1)


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

    if not os.path.exists(settings.base_dir):
        os.makedirs(settings.base_dir)

    if not os.path.exists(settings.cluster_dir):
        os.makedirs(settings.cluster_dir)

    # Run polychord from module library
    _PyPolyChord.run(loglikelihood,
                     prior,
                     nDims,
                     nDerived,
                     settings.nlive,
                     settings.num_repeats,
                     settings.nprior,
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
                     settings.write_prior,
                     settings.update_files,
                     settings.base_dir,
                     settings.file_root,
                     settings.grade_frac,
                     settings.grade_dims)

    return PolyChordOutput(settings.base_dir,settings.file_root)
