from pypolychord.output import PolyChordOutput
import os
import _pypolychord
import numpy as np


def default_prior(cube):
    return cube.copy()


def default_dumper(live, dead, logweights, logZ, logZerr):
    pass


def run_polychord(loglikelihood, nDims, nDerived, settings,
                  prior=default_prior, dumper=default_dumper):
    """
    Runs PolyChord.

    For full details of nested sampling and PolyChord, please refer to:

    * PolyChord paper: http://arxiv.org/abs/1506.00171
    * Nested Sampling paper: http://projecteuclid.org/euclid.ba/1340370944

    To run in mpi mode, just run your script with mpirun as usual.
    Make sure that pypolychord is compiled with MPI:
    $ make veryclean
    $ make pypolychord MPI=1

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
        theta: numpy.array
            physical coordinate. Length nDims.

        Returns
        -------
        (logL, phi): (float, array-like)
            return is a 2-tuple of the log-likelihood (logL) and the derived
            parameters (phi). phi length nDerived.

        Returns
        -------
        logL: float
            log-likelihood

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
        cube: numpy.array
            coordinates in the unit hypercube. Length nDims.

        Returns
        -------
        theta: array-like
            physical coordinates. Length nDims.

    dumper: function
        This function gives run-time access to the posterior and live points.

        Parameters
        ----------
        live: numpy.array
            The live points and their loglikelihood birth and death contours
            Shape (nlive, nDims+nDerived+2)
        dead: numpy.array
            The dead points and their loglikelihood birth and death contours
            Shape (ndead, nDims+nDerived+2)
        logweights: numpy.array
            The posterior weights of the dead points
            Shape (ndead)
        logZ: float
            The current log-evidence estimate
        logZerr: float
            The current log-evidence error estimate

    Returns
    -------
    None. (in Python)

    All output is currently produced in the form of text files in <base_dir>
    directory. If you would like to contribute to pypolychord and improve this,
    please get in touch:

    Will Handley: wh260@cam.ac.uk

    In general the contents of <base_dir> is a set of getdist compatible files.

    <root> = <base_dir>/<file_root>

    <root>.txt                                              (posteriors = True)
        Weighted posteriors. Compatible with getdist. Each line contains:

          weight, -2*loglikelihood, <physical parameters>, <derived parameters>

        Note that here the weights are not normalised, but instead are chosen
        so that the maximum weight is 1.0.

    <root>_equal_weights.txt                                    (equals = True)
        Weighted posteriors. Compatible with getdist. Each line contains:
            1.0, -2*loglikelihood, <physical parameters>, <derived parameters>

    <root>_dead.txt                                         (write_dead = True)
        Dead points. Each line contains:
            loglikelihood, <physical parameters>, <derived parameters>

    <root>_phys_live.txt                                    (write_live = True)
        Live points. Each line contains:
            <physical parameters>, <derived parameters>, loglikelihood
        Mostly useful for run-time monitoring.

    <root>.resume
        Resume file. Human readable.

    <root>.stats
        Final output evidence statistics

    """

    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    except ImportError:
        rank = 0

    try:
        if rank == 0:
            os.makedirs(settings.base_dir)
    except OSError:
        pass
        
    try:
        if rank == 0:
            os.makedirs(settings.cluster_dir)
    except OSError:
        pass

    if settings.cube_samples is not None:
        make_resume_file(settings, loglikelihood, prior)
        read_resume = settings.read_resume
        settings.read_resume=True



    def wrap_loglikelihood(theta, phi):
        logL, phi[:] = loglikelihood(theta)
        return logL

    def wrap_prior(cube, theta):
        theta[:] = prior(cube)

    # Run polychord from module library
    _pypolychord.run(wrap_loglikelihood,
                     wrap_prior,
                     dumper,
                     nDims,
                     nDerived,
                     settings.nlive,
                     settings.num_repeats,
                     settings.nprior,
                     settings.nfail,
                     settings.do_clustering,
                     settings.feedback,
                     settings.precision_criterion,
                     settings.logzero,
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
                     settings.maximise,
                     settings.compression_factor,
                     settings.synchronous,
                     settings.base_dir,
                     settings.file_root,
                     settings.grade_frac,
                     settings.grade_dims,
                     settings.nlives,
                     settings.seed)

    if settings.cube_samples is not None:
        settings.read_resume = read_resume

    return PolyChordOutput(settings.base_dir, settings.file_root)


def make_resume_file(settings, loglikelihood, prior):
    import fortranformat as ff
    resume_filename = os.path.join(settings.base_dir,
                                   settings.file_root)+".resume"

    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    except ImportError:
        rank = 0
        size = 1

    lives = []
    logL_birth = settings.logzero
    for i in np.array_split(np.arange(len(settings.cube_samples)), size)[rank]:
        cube = settings.cube_samples[i]
        theta = prior(cube)
        logL, derived = loglikelihood(theta)
        nDims = len(theta)
        nDerived = len(derived)
        lives.append(np.concatenate([cube,theta,derived,[logL_birth, logL]]))

    try:
        sendbuf = np.array(lives).flatten()
        sendcounts = np.array(comm.gather(len(sendbuf)))
        if rank == 0:
            recvbuf = np.empty(sum(sendcounts), dtype=int)
        else:
            recvbuf = None
        comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)

        lives = np.reshape(sendbuf, (len(settings.cube_samples), len(lives[0])))
    except NameError:
        lives = np.array(lives)
                

    def write(f, var, fmt):
        try:
            i = len(var)
        except TypeError:
            var = [var]
            i = 1
        fmt = '(%i%s)' % (i, fmt)
        writer = ff.FortranRecordWriter(fmt)
        f.write(writer.write(var) + '\n')
    i = 'I12'
    d = 'E24.15E3'

    if rank == 0:
        with open(resume_filename,"w") as f:
            f.write('=== Number of dimensions ===\n')
            write(f, nDims, i)
            f.write('=== Number of derived parameters ===\n')
            write(f, nDerived, i)
            f.write('=== Number of dead points/iterations ===\n')
            write(f, 0, i)
            f.write('=== Number of clusters ===\n')
            write(f, 1, i)
            f.write('=== Number of dead clusters ===\n')
            write(f, 0, i)
            f.write('=== Number of global weighted posterior points ===\n')
            write(f, 0, i)
            f.write('=== Number of global equally weighted posterior points ===\n')
            write(f, 0, i)
            f.write('=== Number of grades ===\n')
            write(f, len(settings.grade_dims), i)
            f.write('=== positions of grades ===\n')
            write(f, settings.grade_dims, i)
            f.write('=== Number of repeats ===\n')
            write(f, settings.num_repeats, i)
            f.write('=== Number of likelihood calls ===\n')
            write(f, len(lives), i)
            f.write('=== Number of live points in each cluster ===\n')
            write(f, len(lives), i)
            f.write('=== Number of phantom points in each cluster ===\n')
            write(f, 0, i)
            f.write('=== Number of weighted posterior points in each cluster ===\n')
            write(f, 0, i)
            f.write('=== Number of equally weighted posterior points in each cluster ===\n')
            write(f, 0, i)
            f.write('=== Minimum loglikelihood positions ===\n')
            write(f, np.argmin(lives[:,-1]), i)
            f.write('=== Number of weighted posterior points in each dead cluster ===\n')
            f.write('=== Number of equally weighted posterior points in each dead cluster ===\n')
            f.write('=== global evidence -- log(<Z>) ===\n')
            write(f, settings.logzero, d)
            f.write('=== global evidence^2 -- log(<Z^2>) ===\n')
            write(f, settings.logzero, d)
            f.write('=== posterior thin factor ===\n')
            write(f, settings.boost_posterior, d)
            f.write('=== local loglikelihood bounds ===\n')
            write(f, lives[:,-1].min(), d)
            f.write('=== local volume -- log(<X_p>) ===\n')
            write(f, 0.0, d)
            f.write('=== last update volume ===\n')
            write(f, 0.0, d)
            f.write('=== global evidence volume cross correlation -- log(<ZX_p>) ===\n')
            write(f, settings.logzero, d)
            f.write('=== local evidence -- log(<Z_p>) ===\n')
            write(f, settings.logzero, d)
            f.write('=== local evidence^2 -- log(<Z_p^2>) ===\n')
            write(f, settings.logzero, d)
            f.write('=== local evidence volume cross correlation -- log(<Z_pX_p>) ===\n')
            write(f, settings.logzero, d)
            f.write('=== local volume cross correlation -- log(<X_pX_q>) ===\n')
            write(f, 0.0, d)
            f.write('=== maximum log weights -- log(w_p) ===\n')
            write(f, settings.logzero, d)
            f.write('=== local dead evidence -- log(<Z_p>) ===\n')
            f.write('=== local dead evidence^2 -- log(<Z_p^2>) ===\n')
            f.write('=== maximum dead log weights -- log(w_p) ===\n')
            f.write('=== covariance matrices ===\n')
            f.write('---------------------------------------\n')
            for x in np.identity(nDims):
                write(f, x, d)
            f.write('=== cholesky decompositions ===\n')
            f.write('---------------------------------------\n')
            for x in np.identity(nDims):
                write(f, x, d)
            f.write('=== live points ===\n')
            f.write('---------------------------------------\n')
            for x in lives:
                write(f, x, d)
            f.write('=== dead points ===\n')
            f.write('=== logweights of dead points ===\n')
            f.write('=== phantom points ===\n')
            f.write('---------------------------------------\n')
            f.write('=== weighted posterior points ===\n')
            f.write('---------------------------------------\n')
            f.write('=== dead weighted posterior points ===\n')
            f.write('=== global weighted posterior points ===\n')
            f.write('=== equally weighted posterior points ===\n')
            f.write('---------------------------------------\n')
            f.write('=== dead equally weighted posterior points ===\n')
            f.write('=== global equally weighted posterior points ===\n')
