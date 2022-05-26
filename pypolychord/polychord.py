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

    settings.grade_dims = [int(d) for d in settings.grade_dims]
    settings.nlives = {float(logL):int(nlive) for logL, nlive in settings.nlives.items()}

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

    if rank == 0:
        with open(resume_filename,"w") as f:
            def write(var):
                var = np.atleast_1d(var)
                if isinstance(var[0], np.integer):
                    fmt = '(%iI12)' % var.size
                elif isinstance(var[0], np.double):
                    fmt = '(%iE24.15E3)' % var.size
                else:
                    fmt = '(A)'
                writer = ff.FortranRecordWriter(fmt)
                f.write(writer.write(var) + '\n')

            write('=== Number of dimensions ===')
            write(nDims)
            write('=== Number of derived parameters ===')
            write(nDerived)
            write('=== Number of dead points/iterations ===')
            write(0)
            write('=== Number of clusters ===')
            write(1)
            write('=== Number of dead clusters ===')
            write(0)
            write('=== Number of global weighted posterior points ===')
            write(0)
            write('=== Number of global equally weighted posterior points ===')
            write(0)
            write('=== Number of grades ===')
            write(len(settings.grade_dims))
            write('=== positions of grades ===')
            write(settings.grade_dims)
            write('=== Number of repeats ===')
            write(settings.num_repeats)
            write('=== Number of likelihood calls ===')
            write(len(lives))
            write('=== Number of live points in each cluster ===')
            write(len(lives))
            write('=== Number of phantom points in each cluster ===')
            write(0)
            write('=== Number of weighted posterior points in each cluster ===')
            write(0)
            write('=== Number of equally weighted posterior points in each cluster ===')
            write(0)
            write('=== Minimum loglikelihood positions ===')
            write(np.argmin(lives[:,-1]))
            write('=== Number of weighted posterior points in each dead cluster ===')
            write('=== Number of equally weighted posterior points in each dead cluster ===')
            write('=== global evidence -- log(<Z>) ===')
            write(settings.logzero)
            write('=== global evidence^2 -- log(<Z^2>) ===')
            write(settings.logzero)
            write('=== posterior thin factor ===')
            write(settings.boost_posterior)
            write('=== local loglikelihood bounds ===')
            write(lives[:,-1].min())
            write('=== local volume -- log(<X_p>) ===')
            write(0.0)
            write('=== last update volume ===')
            write(0.0)
            write('=== global evidence volume cross correlation -- log(<ZX_p>) ===')
            write(settings.logzero)
            write('=== local evidence -- log(<Z_p>) ===')
            write(settings.logzero)
            write('=== local evidence^2 -- log(<Z_p^2>) ===')
            write(settings.logzero)
            write('=== local evidence volume cross correlation -- log(<Z_pX_p>) ===')
            write(settings.logzero)
            write('=== local volume cross correlation -- log(<X_pX_q>) ===')
            write(0.0)
            write('=== maximum log weights -- log(w_p) ===')
            write(settings.logzero)
            write('=== local dead evidence -- log(<Z_p>) ===')
            write('=== local dead evidence^2 -- log(<Z_p^2>) ===')
            write('=== maximum dead log weights -- log(w_p) ===')
            write('=== covariance matrices ===')
            write('---------------------------------------')
            for x in np.identity(nDims):
                write(x)
            write('=== cholesky decompositions ===')
            write('---------------------------------------')
            for x in np.identity(nDims):
                write(x)
            write('=== live points ===')
            write('---------------------------------------')
            for x in lives:
                write(x)
            write('=== dead points ===')
            write('=== logweights of dead points ===')
            write('=== phantom points ===')
            write('---------------------------------------')
            write('=== weighted posterior points ===')
            write('---------------------------------------')
            write('=== dead weighted posterior points ===')
            write('=== global weighted posterior points ===')
            write('=== equally weighted posterior points ===')
            write('---------------------------------------')
            write('=== dead equally weighted posterior points ===')
            write('=== global equally weighted posterior points ===')
