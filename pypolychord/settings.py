import os
import numpy
import warnings

def warn(msg):
    warnings.warn(msg)
    warnings.warn('This will raise an exception in the future.')


class PolyChordSettings:
    """
    PolyChord settings

    For full details of nested sampling and PolyChord, please refer to:

    Parameters
    ----------
    nDims: int
        Dimensionality of the model, i.e. the number of physical parameters.

    nDerived: int
        The number of derived parameters (can be 0).


    Keyword arguments
    -----------------
    nlive: positive int
        (Default: nDims*25)
        The number of live points.
        Increasing nlive increases the accuracy of posteriors and evidences,
        and proportionally increases runtime ~ O(nlive).
        Must be positive.

    num_repeats : int
        (Default: nDims*5)
        The number of slice slice-sampling steps to generate a new point.
        Increasing num_repeats increases the reliability of the algorithm.
        -1 is equivalent to default.
        Typically
        * for reliable evidences need num_repeats ~ O(5*nDims).
        * for reliable posteriors need num_repeats ~ O(nDims)

    nprior : int
        (Default: nlive)
        The number of prior samples to draw before starting compression.
        -1 is equivalent to default.

    nfail : int
        (Default: nlive)
        The number of failed spawns before stopping nested sampling.
        -1 is equivalent to default.


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

    logzero : float
        (Default: -1e30)
        The loglikelihood value at which PolyChord considers points
        'unphysical', and excludes them at the prior level.

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

    write_prior : boolean
        (Default: True)
        Write a prior points file.

    maximise : boolean
        (Default: False)
        Perform maximisation at the end of the run to find the maximum
        likelihood point and value

    compression_factor : double
        (Default: exp(-1))
        How often to update the files and do clustering

    base_dir : string
        (Default: 'chains')
        Where to store output files.

    file_root : string
        (Default: 'test')
        Root name of the files produced.

    grade_frac : List[positive float]
        (Default: [1])
        The amount of time to spend in each speed.
        If any of grade_frac are <= 1, then polychord will time each sub-speed,
        and then choose num_repeats for the number of slowest repeats, and
        spend the proportion of time indicated by grade_frac. Otherwise this
        indicates the number of repeats to spend in each speed.
        Must be positive.

    grade_dims : List[positive int]
        (Default: nDims)
        The number of parameters within each speed.
        Must be positive.

    nlives : dict {double:int}
        (Default: {})
        Variable number of live points option. This dictionary is a mapping
        between loglike contours and nlive.
        You should still set nlive to be a sensible number, as this indicates
        how often to update the clustering, and to define the default value.

    seed : positive int
        (Default: system time in milliseconds)
        Choose the seed to seed the random number generator.
        Note **Positive seeds only**
        a negative seed indicates that you should use the system time in
        milliseconds

    """

    def __init__(self, nDims, nDerived, **kwargs):
        self.nlive = kwargs.pop('nlive', nDims*25)
        self.num_repeats = kwargs.pop('num_repeats', nDims*5)
        self.nprior = kwargs.pop('nprior', -1)
        self.nfail = kwargs.pop('nfail', -1)
        self.do_clustering = kwargs.pop('do_clustering', True)
        self.feedback = kwargs.pop('feedback', 1)
        self.precision_criterion = kwargs.pop('precision_criterion', 0.001)
        self.logzero = kwargs.pop('logzero', -1e30)
        self.max_ndead = kwargs.pop('max_ndead', -1)
        self.boost_posterior = kwargs.pop('boost_posterior', 0.0)
        self.posteriors = kwargs.pop('posteriors', True)
        self.equals = kwargs.pop('equals', True)
        self.cluster_posteriors = kwargs.pop('cluster_posteriors', True)
        self.write_resume = kwargs.pop('write_resume', True)
        self.write_paramnames = kwargs.pop('write_paramnames', False)
        self.read_resume = kwargs.pop('read_resume', True)
        self.write_stats = kwargs.pop('write_stats', True)
        self.write_live = kwargs.pop('write_live', True)
        self.write_dead = kwargs.pop('write_dead', True)
        self.write_prior = kwargs.pop('write_prior', True)
        self.maximise = kwargs.pop('maximise', False)
        self.compression_factor = kwargs.pop('compression_factor',
                                             numpy.exp(-1))
        self.base_dir = kwargs.pop('base_dir', 'chains')
        self.file_root = kwargs.pop('file_root', 'test')
        self.seed = kwargs.pop('seed', -1)
        self.grade_dims = list(kwargs.pop('grade_dims',
                                          [nDims]))
        self.grade_frac = list(kwargs.pop('grade_frac',
                                          [1.0]*len(self.grade_dims)))
        self.nlives = kwargs.pop('nlives', {})
        if kwargs:
            warn('Unexpected **kwargs in Contours constructor: %r' % kwargs)
        if sum(self.grade_dims) != nDims:
            raise ValueError('grade_dims must sum to the total dimensionality:'
                             'sum(grade_dims) = %i /= %i' %
                             (len(self.grade_dims), nDims))
        if not (isinstance(self.grade_dims, list)
                or not all([isinstance(x, int) for x in self.grade_dims])):
            raise ValueError('grade_dims must be a list of integers.')
        if not (isinstance(self.grade_frac, list)
                or not all(
                    [isinstance(x, int) or isinstance(x, float)
                     for x in self.grade_frac]):
                raise ValueError('grade_dims must be a list of doubles.')
        if not (len(self.grade_dims) == len(self.grade_frac)):
                raise ValueError('grade_dims and grade_frac must have the same'
                                 'len.')


    @property
    def cluster_dir(self):
        return os.path.join(self.base_dir, 'clusters')

    @property
    def nDims(self):
        return sum(self.grade_dims)

    @property
    def grade_dims(self):
        return self._grade_dims

    @grade_dims.setter
    def grade_dims(self, value):
        if value is None or value == []:
            warn('invalid grade_dims. Defaulting.')
            self._grade_dims = [self.nDims]
        else:
            self._grade_dims = [__natnum(x) for x in value]
        # try:
        #     self.grade_frac = self.grade_frac[:len(self.grade_dims)]
        # except ValueError:
        #     self.grade_frac = self.grade_frac[:len(self.grade_dims)] + \
        #         [1.0]*(len(self.grade_dims) - len(self.grade_frac))

    @property
    def grade_frac(self):
        if len(self.grade_dims) != self._grade_frac:
            warn('grade_dims doesn\'t match grade_frac.')
        return self._grade_frac

    @grade_frac.setter
    def grade_frac(self, value):
        if len(value) >= len(self.grade_dims):
            self._grade_frac = value[:len(self.grade_dims)]
        else:
            warn('Insufficient values to set grade_frac.'
                 'Need %i, got %i' % (len(value), len(self.grade_dims)))
            warn('Defaulting to set, but this will raise in the future.')


def __natnum(x: int, minimum=1):
    """Helper function to prevent negative or zero numbers being passed in."""
    if x < minimum:
        warn("Expecting a non-zero positive integer: got %i." % x)
    return x
