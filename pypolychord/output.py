try:
    import getdist.mcsamples
except ImportError:
    pass
import re
import os


class PolyChordOutput:
    def __init__(self, base_dir, file_root):
        """
        Returns
        -------
        None. (in Python)

        All output is currently produced in the form of text files in
        <base_dir> directory. If you would like to contribute to pypolychord
        and improve this, please get in touch:

        Will Handley: wh260@cam.ac.uk

        In general the contents of <base_dir> is a set of getdist compatible
        files.

        <root> = <base_dir>/<file_root>

        <root>.txt                                          (posteriors = True)
            Weighted posteriors. Compatible with getdist. Each line contains:

                weight, -2*loglikelihood, <physical params>, <derived params>

            Note that here the weights are not normalised, but instead are
            chosen so that the maximum weight is 1.0.

        <root>_equal_weights.txt                                (equals = True)
            Weighted posteriors. Compatible with getdist. Each line contains:
                1.0, -2*loglikelihood, <physical params>, <derived params>

        <root>.stats
            Final output evidence statistics

        """
        self.base_dir = base_dir
        self.file_root = file_root

        with open('%s.stats' % self.root, 'r') as f:
            for _ in range(9):
                line = f.readline()

            self.logZ = float(line.split()[2])
            self.logZerr = float(line.split()[4])

            for _ in range(6):
                line = f.readline()

            self.logZs = []
            self.logZerrs = []
            while line[:5] == 'log(Z':
                self.logZs.append(float(re.findall(r'=(.*)', line
                                                   )[0].split()[0]))
                self.logZerrs.append(float(re.findall(r'=(.*)', line
                                                      )[0].split()[2]))

                line = f.readline()

            for _ in range(5):
                f.readline()

            self.ncluster = len(self.logZs)
            self.nposterior = int(f.readline().split()[1])
            self.nequals = int(f.readline().split()[1])
            self.ndead = int(f.readline().split()[1])
            self.nlive = int(f.readline().split()[1])
            # Protect against ValueError when .stats file has ******* for nlike
            # (occurs when nlike is too big).
            try:
                self.nlike = int(f.readline().split()[1])
            except ValueError:
                self.nlike = None
            line = f.readline()
            self.avnlike = float(line.split()[1])
            try:
                self.avnlikeslice = float(line.split()[3])
            except ValueError:
                self.avnlikeslice = None

    @property
    def root(self):
        return os.path.join(self.base_dir, self.file_root)

    def cluster_root(self, i):
        return os.path.join(self.base_dir, 'clusters',
                            '%s_%i' % (self.file_root, i))

    @property
    def paramnames_file(self):
        return self.root + '.paramnames'

    @property
    def posterior(self):
        """
        Return a getdist MCSample object.
        Allows access to several chain statistics
        and plotting. 

        NOTE: calling posterior.loglikes will return -2 * loglike 

        :returns: getdist.MCSample object
        :rtype: 

        """

        return getdist.mcsamples.loadMCSamples(self.root)

    @property
    def loglikes(self):
        """
        the log likelihood values of the samples

        :returns: and array of log likelihood values
        :rtype: 

        """

        # grab the getdist mcsample object
        posterior = getdist.mcsamples.loadMCSamples(self.root)

        # note that getdist incorrectly labels
        # -2 loglike as loglike 

        return -0.5 * posterior.loglikes

    def cluster_posterior(self, i):
        return getdist.mcsamples.loadMCSamples(self.cluster_root(i))

    def cluster_paramnames_file(self, i):
        return self.cluster_root(i) + '.paramnames'

    def make_paramnames_files(self, paramnames):
        self.make_paramnames_file(paramnames, self.paramnames_file)
        for i, _ in enumerate(self.logZs):
            self.make_paramnames_file(paramnames,
                                      self.cluster_paramnames_file(i))

    def make_paramnames_file(self, paramnames, filename):
        with open(filename, 'w') as f:
            for name, latex in paramnames:
                f.write('%s   %s\n' % (name, latex))
