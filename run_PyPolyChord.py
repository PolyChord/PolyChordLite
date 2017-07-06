from numpy import pi, log, sqrt
import PyPolyChord
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior

nDims = 3
nDerived = 1


def likelihood(theta):
    """ Simple Gaussian Likelihood"""

    phi = [0.0] * nDerived

    sigma = 0.1
    nDims = len(theta)

    r2 = sum([t**2 for t in theta])
    phi[0] = sqrt(r2)

    logL = -log(2*pi*sigma*sigma)*nDims/2.0
    logL += -r2/2/sigma/sigma

    return logL, phi


def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """

    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = UniformPrior(-1, 1)(x)

    return theta

settings = PolyChordSettings(nDims, nDerived)
settings.file_root = 'gaussian'
settings.do_clustering = True

output = PyPolyChord.run_polychord(likelihood, nDims, nDerived, settings, prior)
paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
paramnames += [('r*', 'r')]
output.make_paramnames_files(paramnames)

try:
    import getdist.plots
    import matplotlib.pyplot as plt
    posterior = output.posterior
    g = getdist.plots.getSubplotPlotter()
    g.triangle_plot(posterior, filled=True)
    plt.show()
except ImportError:
    print("Install matplotlib and getdist for plotting examples")
