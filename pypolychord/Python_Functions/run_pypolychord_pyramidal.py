from numpy import log, exp, sqrt, pi
from scipy.special import loggamma
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior

#EGGBOX 

nDims = 3
nDerived = 1

functionName = 'pyramidal'

# param : array
def likelihood(theta):
    """ Simple pyramidal likelihood"""

    nDims = 0
    factor = 1
    sigma=0.1
    mu=0.5
    logSqrtTwoPi = log(sqrt(2*pi))

    if len(theta) != nDims:
        nDims = len(theta)
        factor = exp(-2/nDims * loggamma(1 + nDims/2)) * (pi/2)
    
    # normalisation
    logL =   - (logSqrtTwoPi+log(sigma))

    # theta dependence
    logL = logL - max(abs(theta-mu)/sigma)**2  / factor

    r2 = 0
    return logL, [r2] # float, array-like

# param : array
def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """
    return UniformPrior(-1, 1)(hypercube) # array

# param : array, array, array, float, float
def dumper(live, dead, logweights, logZ, logZerr):
    print("Last dead point:", dead[-1]) # prints last element of dead (wich is an array)

settings = PolyChordSettings(nDims, nDerived)
settings.file_root = functionName
settings.do_clustering = True
settings.read_resume = False

output = pypolychord.run_polychord(likelihood, nDims, nDerived, settings, prior, dumper)
paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
paramnames += [('r*', 'r')]
output.make_paramnames_files(paramnames)

try:
    import getdist.plots
    import matplotlib.pyplot as plt
    posterior = output.posterior
    g = getdist.plots.getSubplotPlotter()
    g.triangle_plot(posterior, filled=True)
    g.export(functionName + '.pdf')
except ImportError:
    print("Install matplotlib and getdist for plotting examples")
