from numpy import pi, log, sqrt
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior

nDims = 3
nDerived = 1

functionName = 'gaussian'

# param : array
def likelihood(theta):
    """ Simple Gaussian Likelihood"""

    sigma = 0.1
    nDims = len(theta)

    r2 = sum(theta**2)

    logL = -log(2*pi*sigma*sigma)*nDims/2.0
    logL += -r2/2/sigma/sigma

    return logL, [r2] # float, array-like

# param : array
def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """
    return UniformPrior(-1, 1)(hypercube) # array

# param : array, array, array, float, float
def dumper(live, dead, logweights, logZ, logZerr):
    print("Last dead point:", dead[-1]) # prints last element of dead (wich is an array)

settings = PolyChordSettings(nDims, nDerived) #settings is an object
settings.file_root = functionName #string
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
