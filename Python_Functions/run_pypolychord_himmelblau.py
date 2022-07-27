from numpy import pi, log, sqrt
import pypolychord
from pypolychord.priors import UniformPrior

nDims = 3
nDerived = 1

functionName = 'himmelblau'

# param : array
def likelihood(theta):
    """ Himmelblau Likelihood"""

    logL = -log(0.4071069421432255) 

    logL =  logL-(theta[0]**2 + theta[1]- 11)**2 - (theta[0]+theta[1]**2-7)**2 
    r2=0
    return logL, [r2] # float, array-like

# param : array
def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """
    return UniformPrior(-1, 1)(hypercube) # array

# param : array, array, array, float, float
def dumper(live, dead, logweights, logZ, logZerr):
    print("Last dead point:", dead[-1]) # prints last element of dead (wich is an array)

paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
paramnames += [('r*', 'r')]

kwargs = {
    "file_root": functionName, #string
    "do_clustering": True,
    "read_resume": False,
    "paramnames": paramnames,
}

output = pypolychord.run_polychord(likelihood, nDims, nDerived, prior, dumper, **kwargs)

try:
    import getdist.plots
    import matplotlib.pyplot as plt
    posterior = output.posterior
    g = getdist.plots.getSubplotPlotter()
    g.triangle_plot(posterior, filled=True)
    g.export(functionName + '.pdf')
except ImportError:
    print("Install matplotlib and getdist for plotting examples")
