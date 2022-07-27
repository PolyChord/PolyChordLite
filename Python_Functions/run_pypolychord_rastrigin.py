from numpy import pi, log, sqrt, cos
import pypolychord
from pypolychord.priors import UniformPrior

nDims = 3
nDerived = 1

functionName = 'rastrigin'

# param : array
def likelihood(theta):
    """ Rastrigin Likelihood"""

    A = 10
    TwoPi = 2*pi
    r2 = 0

    logL =  -sum( log(4991.21750) + theta**2 - A*cos(TwoPi*theta) )


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
