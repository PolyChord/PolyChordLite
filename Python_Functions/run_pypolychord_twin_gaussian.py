from numpy import pi, log, sqrt, logaddexp, exp
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior

nDims = 3
nDerived = 1

functionName = 'twin_gaussian'

# param : array
def likelihood(theta):
    """ Twin_gaussian Likelihood"""

    sigma = 0.1
    logL = 0

    mu1    = [0] * len(theta)
    mu1[0] = -0.5
    mu1[1] = -0.5

    mu2    = [0] * len(theta)
    mu2[0] = 0.5
    mu2[1] = 0.5
   
    logL1  = - log( sigma ) + log(2*pi)/2 
    logL2  = logL1

    logL1 = logL1 - sum( ( ( theta - mu1 ) / sigma ) ** 2 ) / 2
    logL2 = logL2 - sum( ( ( theta - mu2 ) / sigma ) ** 2 ) / 2

    # def logaddexp(loga, logb):
    #     result = 0
    #     if (loga>logb):
    #         result = loga + log(exp(logb-loga) +1)
    #     else:
    #         result = logb + log(exp(loga-logb) +1)
    #     return result
    
    logL =  logaddexp(logL1, logL2) - log(2)

    if(theta[0] >0.5):
        r2=1
    else:
        r2=-1


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
