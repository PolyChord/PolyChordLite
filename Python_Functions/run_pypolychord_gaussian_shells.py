from numpy import pi, log, sqrt, logaddexp, exp, zeros
from scipy.special import loggamma
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior

nDims = 3
nDerived = 1

functionName = 'gaussian_shells'

# param : array
def likelihood(theta):
    """ Simple Gaussian_shells Likelihood"""
    sigma = 0.1
    radius = 2
    A = -1
    logL ,logL_temp, r0, sigma0, logf0  = 0, 0, 0, 0, 0
    mu = zeros(len(theta))
    logsqrttwopi = log(sqrt(2*pi))
    logpi = log(pi)
    ndims = len(theta)
    i = 0

    if (A==-1):
        r0      = (radius + sqrt(radius**2 + 4*(ndims-1)*sigma**2))/2
        logf0   = -(radius - r0)**2/2/sigma**2 + (ndims-1)*log(r0)
        logf0  += log(ndims+0) + ndims/2*logpi -loggamma(1+ndims/2)
        sigma0  = sigma*sqrt((1+radius/sqrt(radius**2 + 4*(ndims-1)*sigma**2))/2)
        A       = logf0 + logsqrttwopi + log(sigma0)
    
    mu[0] = 3.5        

    #Gaussian normalisation
    logL_temp = -A - ( (sqrt( sum( (mu-theta)**2 ) ) -radius)**2 /(2*sigma*sigma) )
    logL =logL_temp

    mu[1] = 3.5
    logL_temp = -A - ( (sqrt( sum( (mu-theta)**2 ) ) -radius)**2 /(2*sigma*sigma) )

    r2 = 0


    def logincexp(loga,logb,logc=None):
        if (loga>logb):
            loga = loga + log(exp(logb-loga) + 1)
        else:
            loga = logb + log(exp(loga-logb) + 1)
        

        if(logc is not None):
            if (loga>logc):
                loga = loga + log(exp(logc-loga) + 1)
            else:
                loga = logc + log(exp(loga-logc) + 1)
        
        return loga, logb, logc
            
    logL, logL_temp, _ = logincexp(logL, logL_temp)
    logL = logL - log(2)
   
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
