import numpy
from run import run_polychord

def gaussian(theta, nderived, ndims):
    """ A simple gaussian with mean mu and variance sigma """
    mu, sigma = 0.0, 0.1
    rsquared = 0
    for param in theta:
        rsquared += (param-mu)*(param-mu)

    # Normalisation
    loglike = - numpy.log(2*numpy.pi*sigma*sigma)/2*ndims
    # Gaussian
    loglike -= rsquared/2/sigma/sigma

    phi = numpy.zeros(nderived)
    phi[0] = numpy.sqrt(rsquared)

    return loglike, phi

def uniform_prior(cube):
    theta = 2*(cube-0.5)
    return theta

nDims = 5
nDerived = 1

run_polychord(gaussian, nDims, nDerived=nDerived, prior=uniform_prior,  file_root='gaussian')
