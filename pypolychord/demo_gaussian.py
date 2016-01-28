import numpy
from run import run_polychord

def gaussian(theta, nderived, ndims):
    """ A simple gaussian with mean mu and variance sigma """
    mu, sigma = 0.0, 0.1

    loglike = - numpy.log(2*numpy.pi*sigma*sigma)/2*ndims

    for param in theta:
        loglike += -(param-mu)*(param-mu)/2/sigma/sigma

    phi = numpy.zeros(nderived)

    return loglike, phi

def uniform_prior(cube):
    theta = 2*(cube-0.5)
    return theta

run_polychord(gaussian, 2, prior=uniform_prior,  file_root='gaussian')
