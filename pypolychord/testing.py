import numpy
from run import run_polychord

def gaussian(theta, nderived, ndims):
    """ A simple gaussian with mean mu and variance sigma """
    mu, sigma = 0.1, 0.5
    loglike = 0
    for param in theta:
        loglike += -(param-mu)*(param-mu)/2/sigma/sigma
    phi = numpy.zeros(nderived)
    return loglike, phi

run_polychord(gaussian, 20, file_root='gaussian')
