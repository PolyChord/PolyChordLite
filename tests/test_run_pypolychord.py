import pytest
import numpy as np
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior


def gaussian_likelihood(theta):
    """ Simple Gaussian Likelihood"""

    sigma = 0.1

    nDims = len(theta)
    r2 = sum(theta**2)
    logL = -np.log(2*np.pi*sigma*sigma)*nDims/2.0
    logL += -r2/2/sigma/sigma

    return logL, [r2]


default_settings = PolyChordSettings(4, 1)
default_settings.file_root = 'gaussian'
default_settings.nlive = 200
default_settings.do_clustering = True
default_settings.read_resume = False


@pytest.mark.parametrize("settings, likelihood, nDims, nDerived", [(default_settings, gaussian_likelihood, 4, 1)])
def test_run(settings, likelihood, nDims, nDerived):
    # Define a four-dimensional spherical gaussian likelihood,
    # width sigma=0.1, centered on the 0 with one derived parameter.
    # The derived parameter is the squared radius

    # Define a box uniform prior from -1 to 1

    def prior(hypercube):
        """ Uniform prior from [-1,1]^D. """
        return UniformPrior(-1, 1)(hypercube)

    # Optional dumper function giving run-time read access to
    # the live points, dead points, weights and evidences

    def dumper(live, dead, logweights, logZ, logZerr):
        print("Last dead point:", dead[-1])

    # Run PolyChord

    output = pypolychord.run_polychord(likelihood, nDims, nDerived, settings, prior, dumper)

    # Create a paramnames file

    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
    paramnames += [('r*', 'r')]
    output.make_paramnames_files(paramnames)
