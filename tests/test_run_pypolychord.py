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
default_settings.file_root = 'settings'
default_settings.nlive = 200
default_settings.read_resume = False
default_settings.feedback = 0

cube_samples_settings = PolyChordSettings(4, 1)
cube_samples_settings.file_root = 'cube_samples'
cube_samples_settings.nlive = 200
cube_samples_settings.read_resume = False
cube_samples_settings.cube_samples = np.array([[0.1, 0.2, 0.3, 0.4],
                                               [0.5, 0.6, 0.7, 0.8]])
cube_samples_settings.feedback = 0


@pytest.mark.parametrize("settings, likelihood, nDims, nDerived",
                         [(default_settings, gaussian_likelihood, 4, 1),
                          (cube_samples_settings, gaussian_likelihood, 4, 1)])
def test_run(settings, likelihood, nDims, nDerived):
    # Define a box uniform prior from -1 to 1
    def prior(hypercube):
        """ Uniform prior from [-1,1]^D. """
        return UniformPrior(-1, 1)(hypercube)

    # Optional dumper function giving run-time read access to
    # the live points, dead points, weights and evidences

    def dumper(live, dead, logweights, logZ, logZerr):
        print("Last dead point:", dead[-1])

    # Run PolyChord

    print("Running PolyChord")
    output = pypolychord.run_polychord(likelihood, nDims, nDerived, settings, prior, dumper)

    # Create a paramnames file

    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
    paramnames += [('r*', 'r')]
    output.make_paramnames_files(paramnames)


@pytest.mark.parametrize("settings, likelihood, nDims, nDerived",
                         [(default_settings, gaussian_likelihood, 4, 1),
                          (cube_samples_settings, gaussian_likelihood, 4, 1)])
@pytest.mark.mpi
def test_run_mpi(settings, likelihood, nDims, nDerived):
    from mpi4py import MPI

    settings.file_root += '_mpi'

    # Define a box uniform prior from -1 to 1
    def prior(hypercube):
        """ Uniform prior from [-1,1]^D. """
        return UniformPrior(-1, 1)(hypercube)

    # Optional dumper function giving run-time read access to
    # the live points, dead points, weights and evidences

    def dumper(live, dead, logweights, logZ, logZerr):
        print("Last dead point:", dead[-1])

    # Run PolyChord

    print("Running PolyChord")
    output = pypolychord.run_polychord(likelihood, nDims, nDerived, settings, prior, dumper)

    # Create a paramnames file

    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
    paramnames += [('r*', 'r')]
    output.make_paramnames_files(paramnames)
