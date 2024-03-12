import pytest
import numpy as np
from mpi4py import MPI
import anesthetic as ac
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


uniform_prior = UniformPrior(-1, 1)

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
def test_run_polychord(settings, likelihood, nDims, nDerived):
    # Optional dumper function giving run-time read access to
    # the live points, dead points, weights and evidences

    def dumper(live, dead, logweights, logZ, logZerr):
        print("Last dead point:", dead[-1])

    # Run PolyChord

    output = pypolychord.run_polychord(likelihood, nDims, nDerived,
                                       settings, uniform_prior, dumper)

    # Create a paramnames file

    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
    paramnames += [('r*', 'r')]
    output.make_paramnames_files(paramnames)


@pytest.mark.parametrize("likelihood, nDims, nDerived",
                         [(gaussian_likelihood, 4, 1),])
def test_run(likelihood, nDims, nDerived):
    # Create paramnames
    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
    paramnames += [('r*', 'r')]

    # Run PolyChord

    ns = pypolychord.run(likelihood, nDims, nDerived=nDerived,
                         prior=uniform_prior, paramnames=paramnames,
                         read_resume=False)
    assert isinstance(ns, ac.NestedSamples)


@pytest.mark.parametrize("seed", [-1, 0, 1, 2])
def test_seed(seed):
    # Create paramnames
    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(4)]
    paramnames += [('r*', 'r')]

    # Run PolyChord twice
    ns0 = pypolychord.run(gaussian_likelihood, 4, nDerived=1,
                          prior=uniform_prior, paramnames=paramnames,
                          read_resume=False, seed=seed)
    ns1 = pypolychord.run(gaussian_likelihood, 4, nDerived=1,
                          prior=uniform_prior, paramnames=paramnames,
                          read_resume=False, seed=seed)
    assert ns0.equals(ns1) != (seed < 0)


@pytest.mark.filterwarnings("ignore::pandas.errors.PerformanceWarning")
def test_no_derived():
    def no_derived_gaussian_likelihood(theta):
        """ Simple Gaussian Likelihood without derived parameters"""

        sigma = 0.1

        nDims = len(theta)
        r2 = sum(theta**2)
        logL = -np.log(2*np.pi*sigma*sigma)*nDims/2.0
        logL += -r2/2/sigma/sigma

        return logL

    seed = 1
    paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(4)]

    # Run without derived parameters
    ns0 = pypolychord.run(no_derived_gaussian_likelihood, 4,
                          prior=uniform_prior, paramnames=paramnames,
                          read_resume=False, seed=seed)
    # Run with derived parameters
    paramnames += [('r*', 'r')]
    ns1 = pypolychord.run(gaussian_likelihood, 4, nDerived=1,
                          prior=uniform_prior, paramnames=paramnames,
                          read_resume=False, seed=seed)
    assert ns0.equals(ns1.drop(columns='r'))


def test_grade_dims():
    grade_dims = [1, 3]
    pypolychord.run(gaussian_likelihood, 4, nDerived=1,
                    prior=uniform_prior, read_resume=False,
                    grade_dims=grade_dims)
    with pytest.raises(ValueError):
        pypolychord.run(gaussian_likelihood, 5, nDerived=1,
                        prior=uniform_prior, read_resume=False,
                        grade_dims=grade_dims)
