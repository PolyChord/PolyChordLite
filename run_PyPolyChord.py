import matplotlib.pyplot as plt
import getdist.plots
from numpy import pi,log,sqrt,exp
import PyPolyChord.PyPolyChord as PolyChord
from PyPolyChord.settings import PolyChordSettings

def gaussian(theta, mu, sigma):

    nDims = len(theta)
    r2 = sum([(t-m)**2 for t, m in zip(theta, mu)])

    logL = -log(2*pi*sigma**2)*nDims/2.0
    logL -= r2/2/sigma**2

    return logL

nDims = 5
nDerived = 0

def likelihood(theta):

    phi = [0.0] * nDerived

    sigma = 0.1
    mu0 = [0.5] * nDims
    mu0[0] = 0.5 - 3*sigma

    mu1 = [0.5] * nDims
    mu1[0] = 0.5 + 3*sigma

    logL = log((exp(gaussian(theta, mu0, sigma)) + exp(gaussian(theta, mu1, sigma)))/2)

    return logL, phi

def prior(x):

    theta = x

    return theta

settings = PolyChordSettings(nDims, nDerived)
settings.file_root = 'gaussian'
settings.do_clustering = True

PolyChord.mpi_notification()
output = PolyChord.run_polychord(likelihood, nDims, nDerived, settings, prior)
posterior = output.posterior
posterior1 = output.cluster_posterior(1)
posterior2 = output.cluster_posterior(2)

g = getdist.plots.getSubplotPlotter()
g.triangle_plot([posterior,posterior1, posterior2],filled=True)
plt.show()
