from numpy import pi,log,sqrt
import matplotlib.pyplot as plt
import getdist.plots
import PyPolyChord.PyPolyChord as PolyChord
from PyPolyChord.settings import PolyChordSettings

nDims = 5
nDerived = 0

def likelihood(theta):

    phi = [0.0] *nDerived

    sigma = 0.1
    nDims = len(theta)

    r2 = 0.0
    for t in theta:
        r2 += (t-0.5)*(t-0.5)

    logL = -log(2*pi*sigma*sigma)*nDims/2.0
    logL += -r2/2/sigma/sigma

    return logL, phi

PolyChord.mpi_notification()

settings = PolyChordSettings(nDims, nDerived)
settings.file_root = 'gaussian'
settings.do_clustering = True

output = PolyChord.run_polychord(likelihood, nDims, nDerived, settings)

posterior = output.posterior
#posterior1 = output.cluster_posterior(1)
#posterior2 = output.cluster_posterior(2)
#
g = getdist.plots.getSubplotPlotter()
g.triangle_plot(posterior,filled=True)
plt.show()
