from numpy import pi,log,sqrt
import PyPolyChord.PyPolyChord as PolyChord

nDims = 5
nDerived = 0

def gaussian(theta):

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
PolyChord.run_nested_sampling(gaussian, nDims, nDerived, file_root='gaussian')
