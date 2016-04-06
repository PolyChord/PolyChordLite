#!/usr/bin/python
from numpy import pi,log,sqrt
import PyPolyChord.PyPolyChord as PolyChord


def gaussian(theta,phi):
    sigma = 0.1
    nDims = len(theta)

    r2 = 0.0
    for t in theta:
        r2 += (t-0.5)*(t-0.5)

    logL = -log(2*pi*sigma*sigma)*nDims/2.0
    logL += -r2/2/sigma/sigma
    phi[0] = sqrt(r2)

    return logL

nDims = 5
nDerived = 1

PolyChord.run_nested_sampling(gaussian, nDims, nDerived, file_root='gaussian')
