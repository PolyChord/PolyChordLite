#!/usr/bin/python

import _PyPolyChord

nDims = 5
nDerived = 2
nlive = 500
num_repeats = nDims*5
do_clustering = False
feedback = 1
precision_criterion = 0.0001
max_ndead = -1
boost_posterior = 0.0
posteriors = False
equals = False
cluster_posteriors = False
write_resume = False
write_paramnames = False
read_resume = False
write_stats = False
write_live = True
write_dead = False
update_files = nlive
base_dir = 'chains'
file_root = 'PyPolyChord_test'

def gaussian(theta,phi):
    from numpy import pi,log
    sigma = 0.1
    nDims = len(theta)

    logL = -log(2*pi)*nDims/2.0

    for t in theta:
        logL += -(t-0.5)*(t-0.5)/2/sigma/sigma

    return long(logL)

_PyPolyChord.set_loglikelihood(gaussian)

_PyPolyChord.run(nDims, nDerived, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, base_dir, file_root)
