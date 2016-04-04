#!/usr/bin/python

import _PyPolyChord

nDims = 5
nDerived = 1
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
write_live = False
write_dead = False
update_files = nlive

_PyPolyChord.run(nDims, nDerived, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files)
