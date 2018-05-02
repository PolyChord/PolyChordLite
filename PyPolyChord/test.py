import numpy
import _PyPolyChord

nDims = 1 
nDerived = 0
nlive = 1
num_repeats = 1
nprior = 2
do_clustering = False
feedback = -1
precision_criterion = 0
logzero = -1e30
max_ndead = 1
boost_posterior = 0
posteriors = False
equals = False
cluster_posteriors = False
write_resume = False
write_paramnames = False
read_resume = False
write_stats = False
write_live = False
write_dead = False
write_prior = False
compression_factor = 1.
base_dir = ''
file_root = ''
grade_frac = [1.0]
grade_dims = [nDims]
nlives = {}
seed = 0
def loglikelihood(theta):
    return 0., []

def prior(cube):
    return cube

try:
    _PyPolyChord.run(loglikelihood,
                     prior,
                     nDims,
                     nDerived,
                     nlive,
                     num_repeats,
                     nprior,
                     do_clustering,
                     feedback,
                     precision_criterion,
                     logzero,
                     max_ndead,
                     boost_posterior,
                     posteriors,
                     equals,
                     cluster_posteriors,
                     write_resume,
                     write_paramnames,
                     read_resume,
                     write_stats,
                     write_live,
                     write_dead,
                     write_prior,
                     compression_factor,
                     base_dir,
                     file_root,
                     grade_frac,
                     grade_dims,
                     nlives,
                     seed)
except FileNotFoundError:
    pass

