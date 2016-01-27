import numpy
import pypolychord

def gaussian(theta,nderived,ndims):
    loglike = 0
    sigma = 0.1
    for t in theta:
        loglike+= -(t-0.5)*(t-0.5)
    loglike /= sigma*sigma
    loglike /= 2
    phi = numpy.zeros(nderived)
    return loglike, phi

def run_polychord(
        loglikelihood,
        nlive = 50,
        num_repeats = 20,
        do_clustering = True,
        feedback = 2,
        precision_criterion = 0.0001,
        max_ndead = -1,
        boost_posterior = 0.0,
        posteriors = True,
        equals = True,
        cluster_posteriors = True,
        write_resume = True,
        write_paramnames = True,
        read_resume = False,
        write_stats = True,
        write_live = True,
        write_dead = True,
        update_files = True,
        nDims = 20,
        nDerived = 1,
        grade_dims = [0],
        grade_frac = [1.0]
        ):

    pypolychord.interfaces_module.simple_interface(loglikelihood, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived, grade_dims, grade_frac)

run_polychord(gaussian)
