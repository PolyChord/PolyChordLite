import pypolychord


def run_polychord(loglikelihood, nDims, 
                  nDerived=1,
                  nlive=50,
                  num_repeats=None,
                  do_clustering=True,
                  feedback=1,
                  precision_criterion=0.001,
                  max_ndead=-1,
                  boost_posterior=0.0,
                  posteriors=True,
                  equals=True,
                  cluster_posteriors=True,
                  write_resume=True,
                  write_paramnames=True,
                  read_resume=False,
                  write_stats=True,
                  write_live=True,
                  write_dead=True,
                  update_files=True,
                  file_root='test',
                  base_dir='chains',
                  grade_dims=[0],
                  grade_frac=[1.0]):
    """ Interface for running pypolychord. Default parameters may be adjusted as required """

    if num_repeats is None:
        num_repeats = nDims*5

    pypolychord.interfaces_module.simple_interface(loglikelihood, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived, base_dir, file_root, grade_dims, grade_frac)

