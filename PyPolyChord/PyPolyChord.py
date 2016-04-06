try:
    import _PyPolyChord
except ImportError as e:
    if e.message == 'No module named _PyPolyChord':
        print('')
        print('   Could not load Python Extension _PyPolyChord.so')
        print('')
        print('   You have to build it first with:')
        print('')
        print('   $   make PyPolyChord')                  
        print('')
        print('   In the base PolyChord directory.')
        print('')

    elif e.message == 'libchord.so: cannot open shared object file: No such file or directory':
        print('')
        print('   Could not load PolyChord library "libchord.so"')
        print('')
        print('   You have to build it first,')
        print('   and point the LD_LIBRARY_PATH environment variable to it:')
        print('')
        print('   /-- BASH: --------------------------------------------\\')
        print('   |                                                     |')
        print('   |$   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib |')
        print('   |                                                     |')
        print('   +-- CSH: ---------------------------------------------+')
        print('   |                                                     |')
        print('   |$   setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH:$PWD/lib |')
        print('   |                                                     |')
        print('   \\-----------------------------------------------------/')
        print('')
    else:
        raise e

    import sys
    sys.exit(1)

def notification():
    print('/.===============================================================================.\\')
    print('||                                                                               ||')
    print('||    IMPORTANT: PyPolyChord settings                                            ||')
    print('||    -------------------------------                                            ||')
    print('||                                                                               ||')
    print('||    If you get MPI errors, try setting LD_PRELOAD to your mpi installation     ||')
    print('||                                                                               ||')
    print('||   /-- BASH: --------------------------------------------------------\\         ||')
    print('||   |                                                                 |         ||')
    print('||   |$   export LD_PRELOAD=/usr/lib/openmpi/lib/libmpi.so:$LD_PRELOAD |         ||')
    print('||   |                                                                 |         ||')
    print('||   +-- CSH: ---------------------------------------------------------+         ||')
    print('||   |                                                                 |         ||')
    print('||   |$   setenv LD_PRELOAD /usr/lib/openmpi/lib/libmpi.so:$LD_PRELOAD |         ||')
    print('||   |                                                                 |         ||')
    print('||   \\--------------------------------------------------------------- /          ||')
    print('||                                                                               ||')
    print('||   Where /usr/lib/openmpi/lib/libmpi.so should be replaced with the            ||')
    print('||   appropriate loctaion of libmpi.so on your system.                           ||')
    print('||                                                                               ||')
    print('\\.===============================================================================./  ')
    print('')

def run_nested_sampling(
        likelihood,
        nDims,
        nDerived,
        prior=None,
        nlive=None,
        num_repeats = None,
        do_clustering = True,
        feedback = 1,
        precision_criterion = 0.001,
        max_ndead = -1,
        boost_posterior = 0.0,
        posteriors = True,
        equals = True,
        cluster_posteriors = True,
        write_resume = True,
        write_paramnames = False,
        read_resume = True,
        write_stats = True,
        write_live = True,
        write_dead = True,
        update_files = None,
        base_dir = 'chains',
        file_root = 'test'):

    if prior is None:
        def prior(cube):
            theta = cube
            return theta
    if nlive is None:
        nlive = nDims*25
    if num_repeats is None:
        num_repeats = nDims*5
    if update_files is None:
        update_files = nlive

    notification()
    _PyPolyChord.run(likelihood, prior, nDims, nDerived, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, base_dir, file_root)

