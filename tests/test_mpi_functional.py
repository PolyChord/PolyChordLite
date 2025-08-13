#!/usr/bin/env python3
"""
Minimal functional test for custom MPI communicator support.
This actually runs PolyChord with different communicators.
"""

try:
    from mpi4py import MPI
    import pypolychord
    from pypolychord.priors import UniformPrior
    HAS_MPI = True
except ImportError:
    HAS_MPI = False
    exit(0)

def gaussian_likelihood(theta):
    """Simple Gaussian likelihood centered at origin."""
    # 2D Gaussian with unit variance
    logL = -0.5 * sum(theta**2)
    r2 = sum(theta**2)
    return logL, [r2]  # return log likelihood and derived parameter

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    print(f"Rank {rank}/{size}: Testing custom MPI communicator functionality")
    
    # Test 1: Default behavior (should use COMM_WORLD internally)
    try:
        output1 = pypolychord.run(
            gaussian_likelihood, 2,  # 2D Gaussian
            nlive=10, num_repeats=5, max_ndead=50,
            file_root=f'test_default_{rank}',
            write_resume=False, write_dead=False,
            write_stats=False,
            feedback=0, do_clustering=False
        )
        # run() returns None when anesthetic not installed, but that means it completed
        if output1 is None:
            print(f"Rank {rank}: Default comm test - SUCCESS (completed, anesthetic not available)")
        else:
            print(f"Rank {rank}: Default comm test - SUCCESS (logZ ≈ {output1.logZ:.2f})")
    except Exception as e:
        print(f"Rank {rank}: Default comm test - FAILED: {e}")
    
    # Test 2: Explicit COMM_WORLD
    try:
        output2 = pypolychord.run(
            gaussian_likelihood, 2,  # 2D Gaussian
            nlive=10, num_repeats=5, max_ndead=50,
            file_root=f'test_explicit_{rank}',
            write_resume=False, write_dead=False,
            write_stats=False,
            feedback=0, do_clustering=False,
            comm=comm  # Explicit COMM_WORLD
        )
        if output2 is None:
            print(f"Rank {rank}: Explicit COMM_WORLD test - SUCCESS (completed, anesthetic not available)")
        else:
            print(f"Rank {rank}: Explicit COMM_WORLD test - SUCCESS (logZ ≈ {output2.logZ:.2f})")
    except Exception as e:
        print(f"Rank {rank}: Explicit COMM_WORLD test - FAILED: {e}")
    
    # Test 3: Duplicated communicator
    try:
        dup_comm = comm.Dup()
        output3 = pypolychord.run(
            gaussian_likelihood, 2,  # 2D Gaussian
            nlive=10, num_repeats=5, max_ndead=50,
            file_root=f'test_dup_{rank}',
            write_resume=False, write_dead=False,
            write_stats=False,
            feedback=0, do_clustering=False,
            comm=dup_comm  # Custom duplicated communicator
        )
        dup_comm.Free()
        if output3 is None:
            print(f"Rank {rank}: Duplicated communicator test - SUCCESS (completed, anesthetic not available)")
        else:
            print(f"Rank {rank}: Duplicated communicator test - SUCCESS (logZ ≈ {output3.logZ:.2f})")
    except Exception as e:
        print(f"Rank {rank}: Duplicated communicator test - FAILED: {e}")
    
    print(f"Rank {rank}: All tests completed")
    
    # Synchronize before exit
    comm.Barrier()