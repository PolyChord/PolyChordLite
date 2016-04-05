import mpi4py
from mpi4py import MPI

from distutils.core import setup, Extension

print MPI.COMM_WORLD.Get_size()

pypolychord_module = Extension(
        name='_PyPolyChord',
        include_dirs = ['./src/'],
        library_dirs = ['./src/'],
        libraries = ['chord','gfortran'],
        sources=['_PyPolyChord.c']
        )

setup(
    name = 'PyPolyChord',
    version = '1.0',
    description = 'Run PolyChord with Python',
    ext_modules=[pypolychord_module]
)
