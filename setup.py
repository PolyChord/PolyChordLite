"""
Python interface to PolyChord

Polychord is a tool to solve high dimensional problems.
"""

from setuptools import setup, Extension, find_packages
import os, sys
import numpy

NAME = 'pypolychord'
DOCLINES = (__doc__ or '').split("\n")

def readme():
    with open('pypolychord_README.rst') as f:
        return f.read()


def get_version(short=False):
    with open('pypolychord/src/feedback.f90') as f:
        for line in f:
            if 'version' in line:
                return line.split(': ')[1].split('"')[0]

pypolychord_module = Extension(
        name='_pypolychord',
        library_dirs=[os.path.join(os.getcwd(), 'pypolychord/lib')],
        include_dirs=[os.path.join(os.getcwd(), 'pypolychord/include/'),
                      numpy.get_include()],
        runtime_library_dirs=[os.path.join(os.getcwd(), 'pypolychord/lib')],
        libraries=['chord', 'gfortran'],
        sources=['pypolychord/_pypolychord.cpp']
        )

if ('mpi' in sys.argv):
    if ('mpi_usempi' not in pypolychord_module.libraries):
        pypolychord_module.libraries.append('mpi_usempi')
        pypolychord_module.libraries.append('mpi_mpifh')
    sys.argv.remove('mpi')
else:
    NAME += '_nompi'
    DOCLINES[1] = DOCLINES[1] + ' (cannot be used with MPI)'

setup(name=NAME,
      version=get_version(),
      description=DOCLINES[1],
      long_description = "\n".join(DOCLINES[3:-1]),
      url='https://ccpforge.cse.rl.ac.uk/gf/project/polychord/',
      author='Will Handley',
      author_email='wh260@cam.ac.uk',
      license='PolyChord',
      packages=find_packages(),
      install_requires=['numpy', 'scipy'],
      extras_require={'plotting': 'getdist'},
      ext_modules=[pypolychord_module],
      zip_safe=False)
