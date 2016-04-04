import os
from numpy.distutils.core import setup, Extension


base = os.path.dirname(os.path.abspath(__file__))


pypolychord_module = Extension(
        name='_PyPolyChord',
        language="c++",
        include_dirs = [base + '/src/'],
        library_dirs = [base + '/src/'],
        libraries = ['chord','gfortran'],
        sources=[base + '/_PyPolyChord.cpp']
        )

setup(
    name = 'PyPolyChord',
    version = '1.0',
    description = 'Run PolyChord with Python',
    ext_modules=[pypolychord_module]
)
