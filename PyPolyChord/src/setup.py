import os
from distutils.core import setup, Extension

lib_dir = os.environ['LIB_DIR']

pypolychord_module = Extension(
        name= '_PyPolyChord',
        include_dirs = ['src'],
        library_dirs = [lib_dir],
        libraries = ['chord','gfortran'],
        sources=['_PyPolyChord.c']
        )

setup(
    name = 'PyPolyChord',
    version = '1.0',
    description = 'Run PolyChord with Python',
    ext_modules=[pypolychord_module]
)
