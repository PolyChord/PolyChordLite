from distutils.core import setup, Extension

pypolychord_module = Extension(
        name='_PyPolyChord',
        language="c++",
        include_dirs = ['src/'],
        library_dirs = ['src/'],
        libraries = ['chord','gfortran'],
        sources=['_PyPolyChord.cpp']
        )

setup(
    name = 'PyPolyChord',
    version = '1.0',
    description = 'Run PolyChord with Python',
    ext_modules=[pypolychord_module]
)
