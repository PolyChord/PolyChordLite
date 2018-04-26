from setuptools import setup, Extension
import os

def readme():
    with open('PyPolyChord_README.rst') as f:
        return f.read()

pypolychord_module = Extension(
        name= '_PyPolyChord',
        library_dirs = [os.path.join(os.getcwd(),'lib')],
        runtime_library_dirs = [os.path.join(os.getcwd(),'lib')],
        libraries = ['chord'],
        sources=[os.path.join(os.getcwd(),'PyPolyChord/_PyPolyChord.cpp')]
        )

setup(name='PyPolyChord',
      version='1.14',
      description='Python interface to PolyChord 1.14',
      url='https://ccpforge.cse.rl.ac.uk/gf/project/polychord/',
      author='Will Handley',
      author_email='wh260@cam.ac.uk',
      license='PolyChord',
      packages=['PyPolyChord'],
      install_requires=['numpy'],
      extras_require={'plotting':'getdist'},
      ext_modules=[pypolychord_module],
      zip_safe=False)
