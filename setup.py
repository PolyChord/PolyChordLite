from setuptools import setup, Extension, find_packages
import os
import numpy
from setuptools.dist import Distribution


def readme():
    with open('pypolychord_README.rst') as f:
        return f.read()


def get_version(short=False):
    with open('pypolychord/src/feedback.f90') as f:
        for line in f:
            if 'version' in line:
                return line[44:48]

class BinaryDistribution(Distribution):
   def is_pure(self):
       return False


pypolychord_module = Extension(
        name='_pypolychord',
        library_dirs=[os.path.join(os.getcwd(), 'pypolychord/lib')],
        include_dirs=[os.path.join(os.getcwd(), 'pypolychord/include/'),
                      numpy.get_include()],
        runtime_library_dirs=[os.path.join(os.getcwd(), 'pypolychord/lib')],
        libraries=['chord'],
        sources=['pypolychord/_pypolychord.cpp']
        )

setup(name='pypolychord',
      version=get_version(),
      description='Python interface to PolyChord ' + get_version(),
      url='https://ccpforge.cse.rl.ac.uk/gf/project/polychord/',
      author='Will Handley',
      author_email='wh260@cam.ac.uk',
      license='PolyChord',
      packages=find_packages(),
      install_requires=[],
      extras_require={'plotting': 'getdist'},
      ext_modules=[pypolychord_module],
      zip_safe=False,
      include_package_data=True,
      package_data={
        # If any package contains these files, include them:
        '': ['*.o', '*.so']
      },
      distclass=BinaryDistribution)
