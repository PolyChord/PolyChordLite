"""
Python interface to PolyChord

Polychord is a tool to solve high dimensional problems.
"""

from setuptools import setup, Extension, find_packages, Distribution
from setuptools.command.build_py import build_py as _build_py
from distutils.command.clean import clean as _clean

import os, sys, subprocess, shutil

import numpy

BASE_PATH = os.path.dirname(os.path.abspath(__file__))

NAME = 'pypolychord'
DOCLINES = (__doc__ or '').split("\n")

def readme():
    with open('pypolychord_README.rst') as f:
        return f.read()


def get_version(short=False):
    with open('src/polychord/feedback.f90') as f:
        for line in f:
            if 'version' in line:
                return line.split(': ')[1].split('"')[0]

class DistributionWithMPIOption(Distribution):
    def __init__(self, *args, **kwargs):
        self.global_options = self.global_options + [("no-mpi", None, "Don't compile with MPI support.")]
        self.no_mpi = None
        super().__init__(*args, **kwargs)

class CustomBuildPy(_build_py):
    def run(self):
        if self.distribution.no_mpi is None:
            env = {"MPI"   : "1",
                   "PATH"  : os.environ["PATH"],
                   "CC"    : "mpicc",
                   "CXX"   : "mpicxx",
                   "FC"    : "mpifort",}
            # These need to be set so that build_ext uses the right compilers
            os.environ["CC"] = "mpicc"
            os.environ["CXX"] = "mpicxx"
        else:
            env = {"PATH"  : os.environ["PATH"],
                   "CC"    : os.environ["CC"] if "CC" in os.environ else "gcc",
                   "CXX"   : os.environ["CXX"] if "CXX" in os.environ else "g++",
                   "FC"    : os.environ["FC"] if "FC" in os.environ else "gfortran",}

        if sys.platform == "darwin":
            os.environ['MACOSX_DEPLOYMENT_TARGET'] = "10.9"
        
        env["PWD"] = BASE_PATH
        subprocess.run(["make", "libchord.a"], check=True, env=env, cwd=BASE_PATH)
        self.run_command("build_ext")
        return super().run()

class CustomClean(_clean):
    def run(self):
        subprocess.run(["make", "veryclean"], check=True, env=os.environ)
        return super().run()

def get_gfortran_libdir():
    r = subprocess.run(f"gfortran -v", shell=True, env=os.environ, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    gfortran_info = r.stderr.decode("utf-8").split(" ")
    for arg in gfortran_info:
        if "--libdir" in arg:
            libdir = arg[arg.find("=")+1:]
            return libdir
    print("Could not find gfortran library.")
    return ""

if "--no-mpi" in sys.argv:
    NAME += '_nompi'
    DOCLINES[1] = DOCLINES[1] + ' (cannot be used with MPI)'

pypolychord_module = Extension(
        name=f'_pypolychord',
        library_dirs=[os.path.join(BASE_PATH, 'lib'), get_gfortran_libdir()],
        include_dirs=[os.path.join(BASE_PATH, 'src/polychord'),
                      numpy.get_include()],
        libraries=['chord', 'gfortran'],
        # extra_link_args=["-Wl,-rpath,$ORIGIN/lib"],
        # extra_compile_args=["-Wl,-rpath,$ORIGIN/lib"],
        extra_link_args=["-lc++"],
        extra_compile_args=["-std=c++11", "-stdlib=libc++", "-g","-O0"],
        sources=['pypolychord/_pypolychord.cpp']
        )

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
      distclass=DistributionWithMPIOption,
      ext_modules=[pypolychord_module],
      cmdclass={'build_py' : CustomBuildPy,
                'clean' : CustomClean},
    #   package_data={"" : ["lib/libchord.so"]},
    #   include_package_data=True,
      zip_safe=False)
