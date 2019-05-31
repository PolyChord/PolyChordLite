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

class DistributionWithOption(Distribution):
    def __init__(self, *args, **kwargs):
        self.global_options = self.global_options \
                                + [("no-mpi", None, "Don't compile with MPI support."),
                                   ("debug-flags", None, "Compile in debug mode.")]
        self.no_mpi = None
        self.debug_flags = None
        super().__init__(*args, **kwargs)

class CustomBuildPy(_build_py):
    def run(self):
        env = {}
        env["PATH"] = os.environ["PATH"]
        env.update({k : os.environ[k] for k in ["CC", "CXX", "FC"] if k in os.environ})
        if self.distribution.no_mpi is None:
            env["MPI"] = "1"
            # These need to be set so that build_ext uses the right compilers
            if "CC" not in os.environ: os.environ["CC"] = "mpicc"
            if "CXX" not in os.environ: os.environ["CXX"] = "mpicxx"
        else:
            env["MPI"] = "0"

        if self.distribution.debug_flags is not None:
            env["DEBUG"] = "1"

        if sys.platform == "darwin":
            os.environ['MACOSX_DEPLOYMENT_TARGET'] = "10.9"
        
        env["PWD"] = BASE_PATH
        print(env)
        subprocess.run(["make", "libchord.so"], check=True, env=env, cwd=BASE_PATH)
        os.makedirs(os.path.join(BASE_PATH, "pypolychord/lib/"), exist_ok=True)
        shutil.copy(os.path.join(BASE_PATH, "lib/libchord.so"), 
                    os.path.join(BASE_PATH, "pypolychord/lib/"))
        self.run_command("build_ext")
        return super().run()

class CustomClean(_clean):
    def run(self):
        subprocess.run(["make", "veryclean"], check=True, env=os.environ)
        return super().run()

if "--no-mpi" in sys.argv:
    NAME += '_nompi'
    DOCLINES[1] = DOCLINES[1] + ' (cannot be used with MPI)'

pypolychord_module = Extension(
        name='_pypolychord',
        library_dirs=[os.path.join(BASE_PATH, 'lib'),],
        include_dirs=[os.path.join(BASE_PATH, 'src/polychord'),
                      numpy.get_include()],
        libraries=['chord',],
        extra_link_args=["-Wl,-rpath,$ORIGIN/pypolychord/lib"],
        extra_compile_args=["-Wl,-rpath,$ORIGIN/pypolychord/lib", "-std=c++11"],
        runtime_library_dirs=[os.path.join(BASE_PATH, 'lib')],
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
      install_requires=['numpy',],
      extras_require={'plotting': 'getdist'},
      distclass=DistributionWithOption,
      ext_modules=[pypolychord_module],
      cmdclass={'build_py' : CustomBuildPy,
                'clean' : CustomClean},
      package_data={"" : ["lib/libchord.so"]},
      include_package_data=True,
      zip_safe=False)
