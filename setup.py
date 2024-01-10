"""
Python interface to PolyChord

Polychord is a tool to solve high dimensional problems.
"""

from setuptools import setup, Extension, find_packages, Distribution
from setuptools.command.build_py import build_py as _build_py
from distutils.command.clean import clean as _clean

import os, sys, subprocess, shutil

import numpy

def check_compiler(default_CC="gcc"):
    """Checks what compiler is being used (clang, intel, or gcc)."""

    CC = default_CC if "CC" not in os.environ else os.environ["CC"]
    CC_version = subprocess.check_output([CC, "-v"], stderr=subprocess.STDOUT).decode("utf-8").lower()
    
    if "clang" in CC_version:
        CC_family = "clang"
    elif "icc" in CC_version:
        CC_family = "intel"
    elif "gcc" in CC_version:
        CC_family = "gcc"
    else:
        print("Warning: unrecognised compiler: {}".format(CC_version))
        CC_family = ""

    return CC_family


NAME = 'pypolychord'
DOCLINES = (__doc__ or '').split("\n")


# Deal with annoying differences between clang and the other compilers
CC_FAMILY = check_compiler()
CPPRUNTIMELIB_FLAG = []
RPATH_FLAG = []

if CC_FAMILY == "clang":
    CPPRUNTIMELIB_FLAG += ["-stdlib=libc++"]
    if sys.platform == "darwin":
        # macOS idiosyncrasies
        CPPRUNTIMELIB_FLAG += ["-mmacosx-version-min=10.9"]

if sys.platform != "darwin":
    # Set RPATH on Linux machines
    RPATH_FLAG += ["-Wl,-rpath,$ORIGIN/pypolychord/lib"]


def readme():
    with open('pypolychord_README.rst') as f:
        return f.read()


def get_version(short=False):
    with open('src/polychord/feedback.f90') as f:
        for line in f:
            if 'version' in line:
                return line[44:50]


class DistributionWithOption(Distribution, object):
    def __init__(self, *args, **kwargs):
        self.global_options = self.global_options \
                                + [("no-mpi", None, "Don't compile with MPI support."),
                                   ("debug-flags", None, "Compile in debug mode.")]
        self.no_mpi = None
        self.debug_flags = None
        super(DistributionWithOption, self).__init__(*args, **kwargs)

class CustomBuildPy(_build_py, object):
    def run(self):
        env = {}
        env["PATH"] = os.environ["PATH"]
        if self.distribution.no_mpi is None:
            env["MPI"] = "1"
            # These need to be set so that build_ext uses the right compilers
            cc_compiler = subprocess.check_output(["make", "print_CC"]).decode('utf-8').strip()
            os.environ["CC"] = cc_compiler

            cxx_compiler = subprocess.check_output(["make", "print_CXX"]).decode('utf-8').strip()
            os.environ["CXX"] = cxx_compiler
        else:
            env["MPI"] = "0"

        if self.distribution.debug_flags is not None:
            self.distribution.ext_modules[0].extra_compile_args += ["-g", "-O0"]
            env["DEBUG"] = "1"
        
        BASE_PATH = os.path.dirname(os.path.abspath(__file__))
        env["CURDIR"] = BASE_PATH
        env.update({k : os.environ[k] for k in ["CC", "CXX", "FC"] if k in os.environ})
        subprocess.check_call(["make", "-e", "libchord.so"], env=env, cwd=BASE_PATH)
        if not os.path.isdir("pypolychord/lib/"):
            os.makedirs(os.path.join(BASE_PATH, "pypolychord/lib/"))
        shutil.copy(os.path.join(BASE_PATH, "lib/libchord.so"), 
                    os.path.join(BASE_PATH, "pypolychord/lib/"))
        self.run_command("build_ext")
        return super(CustomBuildPy, self).run()

class CustomClean(_clean):
    def run(self):
        subprocess.run(["make", "veryclean"], check=True, env=os.environ)
        return super().run()

if "--no-mpi" in sys.argv:
    NAME += '_nompi'
    DOCLINES[1] = DOCLINES[1] + ' (cannot be used with MPI)'

pypolychord_module = Extension(
        name='_pypolychord',
        library_dirs=['lib'],
        include_dirs=['src/polychord', numpy.get_include()],
        libraries=['chord',],
        extra_link_args=RPATH_FLAG + CPPRUNTIMELIB_FLAG,
        extra_compile_args= ["-std=c++11"] + RPATH_FLAG + CPPRUNTIMELIB_FLAG,
        runtime_library_dirs=['lib'],
        sources=['pypolychord/_pypolychord.cpp']
        )

setup(name=NAME,
      version=get_version(),
      description='Python interface to PolyChord ' + get_version(),
      url='https://github.com/PolyChord/PolyChordLite',
      author='Will Handley',
      author_email='wh260@cam.ac.uk',
      license='PolyChord',
      packages=find_packages(),
      install_requires=['numpy','scipy'],
      extras_require={'plotting': 'getdist'},
      distclass=DistributionWithOption,
      ext_modules=[pypolychord_module],
      cmdclass={'build_py' : CustomBuildPy,
                'clean' : CustomClean},
      package_data={"" : ["lib/libchord.so"]},
      include_package_data=True,
      zip_safe=False)
