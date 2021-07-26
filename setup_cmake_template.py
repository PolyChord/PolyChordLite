"""
Template for setup.py that should be populated by CMake with:
configure_file(${SETUP_PY_IN} ${SETUP_PY})

Used template variables:
    - "target_python_so" : place where _pypolychord.ARCH.so is built (NOT install path)
    - "PACKAGE_VERSION" : version of the package
    - "CMAKE_CURRENT_SOURCE_DIR": root directory, where a standard setup.py would be placed
"""

import os
import shutil
import sys
from distutils.core import setup

from setuptools import Extension
from setuptools.command.build_ext import build_ext


class MyBuildExtension(build_ext):
    """This class 'builds' the extension module, by
    copying it from the place where CMake placed it.
    """
    def build_extension(self, ext):

        # _pypolychord.ARCH.so
        if os.path.exists("${target_python_so}"):
            shutil.copyfile("${target_python_so}", self.get_ext_fullpath(ext.name))
        elif sys.argv[2] == "install":
            # let's warn here, though this should not happen with the current CMake setup
            print("NOT FOUND: ${target_python_so}\nYour installation may be incomplete.")


setup(name='pypolychord',
      version='${PACKAGE_VERSION}',
      package_dir={'pypolychord': '${CMAKE_CURRENT_SOURCE_DIR}/pypolychord'},
      packages=['pypolychord'],
      cmdclass={'build_ext': MyBuildExtension},
      ext_modules=[Extension('_pypolychord', [])],
      )
