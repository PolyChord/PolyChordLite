import os
import shutil
import sys
from distutils.core import setup
from pathlib import Path

from setuptools import Extension
from setuptools.command.build_ext import build_ext


class MyBuildExtension(build_ext):
    def build_extension(self, ext):
        print("MY current working directory:", os.getcwd())
        print("MY args:", sys.argv)

        # _pypolychord.ARCH.so
        if os.path.exists("${target_python_dir}/_pypolychord${python_so_suffix}"):
            print("CP: ${target_python_dir}/_pypolychord${python_so_suffix} ----> ", self.get_ext_fullpath(ext.name))
            shutil.copyfile("${target_python_dir}/_pypolychord${python_so_suffix}", self.get_ext_fullpath(ext.name))
        else:
            print("NOT FOUND: ${target_python_dir}/_pypolychord${python_so_suffix}")

        # lib dir
        if not os.path.isdir("${target_lib_dir}/pypolychord/lib/"):
            print("CREATED: ${target_lib_dir}/pypolychord/lib/")
            Path("${target_lib_dir}/pypolychord/lib/").mkdir(parents=True)

        # libchord.so
        if os.path.exists("${target_lib_dir}/libchord.so"):
            print("CP ${target_lib_dir}/libchord.so -->",
                  os.path.join(os.path.dirname(self.get_ext_fullpath(ext.name)), "pypolychord", "lib", "libchord.so"))
            shutil.copyfile("${target_lib_dir}/libchord.so",
                            os.path.join(os.path.dirname(self.get_ext_fullpath(ext.name)), "pypolychord", "lib", "libchord.so"))

        else:
            print("NOT Found: ${target_lib_dir}/libchord.so")


setup(name='pypolychord',
      version='${PACKAGE_VERSION}',
      package_dir={'pypolychord': '${CMAKE_CURRENT_SOURCE_DIR}/pypolychord'},
      packages=['pypolychord'],
      package_data={'pypolychord': ["${target_lib_dir}/libchord.so"]},
      cmdclass={'build_ext': MyBuildExtension},
      ext_modules=[Extension('_pypolychord', [])],
      )
