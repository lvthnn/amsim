import os
import sys
import subprocess
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        # Get pybind11 cmake directory
        import pybind11
        pybind11_dir = pybind11.get_cmake_dir()
        
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-Dpybind11_DIR={pybind11_dir}",
            "-DBUILD_PYTHON=ON",
        ]

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        # Configure
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)

        # Build
        subprocess.check_call(["cmake", "--build", "."], cwd=self.build_temp)

setup(
    ext_modules=[CMakeExtension("amsim._core", sourcedir="..")],
    cmdclass={"build_ext": CMakeBuild},
)
