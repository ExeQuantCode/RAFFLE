import os
import sys
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import shutil

class CMakeBuild(build_ext):
    """
    Custom build command that uses CMake to build extensions.

    This class extends the `build_ext` command provided by setuptools
    to build C/C++ extensions using CMake.

    Attributes:
        build_temp (str): The directory where the build files will be placed.
        extensions (list): List of extension objects to build.

    Methods:
        run(): Runs the build process.
        build_extension(ext): Builds a specific extension.

    """

    def run(self):
        # Ensure CMake is installed
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        """
        Builds a specific extension using CMake.

        Args:
            ext (Extension): The extension object to build.

        """
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + os.path.abspath(sys.executable)
        ]
        build_args = ['--config', 'Release']

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

        extdir = extdir + "/raffle"
        if not os.path.exists(extdir):
            os.makedirs(extdir)

        # Move the generated .so file to the appropriate location
        for root, _, files in os.walk(self.build_temp):
            for file in files:
                if file.endswith('.so') or file.endswith('.py'):
                    src = os.path.join(root, file)
                    dest = os.path.join(extdir, file) 
                    print(f"Moving {src} to {dest}")
                    shutil.move(src, dest)
                    print(f"Moved {src} to {dest}")


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

minimum_requirements = [
    "f90wrap>=0.2.14,<0.2.15",
    "numpy>=1.26,<2.0.0",
    "ase>=3.23.0",
]

setup(
    name='raffle',
    version='0.3.1',
    author='Ned Thaddeus Taylor',
    author_email='n.t.taylor@exeter.ac.uk',
    description='A Python project with a Fortran library',
    install_requires=minimum_requirements,
    python_requires='>=3.11, <3.12',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    ext_modules=[CMakeExtension('raffle')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    include_package_data=True,
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    package_data={
        'src/raffle': ['*.py'],
        'src': ['*.f90'],
        'src/lib': ['*.f90'],
        'src/wrapper': ['*.f90'],
    },
)
