#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension
import numpy
import os

# Path to the root of the SUNDIALS installation. SUNDIALS must have been
# build for double precision and 32-bit indices.
sundials_root = os.environ['SUNDIALS_ROOT']

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

includes = [numpy_include, f'{sundials_root}/include']

kinsol_module = Extension('_kinsol',
                          include_dirs=includes,
                          library_dirs=[f'{sundials_root}/lib64'],
                          runtime_library_dirs=[f'{sundials_root}/lib64'],
                          libraries=['sundials_kinsol'],
                          sources=['kinsol_wrap.cxx'])

setup(name='py-kinsol',
      version='0.1.0',
      author="Cody J. Balos @ LLNL",
      author_email="balos1@llnl.gov",
      url="computing.llnl.gov/projects/sundials",
      description="""KINSOL python interface""",
      ext_modules=[kinsol_module],
      py_modules=["kinsol"])
