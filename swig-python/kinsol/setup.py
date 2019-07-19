#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension

import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()
includes = [numpy_include, '/home/balos1/Workspace/SUNDIALS/sundials-install-pykinsol/include']

kinsol_module = Extension('_kinsol',
                          include_dirs=includes,
                          library_dirs=['/home/balos1/Workspace/SUNDIALS/sundials-install-pykinsol/lib64'],
                          runtime_library_dirs=['/home/balos1/Workspace/SUNDIALS/sundials-install-pykinsol/lib64'],
                          libraries=['sundials_kinsol'],
                          sources=['kinsol_wrap.c'])

setup (name='kinsol',
       version     = '0.1.0',
       author      = "Cody J. Balos @ LLNL",
       description = """KINSOL python interface""",
       ext_modules = [kinsol_module],
       py_modules  = ["kinsol"],
       )
