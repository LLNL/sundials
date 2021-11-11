#!/bin/bash
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU and Cody J. Balos @ LLNL
# ----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ----------------------------------------------------------------

echo "installing Sphinx"

# Sphinx
pip install sphinx>=4.0.0

# Sphinx fortran domain
pip install sphinx-fortran

# Sphinx readthedocs theme
pip install sphinx_rtd_theme

# Sphinx bibtex extension
pip install sphinxcontrib.bibtex

echo "finished"
