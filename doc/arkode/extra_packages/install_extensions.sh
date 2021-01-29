#!/bin/bash
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
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

# Sphinx
echo "installing Sphinx"
pip install -U Sphinx==1.6.7

# Sphinx fortran domain
echo "installing Sphinx fortran domain"
pip install sphinx-fortran

# Sphinx bootstrap theme
echo "installing Sphinx bootstrap theme"
tar -zxf sphinx-bootstrap-theme.tgz
cd sphinx-bootstrap-theme
python ./setup.py install
cd -

echo "finished"
