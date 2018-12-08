#!/bin/bash
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ----------------------------------------------------------------

# Sphinx
echo "installing Sphinx"
pip install -U Sphinx==1.3.1

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
