# ---------------------------------------------------------------
# Programmer: Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Makefile to generate SUNDIALS python interfaces with SWIG
#
# Must define the environment variable SUNDIALS_ROOT to be the
# path to the sundials installation.
# ---------------------------------------------------------------

# SWIG is only needed to re-generate the interfaces, it is not needed to build them.
SWIG=swig

# The python 3 executable
PYTHON=python3

# Include directories
INCLUDES=-I${SUNDIALS_ROOT}/include

# KINSOL swig interface files
KINSOL=kinsol.i

.PHONY: generate install
all: generate kinsol

generate: $(KINSOL:%:%.i)
	@for i in ${KINSOL} ; do \
		echo "${SWIG} -python -c++ ${INCLUDES} $${i}"; \
		${SWIG} -python -c++ ${INCLUDES} $${i}; \
	done

kinsol: kinsol.py kinsol_wrap.cxx
	python setup.py build_ext

kinsol_inplace: kinsol.py kinsol_wrap.cxx
	python setup.py build_ext --inplace

install:
	python setup.py install

clean:
	rm -rf build
	rm -rf __pycache__
	rm *.cxx
	rm *.so