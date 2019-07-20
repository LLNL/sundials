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
# ---------------------------------------------------------------

SWIG=swig
PYTHON=python

KINSOL=kinsol.i
INCLUDES=-I../../include

.PHONY: kinsol_swig
all: kinsol_swig kinsol

kinsol_swig: $(KINSOL:%:%.i)
	@for i in ${KINSOL} ; do \
		echo "${SWIG} -python -c++ ${INCLUDES} $${i}"; \
		${SWIG} -python -c++ ${INCLUDES} $${i}; \
	done

kinsol: kinsol.py kinsol_wrap.cxx
	python setup.py build_ext --inplace

clean:
	rm *.c
	rm *.so
	rm -rf build
	rm -rf __pycache__