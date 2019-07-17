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

KINSOL=kinsol
GENERIC=sundials
NVECTOR=nvector_serial

INCLUDES=-I../include

.PHONY: generic nvector kinsol

all: generic nvector kinsol

SWIG=swig

kinsol: $(KINSOL:%:kinsol/%.i)
	@for i in ${KINSOL} ; do \
		echo "${SWIG} -python ${INCLUDES} kinsol/$${i}.i" ; \
		${SWIG} -python ${INCLUDES} kinsol/$${i}.i ; \
	done

generic: $(GENERIC:%:sundials/%.i)
	@for i in ${GENERIC} ; do \
		echo "${SWIG} -python ${INCLUDES} sundials/$${i}.i" ; \
		${SWIG} -python ${INCLUDES} sundials/$${i}.i ; \
	done

nvector: $(NVECTOR:%:nvector/%.i)
	@for i in ${NVECTOR} ; do \
		echo "${SWIG} -python ${INCLUDES} nvector/$${i}.i" ; \
		${SWIG} -python ${INCLUDES} nvector/$${i}.i ; \
	done

# sunmatrix: $(SUNMATRIX:%:sunmatrix/fsunmatrix_%_mod.i)
# 	@for i in ${SUNMATRIX} ; do \
# 		echo "swig/swig -fortran -outdir ../src/sunmatrix/$${i}/fmod -o ../src/sunmatrix/$${i}/fmod ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i" ; \
# 		swig/swig -fortran -outdir ../src/sunmatrix/$${i}/fmod -o ../src/sunmatrix/$${i}/fmod/fsunmatrix_$${i}_mod.c ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i ; \
# 	done

# sunlinsol: $(SUNLINSOL:%:sunlinsol/fsunlinsol_%_mod.i)
# 	@for i in ${SUNLINSOL} ; do \
# 		echo "swig/swig -fortran -outdir ../src/sunlinsol/$${i}/fmod -o ../src/sunlinsol/$${i}/fmod ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i" ; \
# 		swig/swig -fortran -outdir ../src/sunlinsol/$${i}/fmod -o ../src/sunlinsol/$${i}/fmod/fsunlinsol_$${i}_mod.c ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i ; \
# 	done

# sunnonlinsol: $(SUNNONLINSOL:%:sunnonlinsol/fsunnonlinsol_%_mod.i)
# 	@for i in ${SUNNONLINSOL} ; do \
# 		echo "swig/swig -fortran -outdir ../src/sunnonlinsol/$${i}/fmod -o ../src/sunnonlinsol/$${i}/fmod ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i" ; \
# 		swig/swig -fortran -outdir ../src/sunnonlinsol/$${i}/fmod -o ../src/sunnonlinsol/$${i}/fmod/fsunnonlinsol_$${i}_mod.c ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i ; \
# 	done

clean:
	rm ./**/*.c
	rm ./**/*.py
