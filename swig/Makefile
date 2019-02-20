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
# Makefile to generate SUNDIALS fortran interfaces with swig
# ---------------------------------------------------------------

ARKODE=arkode/farkode_mod.i arkode/farkode_arkstep_mod.i arkode/farkode_erkstep_mod.i arkode/farkode_mristep_mod.i
CVODE=cvode/fcvode_mod.i
IDA=ida/fida_mod.i
KINSOL=kinsol/fkinsol_mod.i
GENERIC=sundials/fsundials_types.i sundials/fnvector_mod.i sundials/fsunmatrix_mod.i sundials/fsunlinsol_mod.i \
				sundials/fsunnonlinsol_mod.i
NVECTOR=openmp pthreads serial
SUNMATRIX=band dense sparse
SUNLINSOL=band dense spbcgs spfgmr spgmr sptfqmr pcg
SUNNONLINSOL=newton fixedpoint

.PHONY: generic arkode cvode ida kinsol nvector sunmatrix sunlinsol sunnonlinsol

all: generic arkode cvode ida kinsol nvector sunmatrix sunlinsol sunnonlinsol

arkode: ${ARKODE}
	@for i in ${ARKODE} ; do \
		echo "swig -fortran -outdir ../src/arkode/F90 -I../include $${i}" ; \
		swig -fortran -outdir ../src/arkode/F90 -I../include $${i} ; \
	done

cvode: ${CVODE}
	@for i in ${CVODE} ; do \
		echo "swig -fortran -outdir ../src/cvode/F90 -I../include $${i}" ; \
		swig -fortran -outdir ../src/cvode/F90 -I../include $${i} ; \
	done

ida: ${IDA}
	@for i in ${IDA} ; do \
		echo "swig -fortran -outdir ../src/ida/F90 -I../include $${i}" ; \
		swig -fortran -outdir ../src/ida/F90 -I../include $${i} ; \
	done

kinsol: ${KINSOL}
	@for i in ${KINSOL} ; do \
		echo "swig -fortran -outdir ../src/kinsol/F90 -I../include $${i}" ; \
		swig -fortran -outdir ../src/kinsol/F90 -I../include $${i} ; \
	done

generic: ${GENERIC}
	@for i in ${GENERIC} ; do \
		echo "swig -fortran -outdir ../src/sundials/F90 -I../include $${i}" ; \
		swig -fortran -outdir ../src/sundials/F90 -I../include $${i} ; \
	done

nvector: $(NVECTOR:%:nvector/fnvector_%_mod.i)
	@for i in ${NVECTOR} ; do \
		echo "swig -fortran -outdir ../src/nvector/$${i}/F90 -I../include nvector/fnvector_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/nvector/$${i}/F90 -I../include nvector/fnvector_$${i}_mod.i ; \
	done

sunmatrix: $(SUNMATRIX:%:sunmatrix/fsunmatrix_%_mod.i)
	@for i in ${SUNMATRIX} ; do \
		echo "swig -fortran -outdir ../src/sunmatrix/$${i}/F90 -I../include sunmatrix/fsunmatrix_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunmatrix/$${i}/F90 -I../include sunmatrix/fsunmatrix_$${i}_mod.i ; \
	done

sunlinsol: $(SUNLINSOL:%:sunlinsol/fsunlinsol_%_mod.i)
	@for i in ${SUNLINSOL} ; do \
		echo "swig -fortran -outdir ../src/sunlinsol/$${i}/F90 -I../include sunlinsol/fsunlinsol_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunlinsol/$${i}/F90 -I../include sunlinsol/fsunlinsol_$${i}_mod.i ; \
	done

sunnonlinsol: $(SUNNONLINSOL:%:sunnonlinsol/fsunnonlinsol_%_mod.i)
	@for i in ${SUNNONLINSOL} ; do \
		echo "swig -fortran -outdir ../src/sunnonlinsol/$${i}/F90 -I../include sunnonlinsol/fsunnonlinsol_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunnonlinsol/$${i}/F90 -I../include sunnonlinsol/fsunnonlinsol_$${i}_mod.i ; \
	done

clean:
	rm ./**/*.c
