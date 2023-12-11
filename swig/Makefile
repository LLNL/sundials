# ---------------------------------------------------------------
# Programmer: Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
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

ARKODE=farkode_mod farkode_arkstep_mod farkode_erkstep_mod farkode_sprkstep_mod farkode_mristep_mod
CVODE=fcvode_mod
CVODES=fcvodes_mod
IDA=fida_mod
IDAS=fidas_mod
KINSOL=fkinsol_mod
GENERIC=fsundials_types_mod fsundials_nvector_mod fsundials_matrix_mod fsundials_linearsolver_mod fsundials_nonlinearsolver_mod fsundials_futils_mod fsundials_context_mod fsundials_profiler_mod fsundials_logger_mod fsundials_adaptcontroller_mod
NVECTOR=openmp pthreads serial parallel manyvector mpiplusx
SUNMATRIX=band dense sparse
SUNLINSOL=band dense lapackdense klu spbcgs spfgmr spgmr sptfqmr pcg
SUNNONLINSOL=newton fixedpoint
SUNADAPTCONTROLLER=imexgus soderlind

INCLUDES=-I../include

.PHONY: generic arkode cvode cvodes ida idas kinsol nvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller patch

all: generic arkode cvode cvodes ida idas kinsol nvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller patch

arkode: $(ARKODE:%:arkode/%.i)
	@for i in ${ARKODE} ; do \
		echo "swig -fortran -outdir ../src/arkode/fmod -o ../src/arkode/fmod/$${i}.c ${INCLUDES} arkode/$${i}.i" ; \
		swig -fortran -outdir ../src/arkode/fmod -o ../src/arkode/fmod/$${i}.c ${INCLUDES} arkode/$${i}.i ; \
	done

cvode: $(CVODE:%:cvode/%.i)
	@for i in ${CVODE} ; do \
		echo "swig -fortran -outdir ../src/cvode/fmod -o ../src/cvode/fmod/$${i}.c ${INCLUDES} cvode/$${i}.i" ; \
		swig -fortran -outdir ../src/cvode/fmod -o ../src/cvode/fmod/$${i}.c ${INCLUDES} cvode/$${i}.i ; \
	done

cvodes: $(CVODE:%:cvodes/%.i)
	@for i in ${CVODES} ; do \
		echo "swig -fortran -outdir ../src/cvodes/fmod -o ../src/cvodes/fmod/$${i}.c ${INCLUDES} cvodes/$${i}.i" ; \
		swig -fortran -outdir ../src/cvodes/fmod -o ../src/cvodes/fmod/$${i}.c ${INCLUDES} cvodes/$${i}.i ; \
	done

ida: $(IDA:%:ida/%.i)
	@for i in ${IDA} ; do \
		echo "swig -fortran -outdir ../src/ida/fmod -o ../src/ida/fmod/$${i}.c ${INCLUDES} ida/$${i}.i" ; \
		swig -fortran -outdir ../src/ida/fmod -o ../src/ida/fmod/$${i}.c ${INCLUDES} ida/$${i}.i ; \
	done

idas: $(IDAS:%:idas/%.i)
	@for i in ${IDAS} ; do \
		echo "swig -fortran -outdir ../src/idas/fmod -o ../src/idas/fmod/$${i}.c ${INCLUDES} idas/$${i}.i" ; \
		swig -fortran -outdir ../src/idas/fmod -o ../src/idas/fmod/$${i}.c ${INCLUDES} idas/$${i}.i ; \
	done

kinsol: $(KINSOL:%:kinsol/%.i)
	@for i in ${KINSOL} ; do \
		echo "swig -fortran -outdir ../src/kinsol/fmod -o ../src/kinsol/fmod/$${i}.c ${INCLUDES} kinsol/$${i}" ; \
		swig -fortran -outdir ../src/kinsol/fmod -o ../src/kinsol/fmod/$${i}.c ${INCLUDES} kinsol/$${i}.i ; \
	done

generic: $(GENERIC:%:sundials/%.i)
	@for i in ${GENERIC} ; do \
		echo "swig -fortran -outdir ../src/sundials/fmod -o ../src/sundials/fmod/$${i}.c ${INCLUDES} $${i}" ; \
		swig -fortran -outdir ../src/sundials/fmod -o ../src/sundials/fmod/$${i}.c ${INCLUDES} sundials/$${i}.i ; \
	done

nvector: $(NVECTOR:%:nvector/fnvector_%_mod.i) mpimanyvector
	@for i in ${NVECTOR} ; do \
		echo "swig -fortran -outdir ../src/nvector/$${i}/fmod -o ../src/nvector/$${i}/fmod ${INCLUDES} nvector/fnvector_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/nvector/$${i}/fmod -o ../src/nvector/$${i}/fmod/fnvector_$${i}_mod.c ${INCLUDES} nvector/fnvector_$${i}_mod.i ; \
	done

mpimanyvector: nvector/fnvector_mpimanyvector_mod.i
	@echo "swig -fortran -outdir ../src/nvector/manyvector/fmod -o ../src/nvector/manyvector/fmod ${INCLUDES} nvector/fnvector_mpimanyvector_mod.i" ; \
	swig -fortran -outdir ../src/nvector/manyvector/fmod -o ../src/nvector/manyvector/fmod/fnvector_mpimanyvector_mod.c ${INCLUDES} nvector/fnvector_mpimanyvector_mod.i ; \

sunmatrix: $(SUNMATRIX:%:sunmatrix/fsunmatrix_%_mod.i)
	@for i in ${SUNMATRIX} ; do \
		echo "swig -fortran -outdir ../src/sunmatrix/$${i}/fmod -o ../src/sunmatrix/$${i}/fmod ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunmatrix/$${i}/fmod -o ../src/sunmatrix/$${i}/fmod/fsunmatrix_$${i}_mod.c ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i ; \
	done

sunlinsol: $(SUNLINSOL:%:sunlinsol/fsunlinsol_%_mod.i)
	@for i in ${SUNLINSOL} ; do \
		echo "swig -fortran -outdir ../src/sunlinsol/$${i}/fmod -o ../src/sunlinsol/$${i}/fmod ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunlinsol/$${i}/fmod -o ../src/sunlinsol/$${i}/fmod/fsunlinsol_$${i}_mod.c ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i ; \
	done

sunnonlinsol: $(SUNNONLINSOL:%:sunnonlinsol/fsunnonlinsol_%_mod.i)
	@for i in ${SUNNONLINSOL} ; do \
		echo "swig -fortran -outdir ../src/sunnonlinsol/$${i}/fmod -o ../src/sunnonlinsol/$${i}/fmod ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunnonlinsol/$${i}/fmod -o ../src/sunnonlinsol/$${i}/fmod/fsunnonlinsol_$${i}_mod.c ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i ; \
	done

sunadaptcontroller: $(SUNADAPTCONTROLLER:%:sunadaptcontroller/fsunadaptcontroller_%_mod.i)
	@for i in ${SUNADAPTCONTROLLER} ; do \
		echo "swig -fortran -outdir ../src/sunadaptcontroller/$${i}/fmod -o ../src/sunadaptcontroller/$${i}/fmod ${INCLUDES} sunadaptcontroller/fsunadaptcontroller_$${i}_mod.i" ; \
		swig -fortran -outdir ../src/sunadaptcontroller/$${i}/fmod -o ../src/sunadaptcontroller/$${i}/fmod/fsunadaptcontroller_$${i}_mod.c ${INCLUDES} sunadaptcontroller/fsunadaptcontroller_$${i}_mod.i ; \
	done

patch: .FORCE
	cd ..;

clean:
	rm ../src/**/fmod/*.c; rm ../src/**/fmod/*.f90

.FORCE: ;
