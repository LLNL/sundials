# ---------------------------------------------------------------
# Programmer: Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
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
GENERIC=fsundials_core_mod
NVECTOR=openmp pthreads serial parallel manyvector mpiplusx
SUNMATRIX=band dense sparse
SUNLINSOL=band dense lapackdense klu spbcgs spfgmr spgmr sptfqmr pcg
SUNNONLINSOL=newton fixedpoint
SUNADAPTCONTROLLER=imexgus soderlind

INCLUDES=-I../include

.PHONY: generic arkode cvode cvodes ida idas kinsol nvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller patch

.SETINT32:
	$(eval INT_SIZE=32)

.SETINT64:
	$(eval INT_SIZE=64)

all_64: .SETINT64 generic arkode cvode cvodes ida idas kinsol nvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller patch

all_32: .SETINT32 generic arkode cvode cvodes ida idas kinsol nvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller patch

all: all_64 all_32


arkode: $(ARKODE:%:arkode/%.i)
	@for i in ${ARKODE} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/arkode/fmod_int${INT_SIZE} -o ../src/arkode/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} arkode/$${i}.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/arkode/fmod_int${INT_SIZE} -o ../src/arkode/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} arkode/$${i}.i ; \
	done

cvode: $(CVODE:%:cvode/%.i)
	@for i in ${CVODE} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/cvode/fmod_int${INT_SIZE} -o ../src/cvode/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} cvode/$${i}.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/cvode/fmod_int${INT_SIZE} -o ../src/cvode/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} cvode/$${i}.i ; \
	done

cvodes: $(CVODE:%:cvodes/%.i)
	@for i in ${CVODES} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/cvodes/fmod_int${INT_SIZE} -o ../src/cvodes/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} cvodes/$${i}.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/cvodes/fmod_int${INT_SIZE} -o ../src/cvodes/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} cvodes/$${i}.i ; \
	done

ida: $(IDA:%:ida/%.i)
	@for i in ${IDA} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/ida/fmod_int${INT_SIZE} -o ../src/ida/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} ida/$${i}.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/ida/fmod_int${INT_SIZE} -o ../src/ida/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} ida/$${i}.i ; \
	done

idas: $(IDAS:%:idas/%.i)
	@for i in ${IDAS} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/idas/fmod_int${INT_SIZE} -o ../src/idas/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} idas/$${i}.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/idas/fmod_int${INT_SIZE} -o ../src/idas/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} idas/$${i}.i ; \
	done

kinsol: $(KINSOL:%:kinsol/%.i)
	@for i in ${KINSOL} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/kinsol/fmod_int${INT_SIZE} -o ../src/kinsol/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} kinsol/$${i}" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/kinsol/fmod_int${INT_SIZE} -o ../src/kinsol/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} kinsol/$${i}.i ; \
	done

generic: $(GENERIC:%:sundials/%.i)
	@for i in ${GENERIC} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sundials/fmod_int${INT_SIZE} -o ../src/sundials/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} $${i}" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sundials/fmod_int${INT_SIZE} -o ../src/sundials/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} sundials/$${i}.i ; \
	done

nvector: $(NVECTOR:%:nvector/fnvector_%_mod.i) mpimanyvector
	@for i in ${NVECTOR} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/nvector/$${i}/fmod_int${INT_SIZE} -o ../src/nvector/$${i}/fmod_int${INT_SIZE} ${INCLUDES} nvector/fnvector_$${i}_mod.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/nvector/$${i}/fmod_int${INT_SIZE} -o ../src/nvector/$${i}/fmod_int${INT_SIZE}/fnvector_$${i}_mod.c ${INCLUDES} nvector/fnvector_$${i}_mod.i ; \
	done

mpimanyvector: nvector/fnvector_mpimanyvector_mod.i
	@echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/nvector/manyvector/fmod_int${INT_SIZE} -o ../src/nvector/manyvector/fmod_int${INT_SIZE} ${INCLUDES} nvector/fnvector_mpimanyvector_mod.i" ; \
	swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/nvector/manyvector/fmod_int${INT_SIZE} -o ../src/nvector/manyvector/fmod_int${INT_SIZE}/fnvector_mpimanyvector_mod.c ${INCLUDES} nvector/fnvector_mpimanyvector_mod.i ; \

sunmatrix: $(SUNMATRIX:%:sunmatrix/fsunmatrix_%_mod.i)
	@for i in ${SUNMATRIX} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunmatrix/$${i}/fmod_int${INT_SIZE} -o ../src/sunmatrix/$${i}/fmod_int${INT_SIZE} ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunmatrix/$${i}/fmod_int${INT_SIZE} -o ../src/sunmatrix/$${i}/fmod_int${INT_SIZE}/fsunmatrix_$${i}_mod.c ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i ; \
	done

sunlinsol: $(SUNLINSOL:%:sunlinsol/fsunlinsol_%_mod.i)
	@for i in ${SUNLINSOL} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunlinsol/$${i}/fmod_int${INT_SIZE} -o ../src/sunlinsol/$${i}/fmod_int${INT_SIZE} ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunlinsol/$${i}/fmod_int${INT_SIZE} -o ../src/sunlinsol/$${i}/fmod_int${INT_SIZE}/fsunlinsol_$${i}_mod.c ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i ; \
	done

sunnonlinsol: $(SUNNONLINSOL:%:sunnonlinsol/fsunnonlinsol_%_mod.i)
	@for i in ${SUNNONLINSOL} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunnonlinsol/$${i}/fmod_int${INT_SIZE} -o ../src/sunnonlinsol/$${i}/fmod_int${INT_SIZE} ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunnonlinsol/$${i}/fmod_int${INT_SIZE} -o ../src/sunnonlinsol/$${i}/fmod_int${INT_SIZE}/fsunnonlinsol_$${i}_mod.c ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i ; \
	done

sunadaptcontroller: $(SUNADAPTCONTROLLER:%:sunadaptcontroller/fsunadaptcontroller_%_mod.i)
	@for i in ${SUNADAPTCONTROLLER} ; do \
		echo "swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunadaptcontroller/$${i}/fmod_int${INT_SIZE} -o ../src/sunadaptcontroller/$${i}/fmod_int${INT_SIZE} ${INCLUDES} sunadaptcontroller/fsunadaptcontroller_$${i}_mod.i" ; \
		swig -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunadaptcontroller/$${i}/fmod_int${INT_SIZE} -o ../src/sunadaptcontroller/$${i}/fmod_int${INT_SIZE}/fsunadaptcontroller_$${i}_mod.c ${INCLUDES} sunadaptcontroller/fsunadaptcontroller_$${i}_mod.i ; \
	done

patch: .FORCE
	cd ..;

clean:
	rm ../src/**/fmod_int${INT_SIZE}/*.c; rm ../src/**/fmod_int${INT_SIZE}/*.f90

.FORCE: ;
