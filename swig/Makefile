
ARKODE=farkode_mod.i farkode_arkstep_mod.i farkode_erkstep_mod.i farkode_mristep_mod.i

CVODE=fcvode_mod.i

IDA=fida_mod.i

KINSOL=fkinsol_mod.i

GENERIC=fsundials_types.i fnvector_mod.i fsunmatrix_mod.i fsunlinsol_mod.i fsunnonlinsol_mod.i

NVECTOR=fnvector_openmp_mod.i fnvector_pthreads_mod.i fnvector_serial_mod.i

SUNMATRIX=fsunmatrix_band_mod.i fsunmatrix_dense_mod.i fsunmatrix_sparse_mod.i

SUNLINSOL=fsunlinsol_band_mod.i fsunlinsol_dense_mod.i fsunlinsol_spbcgs_mod.i fsunlinsol_spfgmr_mod.i \
					fsunlinsol_spgmr_mod.i fsunlinsol_sptfqmr_mod.i fsunlinsol_pcg_mod.i

SUNNONLINSOL=fsunnonlinsol_newton_mod.i fsunnonlinsol_fixedpoint_mod.i

all: generic arkode cvode ida kinsol nvector sunmatrix sunlinsol sunnonlinsol

arkode: ${ARKODE}
	@for i in ${ARKODE} ; do \
		echo "swig -fortran -outdir ../src/akrode/F90 -I../include $${i}" ; \
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

nvector: ${NVECTOR}
	@for i in ${NVECTOR} ; do \
		echo "swig -fortran -outdir out -I../include $${i}" ; \
		swig -fortran -outdir out -I../include $${i} ; \
	done

sunmatrix: ${SUNMATRIX}
	@for i in ${SUNMATRIX} ; do \
		echo "swig -fortran -outdir out -I../include $${i}" ;\
		swig -fortran -outdir out -I../include $${i} ; \
	done

sunlinsol: ${SUNLINSOL}
	@for i in ${SUNLINSOL} ; do \
		echo "swig -fortran -outdir out -I../include $${i}" ; \
		swig -fortran -outdir out -I../include $${i} ; \
	done

sunnonlinsol: ${SUNNONLINSOL}
	@for i in ${SUNNONLINSOL} ; do \
		echo "swig -fortran -outdir out -I../include $${i}" ; \
		swig -fortran -outdir out -I../include $${i} ; \
	done

clean:
	rm *.c out/*.f90
