
CVODE=fcvode_mod.i
GENERIC=fsundials_types.i fnvector_mod.i fsunmatrix_mod.i fsunlinsol_mod.i fsunnonlinsol_mod.i
NVECTOR=fnvector_openmp_mod.i fnvector_pthreads_mod.i fnvector_serial_mod.i
SUNMATRIX=fsunmatrix_band_mod.i fsunmatrix_dense_mod.i fsunmatrix_sparse_mod.i
SUNLINSOL=
SUNNONLINSOL=

all: generic cvode nvector sunmatrix sunlinsol sunnonlinsol

cvode: ${CVODE}
	@for i in ${CVODE} ; do \
		echo "swig -fortran -outdir ../src/cvode/F90 -I../include $${i}" ; \
		swig -fortran -outdir ../src/cvode/F90 -I../include $${i} ; \
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
