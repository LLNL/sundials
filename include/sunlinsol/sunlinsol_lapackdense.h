/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the LAPACK dense implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the LAPACK dense 
 * implementation of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor 
 * SUNLapackDense as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_LAPDENSE_H
#define _SUNLINSOL_LAPDENSE_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_lapack.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_pthreads.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Interfaces to match 'realtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgetrf_f77    dgetrf_f77
#define xgetrs_f77    dgetrs_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgetrf_f77    sgetrf_f77
#define xgetrs_f77    sgetrs_f77
#else
#error  Incompatible realtype for LAPACK; disable LAPACK and rebuild
#endif

/* Catch to disable LAPACK linear solvers with incompatible sunindextype */
#if defined(SUNDIALS_SIGNED_32BIT_TYPE)
#else  /* incompatible sunindextype for LAPACK */
#error  Incompatible sunindextype for LAPACK; disable LAPACK and rebuild
#endif

/*
 * -----------------------------------------------------------------
 * PART I: LAPACK dense implementation of SUNLinearSolver
 *
 * The LAPACK dense implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     pivots -- index array for partial pivoting in LU factorization
 *     last_flag -- last error return flag from internal setup/solve
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_LapackDense {
  int *pivots;
  long int last_flag;
};

typedef struct _SUNLinearSolverContent_LapackDense *SUNLinearSolverContent_LapackDense;

  
/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_lapackdense
 * 
 * CONSTRUCTOR:
 *    SUNLapackDense creates and allocates memory for a LAPACK dense
 *      matrix solver
 *
 * OTHER:
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNLapackDense(N_Vector y, SUNMatrix A);

/*
 * -----------------------------------------------------------------
 * LAPACK dense implementations of required linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_LapackDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_LapackDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_LapackDense(SUNLinearSolver S, void* A_data,
                                                   ATSetupFn ATSetup, ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_LapackDense(SUNLinearSolver S,
                                                           void* P_data,
                                                           PSetupFn Pset,
                                                           PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_LapackDense(SUNLinearSolver S,
                                                           N_Vector s1,
                                                           N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetup_LapackDense(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_LapackDense(SUNLinearSolver S, SUNMatrix A,
                                               N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_LapackDense(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_LapackDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolNumPSolves_LapackDense(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_LapackDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_LapackDense(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif
