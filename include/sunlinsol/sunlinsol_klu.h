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
 * This is the header file for the KLU implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the KLU implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNLinearSolver type without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor 
 * SLSNew_Band as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_KLU_H
#define _SUNLINSOL_KLU_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "klu.h"
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_pthreads.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: KLU implementation of SUNLinearSolver
 *
 * The KLU implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     last_flag -- last error return flag from internal setup/solve
 *     s_Symbolic -- KLU storage structure for symbolic 
 *       factorization components
 *     s_Numeric -- KLU storage structure for numeric factorization
 *        components
 *     s_Common -- storage structure for common KLU solver 
 *        components
 *     s_ordering -- ordering choice 
 *     sun_klu_solve -- function pointer handling CSR/CSC difference
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_KLU {
  long int      last_flag;
  klu_symbolic *s_Symbolic;
  klu_numeric  *s_Numeric;
  klu_common    s_Common;
  int           s_ordering;
  int          (*sun_klu_solve)(klu_symbolic*, klu_numeric*, int, int, double*, klu_common*);
};

typedef struct _SUNLinearSolverContent_KLU *SUNLinearSolverContent_KLU;

  
/*
 * -----------------------------------------------------------------
 * PART II: macros SLS_CONTENT_K, SLS_LASTFLAG_K, SLS_SYMBOLIC_K, 
 *    SLS_NUMERIC_K, SLS_COMMON_K, SLS_ORDERING_K and SLS_SOLVE_K
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNLinearSolver S;
 * SUNLinearSolverContent_KLU S_cont;
 * long int S_lastflag;
 * klu_symbolic *S_symbolic
 * klu_numeric *S_numeric
 * int S_ordering
 *
 * (1) SLS_CONTENT_K
 *
 *     This macro gives access to the contents of the KLU
 *     SUNLinearSolver
 *
 *     The assignment S_cont = SLS_CONTENT_K(S) sets S_cont to be
 *     a pointer to the KLU SUNLinearSolver content structure.
 *
 * (2) SLS_LASTFLAG_K, SLS_SYMBOLIC_K, SLS_NUMERIC_K, 
 *     SLS_COMMON_K, SLS_ORDERING_K and SLS_SOLVE_K
 *
 *     These macros give access to the individual parts of
 *     the content structure of a KLU SUNLinearSolver.
 *
 *     The assignment S_lastflag = SLS_LASTFLAG_K(S) sets S_lastflag
 *     to be the 'last_flag' entry from the content structure for S.
 *
 *     The assignment S_symbolic = SLS_SYMBOLIC_K(S) sets S_symbolic
 *     to be a pointer to the KLU symbolic solver structure of S.
 *
 *     The assignment S_numeric = SLS_NUMERIC_K(S) sets S_numeric
 *     to be a pointer to the KLU numeric solver structure of S.
 *
 *     The assignment S_ordering = SLS_ORDERING_K(S) sets S_ordering
 *     to be the integer flag for the ordering type of the KLU 
 *     solver interface S.
 *
 *     Access to contents of the 'klu_common' structure may be 
 *     obtained via the macro SLS_COMON_K(S), e.g. 
 *     SLS_COMMON_K(S).condest or SLS_COMMON_K(S).rcond
 *
 *     A function pointer to the specific KLU solver function to be
 *     used (klu_solve vs klu_tsolve) is accessible via the macro 
 *     SLS_SOLVE_K(S), e.g.  (SLS_SOLVE_K(S))(function arguments)
 *
 * -----------------------------------------------------------------
 */

#define SLS_CONTENT_K(S)     ( (SUNLinearSolverContent_KLU)(S->content) )

#define SLS_LASTFLAG_K(S)    ( SLS_CONTENT_K(S)->last_flag )

#define SLS_SYMBOLIC_K(S)    ( SLS_CONTENT_K(S)->s_Symbolic )

#define SLS_NUMERIC_K(S)     ( SLS_CONTENT_K(S)->s_Numeric )

#define SLS_COMMON_K(S)      ( SLS_CONTENT_K(S)->s_Common )

#define SLS_ORDERING_K(S)    ( SLS_CONTENT_K(S)->s_ordering )

#define SLS_SOLVE_K(S)       ( SLS_CONTENT_K(S)->sun_klu_solve )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunlinsol_klu
 * 
 * CONSTRUCTOR:
 *    SUNKLULinearSolver creates and allocates memory for a KLU
 *      sparse-direct linear solver
 *
 * OTHER:
 *    SUNKLUReInit reinitializes memory and flags for a new 
 *      factorization (symbolic and numeric) to be conducted at the 
 *      next solver setup call.  This routine is useful in the 
 *      cases where the number of nonzeroes has changed or if the 
 *      structure of the linear system has changed which would 
 *      require a new symbolic (and numeric factorization).
 *
 *      The reinit_type argumenmt governs the level of 
 *      reinitialization:
 *
 *      reinit_type = 1: The Jacobian matrix will be destroyed and 
 *                       a new one will be allocated based on the 
 *                       nnz value passed to this call. New 
 *                       symbolic and numeric factorizations will 
 *                       be completed at the next solver setup.
 *
 *      reinit_type = 2: Only symbolic and numeric factorizations 
 *                       will be completed.  It is assumed that the 
 *                       Jacobian size has not exceeded the size of 
 *                       nnz given in the sparse matrix provided to
 *                       the original constructor routine.
 *
 *      This routine assumes no other changes to solver use are 
 *      necessary.
 *
 *    SUNKLUSetOrdering sets the ordering used by KLU for reducing 
 *      fill in the linear solve.  Options for ordering_choice are: 
 *          0 for AMD, 
 *          1 for COLAMD, and 
 *          2 for the natural ordering.
 *      The default is 1 for COLAMD.
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNKLULinearSolver(N_Vector y,
                                                   SUNMatrix A);

SUNDIALS_EXPORT int SUNKLUReInit(SUNLinearSolver S, SUNMatrix A,
                                 sunindextype nnz, int reinit_type);

SUNDIALS_EXPORT int SUNKLUSetOrdering(SUNLinearSolver S,
                                      int ordering_choice);

/*
 * -----------------------------------------------------------------
 * KLU implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_KLU(SUNLinearSolver S, void* A_data,
                                           ATSetupFn ATSetup, ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_KLU(SUNLinearSolver S,
                                                   void* P_data,
                                                   PSetupFn Pset,
                                                   PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_KLU(SUNLinearSolver S,
                                                   N_Vector s1,
                                                   N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetup_KLU(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_KLU(SUNLinearSolver S, SUNMatrix A,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolNumPSolves_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_KLU(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif
