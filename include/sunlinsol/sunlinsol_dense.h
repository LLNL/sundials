/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the header file for the dense implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the dense implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNMatrix type without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor 
 * SLSNew_Dense as well as implementation-specific prototypes 
 * for various useful matrix operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can be found
 *     in the header file sundials_linearsolver.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_DENSE_H
#define _SUNLINSOL_DENSE_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_pthreads.h>


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Dense implementation of SUNLinearSolver
 *
 * The dense implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *	 pivots -- index array for partial pivoting in LU factorization
*	 last_flag -- last error return flag from internal setup/solve
 * -----------------------------------------------------------------
 */
  
struct _SUNMatrixContent_Dense {
  long int *pivots;
  long int last_flag;
};

typedef struct _SUNLinearSolverContent_Dense *SUNLinearSolverContent_Dense;

/*
 * -----------------------------------------------------------------
 * PART II: macros SLS_CONTENT_D, SLS_DATA_D, SLS_LASTFLAG_D * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNLinearSolver S;
 * SUNLinearSolverContent_Dense S_cont;
 * long int *S_pivots, S_lastflag;
 *
 * (1) SLS_CONTENT_D
 *
 *     This macro gives access to the contents of the dense
 *     SUNLinearSolver
 *
 *     The assignment S_cont = SLS_CONTENT_B(S) sets S_cont to
 *	  be a pointer to the dense SUNLinearSolver content
 *	  structure.
 *
 * (2) SLS_PIVOTS_D, SLS_LASTFLAG_D
 *
 *     These macros give access to the individual parts of
 *     the content structure of a dense SUNMatrix.
 *
 *     The assignment S_pivots = SLS_PIVOTS_B(S) sets S_pivots
 *     to be a pointer to pivot array of S.
 *
 *	  The assignment S_lastflag = SLS_LASTFLAG_D(S) sets 
 *  	  S_lastflag to be the 'last_flag' eentry from the content 
 * 	  structure for S.
 * ----------------------------------------------------------------- */
  
#define SLS_CONTENT_D(S)     ( (SUNLinearSolverContent_Dense)(S->content) )

#define SLS_PIVOTS_D(S)      ( SLS_CONTENT_D(S)->pivots )

#define SLS_LASTFLAG_D(S)    ( SLS_CONTENT_D(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunlinsol_dense
 * 
 * CONSTRUCTORS:
 *    SUNDenseLinearSolver creates and allocates memory for a 
 * 	 dense matrix solver
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNDenseLinearSolver(N_Vector y,
                                                     SUNMatrix A);

/*
 * -----------------------------------------------------------------
 * dense implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_Dense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_Dense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_Dense(SUNLinearSolver S, void* A_data,
                                             ATSetupFn ATSetup, ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_Dense(SUNLinearSolver S,
                                                     void* P_data,
                                                     PSetupFn Pset,
                                                     PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetup_Dense(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_Dense(SUNLinearSolver S, SUNMatrix A,
                                        N_Vector x, N_Vector b, 
                                        N_Vector w, realtype tol);
SUNDIALS_EXPORT long int SUNLinSolNumIters_Dense(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_Dense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_Dense(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif

