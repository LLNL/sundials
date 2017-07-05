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
 * This is the header file for the band implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the band implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNMatrix type without making explicit
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

#ifndef _SUNLINSOL_BAND_H
#define _SUNLINSOL_BAND_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_band.h>
#include <sunmatrix/sunmatrix_band.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_pthreads.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Band implementation of SUNLinearSolver
 *
 * The band implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     pivots -- index array for partial pivoting in LU factorization
 *     last_flag -- last error return flag from internal setup/solve
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_Band {
  long int *pivots;
  long int last_flag;
};

typedef struct _SUNLinearSolverContent_Band *SUNLinearSolverContent_Band;

  
/*
 * -----------------------------------------------------------------
 * PART II: macros SLS_CONTENT_B, SLS_DATA_B, SLS_LASTFLAG_B
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNLinearSolver S;
 * SUNLinearSolverContent_Band S_cont;
 * long int *S_pivots, S_lastflag;
 *
 * (1) SLS_CONTENT_B
 *
 *     This macro gives access to the contents of the band
 *     SUNLinearSolver
 *
 *     The assignment S_cont = SLS_CONTENT_B(S) sets S_cont to be
 *     a pointer to the band SUNLinearSolver content structure.
 *
 * (2) SLS_PIVOTS_B, SLS_LASTFLAG_B
 *
 *     These macros give access to the individual parts of
 *     the content structure of a band SUNMatrix.
 *
 *     The assignment S_pivots = SLS_PIVOTS_B(S) sets S_pivots
 *     to be a pointer to pivot array of S.
 *
 *     The assignment S_lastflag = SLS_LASTFLAG_B(S) sets S_lastflag
 *     to be the 'last_flag' entry from the content structure for S.
 *
 * -----------------------------------------------------------------
 */

#define SLS_CONTENT_B(S)     ( (SUNLinearSolverContent_Band)(S->content) )

#define SLS_PIVOTS_B(S)      ( SLS_CONTENT_B(S)->pivots )

#define SLS_LASTFLAG_B(S)    ( SLS_CONTENT_B(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunlinsol_band
 * 
 * CONSTRUCTOR:
 *    SUNBandLinearSolver creates and allocates memory for a banded
 *    matrix solver
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNBandLinearSolver(N_Vector y,
                                                    SUNMatrix A);

/*
 * -----------------------------------------------------------------
 * band implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_Band(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_Band(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_Band(SUNLinearSolver S, void* A_data,
                                            ATSetupFn ATSetup, ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_Band(SUNLinearSolver S,
                                                    void* P_data,
                                                    PSetupFn Pset,
                                                    PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetup_Band(SUNLinearSolver S, SUNMatrix A,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3);
SUNDIALS_EXPORT int SUNLinSolSolve_Band(SUNLinearSolver S, SUNMatrix A,
                                        N_Vector x, N_Vector b,
                                        N_Vector w, realtype tol);
SUNDIALS_EXPORT long int SUNLinSolNumIters_Band(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_Band(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_Band(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif
