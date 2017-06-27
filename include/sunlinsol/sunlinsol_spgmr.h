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
 * This is the header file for the SPGMR implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the SPGMR implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNLinearSolver type without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor 
 * SUNSPGMR as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPGMR_H
#define _SUNLINSOL_SPGMR_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_spgmr.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_pthreads.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: SPGMR implementation of SUNLinearSolver
 *
 * The SPGMR implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     maxl -- number of GMRES basis vectors to use
 *     pretype -- flag for type of preconditioning to employ
 *     last_flag -- last error return flag from internal setup/solve
 *     ATSetup -- function pointer to setup routine for ATimes data
 *     ATimes -- function pointer to ATimes routine
 *     Psetup -- function pointer to preconditioner setup routine
 *     Psolve -- function pointer to preconditioner solve routine
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_SPGMR {
  int maxl;
  int pretype;
  long int last_flag;
  ATSetupFn ATSetup;
  ATimesFn ATimes;
  void* ATData;
  PSetupFn Psetup;
  PSolveFn Psolve;
  void* PData;
};

typedef struct _SUNLinearSolverContent_SPGMR *SUNLinearSolverContent_SPGMR;

  
/*
 * -----------------------------------------------------------------
 * PART II: macros SLS_CONTENT_SPGMR
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNLinearSolver S;
 * SUNLinearSolverContent_SPGMR S_cont;
 * int S_maxl, S_pretype;
 * long int S_lastflag;
 *
 * (1) SLS_CONTENT_SPGMR
 *
 *     This macro gives access to the contents of the SPGMR
 *     SUNLinearSolver
 *
 *     The assignment S_cont = SLS_CONTENT_SPGMR(S) sets S_cont to be
 *     a pointer to the SPGMR SUNLinearSolver content structure.
 *
 *     The assignment S_maxl = SLS_CONTENT_SPGMR(S)->maxl sets S_maxl
 *     to be the maxl value from the SPGMR content structure
 *
 *     The assignment S_pretype = SLS_CONTENT_SPGMR(S)->pretype sets
 *     S_pretype to be the pretype value from the SPGMR content
 *
 *     The assignment S_lastflag = SLS_CONTENT_SPGMR(S)->lastflag 
 *     sets S_lastflag to be the lastflag value from the SPGMR content
 *
 * -----------------------------------------------------------------
 */

#define SLS_CONTENT_SPGMR(S)  ( (SUNLinearSolverContent_SPGMR)(S->content) )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunlinsol_spgmr
 * 
 * CONSTRUCTOR:
 *    SUNSPGMR creates and allocates memory for a SPGMR solver
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNSPGMR(N_Vector y, int pretype, int maxl);

/*
 * -----------------------------------------------------------------
 * SPGMR implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_SPGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_SPGMR(SUNLinearSolver S, void* A_data,
                                             ATSetupFn ATSetup, ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_SPGMR(SUNLinearSolver S,
                                                     void* P_data,
                                                     PSetupFn Pset,
                                                     PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetup_SPGMR(SUNLinearSolver S, SUNMatrix A,
                                         N_Vector tmp1, N_Vector tmp2,
                                         N_Vector tmp3);
SUNDIALS_EXPORT int SUNLinSolSolve_SPGMR(SUNLinearSolver S, SUNMatrix A,
                                         N_Vector x, N_Vector b,
                                         N_Vector w, realtype tol);
SUNDIALS_EXPORT int SUNLinSolPerformance_SPGMR(SUNLinearSolver S,
                                               int perftask);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_SPGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_SPGMR(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif
