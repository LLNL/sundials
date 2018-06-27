/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file for the SUNNonlinearSolver module implementation of
 * a Full Newton iteration (for initial testing).
 * 
 * Part I defines the solver-specific content structure.
 * 
 * Part II contains prototypes for the solver constructor and operations.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNNONLINSOL_FULLNEWTON_H
#define _SUNNONLINSOL_FULLNEWTON_H

#include "sundials/sundials_types.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_nonlinearsolver.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------------------
 * I. Content structure
 * ---------------------------------------------------------------------------*/
  
struct _SUNNonlinearSolverContent_FullNewton {

  /* functions provided by the integrator */
  SUNNonlinSolSysFn      Sys;    /* nonlinear system residual function         */
  SUNNonlinSolLSetupFn   LSetup; /* linear solver setup function               */
  SUNNonlinSolLSolveFn   LSolve; /* linear solver solve function               */
  SUNNonlinSolConvTestFn CTest;  /* nonlinear solver convergence test function */

  /* nonlinear solver variables */
  N_Vector         cor;       /* Newton correction vector                  */
  N_Vector         delta;     /* Newton correction vector                  */
  booleantype      callSetup; /* does the setup function need to be called */
  int              maxiters;  /* maximum number of nonlinear iterations    */
  long int         niters;    /* number of nonlinear iterations            */
};

typedef struct _SUNNonlinearSolverContent_FullNewton *SUNNonlinearSolverContent_FullNewton;

/* -----------------------------------------------------------------------------
 * II: Exported functions
 * ---------------------------------------------------------------------------*/
 
/* Constructor to create solver and allocates memory */
SUNDIALS_EXPORT SUNNonlinearSolver SUNFullNewtonSolver(N_Vector y);

/* core functions */
SUNDIALS_EXPORT SUNNonlinearSolver_Type SUNNonlinSolGetType_FullNewton(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT int SUNNonlinSolInit_FullNewton(SUNNonlinearSolver NLS, N_Vector tmpl);

SUNDIALS_EXPORT int SUNNonlinSolSetup_FullNewton(SUNNonlinearSolver NLS, N_Vector y, void* mem);

SUNDIALS_EXPORT int SUNNonlinSolSolve_FullNewton(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y,
                                             N_Vector w, realtype tol, void *mem);

SUNDIALS_EXPORT int SUNNonlinSolFree_FullNewton(SUNNonlinearSolver NLS);

/* set functions */
SUNDIALS_EXPORT int SUNNonlinSolSetSysFn_FullNewton(SUNNonlinearSolver NLS,
                                                SUNNonlinSolSysFn SysFn);

SUNDIALS_EXPORT int SUNNonlinSolSetLSetupFn_FullNewton(SUNNonlinearSolver NLS,
                                                   SUNNonlinSolLSetupFn LSetupFn);

SUNDIALS_EXPORT int SUNNonlinSolSetLSolveFn_FullNewton(SUNNonlinearSolver NLS,
                                                   SUNNonlinSolLSolveFn LSolveFn);

SUNDIALS_EXPORT int SUNNonlinSolSetConvTestFn_FullNewton(SUNNonlinearSolver NLS,
                                                     SUNNonlinSolConvTestFn CTestFn);

SUNDIALS_EXPORT int SUNNonlinSolSetMaxIters_FullNewton(SUNNonlinearSolver NLS, int maxiters);

/* get functions */
SUNDIALS_EXPORT int SUNNonlinSolGetNumIters_FullNewton(SUNNonlinearSolver NLS, long int *niters);

#ifdef __cplusplus
}
#endif

#endif
