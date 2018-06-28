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
 * This is the implementation file for the SUNNonlinearSolver module
 * implementation of a full Newton iteration (for initial testing).
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunnonlinsol/sunnonlinsol_full_newton.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0) /* real 0.0 */
#define ONE  RCONST(1.0) /* real 1.0 */

/* Content structure accessibility macros  */
#define NEWTON_CONTENT(S) ( (SUNNonlinearSolverContent_FullNewton)(S->content) )

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
 * Constructor to create a new Newton solver
 */

SUNNonlinearSolver SUNFullNewtonSolver(N_Vector y)
{
  SUNNonlinearSolver NLS;
  SUNNonlinearSolver_Ops ops;
  SUNNonlinearSolverContent_FullNewton content;
  
  /* Create nonlinear linear solver */
  NLS = NULL;
  NLS = (SUNNonlinearSolver) malloc(sizeof *NLS);
  if (NLS == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNNonlinearSolver_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(NLS); return(NULL); }

  /* Attach operations */
  ops->gettype     = SUNNonlinSolGetType_FullNewton;
  ops->init        = SUNNonlinSolInit_FullNewton;
  ops->setup       = SUNNonlinSolSetup_FullNewton;
  ops->solve       = SUNNonlinSolSolve_FullNewton;
  ops->free        = SUNNonlinSolFree_FullNewton;
  ops->setsysfn    = SUNNonlinSolSetSysFn_FullNewton;
  ops->setlsetupfn = SUNNonlinSolSetLSetupFn_FullNewton;
  ops->setlsolvefn = SUNNonlinSolSetLSolveFn_FullNewton;
  ops->setctestfn  = SUNNonlinSolSetConvTestFn_FullNewton;
  ops->setmaxiters = SUNNonlinSolSetMaxIters_FullNewton;
  ops->getnumiters = SUNNonlinSolGetNumIters_FullNewton;
  /* ops->lastflag          = SUNNonlinSolLastFlag_FullNewton; */

  /* Create content */
  content = NULL;
  content = (SUNNonlinearSolverContent_FullNewton) malloc(sizeof *content);
  if (content == NULL) { free(ops); free(NLS); return(NULL); }

  /* Fill content */
  content->Sys      = NULL;
  content->LSetup   = NULL;
  content->LSolve   = NULL;
  content->CTest    = NULL;
  content->delta    = N_VClone(y);
  content->maxiters = 3;
  content->niters   = 0;
  
  /* Attach content and ops */
  NLS->content = content;
  NLS->ops     = ops;

  return(NLS);
}


/* -----------------------------------------------------------------------------
 * Implementation of nonlinear solver operations
 * ---------------------------------------------------------------------------*/

/*
 * Core functions
 */

SUNNonlinearSolver_Type SUNNonlinSolGetType_FullNewton(SUNNonlinearSolver NLS)
{
  return(SUNNONLINEARSOLVER_ROOTFIND);
}


int SUNNonlinSolInit_FullNewton(SUNNonlinearSolver NLS, N_Vector tmpl)
{
  /* all solver-specific memory has already been allocated */
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetup_FullNewton(SUNNonlinearSolver NLS, N_Vector y, void* mem)
{
  /* check that all necessary function pointer have been set */
  if ( NEWTON_CONTENT(NLS)->Sys    == NULL ||
       NEWTON_CONTENT(NLS)->LSolve == NULL ||
       NEWTON_CONTENT(NLS)->CTest  == NULL    ) {
    return(-1);
  }
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSolve_FullNewton(SUNNonlinearSolver NLS,
                                 N_Vector y0, N_Vector y,
                                 N_Vector w, realtype tol,
                                 booleantype callSetup, void* mem)
{
  int mnewt;
  int retval;
  realtype delnrm;
  N_Vector delta = NEWTON_CONTENT(NLS)->delta;

  /* initialize counter mnewt */
  mnewt = 0;

  /* load prediction into y */
  N_VScale(ONE, y0, y);

  /* looping point for Newton iteration. Break out on any error. */
  for(;;) {

    /* increment number of nonlinear solver iterations */
    NEWTON_CONTENT(NLS)->niters++;

    /* compute the residual */
    retval = NEWTON_CONTENT(NLS)->Sys(y, delta, mem);
    if (retval != SUN_NLS_SUCCESS) break;

    /* setup the linear system */
    if (NEWTON_CONTENT(NLS)->LSetup) {
      retval = NEWTON_CONTENT(NLS)->LSetup(y, delta, mem);
      if (retval != SUN_NLS_SUCCESS) break;
    }

    /* solve the linear system to get correction vector delta */
    retval = NEWTON_CONTENT(NLS)->LSolve(y, delta, mem);
    if (retval != SUN_NLS_SUCCESS) break;

    /* apply delta to y */
    N_VLinearSum(ONE, y, -ONE, delta, y);

    /* compute the norm of the correction */
    delnrm = N_VWrmsNorm(delta, w);

    /* test for convergence, return if successful */
    retval = NEWTON_CONTENT(NLS)->CTest(mnewt, delnrm, tol, mem);
    if (retval == SUN_NLS_SUCCESS)  return(SUN_NLS_SUCCESS);
    if (retval != SUN_NLS_CONTINUE) break;

    /* not yet converged. Increment mnewt and test for max allowed. */
    mnewt++;
    if (mnewt >= NEWTON_CONTENT(NLS)->maxiters) {
      retval = SUN_NLS_CONV_RECVR;
      break;
    }

  } /* end of Newton iteration loop */

  /* all error returns exit here */
  return(retval);
}


int SUNNonlinSolFree_FullNewton(SUNNonlinearSolver NLS)
{
  /* return if NLS is already free */
  if (NLS == NULL)
    return(SUN_NLS_SUCCESS);

  /* free items from contents, then the generic structure */
  if (NLS->content) {

    if (NEWTON_CONTENT(NLS)->delta)
      N_VDestroy(NEWTON_CONTENT(NLS)->delta);

    free(NLS->content);
    NLS->content = NULL;
  }

  /* free the ops structure */
  if (NLS->ops) {
    free(NLS->ops);
    NLS->ops = NULL;
  }

  free(NLS);
  NLS = NULL;

  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetSysFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)
{
  NEWTON_CONTENT(NLS)->Sys = SysFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetLSetupFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn LSetupFn)
{
  NEWTON_CONTENT(NLS)->LSetup = LSetupFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetLSolveFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn LSolveFn)
{
  NEWTON_CONTENT(NLS)->LSolve = LSolveFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetConvTestFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn)
{
  NEWTON_CONTENT(NLS)->CTest = CTestFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetMaxIters_FullNewton(SUNNonlinearSolver NLS, int maxiters)
{
  if (maxiters < 1) return(1);
  NEWTON_CONTENT(NLS)->maxiters = maxiters;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolGetNumIters_FullNewton(SUNNonlinearSolver NLS, long int *niters)
{
  *niters = NEWTON_CONTENT(NLS)->niters;
  return(SUN_NLS_SUCCESS);
}
