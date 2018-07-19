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
 * implementation of a full Newton iteration.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunnonlinsol/sunnonlinsol_full_newton.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0) /* real 0.0 */
#define ONE  RCONST(1.0) /* real 1.0 */

/* Content structure accessibility macros  */
#define NEWTON_CONTENT(S) ( (SUNNonlinearSolverContent_FullNewton)(S->content) )

/* =============================================================================
 * Constructor to create a new Full Newton solver
 * ===========================================================================*/

SUNNonlinearSolver SUNNonlinSol_FullNewton(N_Vector y)
{
  SUNNonlinearSolver NLS;
  SUNNonlinearSolver_Ops ops;
  SUNNonlinearSolverContent_FullNewton content;

  /* Check that the supplied N_Vector is non-NULL */
  if (y == NULL) return(NULL);
  
  /* Check that the supplied N_Vector supports all required operations */
  if ( (y->ops->nvclone     == NULL) ||
       (y->ops->nvdestroy   == NULL) ||
       (y->ops->nvscale     == NULL) ||
       (y->ops->nvlinearsum == NULL) ||
       (y->ops->nvwrmsnorm  == NULL) )
    return(NULL);

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
  ops->initialize  = SUNNonlinSolInitialize_FullNewton;
  ops->setup       = NULL; /* no setup needed */
  ops->solve       = SUNNonlinSolSolve_FullNewton;
  ops->free        = SUNNonlinSolFree_FullNewton;
  ops->setsysfn    = SUNNonlinSolSetSysFn_FullNewton;
  ops->setlsetupfn = SUNNonlinSolSetLSetupFn_FullNewton;
  ops->setlsolvefn = SUNNonlinSolSetLSolveFn_FullNewton;
  ops->setctestfn  = SUNNonlinSolSetConvTestFn_FullNewton;
  ops->setmaxiters = SUNNonlinSolSetMaxIters_FullNewton;
  ops->getnumiters = SUNNonlinSolGetNumIters_FullNewton;
  ops->getcuriter  = SUNNonlinSolGetCurIter_FullNewton;

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
  content->mnewt    = 0;
  content->maxiters = 3;
  content->niters   = 0;
  
  /* check if clone was successful */
  if (content->delta == NULL) { free(ops); free(NLS); return(NULL); }

  /* Attach content and ops */
  NLS->content = content;
  NLS->ops     = ops;

  return(NLS);
}


/* =============================================================================
 * GetType, Initialize, Setup, Solve, and Free operations
 * ===========================================================================*/

SUNNonlinearSolver_Type SUNNonlinSolGetType_FullNewton(SUNNonlinearSolver NLS)
{
  return(SUNNONLINEARSOLVER_ROOTFIND);
}


int SUNNonlinSolInitialize_FullNewton(SUNNonlinearSolver NLS)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* check that all required function pointers have been set */
  if ( (NEWTON_CONTENT(NLS)->Sys    == NULL) ||
       (NEWTON_CONTENT(NLS)->LSolve == NULL) ||
       (NEWTON_CONTENT(NLS)->CTest  == NULL) ) {
    return(SUN_NLS_MEM_NULL);
  }

  /* reset the total number of nonlinear solver iterations */
  NEWTON_CONTENT(NLS)->niters = 0;

  return(SUN_NLS_SUCCESS);
}


/* -----------------------------------------------------------------------------
 * SUNNonlinSolSolve_Newton: Performs the nonlinear solve F(y) = 0
 *
 * Successful solve return code:
 *  SUN_NLS_SUCCESS = 0
 *
 * Recoverable failure return codes (positive):
 *   SUN_NLS_CONV_RECVR
 *   *_RHSFUNC_RECVR (ODEs) or *_RES_RECVR (DAEs)
 *   *_LSETUP_RECVR
 *   *_LSOLVE_RECVR
 *
 * Unrecoverable failure return codes (negative):
 *   *_MEM_NULL
 *   *_RHSFUNC_FAIL (ODEs) or *_RES_FAIL (DAEs)
 *   *_LSETUP_FAIL
 *   *_LSOLVE_FAIL
 *
 * Note return values beginning with * are package specific values returned by
 * the Sys, LSetup, and Solve functions provided to the nonlinear solver.
 * ---------------------------------------------------------------------------*/
int SUNNonlinSolSolve_FullNewton(SUNNonlinearSolver NLS,
                                 N_Vector y0, N_Vector y,
                                 N_Vector w, realtype tol,
                                 booleantype callLSetup, void* mem)
{
  int retval;
  booleantype jcur;
  booleantype jbad;
  N_Vector delta;

  /* check that the inputs are non-null */
  if ( (NLS == NULL) ||
       (y0  == NULL) ||
       (y   == NULL) ||
       (w   == NULL) ||
       (mem == NULL) )
    return(SUN_NLS_MEM_NULL);

  /* initialize jcur */
  jcur = SUNFALSE;

  /* assume the Jacobian is old */
  jbad = SUNTRUE;

  /* shortcut to correction vector */
  delta = NEWTON_CONTENT(NLS)->delta;

  /* initialize counter mnewt */
  NEWTON_CONTENT(NLS)->mnewt = 0;

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
      retval = NEWTON_CONTENT(NLS)->LSetup(y, delta, jbad, &jcur, mem);
      if (retval != SUN_NLS_SUCCESS) break;
    }

    /* compute the negative of the residual for the linear system rhs */
    N_VScale(-ONE, delta, delta);

    /* solve the linear system to get correction vector delta */
    retval = NEWTON_CONTENT(NLS)->LSolve(y, delta, mem);
    if (retval != SUN_NLS_SUCCESS) break;

    /* apply delta to y */
    N_VLinearSum(ONE, y, ONE, delta, y);

    /* test for convergence, return if successful */
    retval = NEWTON_CONTENT(NLS)->CTest(NLS, y, delta, tol, w, mem);
    if (retval == SUN_NLS_SUCCESS)  return(SUN_NLS_SUCCESS);
    if (retval != SUN_NLS_CONTINUE) break;

    /* not yet converged. Increment mnewt and test for max allowed. */
    NEWTON_CONTENT(NLS)->mnewt++;
    if (NEWTON_CONTENT(NLS)->mnewt >= NEWTON_CONTENT(NLS)->maxiters) {
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

  /* free the nonlinear solver */
  free(NLS);
  NLS = NULL;

  return(SUN_NLS_SUCCESS);
}


/* =============================================================================
 * Set functions
 * ===========================================================================*/

int SUNNonlinSolSetSysFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* check that the nonlinear system function is non-null */
  if (SysFn == NULL)
    return(SUN_NLS_ILL_INPUT);

  NEWTON_CONTENT(NLS)->Sys = SysFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetLSetupFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn LSetupFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  NEWTON_CONTENT(NLS)->LSetup = LSetupFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetLSolveFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn LSolveFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* check that the linear solver solve function is non-null */
  if (LSolveFn == NULL)
    return(SUN_NLS_ILL_INPUT);

  NEWTON_CONTENT(NLS)->LSolve = LSolveFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetConvTestFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* check that the convergence test function is non-null */
  if (CTestFn == NULL)
    return(SUN_NLS_ILL_INPUT);

  NEWTON_CONTENT(NLS)->CTest = CTestFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetMaxIters_FullNewton(SUNNonlinearSolver NLS, int maxiters)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* check that maxiters is a vaild */
  if (maxiters < 1)
    return(SUN_NLS_ILL_INPUT);

  NEWTON_CONTENT(NLS)->maxiters = maxiters;
  return(SUN_NLS_SUCCESS);
}


/* =============================================================================
 * Get functions
 * ===========================================================================*/

int SUNNonlinSolGetNumIters_FullNewton(SUNNonlinearSolver NLS, long int *niters)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* return the total number of nonlinear iterations */
  *niters = NEWTON_CONTENT(NLS)->niters;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolGetCurIter_FullNewton(SUNNonlinearSolver NLS, int *iter)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* return the current nonlinear solver iteration count */
  *iter = NEWTON_CONTENT(NLS)->mnewt;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolGetSysFn_FullNewton(SUNNonlinearSolver NLS, SUNNonlinSolSysFn *SysFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* return the nonlinear system defining function */
  *SysFn = NEWTON_CONTENT(NLS)->Sys;
  return(SUN_NLS_SUCCESS);
}
