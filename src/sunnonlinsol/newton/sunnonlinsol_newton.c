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
 * implementation of Newton's method.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0) /* real 0.0 */
#define ONE  RCONST(1.0) /* real 1.0 */

/* Content structure accessibility macros  */
#define NEWTON_CONTENT(S) ( (SUNNonlinearSolverContent_Newton)(S->content) )

/* =============================================================================
 * Constructor to create a new Newton solver
 * ===========================================================================*/

SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y)
{
  SUNNonlinearSolver NLS;
  SUNNonlinearSolver_Ops ops;
  SUNNonlinearSolverContent_Newton content;
  
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
  ops->gettype     = SUNNonlinSolGetType_Newton;
  ops->initialize  = SUNNonlinSolInitialize_Newton;
  ops->setup       = SUNNonlinSolSetup_Newton;
  ops->solve       = SUNNonlinSolSolve_Newton;
  ops->free        = SUNNonlinSolFree_Newton;
  ops->setsysfn    = SUNNonlinSolSetSysFn_Newton;
  ops->setlsetupfn = SUNNonlinSolSetLSetupFn_Newton;
  ops->setlsolvefn = SUNNonlinSolSetLSolveFn_Newton;
  ops->setctestfn  = SUNNonlinSolSetConvTestFn_Newton;
  ops->setmaxiters = SUNNonlinSolSetMaxIters_Newton;
  ops->getnumiters = SUNNonlinSolGetNumIters_Newton;
  /* ops->lastflag          = SUNNonlinSolLastFlag_Newton; */

  /* Create content */
  content = NULL;
  content = (SUNNonlinearSolverContent_Newton) malloc(sizeof *content);
  if (content == NULL) { free(ops); free(NLS); return(NULL); }

  /* Fill content */
  content->Sys      = NULL;
  content->LSetup   = NULL;
  content->LSolve   = NULL;
  content->CTest    = NULL;
  content->delta    = N_VClone(y);
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

SUNNonlinearSolver_Type SUNNonlinSolGetType_Newton(SUNNonlinearSolver NLS)
{
  return(SUNNONLINEARSOLVER_ROOTFIND);
}


int SUNNonlinSolInitialize_Newton(SUNNonlinearSolver NLS)
{
  /* all solver-specific memory has already been allocated */
  if (NLS == NULL)
    return(SUN_NLS_MEM_NULL);

  /* check that all required function pointers have been set */
  if ( NEWTON_CONTENT(NLS)->Sys    == NULL ||
       NEWTON_CONTENT(NLS)->LSolve == NULL ||
       NEWTON_CONTENT(NLS)->CTest  == NULL ) {
    return(SUN_NLS_MEM_NULL);
  }

  /* initialize the total number of nonlinear solver iterations */
  NEWTON_CONTENT(NLS)->niters = 0;

  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetup_Newton(SUNNonlinearSolver NLS, N_Vector y, void* mem)
{
  /* no setup necessary */
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
int SUNNonlinSolSolve_Newton(SUNNonlinearSolver NLS,
                             N_Vector y0, N_Vector y,
                             N_Vector w, realtype tol,
                             booleantype callSetup, void* mem)
{
  int mnewt;
  int retval;
  realtype delnrm;
  N_Vector delta;
  booleantype tryAgain;

  /* shortcut to correction vector */
  delta = NEWTON_CONTENT(NLS)->delta;

  /* looping point for Jacobian/preconditioner setup attempts */
  for(;;) {

      /* compute the residual */
      retval = NEWTON_CONTENT(NLS)->Sys(y0, delta, mem);
      if (retval != SUN_NLS_SUCCESS) break;

      /* if indicated, setup the linear system */
      if (callSetup) {
        retval = NEWTON_CONTENT(NLS)->LSetup(y0, delta, mem);
        if (retval != SUN_NLS_SUCCESS) break;
      }

      /* initialize counter mnewt */
      mnewt = 0;

      /* load prediction into y */
      N_VScale(ONE, y0, y);

      /* looping point for Newton iteration. Break out on any error. */
      for(;;) {

        /* increment number of nonlinear solver iterations */
        NEWTON_CONTENT(NLS)->niters++;

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

        /* call res for new residual and check error flag from res. */
        retval = NEWTON_CONTENT(NLS)->Sys(y, delta, mem);
        if (retval != SUN_NLS_SUCCESS) break;

      } /* end of Newton iteration loop */

      /* all error go here */
      tryAgain = (retval > 0) && (NEWTON_CONTENT(NLS)->LSetup) && (!callSetup);
      
      if (tryAgain)
        callSetup = SUNTRUE;
      else
        break;

  } /* end of setup loop */

  /* all error returns exit here */
  return(retval);
}


int SUNNonlinSolFree_Newton(SUNNonlinearSolver NLS)
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

int SUNNonlinSolSetSysFn_Newton(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)
{
  NEWTON_CONTENT(NLS)->Sys = SysFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetLSetupFn_Newton(SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn LSetupFn)
{
  NEWTON_CONTENT(NLS)->LSetup = LSetupFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetLSolveFn_Newton(SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn LSolveFn)
{
  NEWTON_CONTENT(NLS)->LSolve = LSolveFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetConvTestFn_Newton(SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn)
{
  NEWTON_CONTENT(NLS)->CTest = CTestFn;
  return(SUN_NLS_SUCCESS);
}


int SUNNonlinSolSetMaxIters_Newton(SUNNonlinearSolver NLS, int maxiters)
{
  if (maxiters < 1) return(SUN_NLS_ILL_INPUT);
  NEWTON_CONTENT(NLS)->maxiters = maxiters;
  return(SUN_NLS_SUCCESS);
}


/* =============================================================================
 * Get functions
 * ===========================================================================*/

int SUNNonlinSolGetNumIters_Newton(SUNNonlinearSolver NLS, long int *niters)
{
  *niters = NEWTON_CONTENT(NLS)->niters;
  return(SUN_NLS_SUCCESS);
}
