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
 * This is the implementation file for the SPGMR implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_math.h>

#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPGMR linear solver
 */

SUNLinearSolver SUNSPGMR(N_Vector y, int pretype, int maxl)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_SPGMR content;
  
  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_SPGMR;
  ops->setatimes         = SUNLinSolSetATimes_SPGMR;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_SPGMR;
  ops->initialize        = SUNLinSolInitialize_SPGMR;
  ops->setup             = SUNLinSolSetup_SPGMR;
  ops->solve             = SUNLinSolSolve_SPGMR;
  ops->liniters          = SUNLinSolLinIters_SPGMR;
  ops->lastflag          = SUNLinSolLastFlag_SPGMR;  
  ops->free              = SUNLinSolFree_SPGMR;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPGMR) malloc(sizeof(struct _SUNLinearSolverContent_SPGMR));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  content->maxl = maxl;
  content->pretype = pretype;
  content->ATSetup = NULL;
  content->ATimes = NULL;
  content->ATData = NULL;
  content->Psetup = NULL;
  content->Psolve = NULL;
  content->PData = NULL;
  
  /* Attach content and ops */
  S->content = content;
  S->ops     = ops;

  return(S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPGMR(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_ITERATIVE;
}

int SUNLinSolInitialize_SPGMR(SUNLinearSolver S)
{
  /* allocate solver-specific memory here (Krylov subspace vectors) */


  /* return with success */
  SLS_CONTENT_SPGMR(S)->last_flag = 0;
  return 0;
}

int SUNLinSolSetATimes_SPGMR(SUNLinearSolver S, void* ATData, 
                            ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATSetup 
     and ATimes routines and data, and return with success */
  SLS_CONTENT_SPGMR(S)->ATSetup = ATSetup;
  SLS_CONTENT_SPGMR(S)->ATimes = ATimes;
  SLS_CONTENT_SPGMR(S)->ATData = ATData;
  SLS_CONTENT_SPGMR(S)->last_flag = 0;
  return 0;
}

int SUNLinSolSetPreconditioner_SPGMR(SUNLinearSolver S, void* PData,
                                    PSetupFn PSetup, PSolveFn PSolve)
{
  /* set function pointers to integrator-supplied PSetup and PSolve
     routines and data, and return with success */
  SLS_CONTENT_SPGMR(S)->PSetup = PSetup;
  SLS_CONTENT_SPGMR(S)->PSolve = PSolve;
  SLS_CONTENT_SPGMR(S)->PData = PData;
  SLS_CONTENT_SPGMR(S)->last_flag = 0;
  return 0;
}

int SUNLinSolSetup_SPGMR(SUNLinearSolver S, SUNMatrix A, N_Vector tmp1, 
                         N_Vector tmp2, N_Vector tmp3)
{
  int ierr;
  ATSetupFn ATSetup = SLS_CONTENT_SPGMR(S)->ATSetup;
  PSetupFn PSetup = SLS_CONTENT_SPGMR(S)->PSetup;
  void* ATData = SLS_CONTENT_SPGMR(S)->ATData;
  void* PData = SLS_CONTENT_SPGMR(S)->PData;;
  
  /* no solver-specific setup is required, but if user-supplied 
     ATSetup or PSetup routines exist, call those here */
  if (ATSetup != NULL) {
    SLS_CONTENT_SPGMR(S)->last_flag = ATSetup(ATData);
    if (SLS_CONTENT_SPGMR(S)->last_flag != 0)  return 1;
  }
  if (PSetup != NULL) {
    SLS_CONTENT_SPGMR(S)->last_flag = PSetup(PData);
    if (SLS_CONTENT_SPGMR(S)->last_flag != 0)  return 1;
  }
  
  /* return with success */ 
  return 0;
}

int SUNLinSolSolve_SPGMR(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                         N_Vector b, N_Vector w, realtype tol)
{
   
}

long int SUNLinSolLinIters_SPGMR(SUNLinearSolver S)
{

}

long int SUNLinSolLastFlag_SPGMR(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return (SLS_CONTENT_SPGMR(S)->last_flag);
}

int SUNLinSolFree_SPGMR(SUNLinearSolver S)
{
  /* delete items from the contents structure, then delete generic structures */
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}
