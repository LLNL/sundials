/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This is the implementation file for the diagonal implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_diagonal.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(1.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Diagonal solver structure accessibility macros: 
 * -----------------------------------------------------------------
 */

#define DIAGONAL_CONTENT(S)  ( (SUNLinearSolverContent_Diagonal)(S->content) )
#define LASTFLAG(S)          ( DIAGONAL_CONTENT(S)->last_flag )


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new diagonal linear solver
 */

SUNLinearSolver SUNDiagonalLinearSolver(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_Diagonal content;
  
  /* Check compatibility of supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_DIAGONAL)
    return(NULL);
  if (N_VGetVectorID(y) != N_VGetVectorID(SUNDiagonalMatrix_Diag(A)))
    return(NULL);

  /* verify that supplied N_Vector implements all required operations */
  if (y->ops->nvinvtest == NULL)  return(NULL);
  
  /* It would be great if we could verify that the dimensions of A 
     and y match, but since there is no generic 'length' N_Vector 
     routine we cannot do so. */
  
  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops)
    malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_Diagonal;
  ops->setatimes         = SUNLinSolSetATimes_Diagonal;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_Diagonal;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_Diagonal;
  ops->initialize        = SUNLinSolInitialize_Diagonal;
  ops->setup             = SUNLinSolSetup_Diagonal;
  ops->solve             = SUNLinSolSolve_Diagonal;
  ops->numiters          = SUNLinSolNumIters_Diagonal;
  ops->resnorm           = SUNLinSolResNorm_Diagonal;
  ops->lastflag          = SUNLinSolLastFlag_Diagonal;
  ops->space             = SUNLinSolSpace_Diagonal;
  ops->free              = SUNLinSolFree_Diagonal;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Diagonal)
    malloc(sizeof(struct _SUNLinearSolverContent_Diagonal));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  
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

SUNLinearSolver_Type SUNLinSolGetType_Diagonal(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

int SUNLinSolInitialize_Diagonal(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

int SUNLinSolSetATimes_Diagonal(SUNLinearSolver S, void* A_data,
                                ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* direct solvers do not utilize an 'ATimes' routine, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}

int SUNLinSolSetPreconditioner_Diagonal(SUNLinearSolver S, void* P_data,
                                        PSetupFn Pset, PSolveFn Psol)
{
  /* direct solvers do not utilize preconditioning, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}

int SUNLinSolSetScalingVectors_Diagonal(SUNLinearSolver S, N_Vector s1,
                                        N_Vector s2)
{
  /* direct solvers do not utilize scaling, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}

int SUNLinSolSetup_Diagonal(SUNLinearSolver S, SUNMatrix A)
{
  N_Vector A_diag;
  booleantype retval;
  
  /* check for valid inputs */
  if ( (A == NULL) || (S == NULL) ) 
    return(SUNLS_MEM_NULL);
  
  /* Ensure that A is a diagonal matrix */
  if (SUNMatGetID(A) != SUNMATRIX_DIAGONAL) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(LASTFLAG(S));
  }
  
  /* access diagonal of A (return with failure on NULL) */
  A_diag = NULL;
  A_diag = SUNDiagonalMatrix_Diag(A);
  if (A_diag == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  
  /* perform 'factorization' by inverting diagonal in-place */
  retval = N_VInvTest(A_diag, A_diag);

  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  LASTFLAG(S) = (retval) ? 0 : 1;
  if (!retval)
    return(SUNLS_LUFACT_FAIL);
  else 
    return(SUNLS_SUCCESS);
}

int SUNLinSolSolve_Diagonal(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                            N_Vector b, realtype tol)
{
  /* check for valid inputs and access inverted matrix diagonal */
  N_Vector A_diag_inv;
  if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) ) 
    return(SUNLS_MEM_NULL);
  A_diag_inv = NULL;
  A_diag_inv = SUNDiagonalMatrix_Diag(A);
  if (A_diag_inv == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  
  /* 'solve' is just x = Ainv.*b (MATLAB notation) */
  N_VProd(A_diag_inv, b, x);

  /* return with success */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

int SUNLinSolNumIters_Diagonal(SUNLinearSolver S)
{
  /* direct solvers do not perform 'iterations' */
  return(0);
}

realtype SUNLinSolResNorm_Diagonal(SUNLinearSolver S)
{
  /* direct solvers do not measure the linear residual */
  return(ZERO);
}

long int SUNLinSolLastFlag_Diagonal(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return(LASTFLAG(S));
}

int SUNLinSolSpace_Diagonal(SUNLinearSolver S, 
                            long int *lenrwLS, 
                            long int *leniwLS)
{
  *leniwLS = 1;
  *lenrwLS = 0;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_Diagonal(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL)
    return(SUNLS_SUCCESS);

  /* delete items from contents, then delete generic structure */
  if (S->content) {
    free(S->content);  
    S->content = NULL;
  }
  if (S->ops) {
    free(S->ops);  
    S->ops = NULL;
  }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}

