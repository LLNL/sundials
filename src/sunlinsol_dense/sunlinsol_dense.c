/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the implementation file for the dense implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(1.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Dense solver structure accessibility macros: 
 * -----------------------------------------------------------------
 */

#define DENSE_CONTENT(S)  ( (SUNLinearSolverContent_Dense)(S->content) )
#define PIVOTS(S)         ( DENSE_CONTENT(S)->pivots )
#define LASTFLAG(S)       ( DENSE_CONTENT(S)->last_flag )


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new dense linear solver
 */

SUNLinearSolver SUNDenseLinearSolver(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_Dense content;
  sunindextype MatrixRows;
  
  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
    return(NULL);
  if (SUNDenseMatrix_Rows(A) != SUNDenseMatrix_Columns(A))
    return(NULL);
  MatrixRows = SUNDenseMatrix_Rows(A);
  if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
    return(NULL);

  /* Optimally we would verify that the dimensions of A and y agree, but 
   since there is no generic 'length' routine for N_Vectors we cannot */
  
  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_Dense;
  ops->setatimes         = SUNLinSolSetATimes_Dense;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_Dense;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_Dense;
  ops->initialize        = SUNLinSolInitialize_Dense;
  ops->setup             = SUNLinSolSetup_Dense;
  ops->solve             = SUNLinSolSolve_Dense;
  ops->numiters          = SUNLinSolNumIters_Dense;
  ops->resnorm           = SUNLinSolResNorm_Dense;
  ops->lastflag          = SUNLinSolLastFlag_Dense;
  ops->space             = SUNLinSolSpace_Dense;
  ops->free              = SUNLinSolFree_Dense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Dense) malloc(sizeof(struct _SUNLinearSolverContent_Dense));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->N = MatrixRows;
  content->last_flag = 0;
  content->pivots = NULL;
  content->pivots = (sunindextype *) malloc(MatrixRows * sizeof(sunindextype));
  if (content->pivots == NULL) {
    free(content); free(ops); free(S); return(NULL);
  }
  
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

SUNLinearSolver_Type SUNLinSolGetType_Dense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

int SUNLinSolInitialize_Dense(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

int SUNLinSolSetATimes_Dense(SUNLinearSolver S, void* A_data,
                             ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* direct solvers do not utilize an 'ATimes' routine, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}

int SUNLinSolSetPreconditioner_Dense(SUNLinearSolver S, void* P_data,
                                     PSetupFn Pset, PSolveFn Psol)
{
  /* direct solvers do not utilize preconditioning, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}

int SUNLinSolSetScalingVectors_Dense(SUNLinearSolver S, N_Vector s1,
                                     N_Vector s2)
{
  /* direct solvers do not utilize scaling, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}

int SUNLinSolSetup_Dense(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  sunindextype *pivots;
  
  /* check for valid inputs */
  if ( (A == NULL) || (S == NULL) ) 
    return(SUNLS_MEM_NULL);
  
  /* Ensure that A is a dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(LASTFLAG(S));
  }
  
  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SUNDenseMatrix_Cols(A);
  pivots = PIVOTS(S);
  if ( (A_cols == NULL) || (pivots == NULL) ) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  
  /* perform LU factorization of input matrix */
  LASTFLAG(S) = denseGETRF(A_cols, SUNDenseMatrix_Rows(A),
                           SUNDenseMatrix_Columns(A), pivots);

  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  if (LASTFLAG(S) > 0)
    return(SUNLS_LUFACT_FAIL);
  return(SUNLS_SUCCESS);
}

int SUNLinSolSolve_Dense(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                        N_Vector b, realtype tol)
{
  realtype **A_cols, *xdata;
  sunindextype *pivots;
  
  if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) ) 
    return(SUNLS_MEM_NULL);
  
  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  xdata = NULL;
  pivots = NULL;
  A_cols = SUNDenseMatrix_Cols(A);
  xdata = N_VGetArrayPointer(x);
  pivots = PIVOTS(S);
  if ( (A_cols == NULL) || (xdata == NULL)  || (pivots == NULL) ) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  
  /* solve using LU factors */
  denseGETRS(A_cols, SUNDenseMatrix_Rows(A), pivots, xdata);
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

int SUNLinSolNumIters_Dense(SUNLinearSolver S)
{
  /* direct solvers do not perform 'iterations' */
  return(0);
}

realtype SUNLinSolResNorm_Dense(SUNLinearSolver S)
{
  /* direct solvers do not measure the linear residual */
  return(ZERO);
}

long int SUNLinSolLastFlag_Dense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return(LASTFLAG(S));
}

int SUNLinSolSpace_Dense(SUNLinearSolver S, 
                         long int *lenrwLS, 
                         long int *leniwLS)
{
  *leniwLS = 2 + DENSE_CONTENT(S)->N;
  *lenrwLS = 0;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_Dense(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL)
    return(SUNLS_SUCCESS);

  /* delete items from contents, then delete generic structure */
  if (S->content) {
    if (PIVOTS(S)) {
      free(PIVOTS(S));
      PIVOTS(S) = NULL;
    }
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

