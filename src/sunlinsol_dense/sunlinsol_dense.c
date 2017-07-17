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
  sunindextype MatrixRows, MatrixCols, VecLength;
  
  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
    return NULL;
  MatrixRows = SUNDenseMatrix_Rows(A);
  MatrixCols = SUNDenseMatrix_Columns(A);
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL) {
    VecLength = N_VGetLength_Serial(y);
  }
#ifdef SUNDIALS_OPENMP_ENABLED
  else if (N_VGetVectorID(y) == SUNDIALS_NVEC_OPENMP) {
    VecLength = N_VGetLength_OpenMP(y);
  }
#endif
#ifdef SUNDIALS_PTHREADS_ENABLED
  else if (N_VGetVectorID(y) == SUNDIALS_NVEC_PTHREADS) {
    VecLength = N_VGetLength_Pthreads(y);
  }
#endif
  else
    return NULL;
  if ( (MatrixRows != MatrixCols) || (MatrixRows != VecLength) )
    return NULL;
  
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
  ops->numpsolves        = SUNLinSolNumPSolves_Dense;
  ops->lastflag          = SUNLinSolLastFlag_Dense;
  ops->free              = SUNLinSolFree_Dense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Dense) malloc(sizeof(struct _SUNLinearSolverContent_Dense));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
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
  return SUNLINEARSOLVER_DIRECT;
}

int SUNLinSolInitialize_Dense(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  SLS_LASTFLAG_D(S) = 0;
  return 0;
}

int SUNLinSolSetATimes_Dense(SUNLinearSolver S, void* A_data,
                             ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* direct solvers do not utilize an 'ATimes' routine, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_D(S) = 1;
  return 1;
}

int SUNLinSolSetPreconditioner_Dense(SUNLinearSolver S, void* P_data,
                                     PSetupFn Pset, PSolveFn Psol)
{
  /* direct solvers do not utilize preconditioning, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_D(S) = 1;
  return 1;
}

int SUNLinSolSetScalingVectors_Dense(SUNLinearSolver S, N_Vector s1,
                                     N_Vector s2)
{
  /* direct solvers do not utilize scaling, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_D(S) = 1;
  return 1;
}

int SUNLinSolSetup_Dense(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  sunindextype *pivots;
  
  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SUNDenseMatrix_Cols(A);
  pivots = SLS_PIVOTS_D(S);
  if ( (A_cols == NULL) || (pivots == NULL) )
    return 1;
  
  /* perform LU factorization of input matrix */
  SLS_LASTFLAG_D(S) = denseGETRF(A_cols, SUNDenseMatrix_Rows(A),
                                 SUNDenseMatrix_Columns(A), pivots);

  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  if (SLS_LASTFLAG_D(S) > 0)
    return 1;
  return 0;
}

int SUNLinSolSolve_Dense(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                        N_Vector b, realtype tol)
{
  realtype **A_cols, *xdata;
  sunindextype *pivots;
  
  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  xdata = NULL;
  pivots = NULL;
  A_cols = SUNDenseMatrix_Cols(A);
  xdata = N_VGetArrayPointer(x);
  pivots = SLS_PIVOTS_D(S);
  if ( (A_cols == NULL) || (xdata == NULL)  || (pivots == NULL) )
    return 1;
  
  /* solve using LU factors */
  denseGETRS(A_cols, SUNDenseMatrix_Rows(A), pivots, xdata);
  return 0;
}

int SUNLinSolNumIters_Dense(SUNLinearSolver S)
{
  /* direct solvers do not perform 'iterations' */
  return 0;
}

realtype SUNLinSolResNorm_Dense(SUNLinearSolver S)
{
  /* direct solvers do not measure the linear residual */
  return ZERO;
}

int SUNLinSolNumPSolves_Dense(SUNLinearSolver S)
{
  /* direct solvers do not use preconditioning */
  return 0;
}

long int SUNLinSolLastFlag_Dense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return (SLS_LASTFLAG_D(S));
}

int SUNLinSolFree_Dense(SUNLinearSolver S)
{
  /* delete items from the contents structure, then delete generic structures */
  free(SLS_PIVOTS_D(S));  SLS_PIVOTS_D(S) = NULL;
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}

