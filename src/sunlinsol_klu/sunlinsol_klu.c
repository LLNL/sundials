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
 * This is the implementation file for the KLU implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_klu.h>
#include <sundials/sundials_math.h>

#define ONE RCONST(1.0)
#define TWO RCONST(2.0)
#define TWOTHIRDS RCONST(0.666666666666666666666666666666667)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new KLU linear solver
 */

SUNLinearSolver SUNKLULinearSolver(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_KLU content;
  long int MatrixRows, MatrixCols, VecLength;
  
  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return NULL;
  MatrixRows = SUNSparseMatrix_Rows(A);
  MatrixCols = SUNSparseMatrix_Columns(A);
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
  ops->gettype           = SUNLinSolGetType_KLU;
  ops->setatimes         = SUNLinSolSetATimes_KLU;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_KLU;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_KLU;
  ops->initialize        = SUNLinSolInitialize_KLU;
  ops->setup             = SUNLinSolSetup_KLU;
  ops->solve             = SUNLinSolSolve_KLU;
  ops->numiters          = SUNLinSolNumIters_KLU;
  ops->resnorm           = SUNLinSolResNorm_KLU;
  ops->numpsolves        = SUNLinSolNumPSolves_KLU;
  ops->lastflag          = SUNLinSolLastFlag_KLU;
  ops->free              = SUNLinSolFree_;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_KLU) malloc(sizeof(struct _SUNLinearSolverContent_KLU));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  
  /* Attach content and ops */
  S->content = content;
  S->ops     = ops;

  return(S);
}

int SUNKLUReInit(SUNLinearSolver S, SUNMatrix A,
                 long int nnz, int reinit_type)
{
  
}

int SUNKLUSetOrdering(SUNLinearSolver S, int ordering_choice)
{

}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_KLU(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

int SUNLinSolInitialize_KLU(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  SLS_LASTFLAG_K(S) = 0;
  return 0;
}

int SUNLinSolSetATimes_KLU(SUNLinearSolver S, void* A_data, 
                           ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* direct solvers do not utilize an 'ATimes' routine, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_K(S) = 1;
  return 1;
}

int SUNLinSolSetPreconditioner_KLU(SUNLinearSolver S, void* P_data,
                                    PSetupFn Pset, PSolveFn Psol)
{
  /* direct solvers do not utilize preconditioning, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_K(S) = 1;
  return 1;
}

int SUNLinSolSetScalingVectors_KLU(SUNLinearSolver S, N_Vector s1,
                                   N_Vector s2)
{
  /* direct solvers do not utilize scaling, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_K(S) = 1;
  return 1;
}

int SUNLinSolSetup_KLU(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  long int *pivots;
  
  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SM_COLS_K(A);
  if (A_cols == NULL)
    return 1;

  /* perform LU factorization of input matrix */
  SLS_LASTFLAG_K(S) = bandGBTRF(A_cols, SM_COLUMNS_B(A), SM_UBAND_B(A),
                                SM_LBAND_B(A), SM_SUBAND_B(A), pivots);
  
  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  if (SLS_LASTFLAG_K(S) > 0)
    return 1;
  return 0;
}

int SUNLinSolSolve_KLU(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                       N_Vector b, realtype tol)
{
  realtype **A_cols, *xdata;
  long int *pivots;
  
  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  xdata = NULL;
  pivots = NULL;
  A_cols = SUNSparseMatrix_Cols(A);
  xdata = N_VGetArrayPointer(x);
  pivots = SLS_PIVOTS_B(S);
  if ( (A_cols == NULL) || (xdata == NULL)  || (pivots == NULL) )
    return 1;
  
  /* solve using LU factors */
  bandGBTRS(A_cols, SM_COLUMNS_B(A), SM_SUBAND_B(A), SM_LBAND_B(A), pivots, xdata);
  return 0;
}

int SUNLinSolNumIters_KLU(SUNLinearSolver S)
{
  /* direct solvers do not perform 'iterations' */
  return 0;
}

realtype SUNLinSolResNorm_KLU(SUNLinearSolver S)
{
  /* direct solvers do not measure the linear residual */
  return ZERO;
}

int SUNLinSolNumPSolves_KLU(SUNLinearSolver S)
{
  /* direct solvers do not use preconditioning */
  return 0;
}

long int SUNLinSolLastFlag_KLU(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return (SLS_LASTFLAG_K(S));
}

int SUNLinSolFree_KLU(SUNLinearSolver S)
{
  /* delete items from the contents structure, then delete generic structures */
  /* free(SLS_PIVOTS_K(S));  SLS_PIVOTS_K(S) = NULL; */
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}
