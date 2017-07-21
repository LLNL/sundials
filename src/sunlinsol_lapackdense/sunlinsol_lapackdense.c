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
 * This is the implementation file for the LAPACK dense 
 * implementation of the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_lapackdense.h>
#include <sundials/sundials_math.h>

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Dense solver structure accessibility macros: 
 * -----------------------------------------------------------------
 */

#define LAPACKDENSE_CONTENT(S) ( (SUNLinearSolverContent_LapackDense)(S->content) )
#define PIVOTS(S)              ( LAPACKDENSE_CONTENT(S)->pivots )
#define LASTFLAG(S)            ( LAPACKDENSE_CONTENT(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new LAPACK dense linear solver
 */

SUNLinearSolver SUNLapackDense(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_LapackDense content;
  sunindextype MatrixRows, MatrixCols, VecLength;
  
  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
    return(NULL);
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
    return(NULL);
  if ( (MatrixRows != MatrixCols) || (MatrixRows != VecLength) )
    return(NULL);

  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_LapackDense;
  ops->setatimes         = SUNLinSolSetATimes_LapackDense;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_LapackDense;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_LapackDense;
  ops->initialize        = SUNLinSolInitialize_LapackDense;
  ops->setup             = SUNLinSolSetup_LapackDense;
  ops->solve             = SUNLinSolSolve_LapackDense;
  ops->numiters          = SUNLinSolNumIters_LapackDense;
  ops->resnorm           = SUNLinSolResNorm_LapackDense;
  ops->numpsolves        = SUNLinSolNumPSolves_LapackDense;
  ops->lastflag          = SUNLinSolLastFlag_LapackDense;
  ops->free              = SUNLinSolFree_LapackDense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_LapackDense)
    malloc(sizeof(struct _SUNLinearSolverContent_LapackDense));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  content->pivots = NULL;
  content->pivots = (int *) malloc(MatrixRows * sizeof(int));
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

SUNLinearSolver_Type SUNLinSolGetType_LapackDense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


int SUNLinSolInitialize_LapackDense(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetATimes_LapackDense(SUNLinearSolver S, void* A_data, 
                                   ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* direct solvers do not utilize an 'ATimes' routine, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}


int SUNLinSolSetPreconditioner_LapackDense(SUNLinearSolver S, void* P_data,
                                           PSetupFn Pset, PSolveFn Psol)
{
  /* direct solvers do not utilize preconditioning, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}


int SUNLinSolSetScalingVectors_LapackDense(SUNLinearSolver S, N_Vector s1,
                                           N_Vector s2)
{
  /* direct solvers do not utilize scaling, 
     so return an error is this routine is ever called */
  LASTFLAG(S) = SUNLS_ILL_INPUT;
  return(LASTFLAG(S));
}


int SUNLinSolSetup_LapackDense(SUNLinearSolver S, SUNMatrix A)
{
  int n, ier;

  /* check for valid inputs */
  if ( (A == NULL) || (S == NULL) ) 
    return(SUNLS_MEM_NULL);
  
  /* Ensure that A is a dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(LASTFLAG(S));
  }
  
  /* Call LAPACK to do LU factorization of A */
  n = SUNDenseMatrix_Rows(A);
  xgetrf_f77(&n, &n, SUNDenseMatrix_Data(A), &n, PIVOTS(S), &ier);
  LASTFLAG(S) = (long int) ier;
  if (ier > 0) 
    return(SUNLS_LUFACT_FAIL);
  if (ier < 0) 
    return(SUNLS_PACKAGE_FAIL_UNREC);
  return(SUNLS_SUCCESS);
}


int SUNLinSolSolve_LapackDense(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                              N_Vector b, realtype tol)
{
  int n, one, ier;
  realtype *xdata;
  
  if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) ) 
    return(SUNLS_MEM_NULL);

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  
  /* Call LAPACK to solve the linear system */
  n = SUNDenseMatrix_Rows(A);
  one = 1;
  xgetrs_f77("N", &n, &one, SUNDenseMatrix_Data(A), 
	     &n, PIVOTS(S), xdata, &n, &ier, 1);
  LASTFLAG(S) = (long int) ier;
  if (ier < 0) 
    return(SUNLS_PACKAGE_FAIL_UNREC);

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolNumIters_LapackDense(SUNLinearSolver S)
{
  /* direct solvers do not perform 'iterations' */
  return(0);
}


realtype SUNLinSolResNorm_LapackDense(SUNLinearSolver S)
{
  /* direct solvers do not measure the linear residual */
  return(ZERO);
}


int SUNLinSolNumPSolves_LapackDense(SUNLinearSolver S)
{
  /* direct solvers do not use preconditioning */
  return(0);
}


long int SUNLinSolLastFlag_LapackDense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return(LASTFLAG(S));
}


int SUNLinSolFree_LapackDense(SUNLinearSolver S)
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
