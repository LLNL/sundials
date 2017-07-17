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
 * This is the implementation file for the band implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_band.h>
#include <sundials/sundials_math.h>

#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new band linear solver
 */

SUNLinearSolver SUNBandLinearSolver(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_Band content;
  sunindextype MatrixRows, MatrixCols, VecLength;
  
  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_BAND)
    return NULL;
  MatrixRows = SUNBandMatrix_Rows(A);
  MatrixCols = SUNBandMatrix_Columns(A);
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL)
    VecLength = N_VGetLength_Serial(y);
  else if (N_VGetVectorID(y) == SUNDIALS_NVEC_OPENMP)
    VecLength = N_VGetLength_OpenMP(y);
  else if (N_VGetVectorID(y) == SUNDIALS_NVEC_PTHREADS)
    VecLength = N_VGetLength_Pthreads(y);
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
  ops->gettype           = SUNLinSolGetType_Band;
  ops->setatimes         = SUNLinSolSetATimes_Band;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_Band;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_Band;
  ops->initialize        = SUNLinSolInitialize_Band;
  ops->setup             = SUNLinSolSetup_Band;
  ops->solve             = SUNLinSolSolve_Band;
  ops->numiters          = SUNLinSolNumIters_Band;
  ops->resnorm           = SUNLinSolResNorm_Band;
  ops->numpsolves        = SUNLinSolNumPSolves_Band;
  ops->lastflag          = SUNLinSolLastFlag_Band;
  ops->free              = SUNLinSolFree_Band;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Band) malloc(sizeof(struct _SUNLinearSolverContent_Band));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
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

SUNLinearSolver_Type SUNLinSolGetType_Band(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

int SUNLinSolInitialize_Band(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  SLS_LASTFLAG_B(S) = 0;
  return 0;
}

int SUNLinSolSetATimes_Band(SUNLinearSolver S, void* A_data, 
                            ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* direct solvers do not utilize an 'ATimes' routine, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_B(S) = 1;
  return 1;
}

int SUNLinSolSetPreconditioner_Band(SUNLinearSolver S, void* P_data,
                                    PSetupFn Pset, PSolveFn Psol)
{
  /* direct solvers do not utilize preconditioning, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_B(S) = 1;
  return 1;
}

int SUNLinSolSetScalingVectors_Band(SUNLinearSolver S, N_Vector s1, 
                                    N_Vector s2)
{
  /* direct solvers do not utilize scaling, 
     so return an error is this routine is ever called */
  SLS_LASTFLAG_B(S) = 1;
  return 1;
}

int SUNLinSolSetup_Band(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  sunindextype *pivots;
  
  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SM_COLS_B(A);
  pivots = SLS_PIVOTS_B(S);
  if ( (A_cols == NULL) || (pivots == NULL) )
    return 1;

  /* ensure that storage upper bandwidth is sufficient for fill-in */
  if (SM_SUBAND_B(A) < SUNMIN(SM_COLUMNS_B(A)-1, SM_UBAND_B(A) + SM_LBAND_B(A))) {
    SLS_LASTFLAG_B(S) = -1;
    return 1;
  }
  
  /* perform LU factorization of input matrix */
  SLS_LASTFLAG_B(S) = bandGBTRF(A_cols, SM_COLUMNS_B(A), SM_UBAND_B(A),
                                SM_LBAND_B(A), SM_SUBAND_B(A), pivots);
  
  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  if (SLS_LASTFLAG_B(S) > 0)
    return 1;
  return 0;
}

int SUNLinSolSolve_Band(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
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
  A_cols = SUNBandMatrix_Cols(A);
  xdata = N_VGetArrayPointer(x);
  pivots = SLS_PIVOTS_B(S);
  if ( (A_cols == NULL) || (xdata == NULL)  || (pivots == NULL) )
    return 1;
  
  /* solve using LU factors */
  bandGBTRS(A_cols, SM_COLUMNS_B(A), SM_SUBAND_B(A), SM_LBAND_B(A), pivots, xdata);
  return 0;
}

int SUNLinSolNumIters_Band(SUNLinearSolver S)
{
  /* direct solvers do not perform 'iterations' */
  return 0;
}

realtype SUNLinSolResNorm_Band(SUNLinearSolver S)
{
  /* direct solvers do not measure the linear residual */
  return ZERO;
}

int SUNLinSolNumPSolves_Band(SUNLinearSolver S)
{
  /* direct solvers do not use preconditioning */
  return 0;
}

long int SUNLinSolLastFlag_Band(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return (SLS_LASTFLAG_B(S));
}

int SUNLinSolFree_Band(SUNLinearSolver S)
{
  /* delete items from the contents structure, then delete generic structures */
  free(SLS_PIVOTS_B(S));  SLS_PIVOTS_B(S) = NULL;
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}
