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

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define ROW(i,j,smu) (i-j+smu)

/* private functions */
long int sunlinsolbandGBTRF(realtype **a, long int n, long int mu, 
			    long int ml, long int smu, long int *p);
void sunlinsolbandGBTRS(realtype **a, long int n, long int smu, 
			long int ml, long int *p, realtype *b);


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
  long int MatrixRows, MatrixCols, VecLength;
  
  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_BAND)
    return NULL;
  MatrixRows = SUNBandMatrix_Rows(A);
  MatrixCols = SUNBandMatrix_Columns(A);
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
  ops->gettype           = SUNLinSolGetType_Band;
  ops->setatimes         = SUNLinSolSetATimes_Band;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_Band;
  ops->initialize        = SUNLinSolInitialize_Band;
  ops->setup             = SUNLinSolSetup_Band;
  ops->solve             = SUNLinSolSolve_Band;
  ops->performance       = SUNLinSolPerformance_Band;
  ops->lastflag          = SUNLinSolLastFlag_Band;
  ops->free              = SUNLinSolFree_Band;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Band) malloc(sizeof(struct _SUNLinearSolverContent_Band));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  content->pivots = NULL;
  content->pivots = (long int *) malloc(MatrixRows * sizeof(long int));
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

int SUNLinSolSetup_Band(SUNLinearSolver S, SUNMatrix A)
{
  realtype **A_cols;
  long int *pivots;
  
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
  SLS_LASTFLAG_B(S) = sunlinsolbandGBTRF(A_cols, SM_COLUMNS_B(A), SM_UBAND_B(A),
					 SM_LBAND_B(A), SM_SUBAND_B(A), pivots);
  
  /* store error flag (if nonzero, this row encountered zero-valued pivod) */
  if (SLS_LASTFLAG_B(S) > 0)
    return 1;
  return 0;
}

int SUNLinSolSolve_Band(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                        N_Vector b, N_Vector w, realtype tol)
{
  realtype **A_cols, *xdata;
  long int *pivots;
  
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
  sunlinsolbandGBTRS(A_cols, SM_COLUMNS_B(A), SM_SUBAND_B(A), 
		     SM_LBAND_B(A), pivots, xdata);
  return 0;
}

int SUNLinSolPerformance_Band(SUNLinearSolver S, int perftask)
{
  /* direct solvers do not perform performance monitoring */
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


/*
 Private functions
 */

long int sunlinsolbandGBTRF(realtype **a, long int n, long int mu, 
			    long int ml, long int smu, long int *p)
{
  long int c, r, num_rows;
  long int i, j, k, l, storage_l, storage_k, last_col_k, last_row_k;
  realtype *a_c, *col_k, *diag_k, *sub_diag_k, *col_j, *kptr, *jptr;
  realtype max, temp, mult, a_kj;
  booleantype swap;

  /* zero out the first smu - mu rows of the rectangular array a */

  num_rows = smu - mu;
  if (num_rows > 0) {
    for (c=0; c < n; c++) {
      a_c = a[c];
      for (r=0; r < num_rows; r++) {
	a_c[r] = ZERO;
      }
    }
  }

  /* k = elimination step number */

  for (k=0; k < n-1; k++, p++) {
    
    col_k     = a[k];
    diag_k    = col_k + smu;
    sub_diag_k = diag_k + 1;
    last_row_k = SUNMIN(n-1,k+ml);

    /* find l = pivot row number */

    l=k;
    max = SUNRabs(*diag_k);
    for (i=k+1, kptr=sub_diag_k; i <= last_row_k; i++, kptr++) { 
      if (SUNRabs(*kptr) > max) {
	l=i;
	max = SUNRabs(*kptr);
      }
    }
    storage_l = ROW(l, k, smu);
    *p = l;
    
    /* check for zero pivot element */

    if (col_k[storage_l] == ZERO) return(k+1);
    
    /* swap a(l,k) and a(k,k) if necessary */
    
    if ( (swap = (l != k) )) {
      temp = col_k[storage_l];
      col_k[storage_l] = *diag_k;
      *diag_k = temp;
    }

    /* Scale the elements below the diagonal in         */
    /* column k by -1.0 / a(k,k). After the above swap, */
    /* a(k,k) holds the pivot element. This scaling     */
    /* stores the pivot row multipliers -a(i,k)/a(k,k)  */
    /* in a(i,k), i=k+1, ..., SUNMIN(n-1,k+ml).            */
    
    mult = -ONE / (*diag_k);
    for (i=k+1, kptr = sub_diag_k; i <= last_row_k; i++, kptr++)
      (*kptr) *= mult;

    /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., SUNMIN(n-1,k+ml) */
    /* row k is the pivot row after swapping with row l.                */
    /* The computation is done one column at a time,                    */
    /* column j=k+1, ..., SUNMIN(k+smu,n-1).                               */
    
    last_col_k = SUNMIN(k+smu,n-1);
    for (j=k+1; j <= last_col_k; j++) {
      
      col_j = a[j];
      storage_l = ROW(l,j,smu); 
      storage_k = ROW(k,j,smu); 
      a_kj = col_j[storage_l];

      /* Swap the elements a(k,j) and a(k,l) if l!=k. */
      
      if (swap) {
	col_j[storage_l] = col_j[storage_k];
	col_j[storage_k] = a_kj;
      }

      /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j) */
      /* a_kj = a(k,j), *kptr = - a(i,k)/a(k,k), *jptr = a(i,j) */

      if (a_kj != ZERO) {
	for (i=k+1, kptr=sub_diag_k, jptr=col_j+ROW(k+1,j,smu);
	     i <= last_row_k;
	     i++, kptr++, jptr++)
	  (*jptr) += a_kj * (*kptr);
      }
    }    
  }
  
  /* set the last pivot row to be n-1 and check for a zero pivot */

  *p = n-1; 
  if (a[n-1][smu] == ZERO) return(n);

  /* return 0 to indicate success */

  return(0);
}

void sunlinsolbandGBTRS(realtype **a, long int n, long int smu, 
			long int ml, long int *p, realtype *b)
{
  long int k, l, i, first_row_k, last_row_k;
  realtype mult, *diag_k;
  
  /* Solve Ly = Pb, store solution y in b */
  
  for (k=0; k < n-1; k++) {
    l = p[k];
    mult = b[l];
    if (l != k) {
      b[l] = b[k];
      b[k] = mult;
    }
    diag_k = a[k]+smu;
    last_row_k = SUNMIN(n-1,k+ml);
    for (i=k+1; i <= last_row_k; i++)
      b[i] += mult * diag_k[i-k];
  }
  
  /* Solve Ux = y, store solution x in b */
  
  for (k=n-1; k >= 0; k--) {
    diag_k = a[k]+smu;
    first_row_k = SUNMAX(0,k-smu);
    b[k] /= (*diag_k);
    mult = -b[k];
    for (i=first_row_k; i <= k-1; i++)
      b[i] += mult*diag_k[i-k];
  }
}

