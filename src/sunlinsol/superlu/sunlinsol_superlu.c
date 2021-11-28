/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * Based on codes <solver>_superlu.c, written by
 * Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the SuperLU implementation
 * of the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_superlu.h>
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)

/*
 * -----------------------------------------------------------------
 * SuperLU solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SLU_CONTENT(S)    ( (SUNLinearSolverContent_SuperLU)(S->content) )
#define LASTFLAG(S)         (  SLU_CONTENT(S)->last_flag )
#define FIRSTFACTORIZE(S)   (  SLU_CONTENT(S)->first_factorize )
#define SM_A(S)             (  SLU_CONTENT(S)->A )
#define SM_AC(S)            (  SLU_CONTENT(S)->AC )
#define SM_L(S)             (  SLU_CONTENT(S)->L )
#define SM_U(S)             (  SLU_CONTENT(S)->U )
#define SM_B(S)             (  SLU_CONTENT(S)->B )
#define SLUSTAT(S)          (&(SLU_CONTENT(S)->SLUstat ))
#define PERMR(S)            (  SLU_CONTENT(S)->perm_r )
#define PERMC(S)            (  SLU_CONTENT(S)->perm_c )
#define SIZE(S)             (  SLU_CONTENT(S)->N )
#define DIAGPIVOTTHRESH(S)  (  SLU_CONTENT(S)->diag_pivot_thresh )
#define ORDERING(S)         (  SLU_CONTENT(S)->ordering )
#define OPTIONS(S)          (&(SLU_CONTENT(S)->options) )
#define SM_X(S)             (  SLU_CONTENT(S)->X )
#define EQUED(S)            (  SLU_CONTENT(S)->equed )
#define S_R(S)              (  SLU_CONTENT(S)->R )
#define S_C(S)              (  SLU_CONTENT(S)->C )
#define ETREE(S)            (  SLU_CONTENT(S)->etree )
#define GLU(S)              (&(SLU_CONTENT(S)->Glu) )
#define RPG(S)              (&(SLU_CONTENT(S)->rpg) )
#define RCOND(S)            (&(SLU_CONTENT(S)->rcond) )
#define MEM_USAGE(S)        (&(SLU_CONTENT(S)->mem_usage) )
#define FERR(S)             (&(SLU_CONTENT(S)->ferr) )
#define BERR(S)             (&(SLU_CONTENT(S)->berr) )
#define LWORK(S)            (  SLU_CONTENT(S)->lwork )
#define WORK(S)             (  SLU_CONTENT(S)->work )

/*
 * -----------------------------------------------------------------
 * deprecated wrapper functions
 * -----------------------------------------------------------------
 */

SUNLinearSolver SUNSuperLU(N_Vector y, SUNMatrix A)
{ return(SUNLinSol_SuperLU(y, A)); }

int SUNSuperLUSetOrdering(SUNLinearSolver S, int ordering_choice)
{ return(SUNLinSol_SuperLUSetOrdering(S, ordering_choice)); }


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SuperLU linear solver
 */

SUNLinearSolver SUNLinSol_SuperLU(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_SuperLU content;
  sunindextype MatrixRows;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL);

  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) return(NULL);

  if (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) return(NULL);

  MatrixRows = SUNSparseMatrix_Rows(A);
  if (MatrixRows != N_VGetLength(y)) return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_SuperLU;
  S->ops->getid      = SUNLinSolGetID_SuperLU;
  S->ops->initialize = SUNLinSolInitialize_SuperLU;
  S->ops->setup      = SUNLinSolSetup_SuperLU;
  S->ops->solve      = SUNLinSolSolve_SuperLU;
  S->ops->lastflag   = SUNLinSolLastFlag_SuperLU;
  S->ops->space      = SUNLinSolSpace_SuperLU;
  S->ops->free       = SUNLinSolFree_SuperLU;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SuperLU) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->N                 = MatrixRows;
  content->last_flag         = 0;
  content->diag_pivot_thresh = ONE;
  content->ordering          = SUNSLU_ORDERING_DEFAULT;
  content->perm_r            = NULL;
  content->perm_c            = NULL;
  content->A                 = NULL;
  content->AC                = NULL;
  content->L                 = NULL;
  content->U                 = NULL;
  content->B                 = NULL;
  content->X                 = NULL;
  content->R                 = NULL;
  content->C                 = NULL;
  content->etree             = NULL;
  content->equed[0]          = ' ';
  content->rpg               = 0;
  content->rcond             = 0;
  content->ferr              = 0;
  content->berr              = 0;
  content->lwork             = 0;
  content->work              = NULL;

  /* Allocate content */
  content->perm_r = (sunindextype *) malloc(MatrixRows*sizeof(sunindextype));
  if (content->perm_r == NULL) { SUNLinSolFree(S); return(NULL); }

  content->perm_c = (sunindextype *) malloc(MatrixRows*sizeof(sunindextype));
  if (content->perm_c == NULL) { SUNLinSolFree(S); return(NULL); }

  content->A = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->A == NULL) { SUNLinSolFree(S); return(NULL); }
  content->A->Store = NULL;

  content->AC = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->AC == NULL) { SUNLinSolFree(S); return(NULL); }
  content->AC->Store = NULL;

  content->L = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->L == NULL) { SUNLinSolFree(S); return(NULL); }
  content->L->Store = NULL;

  content->U = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->U == NULL) { SUNLinSolFree(S); return(NULL); }
  content->U->Store = NULL;

  content->B = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->B == NULL) { SUNLinSolFree(S); return(NULL); }
  content->B->Store = NULL;
  xCreate_Dense_Matrix(content->B, MatrixRows, 1, NULL, MatrixRows, SLU_DN,
                       SLU_D, SLU_GE);

  content->X = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->X == NULL) { SUNLinSolFree(S); return(NULL); }
  content->X->Store = NULL;
  xCreate_Dense_Matrix(content->X, MatrixRows, 1, NULL, MatrixRows, SLU_DN,
                       SLU_D, SLU_GE);

  ETREE(S) = intMalloc(MatrixRows);
  content->R = (realtype*)SUPERLU_MALLOC(MatrixRows * sizeof(realtype));
  content->C = (realtype*)SUPERLU_MALLOC(MatrixRows * sizeof(realtype));

  set_default_options(&(content->options));


  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a SuperLU linear solver
 */

int SUNLinSol_SuperLUSetOrdering(SUNLinearSolver S, int ordering_choice)
{
  /* Check for legal ordering_choice */
  if ((ordering_choice < 0) || (ordering_choice > 3))
    return(SUNLS_ILL_INPUT);

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SUNLS_MEM_NULL);

  /* Set ordering_choice */
  ORDERING(S) = ordering_choice;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SuperLU(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


SUNLinearSolver_ID SUNLinSolGetID_SuperLU(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_SUPERLU);
}


int SUNLinSolInitialize_SuperLU(SUNLinearSolver S)
{
  /* force a first factorization */
  FIRSTFACTORIZE(S) = 1;

  /* Initialize statistics variables */
  //StatInit(SLUSTAT(S));

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetup_SuperLU(SUNLinearSolver S, SUNMatrix A)
{
  int_t info;
  int panel_size, relax;
  double drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;

  /* Set option values for SuperLU */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  fact = DOFACT;
  trans = (SUNSparseMatrix_SparseType(A) == CSC_MAT) ? NOTRANS : TRANS;
  usepr = NO;
  drop_tol = ZERO;

  OPTIONS(S)->Trans = trans;

  /* free and reallocate sparse matrix */
  if (SM_A(S)->Store)
    SUPERLU_FREE(SM_A(S)->Store);
  xCreate_CompCol_Matrix(SM_A(S), SUNSparseMatrix_Rows(A),
			 SUNSparseMatrix_Columns(A),
                         SUNSparseMatrix_NNZ(A),
                         SUNSparseMatrix_Data(A),
                         (sunindextype*) SUNSparseMatrix_IndexValues(A),
                         (sunindextype*) SUNSparseMatrix_IndexPointers(A),
			 SLU_NC, SLU_D, SLU_GE);

  /* On first decomposition, set up reusable pieces */
  if (FIRSTFACTORIZE(S)) {
      OPTIONS(S)->Fact = DOFACT;
      refact = NO;
      FIRSTFACTORIZE(S) = 0;

  } else {

      /* Re-initialize statistics variables */
      OPTIONS(S)->Fact = SamePattern_SameRowPerm;
      refact = YES;

  }
  StatInit(SLUSTAT(S));

  /* Initialize the option structure using the user-input parameters.
     Subsequent calls will re-initialize options.  Apply perm_c to
     columns of original A to form AC */

  /* Compute the LU factorization of A.
     The following routine will create num_threads threads. */
     /* ONLY PERFORM THE LU DECOMPOSITION */
  int nrhs =   SM_B(S)->ncol;
  SM_B(S)->ncol = 0;  /* Indicate not to solve the system */
  dgssvx(OPTIONS(S), SM_A(S), PERMC(S), PERMR(S), ETREE(S), EQUED(S), S_R(S), S_C(S),
      SM_L(S), SM_U(S), WORK(S), LWORK(S), SM_B(S), SM_X(S), RPG(S), RCOND(S), FERR(S), BERR(S),
      GLU(S), MEM_USAGE(S), SLUSTAT(S), &info);

  SM_B(S)->ncol = nrhs;  /* Indicate not to solve the system */
  OPTIONS(S)->Fact = FACTORED; /* Indicate the factored form of A is supplied. */

  StatFree(SLUSTAT(S));


  if (info != 0) {
    LASTFLAG(S) = (info < 0) ?
      SUNLS_PACKAGE_FAIL_UNREC : SUNLS_PACKAGE_FAIL_REC;
    return(LASTFLAG(S));
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSolve_SuperLU(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, realtype tol)
{
  int_t info;
  realtype *xdata, *bdata;
  DNformat *Bstore, *Xstore;
  trans_t trans;

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  bdata = N_VGetArrayPointer(b);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  if (bdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }

  Bstore = (DNformat *) (SM_B(S)->Store);
  Xstore = (DNformat *) (SM_X(S)->Store);
  Bstore->nzval = bdata;
  Xstore->nzval = xdata;

  /* Call SuperLU to solve the linear system using L and U */
  trans = (SUNSparseMatrix_SparseType(A) == CSC_MAT) ? NOTRANS : TRANS;

  StatInit(SLUSTAT(S));

  dgssvx(OPTIONS(S), SM_A(S), PERMC(S), PERMR(S), ETREE(S), EQUED(S), S_R(S), S_C(S),
      SM_L(S), SM_U(S), WORK(S), LWORK(S), SM_B(S), SM_X(S), RPG(S), RCOND(S), FERR(S), BERR(S),
      GLU(S), MEM_USAGE(S), SLUSTAT(S), &info);

  StatFree(SLUSTAT(S));

  if (info != 0) {
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    return(LASTFLAG(S));
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


sunindextype SUNLinSolLastFlag_SuperLU(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(-1);
  return(LASTFLAG(S));
}


int SUNLinSolSpace_SuperLU(SUNLinearSolver S,
                             long int *lenrwLS,
                             long int *leniwLS)
{
  /* since the SuperLU structures are opaque objects, we
     omit those from these results */
  *leniwLS = 5 + 2*SIZE(S);
  *lenrwLS = 1;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_SuperLU(SUNLinearSolver S)
{
    int lwork;
    lwork = 0;
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* delete items from the contents structure (if it exists) */
  if (S->content) {
    if (OPTIONS(S) && SM_AC(S)) {
      //pxgstrf_finalize(OPTIONS(S), SM_AC(S));
    }
    if (PERMR(S)) {
      free(PERMR(S));
      PERMR(S) = NULL;
    }
    if (PERMC(S)) {
      free(PERMC(S));
      PERMC(S) = NULL;
    }
    if (SM_X(S)) {
        Destroy_SuperMatrix_Store(SM_X(S));
        free(SM_X(S));
        SM_X(S) = NULL;
    }
    if (SM_B(S)) {
        Destroy_SuperMatrix_Store(SM_B(S));
        free(SM_B(S));
        SM_B(S) = NULL;
    }
    if (lwork == 0) {
        if (SM_L(S)) {
            Destroy_SuperNode_Matrix(SM_L(S));
            free(SM_L(S));
            SM_L(S) = NULL;
        }
        if (SM_U(S)) {
            Destroy_CompCol_Matrix(SM_U(S));
            free(SM_U(S));
            SM_U(S) = NULL;
        }
    }
    else if (lwork > 0) {
        SUPERLU_FREE(WORK(S));
    }

    if (SM_A(S)) {
      if (SM_A(S)->Store) {
        SUPERLU_FREE(SM_A(S)->Store);
        SM_A(S)->Store = NULL;
      }
      free(SM_A(S));
      SM_A(S) = NULL;
    }
    if (SM_AC(S)) {
      free(SM_AC(S));
      SM_AC(S) = NULL;
    }
    if (ETREE(S)) {
        SUPERLU_FREE(ETREE(S));
      ETREE(S) = NULL;
    }
    if (S_R(S)) {
        SUPERLU_FREE(S_R(S));
        S_R(S) = NULL;
    }
    if (S_C(S)) {
        SUPERLU_FREE(S_C(S));
        S_C(S) = NULL;
    }
    free(S->content);
    S->content = NULL;
  }

  /* delete generic structures */
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}
