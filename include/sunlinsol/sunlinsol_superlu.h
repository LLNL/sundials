/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes sundials_superlu_impl.h and <solver>_superlu.h
 *     written by Carol S. Woodward @ LLNL
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
 * This is the header file for the SuperLU implementation of the
 * SUNLINSOL module, SUNLINSOL_SUPERLU.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SLU_H
#define _SUNLINSOL_SLU_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

/* Assume SuperLU library was built with compatible index type */
#if defined(SUNDIALS_INT64_T)
#define _LONGINT
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SuperLU solver parameters */
#define SUNSLU_ORDERING_DEFAULT  3     /* COLAMD */

/* Interfaces to match 'realtype' with the correct SuperLU functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#ifndef _SLU_H
#define _SLU_H
#include "slu_ddefs.h"
#endif
#define xgstrs                  dgstrs
#define pxgstrf                 pdgstrf
#define pxgstrf_init            pdgstrf_init
#define xCreate_Dense_Matrix    dCreate_Dense_Matrix
#define xCreate_CompCol_Matrix  dCreate_CompCol_Matrix
#elif defined(SUNDIALS_SINGLE_PRECISION)
#ifndef _SLU_H
#define _SLU_H
#include "slu_sdefs.h"
#endif
#define xgstrs                  sgstrs
#define pxgstrf                 psgstrf
#define pxgstrf_init            psgstrf_init
#define xCreate_Dense_Matrix    sCreate_Dense_Matrix
#define xCreate_CompCol_Matrix  sCreate_CompCol_Matrix
#else  /* incompatible sunindextype for SuperLU */
#error  Incompatible realtype for SuperLU
#endif


/* --------------------------------------------
 * SuperLU Implementation of SUNLinearSolver
 * -------------------------------------------- */

struct _SUNLinearSolverContent_SuperLU {
  int          last_flag;
  int          first_factorize;
  SuperMatrix  *A, *AC, *L, *U, *B;
  SuperLUStat_t      SLUstat;
  sunindextype *perm_r, *perm_c;
  sunindextype N;
  realtype     diag_pivot_thresh;
  int          ordering;
  superlu_options_t options;
  SuperMatrix  *X;
  char equed[1];
  realtype* R, * C;
  sunindextype *etree;
  GlobalLU_t	   Glu;
  realtype rpg, rcond, ferr, berr;
  mem_usage_t    mem_usage;
  int lwork;
  void* work;
};

typedef struct _SUNLinearSolverContent_SuperLU *SUNLinearSolverContent_SuperLU;


/* -------------------------------------------
 * Exported Functions for SUNLINSOL_SUPERLU
 * ------------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SuperLU(N_Vector y,
                                                    SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSol_SuperLUSetOrdering(SUNLinearSolver S,
                                                   int ordering_choice);

/* deprecated */
SUNDIALS_EXPORT SUNLinearSolver SUNSuperLU(N_Vector y, SUNMatrix A);
/* deprecated */
SUNDIALS_EXPORT int SUNSuperLUSetOrdering(SUNLinearSolver S,
                                            int ordering_choice);

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SuperLU(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SuperLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_SuperLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_SuperLU(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SuperLU(SUNLinearSolver S, SUNMatrix A,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SuperLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_SuperLU(SUNLinearSolver S,
                                             long int *lenrwLS,
                                             long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_SuperLU(SUNLinearSolver S);


#ifdef __cplusplus
}
#endif

#endif
