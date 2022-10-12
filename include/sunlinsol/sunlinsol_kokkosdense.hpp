/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file for a SUNLinearSolver using Kokkoks Kernels
 * ---------------------------------------------------------------------------*/

#ifndef _SUNLINSOL_KOKKOSDENSE_H
#define _SUNLINSOL_KOKKOSDENSE_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

struct _SUNLinearSolverContent_KokkosDense
{
  int last_flag; /* last return flag */
};

typedef struct _SUNLinearSolverContent_KokkosDense *SUNLinearSolverContent_KokkosDense;

/* ----------------------------------
 * Implementation specific functions
 * ----------------------------------*/

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_KokkosDense(N_Vector y, SUNMatrix A,
                                                      SUNContext sunctx);

/* ---------------------------------------
 * Implementation of SUNMatrix operations
 * ---------------------------------------*/

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_KokkosDense(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_KokkosDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_KokkosDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_KokkosDense(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_KokkosDense(SUNLinearSolver S, SUNMatrix A,
                                               N_Vector x, N_Vector b,
                                               sunrealtype tol);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_KokkosDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_KokkosDense(SUNLinearSolver S);

#endif
