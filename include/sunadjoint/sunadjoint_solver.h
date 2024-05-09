/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNAdjointSolver class definition.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_SOLVER_H
#define _SUNADJOINT_SOLVER_H

#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

#include "sundials/sundials_types.h"

struct SUNAdjointSolver_
{
  SUNStepper stepper;
  SUNJacFn Jac;
  SUNJacFn JacP;
  SUNJacTimesFn Jvp;
  SUNJacTimesFn vJp;
  SUNJacTimesFn vJPp;
  SUNAdjointCheckpointScheme checkpoint_scheme;
  SUNContext sunctx;
};

typedef struct SUNAdjointSolver_* SUNAdjointSolver;

#ifdef __cplusplus
extern "C" {
#endif

// IDEA: In lieu of Stepper_ID each package that supports adjoint can have a function that creates the adjoint solver.
// E.g., SUNAdjointSolver ARKStepCreateAdjointSolver();
SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Create(SUNStepper stepper,
                                   sunindextype num_cost_fns, N_Vector sf,
                                   SUNAdjointCheckpointScheme checkpoint_scheme,
                                   SUNContext sunctx,
                                   SUNAdjointSolver* adj_solver);

/*
  Solves the adjoint system.

  :param adj_solver: The adjoint solver object.
  :param tf: The final output time from the forward integration. 
             This is the "starting" time for adjoint solver's backwards integration.
  :param tout: The time at which the adjoint solution is desired.
  :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
  :param tret: On return, the time reached by the adjoint solver.
  :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

  :returns: A SUNErrCode indicating failure or success.
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Solve(SUNAdjointSolver adj_solver, sunrealtype tf,
                                  sunrealtype tout, N_Vector sens,
                                  sunrealtype* tret, int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Step(SUNAdjointSolver, sunrealtype t0, N_Vector sens,
                                 sunrealtype* tret, int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetJacFn(SUNAdjointSolver, SUNJacFn Jac,
                                     SUNJacFn JacP);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetJacTimesVecFn(SUNAdjointSolver, SUNJacTimesFn Jvp,
                                             SUNJacTimesFn JPvp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetVecTimesJacFn(SUNAdjointSolver, SUNJacTimesFn vJp,
                                             SUNJacTimesFn vJPp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Destroy(SUNAdjointSolver*);

#ifdef __cplusplus
}
#endif
#endif /* _SUNADJOINT_SOLVER_H */
