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

struct SUNAdjointSolver_s
{
  SUNStepper stepper;
  SUNJacFn Jac;
  SUNJacFn JacP;
  SUNVecTimesJacFn Jvp;
  SUNVecTimesJacFn vJp;
  SUNVecTimesJacFn vJPp;
};

typedef SUNAdjointSolver_s* SUNAdjointSolver;

#ifdef __cplusplus
extern "C" {

SUNDIALS_EXPORT
SUNAdjointSolver_Create(void* integrator_mem, SUNStepper_ID stepper_id,
                        sunindextype num_cost_fns, N_Vector sf,
                        SUNCheckpointScheme chkpt_scheme, SUNContext sunctx,
                        SUNAdjointSolver* adj_solver);

SUNDIALS_EXPORT
SUNAdjointSolver_Solve(SUNAdjointSolver, sunrealtype t0, N_Vector sens,
                       sunrealtype* tret, int mode);

SUNDIALS_EXPORT
SUNAdjointSolver_SetJacFn(SUNAdjointSolver, SUNJacFn Jac, SUNJacFn JacP);

SUNDIALS_EXPORT
SUNAdjointSolver_SetJacTimesVecFn(SUNAdjointSolver, SUNVecTimesJacFn Jvp,
                                  SUNVecTimesJacFn JPvp);

SUNDIALS_EXPORT
SUNAdjointSolver_SetVecTimesJacFn(SUNAdjointSolver, SUNVecTimesJacFn vJp,
                                  SUNVecTimesJacFn vJPp);

SUNDIALS_EXPORT
SUNAdjointSolver_Destroy(SUNAdjointSolver*);

#ifdef __cplusplus
}
#endif
#endif /* _SUNADJOINT_SOLVER_H */
