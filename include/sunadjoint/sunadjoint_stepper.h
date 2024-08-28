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
 * SUNAdjointStepper class definition.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_STEPPER_H
#define _SUNADJOINT_STEPPER_H

#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

#include "sundials/sundials_types.h"

struct SUNAdjointStepper_
{
  SUNStepper adj_sunstepper;
  SUNStepper fwd_sunstepper;
  SUNAdjointCheckpointScheme checkpoint_scheme;

  sunrealtype tf;
  int64_t step_idx, final_step_idx;

  /* Jacobian-related data */
  SUNMatrix Jac, JacP;
  SUNJacFn JacFn, JacPFn;
  SUNJacTimesFn JvpFn, JPvpFn, vJpFn, vJPpFn;

  /* counters */
  uint64_t nst, njeval, njpeval, njtimesv, njptimesv, nvtimesj, nvtimesjp,
    nrecompute;

  void* user_data;
  void* content;
  SUNContext sunctx;
};

typedef struct SUNAdjointStepper_* SUNAdjointStepper;

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Create(
  SUNStepper fwd_sunstepper, SUNStepper adj_sunstepper, int64_t final_step_idx,
  N_Vector sf, sunrealtype tf, SUNAdjointCheckpointScheme checkpoint_scheme,
  SUNContext sunctx, SUNAdjointStepper* adj_stepper);

/*
  Reinitializes the adjoint stepper to solve a new problem of the same size.

  :param adj_stepper: The SUNAdjointStepper object.
  :param tf: The time to start integrating the adjoint system from.
  :param sf: The terminal condition vector of sensitivity solutions dg/dy0 and dg/dp.

  :return: A SUNErrCode indicating failure or success.
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, N_Vector sf,
                                    sunrealtype tf);

/*
  Integrates the adjoint system.

  :param adj_stepper: The SUNAdjointStepper object.
  :param tout: The time at which the adjoint solution is desired.
  :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
  :param tret: On return, the time reached by the adjoint solver.
  :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

  :return: A SUNErrCode indicating failure or success.
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper,
                                    sunrealtype tout, N_Vector sens,
                                    sunrealtype* tret, int* stop_reason);

/*
  Evolves the adjoint system backwards one step.

  :param adj_stepper: The SUNAdjointStepper object.
  :param tout: The time at which the adjoint solution is desired.
  :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
  :param tret: On return, the time reached by the adjoint solver.
  :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

  :return: A SUNErrCode indicating failure or success.
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper,
                                     sunrealtype tout, N_Vector sens,
                                     sunrealtype* tret, int* stop_reason);

/*
  Evolves the forward system in time from start_idx/t0 to stop_idx/tf with dense checkpointing.

  :param adj_stepper: The SUNAdjointStepper object.
  :param start_idx: the index of the step, w.r.t. the original forward integration, to begin forward integration from.
  :param stop_idx: the index of the step, w.r.t. the original forward integration, to end forward integration.
  :param t0: the initial time, w.r.t. the original forward integration, to start forward integration from.
  :param tf: the final time, w.r.t. the original forward integration, to stop forward integration.
  :param y0: the initial state, w.r.t. the original forward integration, to start forward integration.

  :return: A SUNErrCode indicating failure or success.
*/
SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper,
                                          int64_t start_idx, int64_t stop_idx,
                                          sunrealtype t0, sunrealtype tf,
                                          N_Vector y0);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper, SUNJacFn JacFn,
                                      SUNMatrix Jac, SUNJacFn JacPFn,
                                      SUNMatrix JP);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper,
                                              SUNJacTimesFn Jvp,
                                              SUNJacTimesFn JPvp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetVecTimesJacFn(SUNAdjointStepper,
                                              SUNJacTimesFn vJp,
                                              SUNJacTimesFn vJPp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper, void* user_data);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_stepper,
                                           FILE* outfile, SUNOutputFormat fmt);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Destroy(SUNAdjointStepper*);

#ifdef __cplusplus
}
#endif
#endif /* _SUNADJOINT_STEPPER_H */
