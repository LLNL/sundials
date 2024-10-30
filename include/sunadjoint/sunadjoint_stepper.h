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
  SUNRhsJacFn JacFn, JacPFn;
  SUNRhsJacTimesFn JvpFn, JPvpFn, vJpFn, vJPpFn;

  /* counters */
  int64_t nst, njeval, njpeval, njtimesv, njptimesv, nvtimesj, nvtimesjp,
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

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, N_Vector y0,
                                    sunrealtype t0, N_Vector sf, sunrealtype tf);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper,
                                    sunrealtype tout, N_Vector sens,
                                    sunrealtype* tret, int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper,
                                     sunrealtype tout, N_Vector sens,
                                     sunrealtype* tret, int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper,
                                          int64_t start_idx, sunrealtype t0,
                                          sunrealtype tf, N_Vector y0);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetJacFn(SUNAdjointStepper, SUNRhsJacFn JacFn,
                                      SUNMatrix Jac, SUNRhsJacFn JacPFn,
                                      SUNMatrix JP);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetJacTimesVecFn(SUNAdjointStepper,
                                              SUNRhsJacTimesFn Jvp,
                                              SUNRhsJacTimesFn JPvp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetVecTimesJacFn(SUNAdjointStepper,
                                              SUNRhsJacTimesFn vJp,
                                              SUNRhsJacTimesFn vJPp);

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
