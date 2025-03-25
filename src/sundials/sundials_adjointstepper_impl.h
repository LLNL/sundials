/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_ADJOINTSTEPPER_IMPL_H_
#define SUNDIALS_ADJOINTSTEPPER_IMPL_H_

#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

struct SUNAdjointStepper_
{
  suncountertype njeval, njpeval, njtimesv, njptimesv, nrecompute;
  suncountertype final_step_idx;

  SUNStepper adj_sunstepper;
  SUNStepper fwd_sunstepper;
  SUNAdjointCheckpointScheme checkpoint_scheme;

  SUNMatrix Jac, JacP;
  SUNRhsJacFn JacFn, JacPFn;
  SUNRhsJacTimesFn JvpFn, JPvpFn;

  void* user_data;
  void* content;
  SUNContext sunctx;

  sunrealtype tf;
  int last_flag;
};

#endif
