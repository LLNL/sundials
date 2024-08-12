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
 * -----------------------------------------------------------------*/

#include <sunadjoint/sunadjoint_solver.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sunadjoint/sunadjoint_checkpointscheme.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_stepper.h"
#include "sundials/sundials_types.h"

SUNErrCode SUNAdjointSolver_Create(SUNStepper stepper, int64_t final_step_idx,
                                   N_Vector sf, sunrealtype tf,
                                   SUNAdjointCheckpointScheme checkpoint_scheme,
                                   SUNContext sunctx,
                                   SUNAdjointSolver* adj_solver_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointSolver adj_solver = malloc(sizeof(struct SUNAdjointSolver_));
  SUNAssert(adj_solver, SUN_ERR_MALLOC_FAIL);

  adj_solver->stepper           = stepper;
  adj_solver->Jac               = NULL;
  adj_solver->JacP              = NULL;
  adj_solver->JacFn             = NULL;
  adj_solver->JacPFn            = NULL;
  adj_solver->Jvp               = NULL;
  adj_solver->JPvp              = NULL;
  adj_solver->vJp               = NULL;
  adj_solver->vJPp              = NULL;
  adj_solver->checkpoint_scheme = checkpoint_scheme;
  adj_solver->tf                = tf;
  adj_solver->user_data         = NULL;
  adj_solver->sunctx            = sunctx;
  adj_solver->step_idx          = final_step_idx;

  *adj_solver_ptr = adj_solver;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_Solve(SUNAdjointSolver adj_solver, sunrealtype tout,
                                  N_Vector sens, sunrealtype* tret,
                                  int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper stepper = adj_solver->stepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_solver->tf;
  *stop_reason       = 0;
  while (t > tout)
  {
    SUNCheckCall(SUNAdjointSolver_Step(adj_solver, tout, sens, tret, stop_reason));
    if (*stop_reason < 0)
    {
      retcode = SUN_ERR_ADJOINT_STEPPERFAILED;
      break;
    }
    else { t = *tret; }
  }

  return retcode;
}

SUNErrCode SUNAdjointSolver_Step(SUNAdjointSolver adj_solver, sunrealtype tout,
                                 N_Vector sens, sunrealtype* tret,
                                 int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper stepper = adj_solver->stepper;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = adj_solver->tf;
  *stop_reason       = 0;
  SUNCheckCall(
    SUNStepper_Step(stepper, adj_solver->tf, tout, sens, &t, stop_reason));
  if (*stop_reason < 0) { retcode = SUN_ERR_ADJOINT_STEPPERFAILED; }
  else if (*stop_reason > 1)
  {
    // TODO(CJB): what reasons could this happen, and are they valid?
    // 1==TSTOP_RETURN
    // 2==ROOT_RETURN
    fprintf(stderr, ">>>> HERE, stop_reason = %d\n", *stop_reason);
  }
  adj_solver->step_idx--;
  *tret = t;

  return retcode;
}

SUNErrCode SUNAdjointSolver_Destroy(SUNAdjointSolver* adj_solver_ptr)
{
  SUNAdjointSolver adj_solver = *adj_solver_ptr;
  // SUNAdjointCheckpointScheme_Destroy(adj_solver->checkpoint_scheme);
  SUNStepper_Destroy(&adj_solver->stepper);
  free(adj_solver);
  *adj_solver_ptr = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetJacFn(SUNAdjointSolver adj_solver,
                                     SUNJacFn JacFn, SUNMatrix Jac,
                                     SUNJacFn JacPFn, SUNMatrix JacP)
{
  SUNFunctionBegin(adj_solver->sunctx);

  SUNAssert(JacFn && Jac, SUN_ERR_ARG_CORRUPT);
  SUNAssert(!JacPFn || (JacPFn && JacP), SUN_ERR_ARG_CORRUPT);

  adj_solver->JacFn  = JacFn;
  adj_solver->Jac    = Jac;
  adj_solver->JacPFn = JacPFn;
  adj_solver->JacP   = JacP;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetJacTimesVecFn(SUNAdjointSolver adj_solver,
                                             SUNJacTimesFn Jvp,
                                             SUNJacTimesFn JPvp)
{
  SUNFunctionBegin(adj_solver->sunctx);

  SUNAssert(Jvp, SUN_ERR_ARG_CORRUPT);
  SUNAssert(!JPvp || (JPvp && Jvp), SUN_ERR_ARG_CORRUPT);

  adj_solver->Jvp  = Jvp;
  adj_solver->JPvp = JPvp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetVecTimesJacFn(SUNAdjointSolver adj_solver,
                                             SUNJacTimesFn vJp,
                                             SUNJacTimesFn vJPp)
{
  SUNFunctionBegin(adj_solver->sunctx);

  SUNAssert(vJp, SUN_ERR_ARG_CORRUPT);
  SUNAssert(!vJPp || (vJp && vJPp), SUN_ERR_ARG_CORRUPT);

  adj_solver->vJp  = vJp;
  adj_solver->vJPp = vJPp;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_SetUserData(SUNAdjointSolver adj_solver,
                                        void* user_data)
{
  SUNFunctionBegin(adj_solver->sunctx);

  adj_solver->user_data = user_data;

  return SUN_SUCCESS;
}
