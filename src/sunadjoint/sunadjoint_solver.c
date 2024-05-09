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

SUNErrCode SUNAdjointSolver_Create(SUNStepper stepper,
                                   sunindextype num_cost_fns, N_Vector sf,
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
  adj_solver->Jvp               = NULL;
  adj_solver->vJp               = NULL;
  adj_solver->vJPp              = NULL;
  adj_solver->checkpoint_scheme = checkpoint_scheme;
  adj_solver->sunctx            = sunctx;

  *adj_solver_ptr = adj_solver;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointSolver_Solve(SUNAdjointSolver adj_solver, sunrealtype tf,
                                  sunrealtype tout, N_Vector sens,
                                  sunrealtype* tret, int* stop_reason)

{
  SUNFunctionBegin(adj_solver->sunctx);
  SUNStepper stepper                = adj_solver->stepper;
  SUNAdjointCheckpointScheme scheme = adj_solver->checkpoint_scheme;

  SUNErrCode retcode = SUN_SUCCESS;
  sunrealtype t      = tf;
  *stop_reason       = 0;
  while (t > tout)
  {
    SUNCheckCall(SUNStepper_Advance(stepper, tf, tout, sens, &t, stop_reason));
    if (*stop_reason < 0)
    {
      retcode = SUN_ERR_ADJOINT_STEPPERFAILED;
      break;
    }
    else if (*stop_reason > 0)
    {
      // TODO(CJB): what reasons could this happen, and are they valid?
      // (1) TSTOP_RETURN
      // (2) ROOT_RETURN
      fprintf(stderr, ">>>> HERE, stop_reason = %d\n", *stop_reason);
    }
    else { break; }
  }

  *tret = tf;
  return retcode;
}

SUNErrCode SUNAdjointSolver_Destroy(SUNAdjointSolver* adj_solver_ptr)
{
  SUNAdjointSolver adj_solver = *adj_solver_ptr;
  // SUNAdjointCheckpointScheme_Destroy(adj_solver->checkpoint_scheme);
  free(adj_solver);
  *adj_solver_ptr = NULL;
  return SUN_SUCCESS;
}