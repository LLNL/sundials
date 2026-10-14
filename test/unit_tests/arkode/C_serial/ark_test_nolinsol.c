/* -----------------------------------------------------------------------------
 * Programmer(s): Montana
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Unit test checking that an implicit method using a Newton-type nonlinear
 * solver without an attached linear solver returns an error rather than
 * crashing (GitHub issue #957), while valid configurations still succeed.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_mristep.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#define HFIXED SUN_RCONST(0.01)
#define TOUT   SUN_RCONST(0.1)

/* Linear solver / nonlinear solver configurations to test */
enum config
{
  NO_LS,        /* Newton without a linear solver (invalid) */
  DENSE_LS,     /* Newton with a dense linear solver        */
  FIXED_POINT,  /* fixed-point nonlinear solver             */
  EXPLICIT_NLS, /* explicit method with a Newton solver (MRIStep only) */
  NUM_CONFIGS
};

static const char* config_names[] = {"no linear solver", "dense linear solver",
                                     "fixed-point solver",
                                     "explicit method + Newton solver"};

/* y' = -y */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VScale(-ONE, y, ydot);
  return 0;
}

/* Fast RHS for MRIStep tests */
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return 0;
}

/* Attach the solvers for a given configuration and evolve */
static int configure_and_evolve(void* arkode_mem, enum config cfg, N_Vector y,
                                SUNMatrix* A, SUNLinearSolver* LS,
                                SUNNonlinearSolver* NLS, SUNContext sunctx)
{
  int retval;
  sunrealtype tret;

  if (cfg == DENSE_LS)
  {
    *A  = SUNDenseMatrix(1, 1, sunctx);
    *LS = SUNLinSol_Dense(y, *A, sunctx);
    if (!*A || !*LS) { return -1000; }
    retval = ARKodeSetLinearSolver(arkode_mem, *LS, *A);
    if (retval) { return -1000; }
  }
  else if (cfg == FIXED_POINT)
  {
    *NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
    if (!*NLS) { return -1000; }
    retval = ARKodeSetNonlinearSolver(arkode_mem, *NLS);
    if (retval) { return -1000; }
  }
  else if (cfg == EXPLICIT_NLS)
  {
    *NLS = SUNNonlinSol_Newton(y, sunctx);
    if (!*NLS) { return -1000; }
    retval = ARKodeSetNonlinearSolver(arkode_mem, *NLS);
    if (retval) { return -1000; }
  }

  retval = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-6),
                              SUN_RCONST(1.0e-10));
  if (retval) { return -1000; }

  retval = ARKodeSetFixedStep(arkode_mem, HFIXED);
  if (retval) { return -1000; }

  return ARKodeEvolve(arkode_mem, TOUT, y, &tret, ARK_NORMAL);
}

/* Check the return value from ARKodeEvolve against the expected outcome */
static int check_result(const char* stepper, enum config cfg, int retval)
{
  int expected = (cfg == NO_LS) ? ARK_NLS_INIT_FAIL : ARK_SUCCESS;

  if (retval != expected)
  {
    fprintf(stderr, "FAIL: %s with %s returned %d, expected %d\n", stepper,
            config_names[cfg], retval, expected);
    return 1;
  }

  printf("PASS: %s with %s returned %d\n", stepper, config_names[cfg], retval);
  return 0;
}

static void free_solvers(SUNMatrix* A, SUNLinearSolver* LS,
                         SUNNonlinearSolver* NLS)
{
  if (*LS) { SUNLinSolFree(*LS); }
  if (*A) { SUNMatDestroy(*A); }
  if (*NLS) { SUNNonlinSolFree(*NLS); }
  *LS  = NULL;
  *A   = NULL;
  *NLS = NULL;
}

/* Main program */
int main(int argc, char* argv[])
{
  int retval, cfg;
  int fails                         = 0;
  SUNContext sunctx                 = NULL;
  N_Vector y                        = NULL;
  void* arkode_mem                  = NULL;
  void* inner_mem                   = NULL;
  MRIStepInnerStepper inner_stepper = NULL;
  SUNMatrix A                       = NULL;
  SUNLinearSolver LS                = NULL;
  SUNNonlinearSolver NLS            = NULL;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", retval);
    return 1;
  }

  y = N_VNew_Serial(1, sunctx);
  if (!y)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }

  /* ------------ *
   * Test ARKStep *
   * ------------ */

  /* ARKStep does not allow attaching a nonlinear solver to an explicit
     method, so the explicit configuration is only tested with MRIStep */
  for (cfg = 0; cfg < EXPLICIT_NLS; cfg++)
  {
    N_VConst(ONE, y);

    arkode_mem = ARKStepCreate(NULL, f, ZERO, y, sunctx);
    if (!arkode_mem)
    {
      fprintf(stderr, "ARKStepCreate returned NULL\n");
      return 1;
    }

    /* Backward Euler has no embedding, as in the issue report */
    retval = ARKStepSetTableNum(arkode_mem, ARKODE_BACKWARD_EULER_1_1,
                                ARKODE_ERK_NONE);
    if (retval)
    {
      fprintf(stderr, "ARKStepSetTableNum returned %i\n", retval);
      return 1;
    }

    retval = configure_and_evolve(arkode_mem, cfg, y, &A, &LS, &NLS, sunctx);
    fails += check_result("ARKStep", cfg, retval);

    ARKodeFree(&arkode_mem);
    free_solvers(&A, &LS, &NLS);
  }

  /* ------------ *
   * Test MRIStep *
   * ------------ */

  for (cfg = 0; cfg < NUM_CONFIGS; cfg++)
  {
    N_VConst(ONE, y);

    inner_mem = ARKStepCreate(f0, NULL, ZERO, y, sunctx);
    if (!inner_mem)
    {
      fprintf(stderr, "ARKStepCreate returned NULL\n");
      return 1;
    }

    retval = ARKodeSetFixedStep(inner_mem, HFIXED);
    if (retval)
    {
      fprintf(stderr, "ARKodeSetFixedStep returned %i\n", retval);
      return 1;
    }

    retval = ARKodeCreateMRIStepInnerStepper(inner_mem, &inner_stepper);
    if (retval)
    {
      fprintf(stderr, "ARKodeCreateMRIStepInnerStepper returned %i\n", retval);
      return 1;
    }

    /* Implicit slow time scale, or explicit for the explicit configuration */
    if (cfg == EXPLICIT_NLS)
    {
      arkode_mem = MRIStepCreate(f, NULL, ZERO, y, inner_stepper, sunctx);
    }
    else
    {
      arkode_mem = MRIStepCreate(NULL, f, ZERO, y, inner_stepper, sunctx);
    }
    if (!arkode_mem)
    {
      fprintf(stderr, "MRIStepCreate returned NULL\n");
      return 1;
    }

    retval = configure_and_evolve(arkode_mem, cfg, y, &A, &LS, &NLS, sunctx);
    fails += check_result("MRIStep", cfg, retval);

    ARKodeFree(&arkode_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    ARKodeFree(&inner_mem);
    free_solvers(&A, &LS, &NLS);
  }

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("%d test(s) failed\n", fails);
    return 1;
  }

  printf("SUCCESS\n");
  return 0;
}
