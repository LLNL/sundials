/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Unit test for GetUserData functions
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"

static void err_fn(const int line, const char* const func, const char* const file,
                   const char* const msg, const SUNErrCode err_code,
                   void* const err_user_data, const SUNContext sunctx)
{
  fprintf(stderr, "Error at line %i of %s in %s: %s\n", line, func, file, msg);
  exit(err_code);
}

// RHS for the simple ODE y' = 0
static int f(const sunrealtype t, const N_Vector y, const N_Vector ydot,
             void* const user_data)
{
  N_VConst(SUN_RCONST(0.0), ydot);
  return 0;
}

static int check_step(void* const arkode_mem, const N_Vector y,
                      const sunrealtype h_expected, const int step)
{
  sunrealtype tret;
  /* The ERK method should be able to take the maximum possible timestep for the
     simple ODE y' = 0 without any rejected steps or local error */
  ARKodeEvolve(arkode_mem, SUN_RCONST(1.0), y, &tret, ARK_ONE_STEP);

  long int local_err_fails;
  ARKodeGetNumErrTestFails(arkode_mem, &local_err_fails);
  if (local_err_fails != 0)
  {
    fprintf(stderr, "Expected 0 local error failures at step %i but is %li\n",
            step, local_err_fails);
  }

  const N_Vector err = N_VClone(y);
  ARKodeGetEstLocalErrors(arkode_mem, err);
  const sunrealtype err_norm = N_VMaxNorm(err);
  N_VDestroy(err);

  if (err_norm != 0)
  {
    fprintf(stderr, "Expected local error at step %i to be 0 but is %g\n", step,
            err_norm);
    return 1;
  }

  sunrealtype h_actual;
  ARKodeGetCurrentStep(arkode_mem, &h_actual);
  if (h_expected != h_actual)
  {
    fprintf(stderr, "Expected h at step %i to be %g but is %g\n", step,
            h_expected, h_actual);
    return 1;
  }

  return 0;
}

int main()
{
  SUNContext sunctx;
  int retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval != 0)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", retval);
    return 1;
  }

  retval = SUNContext_PushErrHandler(sunctx, err_fn, NULL);
  if (retval != 0)
  {
    fprintf(stderr, "SUNContext_PushErrHandler returned %i\n", retval);
    return 1;
  }

  const N_Vector y = N_VNew_Serial(1, sunctx);
  N_VConst(SUN_RCONST(1.0), y);

  void* arkode_mem = ERKStepCreate(f, SUN_RCONST(0.0), y, sunctx);

  const sunrealtype h0           = SUN_RCONST(1.0e-4);
  const sunrealtype first_growth = SUN_RCONST(1234.0);
  const sunrealtype growth       = SUN_RCONST(3.0);

  ARKodeSetInitStep(arkode_mem, h0);
  ARKodeSetMaxFirstGrowth(arkode_mem, first_growth);
  ARKodeSetMaxGrowth(arkode_mem, growth);

  sunrealtype h_expect = first_growth * h0;
  /* Take several steps to see the special behavior at step one then to allow
     the adaptivity controller history fill up */
  const int num_steps = 5;
  for (int step = 1; step <= num_steps; step++)
  {
    retval = check_step(arkode_mem, y, h_expect, step);
    if (retval != 0) { return retval; }
    h_expect *= growth;
  }

  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}
