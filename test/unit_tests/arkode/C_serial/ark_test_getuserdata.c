/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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

#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_mristep.h"
#include "nvector/nvector_serial.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Dummy user-supplied function */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  return 0;
}

/* Main program */
int main(int argc, char* argv[])
{
  int retval                        = 0;
  SUNContext sunctx                 = NULL;
  N_Vector y                        = NULL;
  void* arkode_mem                  = NULL;
  void* arkode_inner_mem            = NULL;
  MRIStepInnerStepper inner_stepper = NULL;
  int udata_in                      = 1;
  void* udata_out                   = NULL;

  /* Create the SUNDIALS context object for this simulation. */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", retval);
    return 1;
  }

  /* Create solution vector and initialize to zero */
  y = N_VNew_Serial(1, sunctx);
  if (!y)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ONE, y);

  /* ------------ *
   * Test ARKStep *
   * ------------ */

  /* Create ARKODE mem structure */
  arkode_mem = ARKStepCreate(f, NULL, ZERO, y, sunctx);
  if (!arkode_mem)
  {
    fprintf(stderr, "ARKStepCreate returned NULL\n");
    return 1;
  }

  /* Set user data */
  retval = ARKStepSetUserData(arkode_mem, &udata_in);
  if (retval)
  {
    fprintf(stderr, "ARKStepSetUserData returned %i\n", retval);
    return 1;
  }

  /* Get user data */
  retval = ARKStepGetUserData(arkode_mem, &udata_out);
  if (retval)
  {
    fprintf(stderr, "ARKStepGetUserData returned %i\n", retval);
    return 1;
  }

  if (!udata_out)
  {
    fprintf(stderr, "udata_out is NULL\n");
    return 1;
  }

  if (&udata_in != (int*)udata_out)
  {
    fprintf(stderr, "udata_in != udata_out\n");
    return 1;
  }

  /* Reset output pointer */
  udata_out = NULL;

  /* Free integrator memory */
  ARKStepFree(&arkode_mem);

  /* ------------ *
   * Test ERKStep *
   * ------------ */

  /* Create ARKODE mem structure */
  arkode_mem = ERKStepCreate(f, ZERO, y, sunctx);
  if (!arkode_mem)
  {
    fprintf(stderr, "ERKStepCreate returned NULL\n");
    return 1;
  }

  /* Set user data */
  retval = ERKStepSetUserData(arkode_mem, &udata_in);
  if (retval)
  {
    fprintf(stderr, "ERKStepSetUserData returned %i\n", retval);
    return 1;
  }

  /* Get user data */
  retval = ERKStepGetUserData(arkode_mem, &udata_out);
  if (retval)
  {
    fprintf(stderr, "ERKStepGetUserData returned %i\n", retval);
    return 1;
  }

  if (!udata_out)
  {
    fprintf(stderr, "udata_out is NULL\n");
    return 1;
  }

  if (&udata_in != (int*)udata_out)
  {
    fprintf(stderr, "udata_in != udata_out\n");
    return 1;
  }

  /* Reset output pointer */
  udata_out = NULL;

  /* Free integrator memory */
  ERKStepFree(&arkode_mem);

  /* ------------ *
   * Test MRIStep *
   * ------------ */

  /* Create inner ARKODE mem structure */
  arkode_inner_mem = ARKStepCreate(f, NULL, ZERO, y, sunctx);
  if (!arkode_inner_mem)
  {
    fprintf(stderr, "ARKStepCreate returned NULL\n");
    return 1;
  }

  /* Create inner stepper */
  retval = ARKStepCreateMRIStepInnerStepper(arkode_inner_mem, &inner_stepper);
  if (retval)
  {
    fprintf(stderr, "ARKStepCreateMRIStepInnerStepper returned %i", retval);
    return 1;
  }

  /* Create ARKODE mem structure */
  arkode_mem = MRIStepCreate(f, NULL, ZERO, y, inner_stepper, sunctx);
  if (!arkode_mem)
  {
    fprintf(stderr, "ARKStepCreate returned NULL\n");
    return 1;
  }

  /* Set user data */
  retval = MRIStepSetUserData(arkode_mem, &udata_in);
  if (retval)
  {
    fprintf(stderr, "ARKStepSetUserData returned %i\n", retval);
    return 1;
  }

  /* Get user data */
  retval = MRIStepGetUserData(arkode_mem, &udata_out);
  if (retval)
  {
    fprintf(stderr, "ARKStepGetUserData returned %i\n", retval);
    return 1;
  }

  if (!udata_out)
  {
    fprintf(stderr, "udata_out is NULL\n");
    return 1;
  }

  if (&udata_in != (int*)udata_out)
  {
    fprintf(stderr, "udata_in != udata_out\n");
    return 1;
  }

  /* Reset output pointer */
  udata_out = NULL;

  /* Free integrator memory */
  MRIStepFree(&arkode_mem);
  ARKStepFree(&arkode_inner_mem);
  MRIStepInnerStepper_Free(&inner_stepper);

  /* Clean up */
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
