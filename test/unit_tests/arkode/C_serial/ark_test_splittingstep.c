/* Things to test
Custom inner stepper with exact solution
Orders 1-4
Resizing
Changing coefficients
~3 partitions
Check n_stepper_evolves
*/

/*
Examples
- Fractional Step RK
*/

/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
 * TODO
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>

static int f1(sunrealtype t, N_Vector y, N_Vector yp, void *data) {
  N_VScale(-1, y, yp);
  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  int retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);

  N_Vector y = N_VNew_Serial(4, sunctx);
  N_VConst(1, y);

  void *arkode_mem = ARKStepCreate(f1, NULL, 0, y, sunctx);

  SUNStepper stepper;
  ARKStepCreateMRIStepInnerStepper(arkode_mem, &stepper);

  void *split_mem = SplittingStepCreate(&stepper, 1, 0, y, sunctx);
  ARKodeSetFixedStep(split_mem, 0.25); // TODO: fix valgrind error is this is not called
  sunrealtype tret;
  ARKodeEvolve(split_mem, 1, y, &tret, ARK_NORMAL);

  N_VPrint(y); // TODO: why is valgrind showing errors for these prints?

  ARKodeFree(&split_mem);
  N_VDestroy(y);
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&arkode_mem);
  SUNContext_Free(&sunctx);
}
