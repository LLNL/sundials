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
  printf("%i\n", retval);

  N_Vector y = N_VNew_Serial(4, sunctx);
  N_VConst(1, y);

  void *arkode_mem = ARKStepCreate(f1, NULL, 0, y, sunctx);

  SUNStepper stepper;
  ARKStepCreateMRIStepInnerStepper(arkode_mem, &stepper);



  N_VDestroy(y);
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&arkode_mem);
  SUNContext_Free(&sunctx);
}
