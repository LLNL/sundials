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
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_stepper.h>

#define DT_1   (0.5 / 8)
#define DT_2   (0.125 / 8)
#define T_END  1.0
#define LAMBDA 2.0

static int f1(sunrealtype t, N_Vector y, N_Vector yp, void* data)
{
  printf("%e %e\n", t, NV_Ith_S(y, 0));
  N_VScale(-LAMBDA, y, yp);
  return 0;
}

static int f2(sunrealtype t, N_Vector y, N_Vector yp, void* data)
{
  printf("    %e %e\n", t, NV_Ith_S(y, 0));
  N_VProd(y, y, yp);
  return 0;
}

static sunrealtype exact_sol(sunrealtype t, sunrealtype y0)
{
  return LAMBDA * y0 / (y0 - (y0 - LAMBDA) * exp(LAMBDA * t));
}

static void test()
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  N_Vector y = N_VNew_Serial(1, sunctx);
  N_VConst(1, y);

  void* mem = ARKStepCreate(f1, NULL, 0, y, sunctx);

  ARKodeSetFixedStep(mem, 0.25);

  sunrealtype tret = 0;
  ARKodeEvolve(mem, 1, y, &tret, ARK_NORMAL);

  ARKodeReset(mem, 1, y);
  ARKodeEvolve(mem, 0.5, y, &tret, ARK_NORMAL);

  N_VPrint(y);
}

int main(int argc, char* argv[])
{
  test();
  return 0;
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  N_Vector y = N_VNew_Serial(1, sunctx);
  N_VConst(1, y);
  const sunrealtype sol = exact_sol(T_END, NV_Ith_S(y, 0));

  void* arkode_mem1 = ARKStepCreate(f1, NULL, 0, y, sunctx);
  ARKodeSetOrder(arkode_mem1, 2);
  ARKodeSetFixedStep(arkode_mem1, DT_1);

  void* arkode_mem2 = ARKStepCreate(f2, NULL, 0, y, sunctx);
  ARKodeSetOrder(arkode_mem2, 2);
  ARKodeSetFixedStep(arkode_mem2, DT_2);

  SUNStepper steppers[2];
  ARKStepCreateSUNStepper(arkode_mem1, &steppers[0]);
  ARKStepCreateSUNStepper(arkode_mem2, &steppers[1]);

  void* split_mem = SplittingStepCreate(steppers, 2, 0, y, sunctx);
  SplittingStepCoefficients coeffs = SplittingStepCoefficients_Parallel(2);
  // SplittingStepCoefficients_Write(coeffs, stdout);
  SplittingStep_SetCoefficients(split_mem, coeffs);
  SplittingStepCoefficients_Free(coeffs);
  ARKodeSetFixedStep(split_mem,
                     DT_1);
  sunrealtype tret;
  ARKodeEvolve(split_mem, 1, y, &tret, ARK_NORMAL);
  printf("Final Solution: %e %e\n", T_END, NV_Ith_S(y, 0));

  N_VPrint(y);
  printf("Error: %e\n", sol - NV_Ith_S(y, 0));

  ARKodeFree(&split_mem);
  N_VDestroy(y);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem1);
  ARKodeFree(&arkode_mem2);
  SUNContext_Free(&sunctx);
}
