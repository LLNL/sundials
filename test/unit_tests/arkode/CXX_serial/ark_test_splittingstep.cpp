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
 * Unit tests on several ODEs with analytical solutions to verify the
 * SplittingStep module.
 * ---------------------------------------------------------------------------*/

#include <sundials/sundials_context.hpp>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_stepper.h>

#include <iostream>
#include <vector>
#include <cmath>

/* Integrates the ODE
 *
 * y' = \sum_{i=0}^{P-1} 2^i / (1 - 2^P) * y,    y(0) = 1
 * 
 * with a 3rd order operator splitting method. Each partition is solved by an
 * ERK method. Using the exact solution y(t) = exp(-t), we confirm the numerical
 * solution is sufficiently accurate.
 */
static int test_forward(sundials::Context &ctx, const int partitions) {
  constexpr auto t0 = SUN_RCONST(0.0);
  constexpr auto tf = SUN_RCONST(1.0);
  constexpr auto dt = SUN_RCONST(8.0e-3);
  constexpr auto local_tol = SUN_RCONST(1.0e-6);
  constexpr auto global_tol = 10 * local_tol;
  const auto y = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  const ARKRhsFn f = [](const sunrealtype, const N_Vector z, const N_Vector zdot, void *const user_data) {
    const auto lambda = *static_cast<sunrealtype*>(user_data);
    N_VScale(lambda, z, zdot);
    return 0;
  };

  std::vector<void *> partition_mem(partitions);
  std::vector<sunrealtype> lambda(partitions);
  std::vector<SUNStepper> steppers(partitions);
  for (int i = 0; i < partitions; i++) {
    partition_mem[i] = ARKStepCreate(f, nullptr, t0, y, ctx);
    /* The lambdas sum up to 1 */
    lambda[i] = std::pow(SUN_RCONST(2.0), i) / (1 - std::pow(SUN_RCONST(2.0), partitions));
    ARKodeSetUserData(partition_mem[i], &lambda[i]);
    ARKodeSStolerances(partition_mem[i], local_tol, local_tol);
    ARKStepCreateSUNStepper(partition_mem[i], &steppers[i]);
  }

  auto arkode_mem = SplittingStepCreate(steppers.data(), partitions, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, 3);
  auto tret = t0;
  ARKodeEvolve(arkode_mem, tf, y, &tret, ARK_NORMAL);

  /* Check that the solution matches the exact solution */
  const auto exact_solution = std::exp(t0 - tf);
  if (SUNRCompareTol(exact_solution, NV_Ith_S(y, 0), global_tol)) {
    std::cerr << "Forward solution with " << partitions << " partitions is not within the tolerance\n";
    return 1;
  }

  N_VDestroy(y);
  for (int i = 0; i <partitions; i++) {
    ARKodeFree(&partition_mem[i]);
    SUNStepper_Destroy(&steppers[i]);
  }
  ARKodeFree(&arkode_mem);

  return 0;
}

/* Integrates the ODE
 * 
 * y_1' = y_2 - t
 * y_2' = t - y_1
 * 
 * with initial condition y(0) = [1, 1]^T backward in time to t = -1, then
 * forward to t = 0.4, and backward to t = 0. We apply a splitting method using
 * a component partitioning and check that the numerical solution is close to
 * the original initial condition.
 */
static int test_mixed_directions(const sundials::Context &ctx) {
  constexpr auto t0 = SUN_RCONST(0.0);
  constexpr auto t1 = SUN_RCONST(-1.0);
  constexpr auto t2 = SUN_RCONST(0.4);
  constexpr auto t3 = t0;
  constexpr auto dt = -SUN_RCONST(1.23e-3);
  constexpr auto local_tol = SUN_RCONST(1.0e-6);
  constexpr auto global_tol = 10 * local_tol;
  const auto y = N_VNew_Serial(2, ctx);
  N_VConst(SUN_RCONST(1.0), y);
  const auto err = N_VClone(y);
  N_VConst(SUN_RCONST(1.0), err);

  const ARKRhsFn f1 = [](const sunrealtype t, const N_Vector z, const N_Vector zdot, void *const) {
    NV_Ith_S(zdot, 0) = NV_Ith_S(z, 1) - t;
    NV_Ith_S(zdot, 1) = SUN_RCONST(0.0);
    return 0;
  };

  const ARKRhsFn f2 = [](const sunrealtype t, const N_Vector z, const N_Vector zdot, void *const) {
    NV_Ith_S(zdot, 0) = SUN_RCONST(0.0);
    NV_Ith_S(zdot, 1) = t - NV_Ith_S(z, 0);
    return 0;
  };

  void* parititon_mem[] = {
    ARKStepCreate(f1, nullptr, t0, y, ctx),
    ARKStepCreate(f2, nullptr, t0, y, ctx)
  };
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);
  
  SUNStepper steppers[] = {nullptr, nullptr};
  ARKStepCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKStepCreateSUNStepper(parititon_mem[1], &steppers[1]);

  auto arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
  ARKodeSetMaxNumSteps(arkode_mem, -1);
  const auto coefficients = SplittingStepCoefficients_LoadCoefficients(ARKODE_SPLITTING_RUTH_3_3_2);
  SplittingStep_SetCoefficients(arkode_mem, coefficients);
  SplittingStepCoefficients_Free(coefficients);

  /* Integrate from 0 to -1 */
  auto tret = t0;
  ARKodeEvolve(arkode_mem, t1, y, &tret, ARK_NORMAL);

  /* Integrate from -1 to 0.4 */
  ARKodeReset(arkode_mem, t1, y);
  ARKodeSetStepDirection(arkode_mem, t2 - t1);
  ARKodeEvolve(arkode_mem, t2, y, &tret, ARK_NORMAL);

  /* Integrate from 0.4 to 0 */
  ARKodeReset(arkode_mem, t2, y);
  ARKodeSetStepDirection(arkode_mem, t3 - t2);
  ARKodeEvolve(arkode_mem, t3, y, &tret, ARK_NORMAL);

  N_VLinearSum(SUN_RCONST(1.0), err, -SUN_RCONST(1.0), y, err);
  const auto max_err = N_VMaxNorm(err);

  if (max_err > global_tol) {
    std::cerr << "Mixed direction solution failed with an error of " << max_err << "\n";
    return 1;
  }

  N_VDestroy(y);
  N_VDestroy(err);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);

  return 0;
}

/* Integrates the ODE
 * 
 * y' = [y^2] - [y*(y - 1)] = -y
 * 
 * with initial condition y(0) = 1 and partitioning specified by the square
 * brackets. Powers and products are done componentwise. At t = 0.5, we resize
 * from 1D to 2D:
 * 
 * y_new = [1, y_old]^T
 * 
 * Then we integrate to t = 1 and check the error against the exact solution
 * y(1) = [exp(-0.5), exp(-1)].
 */
static int test_resize(const sundials::Context &ctx) {
  constexpr auto t0 = SUN_RCONST(0.0);
  constexpr auto t1 = SUN_RCONST(0.5);
  constexpr auto t2 = SUN_RCONST(1.0);
  constexpr auto dt = SUN_RCONST(8.0e-3);
  constexpr auto local_tol = SUN_RCONST(1.0e-5);
  constexpr auto global_tol = local_tol;
  const auto y = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  const ARKRhsFn f1 = [](const sunrealtype t, const N_Vector z, const N_Vector zdot, void *const) {
    N_VProd(z, z, zdot);
    return 0;
  };

  const ARKRhsFn f2 = [](const sunrealtype t, const N_Vector z, const N_Vector zdot, void *const) {
    N_VConst(-SUN_RCONST(1.0), zdot);
    N_VLinearSum(SUN_RCONST(1.0), zdot, -SUN_RCONST(1.0), z, zdot);
    N_VProd(zdot, z, zdot);
    return 0;
  };

  void* parititon_mem[] = {
    ARKStepCreate(f1, nullptr, t0, y, ctx),
    ARKStepCreate(f2, nullptr, t0, y, ctx)
  };
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);
  
  SUNStepper steppers[] = {nullptr, nullptr};
  ARKStepCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKStepCreateSUNStepper(parititon_mem[1], &steppers[1]);

  auto arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  const auto coefficients = SplittingStepCoefficients_SymmetricParallel(2);
  SplittingStep_SetCoefficients(arkode_mem, coefficients);
  SplittingStepCoefficients_Free(coefficients);

  /* Integrate from 0 to 0.5 */
  auto tret = t0;
  ARKodeEvolve(arkode_mem, t1, y, &tret, ARK_NORMAL);

  /* Resize */
  const auto y_new = N_VNew_Serial(2, ctx);
  NV_Ith_S(y_new, 0) = SUN_RCONST(1.0);
  NV_Ith_S(y_new, 1) = NV_Ith_S(y, 0);
  N_VDestroy(y);
  ARKodeResize(arkode_mem, y_new, SUN_RCONST(1.0), t1, nullptr, nullptr);
  ARKodeResize(parititon_mem[0], y_new, SUN_RCONST(1.0), t1, nullptr, nullptr);
  ARKodeResize(parititon_mem[1], y_new, SUN_RCONST(1.0), t1, nullptr, nullptr);

  /* Integrate from 0.5 to 1 */
  ARKodeEvolve(arkode_mem, t2, y_new, &tret, ARK_NORMAL);

  const auto err = N_VClone(y_new);
  NV_Ith_S(err, 0) = std::exp(t1 - t2);
  NV_Ith_S(err, 1) = std::exp(t0 - t2);
  N_VLinearSum(SUN_RCONST(1.0), err, -SUN_RCONST(1.0), y_new, err);
  const auto max_err = N_VMaxNorm(err);

  if (max_err > global_tol) {
    std::cerr << "Resized solution failed with an error of " << max_err << "\n";
    return 1;
  }

  N_VDestroy(y_new);
  N_VDestroy(err);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);

  return 0;
}

/* Creates a custom SUNStepper for the linear, scalar ODE y' = lambda*y */
static SUNStepper create_exp_stepper(const sundials::Context &ctx, const sunrealtype &lambda) {
  SUNStepper stepper = nullptr;
  SUNStepper_Create(ctx, &stepper);

  const auto empty_func = [](auto... args) {
    return 0;
  };
  SUNStepper_SetResetFn(stepper, empty_func);
  SUNStepper_SetStopTimeFn(stepper, empty_func);
  SUNStepper_SetSetStepDirectionFn(stepper, empty_func);
  SUNStepper_SetFullRhsFn(stepper, empty_func);
  SUNStepper_SetEvolveFn(stepper, [](const SUNStepper s, const sunrealtype t0,
                                    const sunrealtype tout, const N_Vector y,
                                    sunrealtype* const tret, int* const stop_reason) {
    void *content = nullptr;
    SUNStepper_GetContent(s, &content);
    const auto lam = *static_cast<sunrealtype*>(content);
    N_VScale(std::exp(lam * (tout - t0)), y, y);
    *tret = tout;
    *stop_reason = 0;
    return 0;
  });
  SUNStepper_SetContent(stepper, static_cast<void*>(const_cast<sunrealtype*>(&lambda)));
  return stepper;
}

/* Integrates the ODE
 * 
 * y' = -0.6*y - 0.4*y
 * 
 * with initial condition y(0) = 1 with a sixth order splitting method and
 * exact, custom SUNSteppers for the two partitions. We integrate to t = 1 and
 * compare the numerical solution to the exact solution y(t) = exp(-t).
 */
static int test_custom_stepper(const sundials::Context &ctx) {
  constexpr auto lambda1 = SUN_RCONST(-0.6);
  constexpr auto lambda2 = SUN_RCONST(-0.4);
  constexpr auto t0 = SUN_RCONST(0.0);
  constexpr auto tf = SUN_RCONST(1.0);
  constexpr auto dt = SUN_RCONST(0.1);
  const auto y = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  SUNStepper steppers[] = {
    create_exp_stepper(ctx, lambda1),
    create_exp_stepper(ctx, lambda2)
  };

  auto arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  const auto coefficients = SplittingStepCoefficients_SuzukiFractal(2, 6);
  SplittingStep_SetCoefficients(arkode_mem, coefficients);
  SplittingStepCoefficients_Free(coefficients);

  auto tret = t0;
  ARKodeEvolve(arkode_mem, tf, y, &tret, ARK_NORMAL);

  /* Check that the solution matches the exact solution */
  const auto exact_solution = std::exp(t0 - tf);
  if (SUNRCompare(exact_solution, NV_Ith_S(y, 0))) {
    std::cerr << "Custom SUNStepper is not within the tolerance\n";
    return 1;
  }

  N_VDestroy(y);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);
  return 0;
}

/* Error handling function which prints the error and exits the program */
static void err_fn(const int line, const char* const func, const char* const file,
                   const char* const msg, const SUNErrCode err_code,
                   void* const, const SUNContext)
{
  std::cerr << "Error at line " << line << " of " << func << " in " << file << ": " << msg << "\n";
  exit(err_code);
}

int main() {
  sundials::Context ctx;
  SUNContext_PushErrHandler(ctx, err_fn, nullptr);

  int errors = 0;

  /* Run the tests */
  constexpr auto min_partitions = 2;
  constexpr auto max_partitions = 7;
  for (auto p = min_partitions; p <= max_partitions; p++) {
    errors += test_forward(ctx, p);
  }

  errors += test_mixed_directions(ctx);
  errors += test_resize(ctx);
  errors += test_custom_stepper(ctx);

  if (errors == 0) {
    std::cout << "Success\n";
  } else {
    std::cout << errors << " Test Failures\n";
  }
  
  return errors;
}
