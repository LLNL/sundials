/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol Dense module
 * implementation.
 * ----------------------------------------------------------------- */

#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>

#include "test_sunlinsol.h"

#if defined(USE_HIP)
#define OMP_OR_HIP_OR_CUDA(a, b, c) b
#elif defined(USE_CUDA)
#define OMP_OR_HIP_OR_CUDA(a, b, c) c
#else
#define OMP_OR_HIP_OR_CUDA(a, b, c) a
#endif

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#else
#include <nvector/nvector_serial.h>
#endif

#if defined(USE_CSR)
using GkoMatrixType       = gko::matrix::Csr<sunrealtype>;
using SUNMatrixType       = sundials::ginkgo::Matrix<GkoMatrixType>;
#else
using GkoMatrixType       = gko::matrix::Dense<sunrealtype>;
using SUNMatrixType       = sundials::ginkgo::Matrix<GkoMatrixType>;
#endif

using GkoSolverType       = gko::solver::Gmres<sunrealtype>;
using SUNLinearSolverType = sundials::ginkgo::LinearSolver<GkoSolverType, SUNMatrixType>;

/* ----------------------------------------------------------------------
 * SUNLinSol_Ginkgo Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;                 /* counter for test failures    */
  sunrealtype matcond;           /* matrix condition number      */
  sunindextype matcols, matrows; /* matrix columns, matrows      */
  sunindextype nblocks;          /* number of matrix blocks      */
  sunindextype max_iters;        /* maximum number of iterations */
  N_Vector x, y, b;              /* test vectors                 */
  int print_timing;
  realtype* xdata;
  sundials::Context sunctx;

  auto gko_exec =
      OMP_OR_HIP_OR_CUDA(gko::OmpExecutor::create(), gko::HipExecutor::create(0, gko::OmpExecutor::create(), true),
                         gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true));

  /* check input and set matrix dimensions */
  if (argc < 6) {
    printf("ERROR: FIVE (5) Inputs required: number of matrix columns, condition number, number of blocks, "
           "max iterations, print timing \n");
    return (-1);
  }

  matcols = (sunindextype)atol(argv[1]);
  if (matcols <= 0) {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return (-1);
  }
  matrows = matcols;

  matcond = (sunrealtype)atof(argv[2]);
  if (matcond <= 0) {
    printf("ERROR: matrix condition number must be positive\n");
    return (-1);
  }

  nblocks = (sunindextype)atol(argv[3]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return (-1);
  }

  max_iters = (sunindextype)atol(argv[4]);
  if (max_iters <= 0) {
    printf("ERROR: max iterations must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[5]);
  SetTiming(print_timing);

  printf("\nGinkgo linear solver test: size = %ld x %ld, blocks = %ld, condition = %g, max iters. = %ld\n\n",
         (long int)matrows, (long int)matcols, (long int)nblocks, matcond, (long int)max_iters);

  /* Create vectors and matrices */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution_real(8, 10);

  x = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(matcols, sunctx), N_VNew_Hip(matcols, sunctx), N_VNew_Cuda(matcols, sunctx));
  b = N_VClone(x);
  y = N_VClone(x);

  auto matrix_dim = gko::dim<2>(matrows, matcols);
  auto gko_matdata = gko::matrix_data<sunrealtype, sunindextype>::cond(matrows, gko::remove_complex<sunrealtype>(matcond),
                                                                       distribution_real, engine);
  auto gko_matrixA = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
  auto gko_matrixB = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
  auto gko_ident   = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));

  /* Fill matrices */
  gko_matdata.remove_zeros();
  gko_matrixA->read(gko_matdata);
  gko_matrixB->read(gko_matdata);
  gko_ident->read(gko::matrix_data<sunrealtype>::diag(matrix_dim, 1.0));

  /* Fill x with random data */
  xdata = N_VGetArrayPointer(x);
  for (sunindextype i = 0; i < matcols; i++) {
    xdata[i] = distribution_real(engine);
  }
  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  /* Wrap ginkgo matrices for SUNDIALS.
     sundials::ginkgo::Matrix is overloaded to a SUNMatrix. */
  SUNMatrixType A{gko_matrixA, sunctx};
  SUNMatrixType B{gko_matrixB, sunctx};

  /* You can also create the SUNMatrix via cloning */
  SUNMatrix I = SUNMatClone(A);

  /* Copy A and x into B and y to print in case of solver failure */
  SUNMatCopy(A, B);
  N_VScale(ONE, x, y);

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A);
  fails += SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

    SUNMatDestroy(I);
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);

    return (1);
  }

  /*
   * Create linear solver.
   */

  /* Use default stopping critieria */
  auto crit = sundials::ginkgo::DefaultStop::build() //
                  .with_max_iters(max_iters)         //
                  .on(gko_exec);

  auto precon = gko::preconditioner::Jacobi<sunrealtype, sunindextype>::build() //
                    .on(gko_exec);

  auto gko_solver_factory = GkoSolverType::build()                       //
                                .with_criteria(gko::share(crit))         //
                                .with_preconditioner(gko::share(precon)) //
                                .on(gko_exec);

  SUNLinearSolverType LS{gko::share(gko_solver_factory), sunctx};
  LS.logResNorm(true); // turn on residual norm logging

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS.get(), SUNLINEARSOLVER_GINKGO, 0);
  fails += Test_SUNLinSolGetType(LS.get(), SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS.get(), 0);
  fails += Test_SUNLinSolSetup(LS.get(), A, 0);
  fails += Test_SUNLinSolSolve(LS.get(), A, x, b, 1000 * std::numeric_limits<sunrealtype>::epsilon(), SUNTRUE, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    /* only print if its a small problem */
    if (matcols < 4) {
      printf("\nx (original) =\n");
      // N_VPrint(y);
      printf("\nx (computed) =\n");
      N_VPrint(x);
      printf("\nb =\n");
      N_VPrint(b);
    }
  } else {
    printf("\nSUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Print solve information */
  printf("Number of linear solver iterations: %ld\n", static_cast<long int>(SUNLinSolNumIters(LS)));
  printf("Final residual norm: %g\n", SUNLinSolResNorm(LS));

  /* Free solver, matrix and vectors */
  SUNMatDestroy(I);
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);

  return (fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector expected, N_Vector actual, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(actual), N_VCopyFromDevice_Cuda(actual));
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(expected), N_VCopyFromDevice_Cuda(expected));

  /* get vector data */
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);

  if (xldata != yldata) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return (1);
  }

  /* check vector data */
  for (i = 0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for (i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int)i, xdata[i], ydata[i], SUNRabs(xdata[i] - ydata[i]));
  }

  return failure > 0;
}

void sync_device() { OMP_OR_HIP_OR_CUDA(, hipDeviceSynchronize(), cudaDeviceSynchronize()); }
