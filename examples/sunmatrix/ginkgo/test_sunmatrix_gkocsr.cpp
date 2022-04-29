/*
 * -----------------------------------------------------------------
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#include "test_sunmatrix.h"

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

using GkoMatrixType = gko::matrix::Csr<sunrealtype>;
using SUNMatrixType = sundials::ginkgo::Matrix<GkoMatrixType>;

/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;                 /* counter for test failures  */
  sunindextype matrows, matcols; /* matrix dimensions          */
  N_Vector x, y;                 /* test vectors               */
  realtype* xdata;               /* pointers to vector data    */
  int square;
  sundials::Context sunctx;

  auto gko_exec =
      OMP_OR_HIP_OR_CUDA(gko::OmpExecutor::create(), gko::HipExecutor::create(0, gko::OmpExecutor::create(), true),
                         gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true));

  /* check input and set vector length */
  if (argc < 3) {
    printf("ERROR: TWO (2) Input required: matrix rows, matrix cols\n");
    return (-1);
  }

  matrows = (sunindextype)atol(argv[1]);
  if (matrows <= 0) {
    printf("ERROR: number of rows must be a positive integer \n");
    return (-1);
  }

  matcols = (sunindextype)atol(argv[2]);
  if (matcols <= 0) {
    printf("ERROR: number of cols must be a positive integer \n");
    return (-1);
  }

  SetTiming(0);

  square = (matrows == matcols) ? 1 : 0;
  printf("\n SUNMATRIX_GINKGOCSR test: size %ld by %ld\n\n", (long int)matrows, (long int)matcols);

  /* Create vectors and matrices */
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, matrows);

  x = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(matcols, sunctx), N_VNew_Hip(matcols, sunctx), N_VNew_Cuda(matcols, sunctx));
  y = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(matrows, sunctx), N_VNew_Hip(matrows, sunctx), N_VNew_Cuda(matrows, sunctx));

  auto matrix_dim  = gko::dim<2>(matrows, matcols);
  auto gko_matdata = gko::matrix_data<sunrealtype>(matrix_dim, distribution, generator);
  auto gko_matrix  = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
  auto gko_ident   = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));

  /* Fill matrices and vectors */
  gko_matrix->read(gko_matdata);
  if (square) {
    gko_ident->read(gko::matrix_data<sunrealtype>::diag(matrix_dim, 1.0));
  }

  /* Wrap ginkgo matrices for SUNDIALS.
     sundials::ginkgo::Matrix is overloaded to a SUNMatrix. */
  SUNMatrixType A{gko_matrix, sunctx};
  SUNMatrixType I{gko_ident, sunctx};

  xdata = N_VGetArrayPointer(x);
  for (sunindextype i = 0; i < matcols; i++) {
    xdata[i] = distribution(generator);
  }
  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  /* Compute true solution */
  SUNMatrix Aref      = SUNDenseMatrix(matrows, matcols, sunctx);
  const auto Arowptrs = gko_matrix->get_const_row_ptrs();
  const auto Acolidxs = gko_matrix->get_const_col_idxs();
  const auto Avalues  = gko_matrix->get_const_values();
  for (sunindextype irow = 0; irow < gko_matrix->get_size()[0]; irow++) {
    for (sunindextype inz = Arowptrs[irow]; inz < Arowptrs[irow + 1]; inz++) {
      SM_ELEMENT_D(Aref, irow, Acolidxs[inz]) = Avalues[inz];
    }
  }
  SUNMatMatvec_Dense(Aref, x, y);
  SUNMatDestroy(Aref);

  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(y), N_VCopyToDevice_Cuda(y));

  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_GINKGOCSR, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  if (square) {
    fails += Test_SUNMatScaleAdd(A, I, 0);
    fails += Test_SUNMatScaleAddI(A, I, 0);
  }
  fails += Test_SUNMatMatvecSetup(A, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  /* Print result */
  if (fails)
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
  else
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");

  /* Free vectors and matrices */
  N_VDestroy(x);
  N_VDestroy(y);

  return (fails);
}

/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure         = 0;
  auto Amat           = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  auto Bmat           = static_cast<SUNMatrixType*>(B->content)->gkomtx();
  const auto Arowptrs = Amat->get_const_row_ptrs();
  const auto Acolidxs = Amat->get_const_col_idxs();
  const auto Avalues  = Amat->get_const_values();
  const auto Browptrs = Bmat->get_const_row_ptrs();
  const auto Bcolidxs = Bmat->get_const_col_idxs();
  const auto Bvalues  = Bmat->get_const_values();

  /* check lengths */
  if (Amat->get_size() != Bmat->get_size()) {
    printf(">>> ERROR: check_matrix: Different data array lengths \n");
    return (1);
  }

  /* compare data */
  for (sunindextype irow = 0; irow < Amat->get_size()[0]; irow++) {
    for (sunindextype inz = Arowptrs[irow]; inz < Arowptrs[irow + 1]; inz++) {
      failure += SUNRCompareTol(Avalues[inz], Bvalues[inz], tol);
    }
  }

  return failure > 0;
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure         = 0;
  auto Amat           = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  const auto Arowptrs = Amat->get_const_row_ptrs();
  const auto Acolidxs = Amat->get_const_col_idxs();
  const auto Avalues  = Amat->get_const_values();

  /* compare data */
  for (sunindextype irow = 0; irow < Amat->get_size()[0]; irow++) {
    for (sunindextype inz = Arowptrs[irow]; inz < Arowptrs[irow + 1]; inz++) {
      int check = SUNRCompareTol(Avalues[inz], val, tol);
      if (check) {
        printf("  actual = %g != %e (err = %g)\n", Avalues[inz], val, SUNRabs(Avalues[inz] - val));
        failure += check;
      }
    }
  }

  return failure > 0;
}

int check_vector(N_Vector expected, N_Vector computed, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(computed), N_VCopyFromDevice_Cuda(computed));
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(expected), N_VCopyFromDevice_Cuda(expected));

  /* get vector data */
  xdata = N_VGetArrayPointer(computed);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(computed);
  yldata = N_VGetLength(expected);

  if (xldata != yldata) {
    printf(">>> [ERROR] check_vector: different vector lengths \n");
    return (1);
  }

  /* check vector data */
  for (i = 0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for (i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  computed[%ld] = %g != %e (err = %g)\n", (long int)i, xdata[i], ydata[i], SUNRabs(xdata[i] - ydata[i]));
  }

  return failure > 0;
}

booleantype has_data(SUNMatrix A)
{
  auto Amat = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  if (Amat->get_values() == NULL || Amat->get_size()[0] == 0 || Amat->get_size()[1] == 0)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  auto Amat = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  if (Amat->get_size()[0] == Amat->get_size()[1])
    return SUNTRUE;
  else
    return SUNFALSE;
  return SUNTRUE;
}

void sync_device(SUNMatrix A) { OMP_OR_HIP_OR_CUDA(, hipDeviceSynchronize(), cudaDeviceSynchronize()); }
