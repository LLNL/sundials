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
 * -----------------------------------------------------------------*/

#include <memory>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>

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

using GkoMatrixType      = gko::matrix::Dense<sunrealtype>;
using GkoBatchMatrixType = gko::matrix::BatchDense<sunrealtype>;
using SUNMatrixType      = sundials::ginkgo::BlockMatrix<GkoBatchMatrixType>;
// using SUNMatrixIdentityType = sundials::ginkgo::BlockMatrix<gko::matrix::Identity<sunrealtype>>;

/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;                 /* counter for test failures  */
  sunindextype matrows, matcols; /* matrix dimensions          */
  sunindextype num_blocks;       /* number of matrix blocks    */
  N_Vector x, y;                 /* test vectors               */
  realtype* xdata;               /* pointers to vector data    */
  int square;
  sundials::Context sunctx;

  auto gko_exec =
      OMP_OR_HIP_OR_CUDA(gko::OmpExecutor::create(), gko::HipExecutor::create(0, gko::OmpExecutor::create(), true),
                         gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true));

  /* check input and set vector length */
  if (argc < 3) {
    printf("ERROR: TWO (2) Input required: matrix rows, number of blocks\n");
    return (-1);
  }

  matrows = (sunindextype)atol(argv[1]);
  if (matrows <= 0) {
    printf("ERROR: number of rows must be a positive integer \n");
    return (-1);
  }
  matcols = matrows;

  num_blocks = (sunindextype)atol(argv[2]);
  if (num_blocks < 1) {
    printf("ERROR: number of blocks must >= 1 \n");
    return (-1);
  }

  SetTiming(0);

  square = (matrows == matcols) ? 1 : 0;
  printf("\n SUNMATRIX_GINKGOBLOCKDENSE test: size %ld by %ld, num_blocks %ld\n\n", (long int)matrows,
         (long int)matcols, (long int)num_blocks);

  // Create vectors and matrices
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, matrows);

  x = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(num_blocks * matcols, sunctx), N_VNew_Hip(num_blocks * matcols, sunctx),
                         N_VNew_Cuda(num_blocks * matcols, sunctx));
  y = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(num_blocks * matrows, sunctx), N_VNew_Hip(num_blocks * matrows, sunctx),
                         N_VNew_Cuda(num_blocks * matrows, sunctx));

  // Create batch_dim object to describe the dimensions of the batch matrix.
  auto batch_mat_size = gko::batch_dim<>(num_blocks, gko::dim<2>(matrows, matcols));
  auto batch_vec_size = gko::batch_dim<>(num_blocks, gko::dim<2>(matrows, 1));

  // Fill matrix
  auto gko_matrix = GkoMatrixType::create(gko_exec, batch_mat_size.at(0));
  gko_matrix->read(gko::matrix_data<sunrealtype>(batch_mat_size.at(0), distribution, generator));
  auto gko_batch_matrix = gko::share(GkoBatchMatrixType::create(gko_exec, num_blocks, gko_matrix.get()));

  // Fill identity matrix
  auto gko_ident = gko::share(GkoMatrixType::create(gko_exec, batch_mat_size.at(0)));
  gko_ident->read(gko::matrix_data<sunrealtype>::diag(batch_mat_size.at(0), 1.0));
  auto gko_batch_ident = gko::share(GkoBatchMatrixType::create(gko_exec, num_blocks, gko_ident.get()));

  // Wrap ginkgo matrices for SUNDIALS.
  // sundials::ginkgo::Matrix is overloaded to a SUNMatrix.
  SUNMatrixType A{gko_batch_matrix, sunctx};
  SUNMatrixType I{gko_batch_ident, sunctx};

  xdata = N_VGetArrayPointer(x);
  for (sunindextype i = 0; i < num_blocks * matcols; i++) {
    xdata[i] = distribution(generator);
  }
  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  // Compute true solution
  SUNMatrix Aref = SUNDenseMatrix(num_blocks * matrows, num_blocks * matcols, sunctx);
  for (sunindextype blocki = 0; blocki < num_blocks; blocki++) {
    for (sunindextype i = 0; i < matrows; i++) {
      for (sunindextype j = 0; j < matcols; j++) {
        SM_ELEMENT_D(Aref, blocki * matrows + i, blocki * matcols + j) = gko_batch_matrix->at(blocki, i, j);
      }
    }
  }
  SUNMatMatvec_Dense(Aref, x, y);
  SUNMatDestroy(Aref);

  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(y), N_VCopyToDevice_Cuda(y));

  // Do the SUNMatrix tests
  fails += Test_SUNMatGetID(A, SUNMATRIX_GINKGOBLOCKDENSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  if (square) {
    if (A.get()->ops->scaleadd)
      fails += Test_SUNMatScaleAdd(A, I, 0);
    fails += Test_SUNMatScaleAddI(A, I, 0);
  }
  fails += Test_SUNMatMatvecSetup(A, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  // Print result
  if (fails)
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
  else
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");

  // Free vectors and matrices
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
  sunindextype blocks = Amat->get_num_batch_entries();
  sunindextype rows   = Amat->get_size().at(0)[0];
  sunindextype cols   = Amat->get_size().at(0)[1];

  /* check dimensions */
  if (Amat->get_size() != Bmat->get_size()) {
    printf(">>> ERROR: check_matrix: Different data array lengths \n");
    return (1);
  }

  /* compare data */
  for (sunindextype blocki = 0; blocki < blocks; blocki++) {
    for (sunindextype i = 0; i < rows; i++) {
      for (sunindextype j = 0; j < cols; j++) {
        int check = SUNRCompareTol(Amat->at(blocki, i, j), Bmat->at(blocki, i, j), tol);
        if (check) {
          printf("  A[%ld,%ld,%ld] = %g != B[%ld,%ld,%ld] = %g (err = %g)\n", (long int)blocki, (long int)i,
                 (long int)j, Amat->at(blocki, i, j), (long int)blocki, (long int)i, (long int)j,
                 Bmat->at(blocki, i, j), SUNRabs(Amat->at(blocki, i, j) - Bmat->at(blocki, i, j)));
          failure += check;
        }
      }
    }
  }

  return failure > 0;
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure         = 0;
  auto Amat           = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  sunindextype blocks = Amat->get_num_batch_entries();
  sunindextype rows   = Amat->get_size().at(0)[0];
  sunindextype cols   = Amat->get_size().at(0)[1];

  /* compare data */
  for (sunindextype blocki = 0; blocki < blocks; blocki++) {
    for (sunindextype i = 0; i < rows; i++) {
      for (sunindextype j = 0; j < cols; j++) {
        int check = SUNRCompareTol(Amat->at(blocki, i, j), val, tol);
        if (check) {
          printf("  A[%ld,%ld,%ld] = %g != %g (err = %g)\n", (long int)blocki, (long int)i, (long int)j,
                 Amat->at(blocki, i, j), val, SUNRabs(Amat->at(blocki, i, j) - val));
          failure += check;
        }
      }
    }
  }

  return failure > 0;
}

int check_vector(N_Vector expected, N_Vector computed, realtype tol)
{
  int failure = 0;
  realtype *compu, *expec;
  sunindextype c_len, e_len;
  sunindextype i;

  /* copy vectors to host */
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(computed), N_VCopyFromDevice_Cuda(computed));
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(expected), N_VCopyFromDevice_Cuda(expected));

  /* get vector data */
  compu = N_VGetArrayPointer(computed);
  expec = N_VGetArrayPointer(expected);

  /* check data lengths */
  c_len = N_VGetLength(computed);
  e_len = N_VGetLength(expected);

  if (c_len != e_len) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return (1);
  }

  /* check vector data */
  for (i = 0; i < c_len; i++)
    failure += SUNRCompareTol(compu[i], expec[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for (i = 0; i < c_len; i++)
      if (SUNRCompareTol(compu[i], expec[i], tol) != 0)
        printf("  computed[%ld] = %g != %e (err = %g)\n", (long int)i, compu[i], expec[i], SUNRabs(compu[i] - expec[i]));
  }

  return failure > 0;
}

booleantype has_data(SUNMatrix A)
{
  auto Amat = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  if (Amat->get_values() == NULL || Amat->get_size().get_num_batch_entries() == 0)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  auto Amat = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  for (sunindextype iblk = 0; iblk < Amat->get_size().get_num_batch_entries(); iblk++) {
    if (Amat->get_size().at(iblk)[0] != Amat->get_size().at(iblk)[1])
      return SUNFALSE;
  }
  return SUNTRUE;
}

void sync_device(SUNMatrix A) { OMP_OR_HIP_OR_CUDA(, hipDeviceSynchronize(), cudaDeviceSynchronize()); }
