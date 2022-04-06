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
 * Part of this code has been derived from
 *    https://github.com/ginkgo-project/ginkgo/blob/batch-develop/examples/batched-solver/batched-solver.cpp
 * under the BSD-3 clause license provided below.
 * Ginkgo Copyright Start.
 * Copyright (c) 2017-2022, the Ginkgo authors
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * Ginkgo Copyright End.
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
#define OMP_OR_HIP_OR_CUDA(a,b,c) b
#elif defined(USE_CUDA)
#define OMP_OR_HIP_OR_CUDA(a,b,c) c
#else
#define OMP_OR_HIP_OR_CUDA(a,b,c) a
#endif

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#else
#include <nvector/nvector_serial.h>
#endif

using GkoMatrixType = gko::matrix::Csr<sunrealtype>;
using GkoBatchMatrixType = gko::matrix::BatchCsr<sunrealtype>;
using SUNMatrixType = sundials::ginkgo::BlockMatrix<GkoBatchMatrixType>;
// using SUNMatrixIdentityType = sundials::ginkgo::BlockMatrix<gko::matrix::Identity<sunrealtype>>;

struct ApplSysData {
  // Number of small systems in the batch.
  sunindextype nsystems;
  // Number of rows in each system.
  int nrows;
  // Number of non-zeros in each system matrix.
  int nnz;
  // Row pointers for one matrix
  sunindextype* row_ptrs;
  // Column indices of non-zeros for one matrix
  sunindextype* col_idxs;
  // Nonzero values for all matrices in the batch, concatenated
  sunrealtype* all_values;
};

// Generates a batch of tridiagonal systems.
ApplSysData appl_generate_system(const int nrows, const sunindextype nsystems,
                                 std::shared_ptr<gko::Executor> exec);

// Deallocate application data.
void appl_clean_up(ApplSysData& appl_data, std::shared_ptr<gko::Executor> exec);

/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int               fails = 0;        /* counter for test failures  */
  sunindextype      matrows, matcols; /* matrix dimensions          */
  sunindextype      num_blocks;       /* number of matrix blocks    */
  N_Vector          x, y;             /* test vectors               */
  realtype          *xdata;           /* pointers to vector data    */
  int               square;
  sundials::Context sunctx;

  auto gko_exec = OMP_OR_HIP_OR_CUDA(gko::OmpExecutor::create(),
    gko::HipExecutor::create(0, gko::OmpExecutor::create(), true),
    gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true));

  /* check input and set vector length */
  if (argc < 3) {
    printf("ERROR: TWO (2) Input required: matrix rows, number of blocks\n");
    return(-1);
  }

  matrows = (sunindextype) atol(argv[1]);
  if (matrows <= 0) {
    printf("ERROR: number of rows must be a positive integer \n");
    return(-1);
  }
  matcols = matrows;

  num_blocks = (sunindextype) atol(argv[2]);
  if (num_blocks < 1) {
    printf("ERROR: number of blocks must >= 1 \n");
    return(-1);
  }

  SetTiming(0);

  square = (matrows == matcols) ? 1 : 0;
  printf("\n SUNMATRIX_GINKGOBLOCKCSR test: size %ld by %ld, num_blocks %ld\n\n",
         (long int) matrows, (long int) matcols, (long int) num_blocks);

  // Create vectors and matrices
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, matrows);

  x = OMP_OR_HIP_OR_CUDA( N_VNew_Serial(num_blocks*matcols, sunctx),
                          N_VNew_Hip(num_blocks*matcols, sunctx),
                          N_VNew_Cuda(num_blocks*matcols, sunctx) );
  y = OMP_OR_HIP_OR_CUDA( N_VNew_Serial(num_blocks*matrows, sunctx),
                          N_VNew_Hip(num_blocks*matrows, sunctx),
                          N_VNew_Cuda(num_blocks*matrows, sunctx) );

  // The "application" generates the batch of linear systems on the device
  auto appl_sys = appl_generate_system(matrows, num_blocks, gko_exec);

  // Create batch_dim object to describe the dimensions of the batch matrix.
  auto batch_mat_size =
      gko::batch_dim<>(num_blocks, gko::dim<2>(matrows, matcols));
  auto batch_vec_size =
      gko::batch_dim<>(num_blocks, gko::dim<2>(matrows, 1));

  // Fill matrix
  auto values_view = gko::Array<sunrealtype>::view(
        gko_exec, num_blocks * appl_sys.nnz, appl_sys.all_values);
  auto rowptrs_view = gko::Array<sunindextype>::view(gko_exec, matrows + 1,
                                                          appl_sys.row_ptrs);
  auto colidxs_view = gko::Array<sunindextype>::view(gko_exec, appl_sys.nnz,
                                                          appl_sys.col_idxs);
  auto gko_batch_matrix = gko::share(GkoBatchMatrixType::create(
      gko_exec, batch_mat_size, std::move(values_view), std::move(colidxs_view),
      std::move(rowptrs_view)));

  // Fill identity matrix
  auto gko_ident = gko::share(GkoMatrixType::create(gko_exec, batch_mat_size.at(0)));
  gko_ident->read(gko::matrix_data<sunrealtype>::diag(batch_mat_size.at(0), 1.0));
  auto gko_batch_ident = gko::share(GkoBatchMatrixType::create(gko_exec, num_blocks, gko_ident.get()));

  // Wrap ginkgo matrices for SUNDIALS.
  // sundials::ginkgo::Matrix is overloaded to a SUNMatrix.
  SUNMatrixType A{gko_batch_matrix, sunctx};
  SUNMatrixType I{gko_batch_ident, sunctx};

  xdata = N_VGetArrayPointer(x);
  for(sunindextype i = 0; i < num_blocks*matcols; i++) {
    xdata[i] = distribution(generator);
  }
  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(x),
                       N_VCopyToDevice_Cuda(x) );

  // Compute true solution
  SUNMatrix Aref = SUNDenseMatrix(num_blocks*matrows, num_blocks*matcols, sunctx);
  sunindextype iblk = 0, nnz_count = 0, nnz_offset = 0;
  for (auto& block : A.gkomtx()->unbatch()) {
    for (sunindextype irow = 0; irow < block->get_size()[0]; irow++) {
      for (sunindextype inz = appl_sys.row_ptrs[irow]; inz < appl_sys.row_ptrs[irow + 1]; inz++) {
        SM_ELEMENT_D(Aref, iblk*matrows+irow, iblk*matrows+appl_sys.col_idxs[inz]) = appl_sys.all_values[nnz_offset+inz];
        nnz_count++;
      }
    }
    nnz_offset = nnz_count;
    iblk++;
  }
  SUNMatMatvec_Dense(Aref, x, y);
  SUNMatDestroy(Aref);

  // Do the SUNMatrix tests
  fails += Test_SUNMatGetID(A, SUNMATRIX_GINKGOBLOCKCSR, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  // fails += Test_SUNMatZero(A, 0);
  if (square) {
    // fails += Test_SUNMatScaleAdd(A, I, 0);
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
  appl_clean_up(appl_sys, gko_exec);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  auto Amat = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  auto Bmat = static_cast<SUNMatrixType*>(B->content)->gkomtx();
  const auto Avalues = Amat->get_const_values();
  const auto Bvalues = Bmat->get_const_values();

  /* check lengths */
  if (Amat->get_size() != Bmat->get_size()) {
    printf(">>> [ERROR] check_matrix: different matrix sizes\n");
    return(1);
  }

  /* compare data */
  for (sunindextype ivalue = 0; ivalue < Amat->get_num_stored_elements() ; ivalue++) {
    failure += SUNRCompareTol(Avalues[ivalue], Bvalues[ivalue], tol);
  }

  return failure > 0;
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  auto Amat = static_cast<SUNMatrixType*>(A->content)->gkomtx();
  const auto Avalues = Amat->get_const_values();

  /* compare data */
  for (sunindextype ivalue = 0; ivalue < Amat->get_num_stored_elements(); ivalue++) {
    int check = SUNRCompareTol(Avalues[ivalue], val, tol);
    if (check) {
      printf("  actual = %g != %e (err = %g)\n", Avalues[ivalue], val,  SUNRabs(Avalues[ivalue]-val));
      failure += check;
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
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(computed),
                       N_VCopyFromDevice_Cuda(computed) );
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(expected),
                       N_VCopyFromDevice_Cuda(expected) );

  /* get vector data */
  xdata = N_VGetArrayPointer(computed);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(computed);
  yldata = N_VGetLength(expected);

  if (xldata != yldata) {
    printf(">>> [ERROR] check_vector: different vector lengths \n");
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("check_vector failures:\n");
    for(i=0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  computed[%ld] = %g != %e (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
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

void sync_device(SUNMatrix A)
{
  OMP_OR_HIP_OR_CUDA(, hipDeviceSynchronize(),
                       cudaDeviceSynchronize() );
}

ApplSysData appl_generate_system(const int nrows, const sunindextype nsystems,
                                 std::shared_ptr<gko::Executor> exec)
{
  const int nnz = nrows * 3 - 2;
  std::ranlux48 rgen(15);
  std::normal_distribution<sunrealtype> distb(0.5, 0.1);
  std::vector<sunrealtype> spacings(nsystems * nrows);
  std::generate(spacings.begin(), spacings.end(),
                [&]() { return distb(rgen); });

  std::vector<sunrealtype> allvalues(nnz * nsystems);
  for (sunindextype isys = 0; isys < nsystems; isys++) {
      allvalues[isys * nnz] = 2.0 / spacings[isys * nrows];
      allvalues[isys * nnz + 1] = -1.0;
      for (int irow = 0; irow < nrows - 2; irow++) {
          allvalues[isys * nnz + 2 + irow * 3] = -1.0;
          allvalues[isys * nnz + 2 + irow * 3 + 1] =
              2.0 / spacings[isys * nrows + irow + 1];
          allvalues[isys * nnz + 2 + irow * 3 + 2] = -1.0;
      }
      allvalues[isys * nnz + 2 + (nrows - 2) * 3] = -1.0;
      allvalues[isys * nnz + 2 + (nrows - 2) * 3 + 1] =
          2.0 / spacings[(isys + 1) * nrows - 1];
      assert(isys * nnz + 2 + (nrows - 2) * 3 + 2 == (isys + 1) * nnz);
  }

  std::vector<sunindextype> rowptrs(nrows + 1);
  rowptrs[0] = 0;
  rowptrs[1] = 2;
  for (int i = 2; i < nrows; i++) {
      rowptrs[i] = rowptrs[i - 1] + 3;
  }
  rowptrs[nrows] = rowptrs[nrows - 1] + 2;
  assert(rowptrs[nrows] == nnz);

  std::vector<sunindextype> colidxs(nnz);
  colidxs[0] = 0;
  colidxs[1] = 1;
  const int nnz_per_row = 3;
  for (int irow = 1; irow < nrows - 1; irow++) {
      colidxs[2 + (irow - 1) * nnz_per_row] = irow - 1;
      colidxs[2 + (irow - 1) * nnz_per_row + 1] = irow;
      colidxs[2 + (irow - 1) * nnz_per_row + 2] = irow + 1;
  }
  colidxs[2 + (nrows - 2) * nnz_per_row] = nrows - 2;
  colidxs[2 + (nrows - 2) * nnz_per_row + 1] = nrows - 1;
  assert(2 + (nrows - 2) * nnz_per_row + 1 == nnz - 1);

  sunindextype* const row_ptrs = exec->alloc<sunindextype>(nrows + 1);
  exec->copy_from(exec->get_master().get(), static_cast<sunindextype>(nrows + 1),
                  rowptrs.data(), row_ptrs);
  sunindextype* const col_idxs = exec->alloc<sunindextype>(nnz);
  exec->copy_from(exec->get_master().get(), static_cast<sunindextype>(nnz),
                  colidxs.data(), col_idxs);
  sunrealtype* const all_values = exec->alloc<sunrealtype>(nsystems * nnz);
  exec->copy_from(exec->get_master().get(), nsystems * nnz, allvalues.data(),
                  all_values);
  return {nsystems, nrows, nnz, row_ptrs, col_idxs, all_values};
}

void appl_clean_up(ApplSysData& appl_data, std::shared_ptr<gko::Executor> exec)
{
  // In general, the application would control non-const pointers;
  //  the const casts below would not be needed.
  exec->free(const_cast<sunindextype*>(appl_data.row_ptrs));
  exec->free(const_cast<sunindextype*>(appl_data.col_idxs));
  exec->free(const_cast<sunrealtype*>(appl_data.all_values));
}