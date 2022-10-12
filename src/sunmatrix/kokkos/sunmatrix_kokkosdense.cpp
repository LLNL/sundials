/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for a dense SUNMatrix using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <nvector/nvector_kokkos.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>

#include <KokkosBlas2_gemv.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>

// Constants
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

// Accessor macros
#define VEC_CONTENT(x) ( static_cast<N_VectorContent_Kokkos>(x->content) )
#define MAT_CONTENT(A) ( static_cast<SUNMatrixContent_KokkosDense>(A->content) )

// Private function prototypes
static sunbooleantype compatible_matrix(SUNMatrix A, SUNMatrix B);
static sunbooleantype compatible_matrix(SUNMatrix A, N_Vector x, N_Vector y);

/* ----------------------------------------------------------------------------
 * Implementation specific functions
 * --------------------------------------------------------------------------*/

/* -------------
 * Constructors
 * -------------*/

SUNMatrix SUNMatrix_KokkosDense(sunindextype M, sunindextype N,
                                SUNContext sunctx)
{
  return SUNMatrix_KokkosDenseBlock(1, M, N, sunctx);
}

SUNMatrix SUNMatrix_KokkosDenseBlock(sunindextype nblocks, sunindextype M,
                                     sunindextype N, SUNContext sunctx)
{
  // Return with NULL matrix on illegal dimension input
  if ((M <= 0) || (N <= 0) || (nblocks <= 0))
    return nullptr;

  // Create an empty matrix object
  SUNMatrix A = SUNMatNewEmpty(sunctx);
  if (!A) return nullptr;

  // Attach operations
  A->ops->getid     = SUNMatGetID_KokkosDense;
  A->ops->clone     = SUNMatClone_KokkosDense;
  A->ops->destroy   = SUNMatDestroy_KokkosDense;
  A->ops->zero      = SUNMatZero_KokkosDense;
  A->ops->copy      = SUNMatCopy_KokkosDense;
  A->ops->scaleadd  = SUNMatScaleAdd_KokkosDense;
  A->ops->scaleaddi = SUNMatScaleAddI_KokkosDense;
  A->ops->matvec    = SUNMatMatvec_KokkosDense;

  // Create content
  SUNMatrixContent_KokkosDense content = (SUNMatrixContent_KokkosDense) malloc(sizeof(*content));
  if (!content) { SUNMatDestroy(A); return nullptr; }

  // Fill content
  content->M         = M;
  content->N         = N;
  content->nblocks   = nblocks;
  content->data_view = Kokkos::View<sunrealtype***>("sunmat_view", nblocks, M, N);

  // Attach content
  A->content = content;

  return A;
}

/* -------------------
 * Accessor functions
 * -------------------*/

sunindextype SUNMatrix_KokkosDense_Rows(SUNMatrix Amat)
{
  auto A = MAT_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_KOKKOSDENSE)
    return A->M * A->nblocks;
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNMatrix_KokkosDense_Columns(SUNMatrix Amat)
{
  auto A = MAT_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_KOKKOSDENSE)
    return A->N * A->nblocks;
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNMatrix_KokkosDense_BlockRows(SUNMatrix Amat)
{
  auto A = MAT_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_KOKKOSDENSE)
    return A->M;
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNMatrix_KokkosDense_BlockColumns(SUNMatrix Amat)
{
  auto A = MAT_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_KOKKOSDENSE)
    return A->N;
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNMatrix_KokkosDense_NumBlocks(SUNMatrix Amat)
{
  auto A = MAT_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_KOKKOSDENSE)
    return A->nblocks;
  else
    return SUNMAT_ILL_INPUT;
}

/* -----------------------------------------------------------------------------
 * Implementation of SUNMatrix operations
 * ---------------------------------------------------------------------------*/

SUNMatrix SUNMatClone_KokkosDense(SUNMatrix Amat)
{
  if (!Amat) return nullptr;

  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE) return nullptr;

  auto A = MAT_CONTENT(Amat);

  return SUNMatrix_KokkosDenseBlock(A->nblocks, A->M, A->N, Amat->sunctx);
}

void SUNMatDestroy_KokkosDense(SUNMatrix Amat)
{
  if (!Amat) return;

  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE) return;

  auto A = MAT_CONTENT(Amat);

  // Free content
  if (A)
  {
    free(A);
    Amat->content = nullptr;
  }

  // Free ops
  if (Amat->ops)
  {
    free(Amat->ops);
    Amat->ops = nullptr;
  }

  // Free matrix
  free(Amat);
  Amat = nullptr;

  return;
}

int SUNMatZero_KokkosDense(SUNMatrix Amat)
{
  if (!Amat)
    return SUNMAT_ILL_INPUT;

  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE)
    return SUNMAT_ILL_INPUT;

  const auto blocks = MAT_CONTENT(Amat)->nblocks;
  const auto rows   = MAT_CONTENT(Amat)->M;
  const auto cols   = MAT_CONTENT(Amat)->N;

  auto A_data = MAT_CONTENT(Amat)->data_view;

  // Zero out matrix
  Kokkos::parallel_for("sunmat_zero",
                       Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         A_data(i, j, k) = ZERO;
                       });

  return SUNMAT_SUCCESS;
}

int SUNMatCopy_KokkosDense(SUNMatrix Amat, SUNMatrix Bmat)
{
  if ((!Amat) || (!Bmat))
    return SUNMAT_ILL_INPUT;

  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE)
    return SUNMAT_ILL_INPUT;

  // Verify that A and B are compatible
  if (!compatible_matrix(Amat, Bmat))
    return SUNMAT_ILL_INPUT;

  const auto blocks = MAT_CONTENT(Amat)->nblocks;
  const auto rows   = MAT_CONTENT(Amat)->M;
  const auto cols   = MAT_CONTENT(Amat)->N;

  auto A_data = MAT_CONTENT(Amat)->data_view;
  auto B_data = MAT_CONTENT(Bmat)->data_view;

  // Copy A into B
  Kokkos::parallel_for("sunmat_copy",
                       Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         B_data(i, j, k) = A_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_KokkosDense(sunrealtype c, SUNMatrix Amat)
{
  if (!Amat)
    return SUNMAT_ILL_INPUT;

  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE)
    return SUNMAT_ILL_INPUT;

  const auto blocks = MAT_CONTENT(Amat)->nblocks;
  const auto rows   = MAT_CONTENT(Amat)->M;
  const auto cols   = MAT_CONTENT(Amat)->N;

  auto A_data = MAT_CONTENT(Amat)->data_view;

  // Scale A by c and add I
  Kokkos::parallel_for("sunmat_scale_add_i",
                       Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         if (j == k) A_data(i, j, k) = c * A_data(i, j, k) + ONE;
                         else A_data(i, j, k) = c * A_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_KokkosDense(sunrealtype c, SUNMatrix Amat, SUNMatrix Bmat)
{
  if ((!Amat) || (!Bmat))
    return SUNMAT_ILL_INPUT;

  if (!compatible_matrix(Amat, Bmat))
    return SUNMAT_ILL_INPUT;

  const auto blocks = MAT_CONTENT(Amat)->nblocks;
  const auto rows   = MAT_CONTENT(Amat)->M;
  const auto cols   = MAT_CONTENT(Amat)->N;

  auto A_data = MAT_CONTENT(Amat)->data_view;
  auto B_data = MAT_CONTENT(Bmat)->data_view;

  // Scale A by c and add B
  Kokkos::parallel_for("sunmat_scale_add",
                       Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         A_data(i, j, k) = c * A_data(i, j, k) + B_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

int SUNMatMatvec_KokkosDense(SUNMatrix Amat, N_Vector x, N_Vector y)
{
  if ((!Amat) || (!x) || (!y))
    return SUNMAT_ILL_INPUT;

  if (!compatible_matrix(Amat, x, y))
    return SUNMAT_ILL_INPUT;

  const auto blocks = MAT_CONTENT(Amat)->nblocks;
  const auto rows   = MAT_CONTENT(Amat)->M;
  const auto cols   = MAT_CONTENT(Amat)->N;

  auto A_data = MAT_CONTENT(Amat)->data_view;
  auto x_data = *(VEC_CONTENT(x)->device_data);
  auto y_data = *(VEC_CONTENT(y)->device_data);

  // Use batched or single gemv to do y = alpha * A * x + beta * y
  if (blocks > 1)
  {
    using team_policy = Kokkos::TeamPolicy<ExecSpace>;
    using member_type = Kokkos::TeamPolicy<ExecSpace>::member_type;

    Kokkos::parallel_for("sunmatvec_batch", team_policy(blocks, Kokkos::AUTO, Kokkos::AUTO),
                         KOKKOS_LAMBDA(const member_type &team_member)
                         {
                           const int idx = team_member.league_rank();
                           auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
                           auto x_subdata = Kokkos::subview(x_data, Kokkos::pair<int, int>(idx * cols, (idx+1) * cols));
                           auto y_subdata = Kokkos::subview(y_data, Kokkos::pair<int, int>(idx * rows, (idx+1) * rows));
                           KokkosBatched::TeamVectorGemv<member_type, KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Gemv::Unblocked>
                             ::invoke(team_member, ONE, A_subdata, x_subdata, ZERO, y_subdata);
                         });
  }
  else
  {
    auto A_subdata = Kokkos::subview(A_data, 0, Kokkos::ALL(), Kokkos::ALL());
    KokkosBlas::gemv("N", ONE, A_subdata, x_data, ZERO, y_data);
  }

  return SUNMAT_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

static sunbooleantype compatible_matrix(SUNMatrix Amat, SUNMatrix Bmat)
{
  // Both matrices must be a Kokkos dense matrix
  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE)
    return SUNFALSE;
  if (SUNMatGetID(Bmat) != SUNMATRIX_KOKKOSDENSE)
    return SUNFALSE;

  auto A = MAT_CONTENT(Amat);
  auto B = MAT_CONTENT(Bmat);

  // Both matrices must have the same shape
  if (A->M != B->M)
    return SUNFALSE;
  if (A->N != B->N)
    return SUNFALSE;
  if (A->nblocks != B->nblocks)
    return SUNFALSE;

  return SUNTRUE;
}

static sunbooleantype compatible_matrix(SUNMatrix Amat, N_Vector x, N_Vector y)
{
  // Matrix must be a Kokkos dense matrix
  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE)
    return SUNFALSE;

  //  Vectors must be Kokkos vector
  if (N_VGetVectorID(x) != SUNDIALS_NVEC_KOKKOS)
    return SUNFALSE;

  auto A = MAT_CONTENT(Amat);

  // Inner dimensions must agree
  if (A->N * A->nblocks != N_VGetLength(x))
    return SUNFALSE;

  // Outer dimensions must agree
  if (A->M * A->nblocks != N_VGetLength(y))
    return SUNFALSE;

  return SUNTRUE;
}
