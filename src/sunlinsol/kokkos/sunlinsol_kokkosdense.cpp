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
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>

#include <nvector/nvector_kokkos.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>
#include <sunlinsol/sunlinsol_kokkosdense.hpp>

#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>

// Constants
#define ONE SUN_RCONST(1.0)

// Accessor macros
#define VEC_CONTENT(x) ( static_cast<N_VectorContent_Kokkos>(x->content) )
#define MAT_CONTENT(A) ( static_cast<SUNMatrixContent_KokkosDense>(A->content) )
#define LS_CONTENT(S)  ( static_cast<SUNLinearSolverContent_KokkosDense>(S->content) )

/* -----------------------------------------------------------------------------
 * Implementation specific routines
 * ---------------------------------------------------------------------------*/

/* -------------
 * Constructors
 * -------------*/

SUNLinearSolver SUNLinSol_KokkosDense(N_Vector y, SUNMatrix Amat,
                                      SUNContext sunctx)
{
  // Check inputs
  if (!y || !Amat)
    return nullptr;

  if (!(y->ops) || !(Amat->ops))
    return nullptr;

  if (!(y->ops->nvgetlength) || !(y->ops->nvgetdevicearraypointer) ||
      !(Amat->ops->getid))
    return nullptr;

  // Check compatibility with supplied SUNMatrix
  if (SUNMatGetID(Amat) != SUNMATRIX_KOKKOSDENSE)
    return nullptr;

  if (!(Amat->content))
    return nullptr;

  auto A = (SUNMatrixContent_KokkosDense) Amat->content;

  // Check that the matrix is square
  if (A->M != A->N)
    return nullptr;

  // Check that the matirx and vector dimensions agree
  if (A->M * A->nblocks != N_VGetLength(y))
    return nullptr;

  // Create the linear solver
  SUNLinearSolver S = SUNLinSolNewEmpty(sunctx);
  if (!S) return nullptr;

  // Attach operations
  S->ops->gettype    = SUNLinSolGetType_KokkosDense;
  S->ops->getid      = SUNLinSolGetID_KokkosDense;
  S->ops->initialize = SUNLinSolInitialize_KokkosDense;
  S->ops->setup      = SUNLinSolSetup_KokkosDense;
  S->ops->solve      = SUNLinSolSolve_KokkosDense;
  S->ops->lastflag   = SUNLinSolLastFlag_KokkosDense;
  S->ops->free       = SUNLinSolFree_KokkosDense;

  // Create content
  SUNLinearSolverContent_KokkosDense content = (SUNLinearSolverContent_KokkosDense) malloc(sizeof(*content));
  if (!content) { SUNLinSolFree(S); return nullptr; }

  // Fill content
  content->last_flag = 0;

  // Attach content
  S->content = content;

  return S;
}

/* -----------------------------------------------------------------------------
 * Implementation of SUNLinearSolver operations
 * ---------------------------------------------------------------------------*/

SUNLinearSolver_Type SUNLinSolGetType_KokkosDense(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID SUNLinSolGetID_KokkosDense(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_KOKKOSDENSE;
}

int SUNLinSolInitialize_KokkosDense(SUNLinearSolver S)
{
  // All solver-specific memory has already been allocated
  if (!S) return SUNLS_MEM_NULL;
  LS_CONTENT(S)->last_flag = SUNLS_SUCCESS;
  return SUNLS_SUCCESS;
}

int SUNLinSolSetup_KokkosDense(SUNLinearSolver S, SUNMatrix A)
{
  // Check for valid inputs
  if (!S)
    return SUNLS_MEM_NULL;

  if (!A)
  {
    LS_CONTENT(S)->last_flag = SUNLS_MEM_NULL;
    return SUNLS_MEM_NULL;
  }

  // Ensure that A is a magma dense matrix
  if (SUNMatGetID(A) != SUNMATRIX_KOKKOSDENSE)
  {
    LS_CONTENT(S)->last_flag = SUNLS_ILL_INPUT;
    return SUNLS_ILL_INPUT;
  }

  // Access matrix data
  auto A_data = MAT_CONTENT(A)->data_view;
  const auto blocks = MAT_CONTENT(A)->nblocks;

  // Do LU factorization (no pivoting) of A
  using team_policy = Kokkos::TeamPolicy<ExecSpace>;
  using member_type = Kokkos::TeamPolicy<ExecSpace>::member_type;

  Kokkos::parallel_for("sunlinsol_lu", team_policy(blocks, Kokkos::AUTO, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const member_type &team_member)
                       {
                         const int idx = team_member.league_rank();
                         auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
                         KokkosBatched::TeamLU<member_type, KokkosBatched::Algo::LU::Unblocked>
                           ::invoke(team_member, A_subdata);
                       });

  LS_CONTENT(S)->last_flag = SUNLS_SUCCESS;
  return SUNLS_SUCCESS;
}

int SUNLinSolSolve_KokkosDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol)
{
  // Check for valid inputs
  if (!S)
    return SUNLS_MEM_NULL;

  if (!A || !x || !b)
  {
    LS_CONTENT(S)->last_flag = SUNLS_MEM_NULL;
    return SUNLS_MEM_NULL;
  }

  // Ensure that A is a magma dense matrix
  if (SUNMatGetID(A) != SUNMATRIX_KOKKOSDENSE)
  {
    LS_CONTENT(S)->last_flag = SUNLS_ILL_INPUT;
    return SUNLS_ILL_INPUT;
  }

  // Copy b into x
  N_VScale(ONE, b, x);

  // Access matrix and vector data
  auto A_data = MAT_CONTENT(A)->data_view;
  auto x_data = *(VEC_CONTENT(x)->device_data);

  const auto blocks = MAT_CONTENT(A)->nblocks;
  const auto rows   = MAT_CONTENT(A)->M;

  // Solve the linear system
  using team_policy = Kokkos::TeamPolicy<ExecSpace>;
  using member_type = Kokkos::TeamPolicy<ExecSpace>::member_type;

  Kokkos::parallel_for("sunlinsol_trsv", team_policy(blocks, Kokkos::AUTO, Kokkos::AUTO),
                       KOKKOS_LAMBDA(const member_type &team_member)
                       {
                         const int idx = team_member.league_rank();
                         auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
                         auto x_subdata = Kokkos::subview(x_data, Kokkos::pair<int, int>(idx * rows, (idx+1) * rows));
                         KokkosBatched::TeamVectorTrsv<member_type,
                                                       KokkosBatched::Uplo::Lower,
                                                       KokkosBatched::Trans::NoTranspose,
                                                       KokkosBatched::Diag::Unit,
                                                       KokkosBatched::Algo::Trsv::Unblocked>
                           ::invoke(team_member, ONE, A_subdata, x_subdata);
                         KokkosBatched::TeamVectorTrsv<member_type,
                                                       KokkosBatched::Uplo::Upper,
                                                       KokkosBatched::Trans::NoTranspose,
                                                       KokkosBatched::Diag::NonUnit,
                                                       KokkosBatched::Algo::Trsv::Unblocked>
                           ::invoke(team_member, ONE, A_subdata, x_subdata);
                       });

  LS_CONTENT(S)->last_flag = SUNLS_SUCCESS;
  return SUNLS_SUCCESS;
}

sunindextype SUNLinSolLastFlag_KokkosDense(SUNLinearSolver S)
{
  // return the stored 'last_flag' value
  if (!S) return -1;
  return LS_CONTENT(S)->last_flag;
}

int SUNLinSolFree_KokkosDense(SUNLinearSolver S)
{
  // return if S is already free
  if (!S)
    return SUNLS_SUCCESS;

  // delete items from contents, then delete generic structure
  if (S->content)
  {
    free(S->content);
    S->content = nullptr;
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = nullptr;
  }
  free(S);
  S = nullptr;
  return SUNLS_SUCCESS;
}
