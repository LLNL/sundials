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
 * This is the header file for a dense SUNMarix implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNMATRIX_KOKKOSDENSE_HPP
#define _SUNMATRIX_KOKKOSDENSE_HPP

#include <Kokkos_Core.hpp>
//#include <nvector/nvector_kokkos.h>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_matrix.hpp>

namespace sundials {
namespace kokkos {

// Forward decalaration of regular Matrix class
template<class ExecSpace = Kokkos::DefaultExecutionSpace,
         class MemSpace = class ExecSpace::memory_space>
class DenseMatrix;

// Get the Kokkos dense matrix wrapped by a SUNMatrix
template<class ExecSpace = Kokkos::DefaultExecutionSpace,
         class MemSpace = class ExecSpace::memory_space>
inline DenseMatrix<ExecSpace, MemSpace>* GetDenseMat(SUNMatrix A)
{
  return static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content);
}

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

SUNMatrix_ID SUNMatGetID_KokkosDense(SUNMatrix A)
{
  return SUNMATRIX_KOKKOSDENSE;
}

template<class ExecSpace, class MemSpace>
SUNMatrix SUNMatClone_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};
  auto new_mat{new DenseMatrix<ExecSpace, MemSpace>(A_mat->blocks(),
                                                    A_mat->block_rows(),
                                                    A_mat->block_cols(),
                                                    A_mat->exec_space(),
                                                    A_mat->sunctx())};
  return new_mat->Convert();
}

template<class ExecSpace, class MemSpace>
void SUNMatDestroy_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};
  delete A_mat; // NOLINT
  return;
}

template<class ExecSpace, class MemSpace>
int SUNMatZero_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();

  // Zero out matrix
  Kokkos::parallel_for("sunmat_zero",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                       (A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         A_data(i, j, k) = 0.0;
                       });

  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatCopy_KokkosDense(SUNMatrix A, SUNMatrix B)
{
  auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};
  auto B_mat{GetDenseMat<ExecSpace, MemSpace>(B)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();
  auto B_data = B_mat->view();

  // Copy A into B
  Kokkos::parallel_for("sunmat_copy",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                       (A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         B_data(i, j, k) = A_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatScaleAdd_KokkosDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};
  auto B_mat{GetDenseMat<ExecSpace, MemSpace>(B)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();
  auto B_data = B_mat->view();

  // Scale A by c and add B
  Kokkos::parallel_for("sunmat_scale_add",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                       (A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         A_data(i, j, k) = c * A_data(i, j, k) + B_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatScaleAddI_KokkosDense(sunrealtype c, SUNMatrix A)
{
  auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();

  // Scale A by c and add I
  Kokkos::parallel_for("sunmat_scale_add_i",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>
                       (A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         if (j == k) A_data(i, j, k) = c * A_data(i, j, k) + 1.0;
                         else A_data(i, j, k) = c * A_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatMatvec_KokkosDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  // auto A_mat{GetDenseMat<ExecSpace, MemSpace>(A)};
  // auto x_vec{GetVec<ExecSpace, MemSpace>(x)};
  // auto y_vec{GetVec<ExecSpace, MemSpace>(y)};

  // const auto blocks = A_mat->blocks();
  // const auto rows   = A_mat->block_rows();
  // const auto cols   = A_mat->block_cols();

  // auto A_exec = A_mat->exec_space();
  // auto A_data = A_mat->view();
  // auto x_data = x_vec->view();
  // auto y_data = y_vec->view();

  // // Use batched or single gemv to do y = alpha * A * x + beta * y
  // if (blocks > 1)
  // {
  //   using team_policy = Kokkos::TeamPolicy<ExecSpace>;
  //   using member_type = Kokkos::TeamPolicy<ExecSpace>::member_type;

  //   Kokkos::parallel_for("sunmatvec_batch",
  //                        team_policy(A_exec, blocks, Kokkos::AUTO, Kokkos::AUTO),
  //                        KOKKOS_LAMBDA(const member_type &team_member)
  //                        {
  //                          const int idx = team_member.league_rank();
  //                          auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
  //                          auto x_subdata = Kokkos::subview(x_data, Kokkos::pair<int, int>(idx * cols, (idx+1) * cols));
  //                          auto y_subdata = Kokkos::subview(y_data, Kokkos::pair<int, int>(idx * rows, (idx+1) * rows));
  //                          KokkosBatched::TeamVectorGemv<member_type, KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Gemv::Unblocked>
  //                            ::invoke(team_member, 1.0, A_subdata, x_subdata, 0.0, y_subdata);
  //                        });
  // }
  // else
  // {
  //   auto A_subdata = Kokkos::subview(A_data, 0, Kokkos::ALL(), Kokkos::ALL());
  //   KokkosBlas::gemv("N", 1.0, A_subdata, x_data, 0.0, y_data);
  // }

  return SUNMAT_SUCCESS;
}

} // namespace impl

// =============================================================================
// Public namespace
// =============================================================================

// -----------------------------------------------------------------------------
// Kokkos dense matrix class, convertible to a SUNMatrix
// -----------------------------------------------------------------------------

template<class ExecSpace, class MemSpace>
class DenseMatrix : public sundials::impl::BaseMatrix,
                    public sundials::ConvertibleTo<SUNMatrix>
{
public:
  // Default constructor - means the matrix must be copied or moved to
  DenseMatrix() = default;

  // Single matrix constructors
  DenseMatrix(sunindextype rows, sunindextype cols, SUNContext sunctx)
    : DenseMatrix(1, rows, cols, ExecSpace(), sunctx)
  { }

  DenseMatrix(sunindextype rows, sunindextype cols, ExecSpace exec_space,
              SUNContext sunctx)
    : DenseMatrix(1, rows, cols, exec_space, sunctx)
  { }

  // Block-diagonal matrix constructors
  DenseMatrix(sunindextype blocks, sunindextype block_rows, sunindextype block_cols,
              SUNContext sunctx)
    : DenseMatrix(blocks, block_rows, block_cols, ExecSpace(), sunctx)
  { }

  // Block-diagonal matrix with user-supplied execution space instance
  DenseMatrix(sunindextype blocks, sunindextype block_rows,
              sunindextype block_cols, ExecSpace exec_space, SUNContext sunctx)
    : sundials::impl::BaseMatrix(sunctx), exec_space_(exec_space)
  {
    initSUNMatrix();
    view_ = Kokkos::View<sunrealtype***, MemSpace>("sunmat_view", blocks,
                                                   block_rows, block_cols);
  }

  // Move constructor
  DenseMatrix(DenseMatrix&& that_matrix) noexcept
    : sundials::impl::BaseMatrix(std::forward<DenseMatrix>(that_matrix)),
      exec_space_(std::move(that_matrix.exec_space_)), view_(std::move(that_matrix.exec_space_))
  { }

  // Copy constructor
  DenseMatrix(const DenseMatrix& that_matrix)
    : sundials::impl::BaseMatrix(that_matrix),
      exec_space_(that_matrix.exec_space_)
  {
    view_ = Kokkos::View<sunrealtype***, MemSpace>("sunmat_view",
                                                   that_matrix.blocks(),
                                                   that_matrix.block_rows(),
                                                   that_matrix.block_cols());
    deep_copy(view_, that_matrix.view_);
  }

  // Move assignment
  DenseMatrix& operator=(DenseMatrix&& rhs) noexcept
  {
    exec_space_ = std::move(rhs.exec_space_);
    view_       = std::move(rhs.view_);

    sundials::impl::BaseMatrix::operator=(std::forward<DenseMatrix>(rhs));

    return *this;
  }

  // Copy assignment
  DenseMatrix& operator=(const DenseMatrix& rhs)
  {
    exec_space_ = rhs.exec_space_;
    view_       = rhs.view_;

    sundials::impl::BaseMatrix::operator=(rhs);

    return *this;
  }

  // Default destructor since all members are RAII
  virtual ~DenseMatrix() = default;

  // Get the Kokkos execution space
  ExecSpace exec_space()
  {
    return exec_space_;
  }

  // Get the Kokkos view
  Kokkos::View<sunrealtype***, MemSpace> view()
  {
    return view_;
  }

  // Get the number of blocks
  sunindextype blocks() const
  {
    return view_.extent(0);
  }

  // Get the number of rows in a block
  sunindextype block_rows() const
  {
    return view_.extent(1);
  }

  // Get the number of columns in a block
  sunindextype block_cols() const
  {
    return view_.extent(2);
  }

  // Get the number of rows
  sunindextype rows() const
  {
    return view_.extent(0) * view_.extent(1);
  }

  // Get the number of columns
  sunindextype cols() const
  {
    return view_.extent(0) * view_.extent(2);
  }

  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods

  // Implicit conversion to a SUNMatrix
  operator SUNMatrix() override
  {
    return object_.get();
  }

  // Implicit conversion to SUNMatrix
  operator SUNMatrix() const override
  {
    return object_.get();
  }

  // Explicit conversion to a SUNMatrix
  SUNMatrix Convert() override
  {
    return object_.get();
  }

  // Explicit conversion to a SUNMatrix
  SUNMatrix Convert() const override
  {
    return object_.get();
  }

private:
  // Kokkos execution space
  ExecSpace exec_space_;

  // Matrix data view [blocks, rows, cols]
  Kokkos::View<sunrealtype***, MemSpace> view_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid     = impl::SUNMatGetID_KokkosDense;
    this->object_->ops->clone     = impl::SUNMatClone_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->destroy   = impl::SUNMatDestroy_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->zero      = impl::SUNMatZero_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->copy      = impl::SUNMatCopy_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->scaleadd  = impl::SUNMatScaleAdd_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->scaleaddi = impl::SUNMatScaleAddI_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->matvec    = impl::SUNMatMatvec_KokkosDense<ExecSpace, MemSpace>;
  }
};

} // namespace kokkos
} // namespace sundials

#endif
