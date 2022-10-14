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
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_matrix.hpp>

#include <KokkosBlas2_gemv.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>

namespace sundials {
namespace kokkos {

// Forward decalaration of Matrix class
template<class ExecSpace = Kokkos::DefaultExecutionSpace,
         class MemSpace = class ExecSpace::memory_space>
class DenseMatrix;

// Get the Kokkos dense matrix wrapped by a SUNMatrix
template<class MatrixType>
inline MatrixType* GetDenseMat(SUNMatrix A)
{
  return static_cast<MatrixType*>(A->content);
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

template<class MatrixType>
SUNMatrix SUNMatClone_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto new_mat{new MatrixType(A_mat->blocks(),
                              A_mat->block_rows(),
                              A_mat->block_cols(),
                              A_mat->exec_space(),
                              A_mat->sunctx())};
  return new_mat->Convert();
}

template<class MatrixType>
void SUNMatDestroy_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  delete A_mat; // NOLINT
  return;
}

template<class MatrixType>
int SUNMatZero_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Zero out matrix
  Kokkos::parallel_for("sunmat_zero",
                       range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k)
                       {
                         A_data(i, j, k) = 0.0;
                       });

  return SUNMAT_SUCCESS;
}

template<class MatrixType>
int SUNMatCopy_KokkosDense(SUNMatrix A, SUNMatrix B)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto B_mat{GetDenseMat<MatrixType>(B)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();
  auto B_data = B_mat->view();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Copy A into B
  Kokkos::parallel_for("sunmat_copy",
                       range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k)
                       {
                         B_data(i, j, k) = A_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

template<class MatrixType>
int SUNMatScaleAdd_KokkosDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto B_mat{GetDenseMat<MatrixType>(B)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();
  auto B_data = B_mat->view();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Scale A by c and add B
  Kokkos::parallel_for("sunmat_scale_add",
                       range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k)
                       {
                         A_data(i, j, k) = c * A_data(i, j, k) + B_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

template<class MatrixType>
int SUNMatScaleAddI_KokkosDense(sunrealtype c, SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Scale A by c and add I
  Kokkos::parallel_for("sunmat_scale_add_i",
                       range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k)
                       {
                         if (j == k) A_data(i, j, k) = c * A_data(i, j, k) + 1.0;
                         else A_data(i, j, k) = c * A_data(i, j, k);
                       });

  return SUNMAT_SUCCESS;
}

template<class VectorType, class MatrixType>
int SUNMatMatvec_KokkosDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto x_vec{GetVec<VectorType>(x)};
  auto y_vec{GetVec<VectorType>(y)};

  const auto blocks = A_mat->blocks();
  const auto rows   = A_mat->block_rows();
  const auto cols   = A_mat->block_cols();

  auto A_exec = A_mat->exec_space();
  auto A_data = A_mat->view();
  auto x_data = x_vec->View();
  auto y_data = y_vec->View();

  // Use batched or single gemv to do y = alpha * A * x + beta * y
  if (blocks > 1)
  {
    using team_policy = typename MatrixType::team_policy;
    using member_type = typename MatrixType::member_type;
    using size_type   = typename MatrixType::size_type;

    Kokkos::parallel_for
      ("sunmatvec_batch",
       team_policy(A_exec, static_cast<int>(blocks), Kokkos::AUTO, Kokkos::AUTO),
       KOKKOS_LAMBDA(const member_type &team_member)
       {
         const int idx = team_member.league_rank();
         auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
         auto x_subdata = Kokkos::subview(x_data, Kokkos::pair<size_type, size_type>(idx * cols, (idx+1) * cols));
         auto y_subdata = Kokkos::subview(y_data, Kokkos::pair<size_type, size_type>(idx * rows, (idx+1) * rows));
         KokkosBatched::TeamVectorGemv<member_type, KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Gemv::Unblocked>
           ::invoke(team_member, 1.0, A_subdata, x_subdata, 0.0, y_subdata);
       });
  }
  else
  {
    auto A_subdata = Kokkos::subview(A_data, 0, Kokkos::ALL(), Kokkos::ALL());
    KokkosBlas::gemv("N", 1.0, A_subdata, x_data, 0.0, y_data);
  }

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
  using view_type    = Kokkos::View<sunrealtype***, MemSpace>;
  using size_type    = typename view_type::size_type;
  using range_policy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>;
  using team_policy  = typename Kokkos::TeamPolicy<ExecSpace>;
  using member_type  = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  // Default constructor - means the matrix must be copied or moved to
  DenseMatrix() = default;

  // Single matrix constructors
  DenseMatrix(size_type rows, size_type cols, SUNContext sunctx)
    : DenseMatrix(1, rows, cols, ExecSpace(), sunctx)
  { }

  DenseMatrix(size_type rows, size_type cols, ExecSpace exec_space,
              SUNContext sunctx)
    : DenseMatrix(1, rows, cols, exec_space, sunctx)
  { }

  // Block-diagonal matrix constructors
  DenseMatrix(size_type blocks, size_type block_rows, size_type block_cols,
              SUNContext sunctx)
    : DenseMatrix(blocks, block_rows, block_cols, ExecSpace(), sunctx)
  { }

  // Block-diagonal matrix with user-supplied execution space instance
  DenseMatrix(size_type blocks, size_type block_rows,
              size_type block_cols, ExecSpace exec_space, SUNContext sunctx)
    : sundials::impl::BaseMatrix(sunctx), exec_space_(exec_space)
  {
    initSUNMatrix();
    view_ = view_type("sunmat_view", blocks, block_rows, block_cols);
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
    view_ = view_type("sunmat_view", that_matrix.blocks(),
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
  view_type view()
  {
    return view_;
  }

  // Get the number of blocks
  size_type blocks() const
  {
    return static_cast<size_type>(view_.extent(0));
  }

  // Get the number of rows in a block
  size_type block_rows() const
  {
    return static_cast<size_type>(view_.extent(1));
  }

  // Get the number of columns in a block
  size_type block_cols() const
  {
    return static_cast<size_type>(view_.extent(2));
  }

  // Get the number of rows
  size_type rows() const
  {
    return static_cast<size_type>(view_.extent(0) * view_.extent(1));
  }

  // Get the number of columns
  size_type cols() const
  {
    return static_cast<size_type>(view_.extent(0) * view_.extent(2));
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
  view_type view_;

  void initSUNMatrix()
  {
    using vec_type = Vector<ExecSpace, MemSpace>;
    using mat_type = DenseMatrix<ExecSpace, MemSpace>;

    this->object_->content = this;

    this->object_->ops->getid     = impl::SUNMatGetID_KokkosDense;
    this->object_->ops->clone     = impl::SUNMatClone_KokkosDense<mat_type>;
    this->object_->ops->destroy   = impl::SUNMatDestroy_KokkosDense<mat_type>;
    this->object_->ops->zero      = impl::SUNMatZero_KokkosDense<mat_type>;
    this->object_->ops->copy      = impl::SUNMatCopy_KokkosDense<mat_type>;
    this->object_->ops->scaleadd  = impl::SUNMatScaleAdd_KokkosDense<mat_type>;
    this->object_->ops->scaleaddi = impl::SUNMatScaleAddI_KokkosDense<mat_type>;
    this->object_->ops->matvec    = impl::SUNMatMatvec_KokkosDense<vec_type, mat_type>;
  }
};

} // namespace kokkos
} // namespace sundials

#endif
