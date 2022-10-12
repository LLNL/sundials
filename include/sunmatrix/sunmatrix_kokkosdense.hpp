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

// #define VEC_CONTENT(x) ( static_cast<N_VectorContent_Kokkos>(x->content) )
// #define VEC_VIEW(x)    ( *VEC_CONTENT(x)->device_data )

namespace sundials {
namespace kokkos {

// Forward decalaration of regular Matrix class
template<class ExecSpace = Kokkos::DefaultExecutionSpace,
         class MemSpace = class ExecSpace::memory_space>
class DenseMatrix;

// -----------------------------------------------------------------------------
// Functions that operate on a DenseMatrix
// -----------------------------------------------------------------------------

namespace impl {

template<class ExecSpace, class MemSpace>
void Zero(DenseMatrix<ExecSpace, MemSpace>& A)
{
  const auto blocks = A.blocks();
  const auto rows   = A.block_rows();
  const auto cols   = A.block_cols();

  auto A_data = A.view();

  // Zero out matrix
  Kokkos::parallel_for("sunmat_zero",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         A_data(i, j, k) = 0.0;
                       });
}

template<class ExecSpace, class MemSpace>
void Copy(DenseMatrix<ExecSpace, MemSpace>& A, DenseMatrix<ExecSpace, MemSpace>& B)
{
  const auto blocks = A.blocks();
  const auto rows   = A.block_rows();
  const auto cols   = A.block_cols();

  auto A_data = A.view();
  auto B_data = B.view();

  // Copy A into B
  Kokkos::parallel_for("sunmat_copy",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         B_data(i, j, k) = A_data(i, j, k);
                       });
}

template<class ExecSpace, class MemSpace>
void ScaleAdd(const sunrealtype c, DenseMatrix<ExecSpace, MemSpace>& A,
              DenseMatrix<ExecSpace, MemSpace>& B)
{
  const auto blocks = A.blocks();
  const auto rows   = A.block_rows();
  const auto cols   = A.block_cols();

  auto A_data = A.view();
  auto B_data = B.view();

  // Scale A by c and add B
  Kokkos::parallel_for("sunmat_scale_add",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         A_data(i, j, k) = c * A_data(i, j, k) + B_data(i, j, k);
                       });
}

template<class ExecSpace, class MemSpace>
void ScaleAddI(const sunrealtype c, DenseMatrix<ExecSpace, MemSpace>& A)
{
  const auto blocks = A.blocks();
  const auto rows   = A.block_rows();
  const auto cols   = A.block_cols();

  auto A_data = A.view();

  // Scale A by c and add I
  Kokkos::parallel_for("sunmat_scale_add_i",
                       Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0}, {blocks, rows, cols}),
                       KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k)
                       {
                         if (j == k) A_data(i, j, k) = c * A_data(i, j, k) + 1.0;
                         else A_data(i, j, k) = c * A_data(i, j, k);
                       });
}

template<class ExecSpace, class MemSpace>
void Matvec(DenseMatrix<ExecSpace, MemSpace>& A, N_Vector x, N_Vector y)
{
  // const auto blocks = A.blocks();
  // const auto rows   = A.block_rows();
  // const auto cols   = A.block_cols();

  // auto A_data = A.view();
  // // auto x_data = x.view();
  // // auto y_data = y.view();
  // auto x_data = VEC_VIEW(x);
  // auto y_data = VEC_VIEW(y);

  // // Use batched or single gemv to do y = alpha * A * x + beta * y
  // if (blocks > 1)
  // {
  //   using team_policy = Kokkos::TeamPolicy<ExecSpace>;
  //   using member_type = Kokkos::TeamPolicy<ExecSpace>::member_type;

  //   Kokkos::parallel_for("sunmatvec_batch", team_policy(blocks, Kokkos::AUTO, Kokkos::AUTO),
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
}

} // namespace impl

// -----------------------------------------------------------------------------
// Methods that operate on a SUNMatrix
// -----------------------------------------------------------------------------

template<class ExecSpace = Kokkos::DefaultExecutionSpace,
         class MemSpace = class ExecSpace::memory_space>
DenseMatrix<ExecSpace, MemSpace>* GetDenseMat(SUNMatrix A)
{
  return static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content);
}

template<class ExecSpace, class MemSpace>
SUNMatrix_ID SUNMatGetID_KokkosDense(SUNMatrix A)
{
  return SUNMATRIX_KOKKOSDENSE;
}

template<class ExecSpace, class MemSpace>
SUNMatrix SUNMatClone_KokkosDense(SUNMatrix A)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  auto new_mat{new DenseMatrix<ExecSpace, MemSpace>(*A_mat)}; // NOLINT
  return new_mat->Convert();
}

template<class ExecSpace, class MemSpace>
void SUNMatDestroy_KokkosDense(SUNMatrix A)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  delete A_mat; // NOLINT
  return;
}

template<class ExecSpace, class MemSpace>
int SUNMatZero_KokkosDense(SUNMatrix A)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  impl::Zero(*A_mat);
  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatCopy_KokkosDense(SUNMatrix A, SUNMatrix B)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  auto B_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(B->content)};
  impl::Copy(*A_mat, *B_mat);
  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatScaleAdd_KokkosDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  auto B_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(B->content)};
  impl::ScaleAdd(c, *A_mat, *B_mat);
  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatScaleAddI_KokkosDense(sunrealtype c, SUNMatrix A)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  impl::ScaleAddI(c, *A_mat);
  return SUNMAT_SUCCESS;
}

template<class ExecSpace, class MemSpace>
int SUNMatMatvec_KokkosDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto A_mat{static_cast<DenseMatrix<ExecSpace, MemSpace>*>(A->content)};
  impl::Matvec(*A_mat, x, y);
  return SUNMAT_SUCCESS;
}

// -----------------------------------------------------------------------------
// Class that wraps a matrix and allows it to convert to a SUNMatrix
// -----------------------------------------------------------------------------

template<class ExecSpace, class MemSpace>
class DenseMatrix : public sundials::impl::BaseMatrix,
                    public sundials::ConvertibleTo<SUNMatrix>
{
public:
  // Single matrix constructors
  DenseMatrix(sunindextype rows, sunindextype cols, SUNContext sunctx)
    : DenseMatrix(1, rows, cols, ExecSpace(), sunctx)
  { }

  DenseMatrix(sunindextype rows, sunindextype cols, ExecSpace ex,
              SUNContext sunctx)
    : DenseMatrix(1, rows, cols, ex, sunctx)
  { }

  // Block-diagonal matrix constructors
  DenseMatrix(sunindextype blocks, sunindextype rows, sunindextype cols,
              SUNContext sunctx)
    : DenseMatrix(blocks, rows, cols, ExecSpace(), sunctx)
  { }

  // Block-diagonal matrix with user-supplied execution space instance
  DenseMatrix(sunindextype blocks, sunindextype rows, sunindextype cols,
              ExecSpace ex, SUNContext sunctx)
    : sundials::impl::BaseMatrix(sunctx), ex_(ex)
  {
    initSUNMatrix();
    view_ = Kokkos::View<sunrealtype***, MemSpace>("sunmat_view", blocks, rows,
                                                   cols);
  }

  // Move constructor
  DenseMatrix(DenseMatrix&& that_matrix) noexcept
    : sundials::impl::BaseMatrix(std::forward<DenseMatrix>(that_matrix)),
      ex_(std::move(that_matrix.ex_)), view_(std::move(that_matrix.ex_))
  { }

  // Copy constructor (does not copy view data)
  DenseMatrix(const DenseMatrix& that_matrix)
    : sundials::impl::BaseMatrix(that_matrix),
      ex_(that_matrix.ex_)
  {
    view_ = Kokkos::View<sunrealtype***, MemSpace>("sunmat_view",
                                                   that_matrix.blocks(),
                                                   that_matrix.block_rows(),
                                                   that_matrix.block_cols());
  }

  // Move assignment
  DenseMatrix& operator=(DenseMatrix&& rhs) noexcept
  {
    ex_   = std::move(rhs.ex_);
    view_ = std::move(rhs.view_);

    sundials::impl::BaseMatrix::operator=(std::forward<DenseMatrix>(rhs));

    return *this;
  }

  // Copy assignment
  DenseMatrix& operator=(const DenseMatrix& rhs)
  {
    ex_   = rhs.ex_;
    view_ = rhs.view_;

    sundials::impl::BaseMatrix::operator=(rhs);

    return *this;
  }

  // Default destructor since all members are RAII
  virtual ~DenseMatrix() = default;

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
  // matrix data view [blocks, rows, cols]
  ExecSpace ex_;
  Kokkos::View<sunrealtype***, MemSpace> view_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid     = SUNMatGetID_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->clone     = SUNMatClone_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->destroy   = SUNMatDestroy_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->zero      = SUNMatZero_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->copy      = SUNMatCopy_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->scaleadd  = SUNMatScaleAdd_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->scaleaddi = SUNMatScaleAddI_KokkosDense<ExecSpace, MemSpace>;
    this->object_->ops->matvec    = SUNMatMatvec_KokkosDense<ExecSpace, MemSpace>;
  }
};

} // namespace kokkos
} // namespace sundials

#endif
