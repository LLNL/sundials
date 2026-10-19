/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file for a dense SUNMarix implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_SUNMATRIX_KOKKOSDENSE_HPP
#define SUNDIALS_SUNMATRIX_KOKKOSDENSE_HPP

#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <Kokkos_Core.hpp>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_matrix.hpp>

namespace sundials {
namespace kokkos {

// Forward declaration of Matrix class
template<class ExecutionSpace = Kokkos::DefaultExecutionSpace,
         class MemorySpace    = typename ExecutionSpace::memory_space>
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

static SUNMatrix_ID SUNMatGetID_KokkosDense(SUNMatrix A)
{
  return SUNMATRIX_KOKKOSDENSE;
}

template<class MatrixType>
SUNMatrix SUNMatClone_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto new_mat{new MatrixType(*A_mat)};
  return new_mat->get();
}

template<class MatrixType>
void SUNMatDestroy_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  delete A_mat; // NOLINT
  return;
}

template<class MatrixType>
SUNErrCode SUNMatZero_KokkosDense(SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};

  const auto blocks = A_mat->Blocks();
  const auto rows   = A_mat->BlockRows();
  const auto cols   = A_mat->BlockCols();

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Zero out matrix
  Kokkos::parallel_for(
    "sunmat_zero", range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
    KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k) {
      A_data(i, j, k) = SUN_RCONST(0.0);
    });

  return SUN_SUCCESS;
}

template<class MatrixType>
SUNErrCode SUNMatCopy_KokkosDense(SUNMatrix A, SUNMatrix B)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto B_mat{GetDenseMat<MatrixType>(B)};

  const auto blocks = A_mat->Blocks();
  const auto rows   = A_mat->BlockRows();
  const auto cols   = A_mat->BlockCols();

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();
  auto B_data = B_mat->View();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Copy A into B
  Kokkos::parallel_for(
    "sunmat_copy", range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
    KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k) {
      B_data(i, j, k) = A_data(i, j, k);
    });

  return SUN_SUCCESS;
}

template<class MatrixType>
SUNErrCode SUNMatScaleAdd_KokkosDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto B_mat{GetDenseMat<MatrixType>(B)};

  const auto blocks = A_mat->Blocks();
  const auto rows   = A_mat->BlockRows();
  const auto cols   = A_mat->BlockCols();

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();
  auto B_data = B_mat->View();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Scale A by c and add B
  Kokkos::parallel_for(
    "sunmat_scale_add", range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
    KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k) {
      A_data(i, j, k) = c * A_data(i, j, k) + B_data(i, j, k);
    });

  return SUN_SUCCESS;
}

template<class MatrixType>
SUNErrCode SUNMatScaleAddI_KokkosDense(sunrealtype c, SUNMatrix A)
{
  auto A_mat{GetDenseMat<MatrixType>(A)};

  const auto blocks = A_mat->Blocks();
  const auto rows   = A_mat->BlockRows();
  const auto cols   = A_mat->BlockCols();

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();

  using range_policy = typename MatrixType::range_policy;
  using size_type    = typename MatrixType::size_type;

  // Scale A by c and add I
  Kokkos::parallel_for(
    "sunmat_scale_add_i", range_policy(A_exec, {0, 0, 0}, {blocks, rows, cols}),
    KOKKOS_LAMBDA(const size_type i, const size_type j, const size_type k) {
      if (j == k) A_data(i, j, k) = c * A_data(i, j, k) + SUN_RCONST(1.0);
      else A_data(i, j, k) = c * A_data(i, j, k);
    });

  return SUN_SUCCESS;
}

template<class VectorType, class MatrixType>
SUNErrCode SUNMatMatvec_KokkosDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  if (!x || !y || !VectorType::IsNVector(x) || !VectorType::IsNVector(y) ||
      N_VGetVectorID(x) != SUNDIALS_NVEC_KOKKOS ||
      N_VGetVectorID(y) != SUNDIALS_NVEC_KOKKOS)
    return SUN_ERR_ARG_WRONGTYPE;
  if (N_VGetLength(x) != static_cast<sunindextype>(GetDenseMat<MatrixType>(A)->BlockCols() * GetDenseMat<MatrixType>(A)->Blocks()) ||
      N_VGetLength(y) != static_cast<sunindextype>(GetDenseMat<MatrixType>(A)->BlockRows() * GetDenseMat<MatrixType>(A)->Blocks()))
    return SUN_ERR_ARG_DIMSMISMATCH;
  // Keep this a runtime error: the callback is part of the C API and must
  // reject an accidentally installed rank-N vector rather than making the
  // integration ill-formed when its concrete vector type changes.
  if constexpr (VectorType::view_type::rank != 1)
    return SUN_ERR_ARG_WRONGTYPE;
  static_assert(std::is_same_v<typename VectorType::execution_space,
                               typename MatrixType::exec_space>,
                "Kokkos dense matrix/vector execution-space mismatch");
  static_assert(std::is_same_v<typename VectorType::memory_space,
                               typename MatrixType::memory_space>,
                "Kokkos dense matrix/vector memory-space mismatch");
  auto A_mat{GetDenseMat<MatrixType>(A)};
  auto x_vec{VectorType::FromNVectorConst(x)};
  auto y_vec{VectorType::FromNVector(y)};

  const auto blocks = A_mat->Blocks();
  const auto rows   = A_mat->BlockRows();
  const auto cols   = A_mat->BlockCols();

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();
  auto x_data = x_vec->View();
  auto y_data = y_vec->View();

  // Use batched or single gemv to do y = alpha * A * x + beta * y
  if (blocks > 1)
  {
    using team_policy = typename MatrixType::team_policy;
    using member_type = typename MatrixType::member_type;
    using size_type   = typename MatrixType::size_type;

    Kokkos::parallel_for(
      "sunmatvec_batch",
      team_policy(A_exec, static_cast<int>(blocks), Kokkos::AUTO, Kokkos::AUTO),
      KOKKOS_LAMBDA(const member_type& team_member) {
        const int idx  = team_member.league_rank();
        auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(),
                                         Kokkos::ALL());
        auto x_subdata =
          Kokkos::subview(x_data,
                          Kokkos::pair<size_type, size_type>(idx * cols,
                                                             (idx + 1) * cols));
        auto y_subdata =
          Kokkos::subview(y_data,
                          Kokkos::pair<size_type, size_type>(idx * rows,
                                                             (idx + 1) * rows));
        KokkosBlas::TeamVectorGemv<
          member_type, KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Algo::Gemv::Unblocked>::invoke(team_member,
                                                        SUN_RCONST(1.0),
                                                        A_subdata, x_subdata,
                                                        SUN_RCONST(0.0),
                                                        y_subdata);
      });
  }
  else
  {
    auto A_subdata = Kokkos::subview(A_data, 0, Kokkos::ALL(), Kokkos::ALL());
    KokkosBlas::gemv("N", SUN_RCONST(1.0), A_subdata, x_data, SUN_RCONST(0.0),
                     y_data);
  }

  return SUN_SUCCESS;
}

} // namespace impl

// =============================================================================
// Public namespace
// =============================================================================

// -----------------------------------------------------------------------------
// Kokkos dense matrix class, convertible to a SUNMatrix
// -----------------------------------------------------------------------------

template<class ExecutionSpace, class MemorySpace>
class DenseMatrix : public sundials::impl::BaseMatrix,
                    public sundials::ConvertibleTo<SUNMatrix>
{
public:
  using exec_space   = ExecutionSpace;
  using memory_space = MemorySpace;
  using view_type    = Kokkos::View<sunrealtype***, memory_space>;
  using size_type    = typename view_type::size_type;
  using range_policy = Kokkos::MDRangePolicy<exec_space, Kokkos::Rank<3>>;
  using team_policy  = typename Kokkos::TeamPolicy<exec_space>;
  using member_type  = typename Kokkos::TeamPolicy<exec_space>::member_type;

  // Default constructor - means the matrix must be copied or moved to
  DenseMatrix() = default;

  // Single matrix constructors
  DenseMatrix(size_type rows, size_type cols, SUNContext sunctx)
    : DenseMatrix(1, rows, cols, exec_space(), sunctx)
  {}

  DenseMatrix(size_type rows, size_type cols, exec_space ex, SUNContext sunctx)
    : DenseMatrix(1, rows, cols, ex, sunctx)
  {}

  // Block-diagonal matrix constructors
  DenseMatrix(size_type blocks, size_type block_rows, size_type block_cols,
              SUNContext sunctx)
    : DenseMatrix(blocks, block_rows, block_cols, exec_space(), sunctx)
  {}

  // Block-diagonal matrix with user-supplied execution space instance
  DenseMatrix(size_type blocks, size_type block_rows, size_type block_cols,
              exec_space ex, SUNContext sunctx)
    : sundials::impl::BaseMatrix(sunctx),
      exec_space_(ex),
      view_("sunmat_view", blocks, block_rows, block_cols)
  {
    initSUNMatrix();
  }

  // Move constructor
  DenseMatrix(DenseMatrix&& that_matrix) noexcept
    : sundials::impl::BaseMatrix(std::forward<DenseMatrix>(that_matrix)),
      exec_space_(std::move(that_matrix.exec_space_)),
      view_(std::move(that_matrix.exec_space_))
  {
    initSUNMatrix();
  }

  // Copy constructor
  DenseMatrix(const DenseMatrix& that_matrix)
    : sundials::impl::BaseMatrix(that_matrix),
      exec_space_(that_matrix.exec_space_),
      view_("sunmat_view", that_matrix.Blocks(), that_matrix.BlockRows(),
            that_matrix.BlockCols())
  {
    initSUNMatrix();
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
    view_       = view_type("sunmat_view", rhs.Blocks(), rhs.BlockRows(),
                            rhs.BlockCols());

    sundials::impl::BaseMatrix::operator=(rhs);

    return *this;
  }

  // Default destructor since all members are RAII
  virtual ~DenseMatrix() = default;

  // Get the Kokkos execution space
  exec_space ExecSpace() { return exec_space_; }

  // Get the Kokkos view
  view_type View() { return view_; }

  // Get the number of blocks
  size_type Blocks() const { return static_cast<size_type>(view_.extent(0)); }

  // Get the number of rows in a block
  size_type BlockRows() const
  {
    return static_cast<size_type>(view_.extent(1));
  }

  // Get the number of columns in a block
  size_type BlockCols() const
  {
    return static_cast<size_type>(view_.extent(2));
  }

  // Get the number of rows
  size_type Rows() const
  {
    return static_cast<size_type>(view_.extent(0) * view_.extent(1));
  }

  // Get the number of columns
  size_type Cols() const
  {
    return static_cast<size_type>(view_.extent(0) * view_.extent(2));
  }

  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods

  // Implicit conversion to a SUNMatrix
  operator SUNMatrix() noexcept override { return object_.get(); }

  // Implicit conversion to SUNMatrix
  operator SUNMatrix() const noexcept override { return object_.get(); }

  // Explicit conversion to a SUNMatrix
  SUNMatrix get() noexcept override { return object_.get(); }

  // Explicit conversion to a SUNMatrix
  SUNMatrix get() const noexcept override { return object_.get(); }

private:
  exec_space exec_space_; // Kokkos execution space
  view_type view_;        // Matrix data view [blocks, rows, cols]

  void initSUNMatrix()
  {
    // Use Kokkos' default rank-1 layout for this execution space (LayoutRight
    // on the host, and the corresponding supported default on accelerators).
    // The callback's concrete-handle check below rejects other layouts rather
    // than interpreting their mapping as a flattened vector.
    using vec_view = Kokkos::View<sunrealtype*, exec_space,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using vec_type = VectorFor<vec_view>;
    using mat_type = DenseMatrix<exec_space, memory_space>;

    this->object_->content = this;

    this->object_->ops->getid     = impl::SUNMatGetID_KokkosDense;
    this->object_->ops->clone     = impl::SUNMatClone_KokkosDense<mat_type>;
    this->object_->ops->destroy   = impl::SUNMatDestroy_KokkosDense<mat_type>;
    this->object_->ops->zero      = impl::SUNMatZero_KokkosDense<mat_type>;
    this->object_->ops->copy      = impl::SUNMatCopy_KokkosDense<mat_type>;
    this->object_->ops->scaleadd  = impl::SUNMatScaleAdd_KokkosDense<mat_type>;
    this->object_->ops->scaleaddi = impl::SUNMatScaleAddI_KokkosDense<mat_type>;
    this->object_->ops->matvec =
      impl::SUNMatMatvec_KokkosDense<vec_type, mat_type>;
  }
};

} // namespace kokkos
} // namespace sundials

#endif
