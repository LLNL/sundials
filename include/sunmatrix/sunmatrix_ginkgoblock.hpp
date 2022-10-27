/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------*/

#ifndef _SUNMATRIX_GINKGOBLOCK_HPP
#define _SUNMATRIX_GINKGOBLOCK_HPP

#include <cstring>
#include <ginkgo/ginkgo.hpp>
#include <memory>
#include <sundials/sundials_matrix.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

namespace sundials {
namespace ginkgo {

using GkoBatchDenseMat = gko::matrix::BatchDense<sunrealtype>;
using GkoBatchCsrMat   = gko::matrix::BatchCsr<sunrealtype, sunindextype>;
using GkoBatchVecType  = GkoBatchDenseMat;

// Forward declare BlockMatrix class
template<class GkoBatchMatType>
class BlockMatrix;

//
// Methods that operate on SUNMatrix
//
template<class GkoBatchMatType>
SUNMatrix_ID SUNMatGetID_GinkgoBlock(SUNMatrix A)
{
  return SUNMATRIX_GINKGOBLOCK;
}

template<class GkoBatchMatType>
SUNMatrix SUNMatClone_GinkgoBlock(SUNMatrix A)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  auto new_mat{new BlockMatrix<GkoBatchMatType>(*Amat)};
  return new_mat->get();
}

template<class GkoBatchMatType>
void SUNMatDestroy_GinkgoBlock(SUNMatrix A)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  delete Amat; // NOLINT
  return;
}

// template<class GkoBatchMatType>
// int SUNMatZero_GinkgoBlock(SUNMatrix A)
// {
//   auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
//   Zero(*Amat);
//   return SUNMAT_SUCCESS;
// }

template<class GkoBatchMatType>
int SUNMatCopy_GinkgoBlock(SUNMatrix A, SUNMatrix B)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  Copy(*Amat, *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<class GkoBatchMatType>
int SUNMatScaleAdd_GinkgoBlock(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  ScaleAdd(c, *Amat, *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<class GkoBatchMatType>
int SUNMatScaleAddI_GinkgoBlock(sunrealtype c, SUNMatrix A)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  ScaleAddI(c, *Amat);
  return SUNMAT_SUCCESS;
}

template<class GkoBatchMatType>
int SUNMatMatvec_GinkgoBlock(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  Matvec(*Amat, x, y);
  return SUNMAT_SUCCESS;
}

template<class GkoBatchMatType>
class BlockMatrix : public sundials::impl::BaseMatrix, public sundials::ConvertibleTo<SUNMatrix>
{
public:
  // Default constructor means the matrix must be copied or moved to
  BlockMatrix() = default;

  // We do not have implementations of these two constructors for general GkoBatchMatType
  BlockMatrix(gko::size_type num_blocks, sunindextype M, sunindextype N, std::shared_ptr<const gko::Executor> gko_exec,
              SUNContext sunctx);
  BlockMatrix(gko::size_type num_blocks, sunindextype M, sunindextype N, sunindextype num_nonzeros,
              std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx);

  BlockMatrix(std::shared_ptr<GkoBatchMatType> gko_mat, SUNContext sunctx)
      : gkomtx_(gko_mat), sundials::impl::BaseMatrix(sunctx)
  {
    initSUNMatrix();
  }

  // Move constructor
  BlockMatrix(BlockMatrix&& that_matrix) noexcept
      : gkomtx_(std::move(that_matrix.gkomtx_)), sundials::impl::BaseMatrix(std::forward<BlockMatrix>(that_matrix))
  {}

  // Copy constructor clones the gko::matrix and SUNMatrix
  BlockMatrix(const BlockMatrix& that_matrix)
      : gkomtx_(gko::clone(that_matrix.gkomtx_)), sundials::impl::BaseMatrix(that_matrix)
  {}

  // Move assignment
  BlockMatrix& operator=(BlockMatrix&& rhs) noexcept
  {
    gkomtx_ = std::move(rhs.gkomtx_);
    sundials::impl::BaseMatrix::operator=(std::forward<BlockMatrix>(rhs));
    return *this;
  }

  // Copy assignment clones the gko::matrix and SUNMatrix
  BlockMatrix& operator=(const BlockMatrix& rhs)
  {
    gkomtx_ = gko::clone(rhs.gkomtx_);
    sundials::impl::BaseMatrix::operator=(rhs);
    return *this;
  }

  // Default destructor is fine since all members are RAII
  virtual ~BlockMatrix() = default;

  // Getters
  std::shared_ptr<GkoBatchMatType> gkomtx() const { return gkomtx_; }
  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx()->get_executor(); }
  const gko::batch_dim<2>& gkodim() const { return gkomtx()->get_size(); }
  const gko::dim<2>& blockSize(gko::size_type block = 0) const { return gkodim().at(block); }
  sunindextype blockDim(gko::size_type block = 0, sunindextype dim = 0) const { return gkodim().at(block)[dim]; }
  sunindextype blockNNZ(gko::size_type block = 0) const
  {
    return gkomtx()->get_num_stored_elements() / gkomtx()->get_num_batch_entries();
  }
  sunindextype numBlocks() const { return gkodim().get_num_batch_entries(); }
  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods
  operator SUNMatrix() override { return object_.get(); }
  operator SUNMatrix() const override { return object_.get(); }
  SUNMatrix get() override { return object_.get(); }
  SUNMatrix get() const override { return object_.get(); }

private:
  std::shared_ptr<GkoBatchMatType> gkomtx_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid = SUNMatGetID_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->clone = SUNMatClone_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->copy  = SUNMatCopy_GinkgoBlock<GkoBatchMatType>;
    // Ginkgo does not provide what we need for ScaleAdd except with BatchDense
    this->object_->ops->scaleadd =
        std::is_same<GkoBatchDenseMat, GkoBatchMatType>::value ? SUNMatScaleAdd_GinkgoBlock<GkoBatchMatType> : nullptr;
    this->object_->ops->scaleaddi = SUNMatScaleAddI_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->matvec    = SUNMatMatvec_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->destroy   = SUNMatDestroy_GinkgoBlock<GkoBatchMatType>;
  }
};

//
// Class method specializations.
//

template<>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(gko::size_type num_blocks, sunindextype M, sunindextype N,
                                                  std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(GkoBatchDenseMat::create(gko_exec, gko::batch_dim<2>(num_blocks, gko::dim<2>(M, N)))),
      sundials::impl::BaseMatrix(sunctx)
{
  initSUNMatrix();
}

template<>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(gko::size_type num_blocks, sunindextype M, sunindextype N,
                                                sunindextype num_nonzeros,
                                                std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(GkoBatchCsrMat::create(gko_exec, gko::batch_dim<>(num_blocks, gko::dim<2>(M, N)), num_nonzeros)),
      sundials::impl::BaseMatrix(sunctx)
{
  initSUNMatrix();
}

//
// Non-class methods
//

inline std::unique_ptr<GkoBatchVecType> WrapBatchVector(std::shared_ptr<const gko::Executor> gko_exec,
                                                        gko::size_type num_blocks, N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x)};
  const auto xvec_len{N_VGetLength(x)};
  auto batch_vec_stride{gko::batch_stride(num_blocks, 1)};
  auto batch_xvec_size{gko::batch_dim<2>(num_blocks, gko::dim<2>(xvec_len / num_blocks, 1))};
  auto xvec_view{gko::Array<sunrealtype>::view(gko_exec, xvec_len, x_arr)};
  return GkoBatchVecType::create(gko_exec, batch_xvec_size, std::move(xvec_view), batch_vec_stride);
}

inline std::unique_ptr<const GkoBatchVecType> WrapConstBatchVector(std::shared_ptr<const gko::Executor> gko_exec,
                                                                   gko::size_type num_blocks, N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x)};
  const auto xvec_len{N_VGetLength(x)};
  auto batch_vec_stride{gko::batch_stride(num_blocks, 1)};
  auto batch_xvec_size{gko::batch_dim<2>(num_blocks, gko::dim<2>(xvec_len / num_blocks, 1))};
  auto xvec_view{gko::Array<sunrealtype>::const_view(gko_exec, xvec_len, x_arr)};
  return GkoBatchVecType::create_const(gko_exec, batch_xvec_size, std::move(xvec_view), batch_vec_stride);
}

template<class GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, GkoBatchVecType* x, GkoBatchVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template<class GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y) {
    auto x_vec = WrapConstBatchVector(A.gkoexec(), A.numBlocks(), x);
    auto y_vec = WrapBatchVector(A.gkoexec(), A.numBlocks(), y);

    // y = Ax
    A.gkomtx()->apply(x_vec.get(), y_vec.get());
  }
  else {
    auto x_vec = WrapBatchVector(A.gkoexec(), A.numBlocks(), x);

    // x = Ax
    A.gkomtx()->apply(x_vec.get(), x_vec.get());
  }
}

template<typename GkoBatchMatrixType>
void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchMatrixType>& A, BlockMatrix<GkoBatchMatrixType>& B)
{
  // NOTE: This is not implemented by Ginkgo for BatchCsr yet
  const auto I{gko::matrix::BatchIdentity<sunrealtype>::create(A.gkoexec(), A.gkodim())};
  const auto one{gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec())};
  const auto cmat{gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec())};
  // A = B + cA
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<class GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A)
{
  const auto one{gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec())};
  const auto cmat{gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec())};
  // A = 1*I + c*A = cA + I
  A.gkomtx()->add_scaled_identity(one.get(), cmat.get());
}

template<class GkoBatchMatType>
void Copy(BlockMatrix<GkoBatchMatType>& A, BlockMatrix<GkoBatchMatType>& B)
{
  B.gkomtx()->copy_from(A.gkomtx().get());
}

template<class GkoBatchMatType>
void Print(BlockMatrix<GkoBatchMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.gkomtx().get());
}

} // namespace ginkgo
} // namespace sundials

#endif