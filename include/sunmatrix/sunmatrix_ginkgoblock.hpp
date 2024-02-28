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
#include <ginkgo/core/base/batch_multi_vector.hpp>
#include <ginkgo/ginkgo.hpp>
#include <memory>
#include <sundials/sundials_matrix.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#include "sundials/sundials_types.h"

namespace sundials {
namespace ginkgo {

using gkoblock_indextype = int;

using GkoBatchDenseMat = gko::batch::matrix::Dense<sunrealtype>;
using GkoBatchCsrMat = gko::batch::matrix::Csr<sunrealtype, gkoblock_indextype>;
using GkoBatchVecType = gko::batch::MultiVector<sunrealtype>;

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
  return new_mat->Convert();
}

template<class GkoBatchMatType>
void SUNMatDestroy_GinkgoBlock(SUNMatrix A)
{
  auto Amat{static_cast<BlockMatrix<GkoBatchMatType>*>(A->content)};
  delete Amat; // NOLINT
  return;
}

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
  ScaleAdd(c, *Amat,
           *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(B->content));
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
class BlockMatrix : public sundials::impl::BaseMatrix,
                    public sundials::ConvertibleTo<SUNMatrix>
{
public:
  // Default constructor means the matrix must be copied or moved to
  BlockMatrix() = default;

  // We do not have implementations of these two constructors for general
  // GkoBatchMatType
  BlockMatrix(gko::size_type num_blocks, sunindextype M, sunindextype N,
              std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx);
  BlockMatrix(gko::size_type num_blocks, sunindextype M, sunindextype N,
              sunindextype num_nonzeros,
              std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx);

  BlockMatrix(std::shared_ptr<GkoBatchMatType> gko_mat, SUNContext sunctx)
    : gkomtx_(gko_mat), sundials::impl::BaseMatrix(sunctx)
  {
    initSUNMatrix();
  }

  // Move constructor
  BlockMatrix(BlockMatrix&& that_matrix) noexcept
    : gkomtx_(std::move(that_matrix.gkomtx_)),
      sundials::impl::BaseMatrix(std::forward<BlockMatrix>(that_matrix))
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
  std::shared_ptr<GkoBatchMatType> GkoMtx() const { return gkomtx_; }

  std::shared_ptr<const gko::Executor> GkoExec() const
  {
    return GkoMtx()->get_executor();
  }

  const gko::batch_dim<2>& GkoSize() const { return GkoMtx()->get_size(); }

  // const gko::dim<2>& blockSize(gko::size_type block = 0) const { return
  // GkoSize().at(block); } sunindextype blockDim(gko::size_type block = 0,
  // sunindextype dim = 0) const { return GkoSize().at(block)[dim]; }
  // sunindextype blockNNZ(gko::size_type block = 0) const
  // {
  // return GkoMtx()->get_num_stored_elements() /
  // GkoMtx()->get_num_batch_entries();
  // }

  sunindextype NumBlocks() const { return GkoSize().get_num_batch_items(); }

  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods
  operator SUNMatrix() override { return object_.get(); }

  operator SUNMatrix() const override { return object_.get(); }

  SUNMatrix Convert() override { return object_.get(); }

  SUNMatrix Convert() const override { return object_.get(); }

private:
  std::shared_ptr<GkoBatchMatType> gkomtx_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid = SUNMatGetID_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->clone = SUNMatClone_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->copy  = SUNMatCopy_GinkgoBlock<GkoBatchMatType>;
    // Ginkgo does not provide what we need for ScaleAdd except with BatchDense
    this->object_->ops->scaleadd  = SUNMatScaleAdd_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->scaleaddi = SUNMatScaleAddI_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->matvec    = SUNMatMatvec_GinkgoBlock<GkoBatchMatType>;
    this->object_->ops->destroy   = SUNMatDestroy_GinkgoBlock<GkoBatchMatType>;
  }
};

//
// Class method specializations.
//

template<>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(
  gko::size_type num_blocks, sunindextype M, sunindextype N,
  std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
  : gkomtx_(
      GkoBatchDenseMat::create(gko_exec, gko::batch_dim<2>(num_blocks,
                                                           gko::dim<2>(M, N)))),
    sundials::impl::BaseMatrix(sunctx)
{
  initSUNMatrix();
}

template<>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(
  gko::size_type num_blocks, sunindextype M, sunindextype N,
  sunindextype num_nonzeros, std::shared_ptr<const gko::Executor> gko_exec,
  SUNContext sunctx)
  : gkomtx_(
      GkoBatchCsrMat::create(gko_exec,
                             gko::batch_dim<2>(num_blocks, gko::dim<2>(M, N)),
                             num_nonzeros)),
    sundials::impl::BaseMatrix(sunctx)
{
  initSUNMatrix();
}

//
// Non-class methods
//

inline std::unique_ptr<GkoBatchVecType> WrapBatchVector(
  std::shared_ptr<const gko::Executor> gko_exec, gko::size_type num_blocks,
  N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x)
                                               : N_VGetArrayPointer(x)};
  const auto xvec_len{N_VGetLength(x)};
  auto batch_xvec_size{
    gko::batch_dim<2>(num_blocks, gko::dim<2>(xvec_len / num_blocks, 1))};
  auto xvec_view{gko::array<sunrealtype>::view(gko_exec, xvec_len, x_arr)};
  return GkoBatchVecType::create(gko_exec, batch_xvec_size, std::move(xvec_view));
}

inline std::unique_ptr<const GkoBatchVecType> WrapConstBatchVector(
  std::shared_ptr<const gko::Executor> gko_exec, gko::size_type num_blocks,
  N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x)
                                               : N_VGetArrayPointer(x)};
  const auto xvec_len{N_VGetLength(x)};
  auto batch_xvec_size{
    gko::batch_dim<2>(num_blocks, gko::dim<2>(xvec_len / num_blocks, 1))};
  auto xvec_view{gko::array<sunrealtype>::const_view(gko_exec, xvec_len, x_arr)};
  return GkoBatchVecType::create_const(gko_exec, batch_xvec_size,
                                       std::move(xvec_view));
}

template<class GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, GkoBatchVecType* x,
            GkoBatchVecType* y)
{
  A.GkoMtx()->apply(x, y);
}

template<class GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y)
  {
    auto x_vec = WrapConstBatchVector(A.GkoExec(), A.NumBlocks(), x);
    auto y_vec = WrapBatchVector(A.GkoExec(), A.NumBlocks(), y);

    // y = Ax
    A.GkoMtx()->apply(x_vec.get(), y_vec.get());
  }
  else
  {
    auto x_vec = WrapBatchVector(A.GkoExec(), A.NumBlocks(), x);

    // x = Ax
    A.GkoMtx()->apply(x_vec.get(), x_vec.get());
  }
}

void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchDenseMat>& A,
              BlockMatrix<GkoBatchDenseMat>& B)
{
  auto cmat =
    GkoBatchVecType::create(A.GkoExec(),
                            gko::batch_dim<2>(A.NumBlocks(), gko::dim<2>(1, 1)));
  cmat->fill(c);
  // A = B + cA
  A.GkoMtx()->scale_add(cmat.get(), B.GkoMtx().get());
}

void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchCsrMat>& A,
              BlockMatrix<GkoBatchCsrMat>& B)
{
  // NOTE: This is not implemented by Ginkgo for BatchCsr yet
  throw("scale_add not implemented for gko::batch::matrix::Csr");
}

template<class GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A)
{
  auto one =
    GkoBatchVecType::create(A.GkoExec(),
                            gko::batch_dim<2>(A.NumBlocks(), gko::dim<2>(1, 1)));
  one->fill(gko::one<sunrealtype>());
  auto cmat =
    GkoBatchVecType::create(A.GkoExec(),
                            gko::batch_dim<2>(A.NumBlocks(), gko::dim<2>(1, 1)));
  cmat->fill(c);
  // A = 1*I + c*A = cA + I
  A.GkoMtx()->add_scaled_identity(one.get(), cmat.get());
}

template<class GkoBatchMatType>
void Copy(BlockMatrix<GkoBatchMatType>& A, BlockMatrix<GkoBatchMatType>& B)
{
  B.GkoMtx()->copy_from(A.GkoMtx().get());
}

template<class GkoBatchMatType>
void Print(BlockMatrix<GkoBatchMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.GkoMtx().get());
}

} // namespace ginkgo
} // namespace sundials

#endif
