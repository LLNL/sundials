#include <cstring>
#include <ginkgo/ginkgo.hpp>
#include <memory>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#ifndef _SUNMATRIX_GINKGOBLOCK_HPP
#define _SUNMATRIX_GINKGOBLOCK_HPP

namespace sundials {
namespace ginkgo {

using GkoBatchDenseMat = gko::matrix::BatchDense<sunrealtype>;
using GkoBatchCsrMat   = gko::matrix::BatchCsr<sunrealtype, sunindextype>;
using GkoBatchVecType  = GkoBatchDenseMat;

template<typename GkoBatchMatType>
class BlockMatrix {
public:
  BlockMatrix(sunindextype num_blocks, sunindextype M, sunindextype N, std::shared_ptr<const gko::Executor> gko_exec,
              SUNContext sunctx)
  {
    std::runtime_error("Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  BlockMatrix(sunindextype num_blocks, sunindextype M, sunindextype N, sunindextype num_nonzeros,
              std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
  {
    std::runtime_error("Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  BlockMatrix(std::shared_ptr<GkoBatchMatType> gko_mat, SUNContext sunctx)
  {
    std::runtime_error("Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  operator SUNMatrix() { return sunmtx_.get(); }

  operator SUNMatrix() const { return sunmtx_.get(); }

  SUNMatrix get() { return sunmtx_.get(); }

  SUNMatrix get() const { return sunmtx_.get(); }

  SUNContext sunctx() const { return sunmtx_->sunctx; }

  std::shared_ptr<GkoBatchMatType> gkomtx() const { return gkomtx_; }

  const gko::dim<2>& blockSize(sunindextype block = 0) const { return gkodim().at(block); }

  sunindextype blockDim(sunindextype block = 0, sunindextype dim = 0) const { return gkodim().at(block)[dim]; }

  const gko::batch_dim<>& gkodim() const { return gkomtx_->get_size(); }

  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx_->get_executor(); }

  std::shared_ptr<GkoBatchMatType> gkomtx() { return gkomtx_; }

  sunindextype numBlocks() const { return gkodim().get_num_batch_entries(); }

  bool isBlockDiagonal() const { return true; }

  long int workspaceSize() const { return blockDim(0, 0) * blockDim(0, 1) * numBlocks(); }

  gko::LinOp* gkolinop() { return static_cast<gko::LinOp*>(gkomtx().get()); }

private:
  std::shared_ptr<GkoBatchMatType> gkomtx_;
  std::unique_ptr<struct _generic_SUNMatrix> sunmtx_;
  std::unique_ptr<struct _generic_SUNMatrix_Ops> sunmtx_ops_;

  void initSUNMatrix(SUNContext sunctx);
};

//
// Methods that operate on SUNMatrix
//

template<typename GkoBatchMatType>
SUNMatrix_ID SUNMatGetID_GinkgoBlock(SUNMatrix A)
{
  return std::is_same<GkoBatchCsrMat, GkoBatchMatType>::value ? SUNMATRIX_GINKGOBLOCKCSR : SUNMATRIX_GINKGOBLOCKDENSE;
}

template<typename GkoBatchMatType>
SUNMatrix SUNMatClone_GinkgoBlock(SUNMatrix A)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);

  auto new_mat =
      std::is_same<GkoBatchCsrMat, GkoBatchMatType>::value
          ? new ginkgo::BlockMatrix<GkoBatchMatType>(Amat->numBlocks(), Amat->blockDim(0, 0), Amat->blockDim(0, 1),
                                                     Amat->gkomtx()->get_num_stored_elements(), Amat->gkoexec(),
                                                     Amat->sunctx())
          : new ginkgo::BlockMatrix<GkoBatchMatType>(Amat->numBlocks(), Amat->blockDim(0, 0), Amat->blockDim(0, 1),
                                                     Amat->gkoexec(), Amat->sunctx());

  return new_mat->get();
}

template<typename GkoBatchMatType>
void SUNMatDestroy_GinkgoBlock(SUNMatrix A)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);
  delete Amat;
  return;
}

template<typename GkoBatchMatType>
int SUNMatZero_GinkgoBlock(SUNMatrix A)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);
  Zero(*Amat);
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchMatType>
int SUNMatCopy_GinkgoBlock(SUNMatrix A, SUNMatrix B)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);
  Copy(*Amat, *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchMatType>
int SUNMatScaleAdd_GinkgoBlock(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);
  ScaleAdd(c, *Amat, *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchMatType>
int SUNMatScaleAddI_GinkgoBlock(sunrealtype c, SUNMatrix A)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);
  ScaleAddI(c, *Amat);
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchMatType>
int SUNMatMatvecSetup_GinkgoBlock(SUNMatrix A)
{
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchMatType>
int SUNMatMatvec_GinkgoBlock(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto Amat = static_cast<BlockMatrix<GkoBatchMatType>*>(A->content);
  Matvec(*Amat, x, y);
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchMatType>
int SUNMatSpace_GinkgoBlock(SUNMatrix A, long int* lenrw, long int* leniw)
{
  auto Amat = static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(A->content);
  *lenrw    = Amat->workspaceSize();
  *leniw    = 0;
  return SUNMAT_SUCCESS;
}

//
// Class method specializations.
//

template<typename GkoBatchMatType>
void BlockMatrix<GkoBatchMatType>::initSUNMatrix(SUNContext sunctx)
{
  sunmtx_->content = this;
  sunmtx_->ops     = sunmtx_ops_.get();
  sunmtx_->sunctx  = sunctx;

  std::memset(sunmtx_->ops, 0, sizeof(_generic_SUNMatrix_Ops));
  sunmtx_->ops->getid = SUNMatGetID_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->clone = SUNMatClone_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->zero  = SUNMatZero_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->copy  = SUNMatCopy_GinkgoBlock<GkoBatchMatType>;
  // Ginkgo does not provide what we need for ScaleAdd with BatchCsr yet
  sunmtx_->ops->scaleadd =
      std::is_same<GkoBatchCsrMat, GkoBatchMatType>::value ? NULL : SUNMatScaleAdd_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->matvec      = SUNMatMatvec_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->destroy     = SUNMatDestroy_GinkgoBlock<GkoBatchMatType>;
  sunmtx_->ops->space       = SUNMatSpace_GinkgoBlock<GkoBatchMatType>;
}

template<>
BlockMatrix<GkoBatchDenseMat>::BlockMatrix(sunindextype num_blocks, sunindextype M, sunindextype N,
                                           std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(gko::share(GkoBatchDenseMat::create(gko_exec, gko::batch_dim<2>(num_blocks, gko::dim<2>(M, N))))),
      sunmtx_(std::make_unique<_generic_SUNMatrix>()), sunmtx_ops_(std::make_unique<_generic_SUNMatrix_Ops>())
{
  initSUNMatrix(sunctx);
}

template<>
BlockMatrix<GkoBatchDenseMat>::BlockMatrix(std::shared_ptr<GkoBatchDenseMat> gko_mat, SUNContext sunctx)
    : gkomtx_(gko_mat), sunmtx_(std::make_unique<_generic_SUNMatrix>()),
      sunmtx_ops_(std::make_unique<_generic_SUNMatrix_Ops>())
{
  initSUNMatrix(sunctx);
}

template<>
BlockMatrix<GkoBatchCsrMat>::BlockMatrix(sunindextype num_blocks, sunindextype M, sunindextype N,
                                         std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(gko::share(GkoBatchCsrMat::create(gko_exec, num_blocks, gko::dim<2>(M, N)))),
      sunmtx_(std::make_unique<_generic_SUNMatrix>()), sunmtx_ops_(std::make_unique<_generic_SUNMatrix_Ops>())
{
  initSUNMatrix(sunctx);
}

template<>
BlockMatrix<GkoBatchCsrMat>::BlockMatrix(sunindextype num_blocks, sunindextype M, sunindextype N,
                                         sunindextype num_nonzeros, std::shared_ptr<const gko::Executor> gko_exec,
                                         SUNContext sunctx)
    : gkomtx_(gko::share(GkoBatchCsrMat::create(gko_exec, num_blocks, gko::dim<2>(M, N), num_nonzeros))),
      sunmtx_(std::make_unique<_generic_SUNMatrix>()), sunmtx_ops_(std::make_unique<_generic_SUNMatrix_Ops>())
{
  initSUNMatrix(sunctx);
}

template<>
BlockMatrix<GkoBatchCsrMat>::BlockMatrix(std::shared_ptr<GkoBatchCsrMat> gko_mat, SUNContext sunctx)
    : gkomtx_(gko_mat), sunmtx_(std::make_unique<_generic_SUNMatrix>()),
      sunmtx_ops_(std::make_unique<_generic_SUNMatrix_Ops>())
{
  initSUNMatrix(sunctx);
}

//
// Non-class methods
//

inline std::unique_ptr<GkoBatchVecType> WrapVector(std::shared_ptr<const gko::Executor> gko_exec,
                                                   sunindextype num_blocks, N_Vector x)
{
  sunrealtype* x_arr          = (x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x);
  const sunindextype xvec_len = N_VGetLength(x);
  auto batch_vec_stride       = gko::batch_stride(num_blocks, 1);
  auto batch_xvec_size        = gko::batch_dim<>(num_blocks, gko::dim<2>(xvec_len / num_blocks, 1));
  auto xvec_view              = gko::Array<sunrealtype>::view(gko_exec, xvec_len, x_arr);
  return GkoBatchVecType::create(gko_exec, batch_xvec_size, std::move(xvec_view), batch_vec_stride);
}

inline std::unique_ptr<const GkoBatchVecType> WrapConstVector(std::shared_ptr<const gko::Executor> gko_exec,
                                                              sunindextype num_blocks, N_Vector x)
{
  sunrealtype* x_arr          = (x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x);
  const sunindextype xvec_len = N_VGetLength(x);
  auto batch_vec_stride       = gko::batch_stride(num_blocks, 1);
  auto batch_xvec_size        = gko::batch_dim<>(num_blocks, gko::dim<2>(xvec_len / num_blocks, 1));
  auto xvec_view              = gko::Array<sunrealtype>::const_view(gko_exec, xvec_len, x_arr);
  return GkoBatchVecType::create_const(gko_exec, batch_xvec_size, std::move(xvec_view), batch_vec_stride);
}

template<typename GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, GkoBatchVecType* x, GkoBatchVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template<typename GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y) {
    auto x_vec = WrapConstVector(A.gkoexec(), A.numBlocks(), x);
    auto y_vec = WrapVector(A.gkoexec(), A.numBlocks(), y);

    // y = Ax
    A.gkomtx()->apply(x_vec.get(), y_vec.get());
  } else {
    auto x_vec = WrapVector(A.gkoexec(), A.numBlocks(), x);

    // x = Ax
    A.gkomtx()->apply(x_vec.get(), x_vec.get());
  }
}

template<typename GkoBatchMatrixType>
void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchMatrixType>& A, BlockMatrix<GkoBatchMatrixType>& B)
{
  // NOTE: This is not implemented by Ginkgo for BatchCsr yet
  const auto I    = gko::matrix::BatchIdentity<sunrealtype>::create(A.gkoexec(), A.gkodim());
  const auto one  = gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec());
  const auto cmat = gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec());
  // A = B + cA
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<typename GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A)
{
  const auto one  = gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec());
  const auto cmat = gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec());
  // A = 1*I + c*A = cA + I
  A.gkomtx()->add_scaled_identity(one.get(), cmat.get());
}

template<typename GkoBatchMatType>
void Zero(BlockMatrix<GkoBatchMatType>& A)
{
  // TODO(CJB): This seems really inefficient, need to find a better way.
  for (auto& block : A.gkomtx()->unbatch()) {
    block->scale(gko::initialize<GkoDenseMat>({0.0}, A.gkoexec()).get());
  }
}

template<>
inline void Zero(BlockMatrix<GkoBatchDenseMat>& A)
{
  A.gkomtx()->scale(gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {0.0}, A.gkoexec()).get());
}

template<typename GkoBatchMatType>
void Copy(BlockMatrix<GkoBatchMatType>& A, BlockMatrix<GkoBatchMatType>& B)
{
  B.gkomtx()->copy_from(A.gkomtx().get());
}

template<typename GkoBatchMatType>
void Print(BlockMatrix<GkoBatchMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.gkomtx().get());
}

} // namespace ginkgo
} // namespace sundials

#endif