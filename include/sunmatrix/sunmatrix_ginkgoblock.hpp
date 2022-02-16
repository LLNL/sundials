#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#ifndef _SUNMATRIX_GINKGOBLOCK_HPP
#define _SUNMATRIX_GINKGOBLOCK_HPP

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

constexpr static int debug_print_off = 0;
#define debugstream    \
  if (debug_print_off) \
  {}                   \
  else                 \
    std::cerr

//
// Dense
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoDenseBlock(sunindextype nblocks, sunindextype M,
                                     sunindextype N,
                                     std::shared_ptr<gko::Executor> gko_exec,
                                     SUNContext sunctx);

//
// CSR
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoCsrBlock(sunindextype nblocks, sunindextype M,
                                   sunindextype N, sunindextype NNZ,
                                   std::shared_ptr<gko::Executor> gko_exec,
                                   SUNContext sunctx);

#ifdef __cplusplus
}
#endif

namespace sundials {
namespace ginkgo {

using GkoBatchDenseMat = gko::matrix::BatchDense<sunrealtype>;
using GkoBatchCsrMat   = gko::matrix::BatchCsr<sunrealtype, sunindextype>;
using GkoBatchVecType  = GkoBatchDenseMat;

template <typename GkoBatchMatType> class BlockMatrix : public BaseMatrix {
public:
  BlockMatrix(sunindextype nblocks, sunindextype M, sunindextype N,
              std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  {
    std::runtime_error(
        "Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  BlockMatrix(sunindextype nblocks, sunindextype M, sunindextype N,
              sunindextype NNZ, std::shared_ptr<gko::Executor> gko_exec,
              SUNContext sunctx)
  {
    std::runtime_error(
        "Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  BlockMatrix(std::shared_ptr<GkoBatchMatType> gko_mat, SUNContext sunctx)
  {
    std::runtime_error(
        "Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  BlockMatrix(const BlockMatrix<GkoBatchMatType>& A) {}

  BlockMatrix(BlockMatrix<GkoBatchMatType>&& A)
      : gkomtx_(A.gkomtx_), sunmtx_(std::move(A.sunmtx_))
  {}

  BlockMatrix& operator=(const BlockMatrix& A)
  {
    return *this = BlockMatrix<GkoBatchMatType>(A);
  }

  BlockMatrix& operator=(BlockMatrix&& A)
  {
    sunmtx_ = std::move(A.sunmtx_);
    return *this;
  }

  ~BlockMatrix() {}

  operator SUNMatrix() { return sunmtx_.get(); }
  operator SUNMatrix() const { return sunmtx_.get(); }

  const gko::dim<2>& blockDim() const { return gkodim().at(0); }
  const gko::batch_dim<>& gkodim() const { return gkomtx_->get_size(); }
  std::shared_ptr<const gko::Executor> gkoexec() const
  {
    return gkomtx_->get_executor();
  }
  std::shared_ptr<GkoBatchMatType> gkomtx() { return gkomtx_; }

  sunindextype numBlocks() const { return gkodim().get_num_batch_entries(); }

  bool isBlockDiagonal() const override { return true; }
  long int workspaceSize() const override
  {
    return blockDim()[0] * blockDim()[1] * numBlocks();
  }

  BlockMatrix reset()
  {
    gkomtx_ = gko::share(gkomtx_->clone());
    return *this;
  }

  gko::LinOp* gkolinop() override { return nullptr; }

private:
  std::shared_ptr<GkoBatchMatType> gkomtx_;
  std::unique_ptr<struct _generic_SUNMatrix> sunmtx_;
};

//
// Class method specializations.
//

template <>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(
    sunindextype nblocks, sunindextype M, sunindextype N,
    std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(gko::share(
          GkoBatchDenseMat::create(gko_exec,
                                   gko::batch_dim<>(nblocks, gko::dim<2>(M, N))))),
      sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->clone       = SUNMatClone_GinkgoDense;
  sunmtx_->ops->zero        = SUNMatZero_GinkgoDense;
  sunmtx_->ops->copy        = SUNMatCopy_GinkgoDense;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_GinkgoDense;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_GinkgoDense;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_GinkgoDense;
  sunmtx_->ops->matvec      = SUNMatMatvec_GinkgoDense;
  sunmtx_->ops->destroy     = SUNMatDestroy_GinkgoDense;
  sunmtx_->ops->space       = SUNMatSpace_GinkgoDense;
}

template <>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(
    sunindextype nblocks, sunindextype M, sunindextype N,
    std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(gko::share(
          GkoBatchCsrMat::create(gko_exec, nblocks, gko::dim<2>(M, N)))),
      sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->clone       = SUNMatClone_GinkgoCsr;
  sunmtx_->ops->zero        = SUNMatZero_GinkgoCsr;
  sunmtx_->ops->copy        = SUNMatCopy_GinkgoCsr;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_GinkgoCsr;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_GinkgoCsr;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_GinkgoCsr;
  sunmtx_->ops->matvec      = SUNMatMatvec_GinkgoCsr;
  sunmtx_->ops->destroy     = SUNMatDestroy_GinkgoCsr;
  sunmtx_->ops->space       = SUNMatSpace_GinkgoCsr;
}

template <>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(
    sunindextype nblocks, sunindextype M, sunindextype N, sunindextype NNZ,
    std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(gko::share(
          GkoBatchCsrMat::create(gko_exec, nblocks, gko::dim<2>(M, N), NNZ))),
      sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->clone       = SUNMatClone_GinkgoCsr;
  sunmtx_->ops->zero        = SUNMatZero_GinkgoCsr;
  sunmtx_->ops->copy        = SUNMatCopy_GinkgoCsr;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_GinkgoCsr;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_GinkgoCsr;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_GinkgoCsr;
  sunmtx_->ops->matvec      = SUNMatMatvec_GinkgoCsr;
  sunmtx_->ops->destroy     = SUNMatDestroy_GinkgoCsr;
  sunmtx_->ops->space       = SUNMatSpace_GinkgoCsr;
}

template <>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(
    std::shared_ptr<GkoBatchDenseMat> gko_mat, SUNContext sunctx)
    : gkomtx_(gko_mat), sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->clone       = SUNMatClone_GinkgoDense;
  sunmtx_->ops->zero        = SUNMatZero_GinkgoDense;
  sunmtx_->ops->copy        = SUNMatCopy_GinkgoDense;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_GinkgoDense;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_GinkgoDense;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_GinkgoDense;
  sunmtx_->ops->matvec      = SUNMatMatvec_GinkgoDense;
  sunmtx_->ops->destroy     = SUNMatDestroy_GinkgoDense;
  sunmtx_->ops->space       = SUNMatSpace_GinkgoDense;
}

template <>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(
    std::shared_ptr<GkoBatchCsrMat> gko_mat, SUNContext sunctx)
    : gkomtx_(gko_mat), sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->clone       = SUNMatClone_GinkgoCsr;
  sunmtx_->ops->zero        = SUNMatZero_GinkgoCsr;
  sunmtx_->ops->copy        = SUNMatCopy_GinkgoCsr;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_GinkgoCsr;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_GinkgoCsr;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_GinkgoCsr;
  sunmtx_->ops->matvec      = SUNMatMatvec_GinkgoCsr;
  sunmtx_->ops->destroy     = SUNMatDestroy_GinkgoCsr;
  sunmtx_->ops->space       = SUNMatSpace_GinkgoCsr;
}

template <>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(
    const BlockMatrix<GkoBatchDenseMat>& A)
    : gkomtx_(gko::share(GkoBatchDenseMat::create(A.gkoexec(), A.gkodim())))
{
  const SUNMatrix Asun = A;
  sunmtx_.reset(SUNMatNewEmpty(Asun->sunctx));
  sunmtx_->content = this;
  SUNMatCopyOps(Asun, sunmtx_.get());
}

template <>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(const BlockMatrix<GkoBatchCsrMat>& A)
    : gkomtx_(gko::share(
          GkoBatchCsrMat::create(A.gkoexec(), A.numBlocks(), A.blockDim())))
{
  const SUNMatrix Asun = A;
  sunmtx_.reset(SUNMatNewEmpty(Asun->sunctx));
  sunmtx_->content = this;
  SUNMatCopyOps(Asun, sunmtx_.get());
}

//
// Non-class methods
//

inline std::unique_ptr<GkoBatchVecType> WrapVector(
    std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_blocks,
    N_Vector x)
{
  sunrealtype* x_arr          = (x->ops->nvgetdevicearraypointer)
                                    ? N_VGetDeviceArrayPointer(x)
                                    : N_VGetArrayPointer(x);
  const sunindextype xvec_len = N_VGetLength(x);
  auto batch_vec_stride       = gko::batch_stride(num_blocks, 1);
  auto batch_xvec_size = gko::batch_dim<>(num_blocks, gko::dim<2>(xvec_len, 1));
  auto xvec_view = gko::Array<sunrealtype>::view(gko_exec, xvec_len, x_arr);
  return GkoBatchVecType::create(gko_exec, batch_xvec_size, xvec_view,
                                 batch_vec_stride);
}

template <typename GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, GkoBatchVecType* x,
            GkoBatchVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template <typename GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y)
  {
    auto x_vec = WrapVector(A.gkoexec(), A.numBlocks(), x);
    auto y_vec = WrapVector(A.gkoexec(), A.numBlocks(), y);

    // y = Ax
    A.gkomtx()->apply(gko::lend(x_vec), gko::lend(y_vec));
  }
  else
  {
    auto x_vec = WrapVector(A.gkoexec(), A.numBlocks(), x);

    // x = Ax
    A.gkomtx()->apply(gko::lend(x_vec), gko::lend(x_vec));
  }
}

template <typename GkoMatType, typename GkoBatchMatType>
std::shared_ptr<GkoBatchMatType> CreateBatchIdentity(
    BlockMatrix<GkoBatchMatType>& tmpl)
{
  auto Iblock = GkoMatType::create(tmpl.gkoexec(), tmpl.blockDim());
  Iblock->read(
      gko::matrix_data<sunrealtype, sunindextype>::diag(Iblock->get_size(), 1.0));
  auto I = gko::share(
      GkoBatchMatType::create(tmpl.gkoexec(), tmpl.numBlocks(), Iblock.get()));
  return I;
}

template <typename GkoBatchMatType>
void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A,
              BlockMatrix<GkoBatchMatType>& B)
{
  throw std::runtime_error("ScaleAdd not implemented for matrix type\n");
}

template <>
inline void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchDenseMat>& A,
                     BlockMatrix<GkoBatchDenseMat>& B)
{
  debugstream << "\n>>> Called " << __PRETTY_FUNCTION__ << "\n";
  const auto I = CreateBatchIdentity<GkoDenseMat, GkoBatchDenseMat>(A);
  const auto one =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec());
  const auto constant =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec());
  // TODO: this is not implemented for CUDA (dense or csr) and OMP (csr)
  // A = B + cA
  I->apply(one.get(), B.gkomtx().get(), constant.get(), A.gkomtx().get());
}

template <>
inline void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchCsrMat>& A,
                     BlockMatrix<GkoBatchCsrMat>& B)
{
  debugstream << "\n>>> Called " << __PRETTY_FUNCTION__ << "\n";
  const auto I = CreateBatchIdentity<GkoCsrMat, GkoBatchCsrMat>(A);
  const auto one =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec());
  const auto constant =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec());
  // TODO: this is not implemented for CUDA (dense or csr) and OMP (csr)
  // A = B + cA
  I->apply(one.get(), B.gkomtx().get(), constant.get(), A.gkomtx().get());
}

template <typename GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A)
{
  throw std::runtime_error("ScaleAddI not implemented for matrix type\n");
}

template <>
inline void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchDenseMat>& A)
{
  debugstream << "\n>>> Called " << __PRETTY_FUNCTION__ << "\n";
  const auto I = CreateBatchIdentity<GkoDenseMat, GkoBatchDenseMat>(A);
  const auto one =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec());
  const auto constant =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec());
  // TODO: this is not implemented for CUDA (dense or csr) and OMP (csr)
  // A = I + cA
  I->apply(one.get(), I.get(), constant.get(), A.gkomtx().get());
}

template <>
inline void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchCsrMat>& A)
{
  debugstream << "\n>>> Called " << __PRETTY_FUNCTION__ << "\n";
  const auto I = CreateBatchIdentity<GkoCsrMat, GkoBatchCsrMat>(A);
  const auto one =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {1.0}, A.gkoexec());
  const auto constant =
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {c}, A.gkoexec());
  // TODO: this is not implemented for CUDA (dense or csr) and OMP (csr)
  // A = I + cA
  I->apply(one.get(), I.get(), constant.get(), A.gkomtx().get());
}

template <typename GkoBatchMatType> void Zero(BlockMatrix<GkoBatchMatType>& A)
{
  A.reset();
}

template <> inline void Zero(BlockMatrix<GkoBatchDenseMat>& A)
{
  A.gkomtx()->scale(
      gko::batch_initialize<GkoBatchDenseMat>(A.numBlocks(), {0.0}, A.gkoexec())
          .get());
}

template <typename GkoBatchMatType>
void Copy(BlockMatrix<GkoBatchMatType>& A, BlockMatrix<GkoBatchMatType>& B)
{
  B.gkomtx()->copy_from(A.gkomtx().get());
}

template <typename GkoBatchMatType>
void Print(BlockMatrix<GkoBatchMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.gkomtx().get());
}

} // namespace ginkgo
} // namespace sundials

#endif