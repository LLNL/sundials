#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#ifndef _SUNMATRIX_GINKGOBLOCK_HPP
#define _SUNMATRIX_GINKGOBLOCK_HPP

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

//
// Dense
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoDenseBlock(std::shared_ptr<gko::Executor> gko_exec,
                                     sunindextype nblocks, sunindextype M, sunindextype N, SUNContext sunctx);

//
// CSR
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoCsrBlock(std::shared_ptr<gko::Executor> gko_exec,
                                   sunindextype nblocks, sunindextype M, sunindextype N, SUNContext sunctx);

#ifdef __cplusplus
}
#endif


namespace sundials
{
namespace ginkgo
{

using GkoBatchDenseMat = gko::matrix::BatchDense<sunrealtype>;
using GkoBatchCsrMat = gko::matrix::BatchCsr<sunrealtype, sunindextype>;
using GkoBatchVecType = GkoBatchDenseMat;

template<typename GkoBatchMatType>
class BlockMatrix : public BaseMatrix
{

public:
  BlockMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype nblocks, sunindextype M, sunindextype N, SUNContext sunctx)
  { }

  BlockMatrix(const BlockMatrix<GkoBatchMatType>& A)
  { }

  BlockMatrix& operator=(const BlockMatrix& A)
  {
    return *this = BlockMatrix<GkoBatchMatType>(A);
  }

  ~BlockMatrix()
  {
    if (sunmtx_) free(sunmtx_);
  }

  operator SUNMatrix() { return sunmtx_; }
  operator SUNMatrix() const { return sunmtx_; }

  const gko::dim<2>& blockDim() const { return gkodim().at(0); }

  const gko::batch_dim<>& gkodim() const { return gkomtx_->get_size(); }

  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx_->get_executor(); }

  std::shared_ptr<GkoBatchMatType> gkomtx() { return gkomtx_; }

  bool isBlockDiagonal() const override { return true; }

  sunindextype numBlocks() const { return gkodim().get_num_batch_entries(); }

  long int workspaceSize() const override
  {
    return blockDim()[0] * blockDim()[1] * numBlocks();
  }

  void Zero()
  {
    gkomtx_ = gko::share(GkoBatchMatType::create(gkoexec(), gkodim()));
  }

private:
  std::shared_ptr<GkoBatchMatType> gkomtx_;
  SUNMatrix sunmtx_;

};

//
// Class method specializations.
//

template<>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype nblocks, sunindextype M, sunindextype N, SUNContext sunctx)
  : gkomtx_(gko::share(GkoBatchDenseMat::create(
      gko_exec, gko::batch_dim<>(nblocks, gko::dim<2>(M, N))
    ))),
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

template<>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype nblocks, sunindextype M, sunindextype N, SUNContext sunctx)
  : gkomtx_(gko::share(GkoBatchCsrMat::create(
       gko_exec, nblocks, gko::dim<2>(M, N)
     ))),
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

template<>
inline BlockMatrix<GkoBatchDenseMat>::BlockMatrix(const BlockMatrix<GkoBatchDenseMat>& A)
  : gkomtx_(gko::share(GkoBatchDenseMat::create(A.gkoexec(), A.gkodim())))
{
  const SUNMatrix Asun = A;
  sunmtx_ = SUNMatNewEmpty(Asun->sunctx);
  sunmtx_->content = this;
  SUNMatCopyOps(Asun, sunmtx_);
}

template<>
inline BlockMatrix<GkoBatchCsrMat>::BlockMatrix(const BlockMatrix<GkoBatchCsrMat>& A)
  : gkomtx_(gko::share(GkoBatchCsrMat::create(
      A.gkoexec(), A.numBlocks(),  A.blockDim())
    ))
{
  const SUNMatrix Asun = A;
  sunmtx_ = SUNMatNewEmpty(Asun->sunctx);
  sunmtx_->content = this;
  SUNMatCopyOps(Asun, sunmtx_);
}

template<>
inline void BlockMatrix<GkoBatchCsrMat>::Zero()
{
  gkomtx_ = gko::share(GkoBatchCsrMat::create(
    gkoexec(), numBlocks(), blockDim()
  ));
}


//
// Non-class methods
//

template<typename GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, GkoBatchVecType* x, GkoBatchVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template<typename GkoBatchMatType>
void Matvec(BlockMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y)
{
  auto batch_vec_stride = gko::batch_stride(A.numBlocks(), 1);

  if (x != y)
  {
    sunrealtype* x_arr =
      (x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x);
    const sunindextype xvec_len = N_VGetLength(x);
    auto batch_xvec_size = gko::batch_dim<>(A.numBlocks(), gko::dim<2>(xvec_len, 1));
    auto xvec_view = gko::Array<sunrealtype>::view(A.gkoexec(), xvec_len, x_arr);
    auto xvec = GkoBatchVecType::create(A.gkoexec(), batch_xvec_size, xvec_view, batch_vec_stride);

    sunrealtype* y_arr =
      (y->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(y) : N_VGetArrayPointer(y);
    const sunindextype yvec_len = N_VGetLength(y);
    auto batch_yvec_size = gko::batch_dim<>(A.numBlocks(), gko::dim<2>(yvec_len, 1));
    auto yvec_view = gko::Array<sunrealtype>::view(A.gkoexec(), yvec_len, y_arr);
    auto yvec = GkoBatchVecType::create(A.gkoexec(), batch_yvec_size, yvec_view, batch_vec_stride);

    // y = Ax
    A.gkomtx()->apply(xvec.get(), yvec.get());
  }
  else
  {
    sunrealtype* x_arr =
      (x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x);
    const sunindextype xvec_len = N_VGetLength(x);
    auto batch_xvec_size = gko::batch_dim<>(A.numBlocks(), gko::dim<2>(xvec_len, 1));
    auto xvec_view = gko::Array<sunrealtype>::view(A.gkoexec(), xvec_len, x_arr);
    auto xvec = GkoBatchVecType::create(A.gkoexec(), batch_xvec_size, xvec_view, batch_vec_stride);

    // x = Ax
    A.gkomtx()->apply(xvec.get(), xvec.get());
  }

}

template<typename GkoMatType, typename GkoBatchMatType>
std::shared_ptr<GkoBatchMatType> CreateBatchIdentity(BlockMatrix<GkoBatchMatType>& tmpl)
{
  auto Iblock = GkoMatType::create(tmpl.gkoexec(), tmpl.blockDim());
  Iblock->read(gko::matrix_data<sunrealtype, sunindextype>::diag(Iblock->get_size(), 1.0));
  auto I = gko::share(GkoBatchMatType::create(tmpl.gkoexec(), tmpl.numBlocks(), Iblock.get()));
  return I;
}

template<typename GkoBatchMatType>
void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A, BlockMatrix<GkoBatchMatType>& B)
{
  throw std::runtime_error("ScaleAdd not implemented for matrix type\n");
}

template<>
inline void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchDenseMat>& A, BlockMatrix<GkoBatchDenseMat>& B)
{
  auto I = CreateBatchIdentity<GkoDenseMat, GkoBatchDenseMat>(A);
  auto one = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  auto constant = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  I->apply(
    GkoBatchDenseMat::create(A.gkoexec(), A.numBlocks(), one.get()).get(),
    B.gkomtx().get(),
    GkoBatchDenseMat::create(A.gkoexec(), A.numBlocks(), constant.get()).get(),
    A.gkomtx().get()
  );
}

template<>
inline void ScaleAdd(const sunrealtype c, BlockMatrix<GkoBatchCsrMat>& A, BlockMatrix<GkoBatchCsrMat>& B)
{
  auto I = CreateBatchIdentity<GkoCsrMat, GkoBatchCsrMat>(A);
  auto one = gko::initialize<GkoCsrMat>({1.0}, A.gkoexec());
  auto constant = gko::initialize<GkoCsrMat>({c}, A.gkoexec());
  I->apply(
    GkoBatchCsrMat::create(A.gkoexec(), A.numBlocks(), one.get()).get(),
    B.gkomtx().get(),
    GkoBatchCsrMat::create(A.gkoexec(), A.numBlocks(), constant.get()).get(),
    A.gkomtx().get()
  );
}

template<typename GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchMatType>& A)
{
  throw std::runtime_error("ScaleAddI not implemented for matrix type\n");
}

template<>
inline void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchDenseMat>& A)
{
  auto I = CreateBatchIdentity<GkoDenseMat, GkoBatchDenseMat>(A);
  auto one = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  auto constant = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  I->apply(
    GkoBatchDenseMat::create(A.gkoexec(), A.numBlocks(), one.get()).get(),
    I.get(),
    GkoBatchDenseMat::create(A.gkoexec(), A.numBlocks(), constant.get()).get(),
    A.gkomtx().get()
  );
}

template<>
inline void ScaleAddI(const sunrealtype c, BlockMatrix<GkoBatchCsrMat>& A)
{
  auto I = CreateBatchIdentity<GkoCsrMat, GkoBatchCsrMat>(A);
  auto one = gko::initialize<GkoCsrMat>({1.0}, A.gkoexec());
  auto constant = gko::initialize<GkoCsrMat>({c}, A.gkoexec());
  I->apply(
    GkoBatchCsrMat::create(A.gkoexec(), A.numBlocks(), one.get()).get(),
    I.get(),
    GkoBatchCsrMat::create(A.gkoexec(), A.numBlocks(), constant.get()).get(),
    A.gkomtx().get()
  );
}

// template<typename GkoBatchMatType>
// void Zero(BlockMatrix<GkoBatchMatType>& A)
// {
//   // TODO: Do we need to actually zero the matrix or can we just create a new empty matrix?
//   Fill(A, 0.0);
// }

template<typename GkoBatchMatType>
void Fill(BlockMatrix<GkoBatchMatType>& A, const sunrealtype c)
{
  // This will trigger a host-device transfer.
  // TODO: Look at using gko::device_matrix_data if the executor is a device executor.
  // A.gkomtx()->read(gko::matrix_data<sunrealtype, sunindextype>(A.gkodim(), c));
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

}//ginkgo
}//sundials

#endif