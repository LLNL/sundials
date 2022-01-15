
#include <memory>

#include <sundials/sundials_matrix.h>
#include <ginkgo/ginkgo.hpp>

#ifndef _SUNMATRIX_GINKGO_HPP
#define _SUNMATRIX_GINKGO_HPP

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

//
// Dense
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoDense(sunindextype M, sunindextype N,
                                std::shared_ptr<gko::Executor> gko_exec,
                                SUNContext sunctx);

// SUNMatrix overrides
SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_GinkgoDense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_GinkgoDense(sunrealtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_GinkgoDense(sunrealtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvecSetup_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_GinkgoDense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_GinkgoDense(SUNMatrix A, long int *lenrw, long int *leniw);

// Additional functions
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_LData(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_BlockRows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_BlockColumns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_BlockLData(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_GinkgoDense_NumBlocks(SUNMatrix A);
SUNDIALS_EXPORT gko::matrix::Dense<sunrealtype>* SUNMatrix_GinkgoDense_GetGkoMat(SUNMatrix A);
SUNDIALS_EXPORT gko::matrix::BatchDense<sunrealtype>* SUNMatrix_GinkgoDense_GetGkoBatchMat(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatPrint_GinkgoDense(SUNMatrix A);

//
// CSR
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoCsr(sunindextype M, sunindextype N, sunindextype NNZ, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx);

// SUNMatrix overrides
SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_GinkgoCsr(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_GinkgoCsr(sunrealtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_GinkgoCsr(sunrealtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvecSetup_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_GinkgoCsr(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_GinkgoCsr(SUNMatrix A, long int *lenrw, long int *leniw);

// Additional functions
SUNDIALS_EXPORT int SUNMatPrint_GinkgoCsr(SUNMatrix A);

#ifdef __cplusplus
}
#endif


namespace sundials
{
namespace ginkgo
{

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;
using GkoCsrMat = gko::matrix::Csr<sunrealtype, sunindextype>;
using GkoVecType = GkoDenseMat;

class BaseMatrix
{
public:
  virtual bool isBlockDiagonal() const = 0;
  virtual long int workspaceSize() const = 0;
};

template<typename GkoMatType>
class Matrix : public BaseMatrix
{

public:
  Matrix(sunindextype M, sunindextype N, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  { std::runtime_error("Constructor is not implemented for the Ginkgo matrix type provided\n"); }

  Matrix(sunindextype M, sunindextype N, sunindextype NNZ, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  { std::runtime_error("Constructor is not implemented for the Ginkgo matrix type provided\n"); }

  Matrix(std::shared_ptr<GkoMatType> gko_mat, SUNContext sunctx)
  { std::runtime_error("Constructor is not implemented for the Ginkgo matrix type provided\n"); }

  Matrix(const Matrix<GkoMatType>& A)
    : gkomtx_(gko::share(GkoMatType::create(A.gkoexec(), A.gkodim())))
  {
    const SUNMatrix Asun = A;
    sunmtx_.reset(SUNMatNewEmpty(Asun->sunctx));
    sunmtx_->content = this;
    SUNMatCopyOps(Asun, sunmtx_.get());
  }

  Matrix(Matrix<GkoMatType>&& A)
    : gkomtx_(A.gkomtx_),
      sunmtx_(std::move(A.sunmtx_))
  { }

  Matrix& operator=(const Matrix& A)
  {
    return *this = Matrix<GkoMatType>(A);
  }

  Matrix& operator=(Matrix&& A)
  {
    sunmtx_ = std::move(A.sunmtx_);
    return *this;
  }

  ~Matrix()
  { }

  operator SUNMatrix() { return sunmtx_.get(); }
  operator SUNMatrix() const { return sunmtx_.get(); }

  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx_->get_executor(); }
  const gko::dim<2>& gkodim() const { return gkomtx_->get_size(); }
  std::shared_ptr<GkoMatType> gkomtx() { return gkomtx_; }

  bool isBlockDiagonal() const override { return false; }
  long int workspaceSize() const override { return gkodim()[0] * gkodim()[1]; }

private:
  std::shared_ptr<GkoMatType> gkomtx_;
  std::unique_ptr<struct _generic_SUNMatrix> sunmtx_;

};

//
// Class method specializations.
//

template<>
inline Matrix<GkoDenseMat>::Matrix(sunindextype M, sunindextype N, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  : gkomtx_(gko::share(GkoDenseMat::create(gko_exec, gko::dim<2>(M, N)))),
    sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_GinkgoDense;
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
inline Matrix<GkoCsrMat>::Matrix(sunindextype M, sunindextype N, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  : gkomtx_(gko::share(GkoCsrMat::create(gko_exec, gko::dim<2>(M, N)))),
    sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_GinkgoCsr;
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
inline Matrix<GkoCsrMat>::Matrix(sunindextype M, sunindextype N, sunindextype NNZ, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  : gkomtx_(gko::share(GkoCsrMat::create(gko_exec, gko::dim<2>(M, N), NNZ))),
    sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_GinkgoCsr;
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
inline Matrix<GkoDenseMat>::Matrix(std::shared_ptr<GkoDenseMat> gko_mat, SUNContext sunctx)
  : gkomtx_(gko_mat),
    sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_GinkgoDense;
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
inline Matrix<GkoCsrMat>::Matrix(std::shared_ptr<GkoCsrMat> gko_mat, SUNContext sunctx)
  : gkomtx_(gko_mat),
    sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_GinkgoCsr;
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

//
// Non-class methods
//

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, GkoVecType* x, GkoVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, N_Vector x, N_Vector y)
{
  sunrealtype* x_arr = (x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x);

  if (x != y)
  {
    sunrealtype* y_arr = (y->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(y) : N_VGetArrayPointer(y);

    const sunindextype x_len = N_VGetLength(x);
    auto x_vec = GkoVecType::create_const(A.gkoexec(), gko::dim<2>(x_len, 1),
      gko::Array<sunrealtype>::const_view(A.gkoexec(), x_len, x_arr), 1);

    const sunindextype y_len = N_VGetLength(y);
    auto y_vec = GkoVecType::create(A.gkoexec(), gko::dim<2>(y_len, 1),
      gko::Array<sunrealtype>::view(A.gkoexec(), y_len, y_arr), 1);

    // y = Ax
    A.gkomtx()->apply(x_vec.get(), y_vec.get());
  }
  else
  {
    const sunindextype x_len = N_VGetLength(x);
    auto x_vec = GkoVecType::create(A.gkoexec(), gko::dim<2>(x_len, 1),
      gko::Array<sunrealtype>::view(A.gkoexec(), x_len, x_arr), 1);

    // x = Ax
    A.gkomtx()->apply(x_vec.get(), x_vec.get());
  }

}

template<typename GkoMatType>
std::shared_ptr<GkoMatType> CreateIdentity(std::shared_ptr<const gko::Executor> gko_exec, const gko::dim<2>& gko_dim)
{
  auto I = gko::share(GkoMatType::create(gko_exec, gko_dim));
  I->read(gko::matrix_data<sunrealtype, sunindextype>(gko_dim, 1.0));
  auto Idia = I->extract_diagonal();
  auto Icsr = GkoCsrMat::create(gko_exec, gko_dim);
  Idia->move_to(Icsr.get());
  Icsr->move_to(I.get());
  return I;
}

template<typename GkoMatType>
void ScaleAdd(const sunrealtype c, Matrix<GkoMatType>& A, Matrix<GkoMatType>& B)
{
  throw std::runtime_error("ScaleAdd not implemented for matrix type\n");
}

template<>
inline void ScaleAdd(const sunrealtype c, Matrix<GkoDenseMat>& A, Matrix<GkoDenseMat>& B)
{
  const auto I = CreateIdentity<GkoDenseMat>(A.gkoexec(), A.gkodim());
  const auto one = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  const auto cmat = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<>
inline void ScaleAdd(const sunrealtype c, Matrix<gko::matrix::Csr<sunrealtype>>& A, Matrix<gko::matrix::Csr<sunrealtype>>& B)
{
  const auto I = gko::matrix::Identity<sunrealtype>::create(A.gkoexec(), A.gkodim());
  const auto one = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  const auto cmat = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  // A = B + cA
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<typename GkoMatType>
void ScaleAddI(const sunrealtype c, Matrix<GkoMatType>& A)
{
  auto I = CreateIdentity<GkoMatType>(A.gkoexec(), A.gkodim());
  // A = I + cA
  I->apply(
    gko::initialize<GkoDenseMat>({1.0}, A.gkoexec()).get(),
    I.get(),
    gko::initialize<GkoDenseMat>({c}, A.gkoexec()).get(),
    A.gkomtx().get()
  );
}

template<typename GkoMatType>
void Zero(Matrix<GkoMatType>& A)
{
  A.gkomtx()->scale(gko::initialize<GkoDenseMat>({0.0}, A.gkoexec()).get());
}

template<>
inline void Zero(Matrix<GkoDenseMat>& A)
{
  A.gkomtx()->fill(0.0);
}

template<typename GkoMatType>
void Copy(Matrix<GkoMatType>& A, Matrix<GkoMatType>& B)
{
  B.gkomtx()->copy_from(A.gkomtx().get());
}

template<typename GkoMatType>
void Print(Matrix<GkoMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.gkomtx().get());
}

}//ginkgo
}//sundials

#endif