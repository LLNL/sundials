
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
SUNMatrix SUNMatrix_GinkgoDense(std::shared_ptr<gko::Executor> gko_exec,
                                sunindextype M, sunindextype N, SUNContext sunctx);

// SUNMatrix overrides
SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_GinkgoDense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_GinkgoDense(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_GinkgoDense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvecSetup_GinkgoDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_GinkgoDense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_GinkgoDense(SUNMatrix A, long int *lenrw, long int *leniw);

// Additional functions
SUNDIALS_EXPORT int SUNMatFill_GinkgoDense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatPrint_GinkgoDense(SUNMatrix A);

//
// CSR
//

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_GinkgoCsr(std::shared_ptr<gko::Executor> gko_exec,
                              sunindextype M, sunindextype N, SUNContext sunctx);

// SUNMatrix overrides
SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_GinkgoCsr(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_GinkgoCsr(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_GinkgoCsr(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvecSetup_GinkgoCsr(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_GinkgoCsr(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_GinkgoCsr(SUNMatrix A, long int *lenrw, long int *leniw);

// Additional functions
SUNDIALS_EXPORT int SUNMatFill_GinkgoCsr(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatPrint_GinkgoCsr(SUNMatrix A);

#ifdef __cplusplus
}
#endif


namespace sundials
{
namespace ginkgo
{
namespace matrix
{

using GkoVecType = gko::matrix::Dense<realtype>;


template<typename GkoMatType>
class GinkgoMatrix
{

public:
  GinkgoMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype M, sunindextype N, SUNContext sunctx)
  { }

  GinkgoMatrix(const GinkgoMatrix<GkoMatType>& A)
    : gkomtx_(gko::share(GkoMatType::create(A.gkoexec(), A.gkodim())))
  {
    const SUNMatrix Asun = A;
    sunmtx_ = SUNMatNewEmpty(Asun->sunctx);
    sunmtx_->content = this;
    SUNMatCopyOps(Asun, sunmtx_);
  }

  GinkgoMatrix& operator=(const GinkgoMatrix& A)
  {
    return *this = GinkgoMatrix<GkoMatType>(A);
  }

  ~GinkgoMatrix()
  {
    if (sunmtx_) free(sunmtx_);
  }

  operator SUNMatrix() { return sunmtx_; }
  operator SUNMatrix() const { return sunmtx_; }

  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx_->get_executor(); }
  const gko::dim<2>& gkodim() const { return gkomtx_->get_size(); }
  std::shared_ptr<GkoMatType> gkomtx() { return gkomtx_; }

  long int WorkspaceSize() { return gkodim()[0] * gkodim()[1]; }

private:
  SUNMatrix sunmtx_;
  std::shared_ptr<GkoMatType> gkomtx_;

};

//
// Class methods
//

template<>
inline GinkgoMatrix<gko::matrix::Dense<realtype>>::GinkgoMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype M, sunindextype N, SUNContext sunctx)
  : sunmtx_(SUNMatNewEmpty(sunctx)),
    gkomtx_(gko::share(gko::matrix::Dense<realtype>::create(gko_exec, gko::dim<2>(M, N))))
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
inline GinkgoMatrix<gko::matrix::Csr<realtype>>::GinkgoMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype M, sunindextype N, SUNContext sunctx)
  : sunmtx_(SUNMatNewEmpty(sunctx)),
    gkomtx_(gko::share(gko::matrix::Csr<realtype>::create(gko_exec, gko::dim<2>(M, N))))
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

//
// Non-class methods
//

template<typename GkoMatType>
void Matvec(GinkgoMatrix<GkoMatType>& A, GkoVecType* x, GkoVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template<typename GkoMatType>
void Matvec(GinkgoMatrix<GkoMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y)
  {
    const sunindextype x_len = N_VGetLength(x);
    auto x_vec = GkoVecType::create_const(A.gkoexec(), gko::dim<2>(x_len, 1),
      gko::Array<realtype>::const_view(A.gkoexec(), x_len, N_VGetArrayPointer(x)), 1);

    const sunindextype y_len = N_VGetLength(y);
    auto y_vec = GkoVecType::create(A.gkoexec(), gko::dim<2>(y_len, 1),
      gko::Array<realtype>::view(A.gkoexec(), y_len, N_VGetArrayPointer(y)), 1);

    A.gkomtx()->apply(x_vec.get(), y_vec.get());
  }
  else
  {
    const sunindextype x_len = N_VGetLength(x);
    auto x_vec = GkoVecType::create(A.gkoexec(), gko::dim<2>(x_len, 1),
      gko::Array<realtype>::view(A.gkoexec(), x_len, N_VGetArrayPointer(x)), 1);
    A.gkomtx()->apply(x_vec.get(), x_vec.get());
  }

}

template<typename GkoMatType>
std::shared_ptr<GkoMatType> CreateIdentity(std::shared_ptr<const gko::Executor> gko_exec, const gko::dim<2>& gko_dim)
{
  auto I = gko::share(GkoMatType::create(gko_exec, gko_dim));
  // I->fill(1.0);
  I->read(gko::matrix_data<realtype, sunindextype>(gko_dim, 1.0));
  auto Idia = I->extract_diagonal();
  auto Icsr = gko::matrix::Csr<realtype>::create(gko_exec, gko_dim);
  Idia->move_to(Icsr.get());
  Icsr->move_to(I.get());
  return I;
}

template<typename GkoMatType>
void ScaleAdd(const realtype c, GinkgoMatrix<GkoMatType>& A, GinkgoMatrix<GkoMatType>& B)
{
  // B.gkolinop_ =
  //   gko::share(
  //     gko::Combination<realtype>::create(
  //       gko::initialize<GkoMatType>({1.0}, gkoexec()),
  //       B.gkolinop_,
  //       gko::initialize<GkoMatType>({c}, gkoexec()),
  //       gkolinop_
  //     )
  //   );
  // B.linop_was_updated_ = true;
  auto I = CreateIdentity<GkoMatType>(A.gkoexec(), A.gkodim());
  I->apply(
    gko::initialize<gko::matrix::Dense<realtype>>({1.0}, A.gkoexec()).get(),
    B.gkomtx().get(),
    gko::initialize<gko::matrix::Dense<realtype>>({c}, A.gkoexec()).get(),
    A.gkomtx().get()
  );
}

template<typename GkoMatType>
void ScaleAddI(const realtype c, GinkgoMatrix<GkoMatType>& A)
{
  // A is the linear operator c*A+I, it is not computed upon creation.
  // gkolinop_ =
  //   gko::share(
  //     gko::Combination<realtype>::create(
  //       gko::initialize<GkoMatType>({c}, gkoexec()),
  //       gkolinop_,
  //       gko::initialize<GkoMatType>({1.0}, gkoexec()),
  //       gko::matrix::Identity<realtype>::create(gkoexec(), gkodim()[0])
  //     )
  //   );
  // linop_was_updated_ = true;
  auto I = CreateIdentity<GkoMatType>(A.gkoexec(), A.gkodim());
  I->apply(
    gko::initialize<gko::matrix::Dense<realtype>>({1.0}, A.gkoexec()).get(),
    I.get(),
    gko::initialize<gko::matrix::Dense<realtype>>({c}, A.gkoexec()).get(),
    A.gkomtx().get()
  );
}

template<typename GkoMatType>
void Zero(GinkgoMatrix<GkoMatType>& A)
{
  Fill(A, 0.0);
}

template<typename GkoMatType>
void Fill(GinkgoMatrix<GkoMatType>& A, const realtype c)
{
  A.gkomtx()->read(gko::matrix_data<realtype, sunindextype>(A.gkodim(), c));
}

template<typename GkoMatType>
void Copy(GinkgoMatrix<GkoMatType>& A, GinkgoMatrix<GkoMatType>& B)
{
  B.gkomtx()->copy_from(A.gkomtx().get());
}

template<typename GkoMatType>
void Print(GinkgoMatrix<GkoMatType>& A)
{
  gko::write(std::cout, A.gkomtx().get());
}

}//matrix
}//ginkgo
}//sundials

#endif