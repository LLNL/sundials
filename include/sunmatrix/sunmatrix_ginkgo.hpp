
#include <sundials/sundials_matrix.h>
#include <ginkgo/ginkgo.hpp>


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

// Implementation specific methods
SUNDIALS_EXPORT
SUNMatrix SUNMatrix_Ginkgo(std::shared_ptr<gko::Executor> gko_exec, sunindextype M,
                           sunindextype N, SUNContext sunctx);

// SUNMatrix overrides
SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_Ginkgo(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_Ginkgo(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_Ginkgo(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_Ginkgo(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_Ginkgo(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_Ginkgo(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_Ginkgo(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvecSetup_Ginkgo(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_Ginkgo(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_Ginkgo(SUNMatrix A, long int *lenrw, long int *leniw);

// Additional functions
SUNDIALS_EXPORT int SUNMatFill_Ginkgo(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatPrint_Ginkgo(SUNMatrix A);

#ifdef __cplusplus
}
#endif


namespace sundials
{
namespace ginkgo
{

using GkoVecType = gko::matrix::Dense<realtype>;


template<typename GkoMatType>
class GinkgoMatrix
{

public:
  static std::shared_ptr<GkoMatType> CreateIdentity(std::shared_ptr<const gko::Executor> exec, const gko::dim<2>& size);

  GinkgoMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype M, sunindextype N, SUNContext sunctx);
  ~GinkgoMatrix();

  operator std::shared_ptr<GkoMatType>() { return gkomtx(); }
  operator SUNMatrix() { return sunmtx_; }

  std::shared_ptr<gko::LinOp> gkolinop() { return gkolinop_; }
  std::shared_ptr<GkoMatType> gkomtx() { return gkomtx_; }
  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx_->get_executor(); }

  const gko::dim<2>& Size() const;
  void Reset();
  void Evaluate();

  void ScaleAdd(const realtype c, GinkgoMatrix& B);
  void ScaleAddI(const realtype c);
  long int WorkspaceSize();

private:
  SUNMatrix sunmtx_;
  std::shared_ptr<GkoMatType> gkomtx_;
  std::shared_ptr<gko::LinOp> gkolinop_;
  bool linop_was_updated_;

};

//
// These functions do not modify the matrix data or linear operator.
//

template<typename GkoMatType>
void Matvec(GinkgoMatrix<GkoMatType>& A, GkoVecType* x, GkoVecType* y)
{
  A.gkolinop()->apply(x, y);
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

    A.gkolinop()->apply(x_vec.get(), y_vec.get());
  }
  else
  {
    const sunindextype x_len = N_VGetLength(x);
    auto x_vec = GkoVecType::create(A.gkoexec(), gko::dim<2>(x_len, 1),
      gko::Array<realtype>::view(A.gkoexec(), x_len, N_VGetArrayPointer(x)), 1);
    A.gkolinop()->apply(x_vec.get(), x_vec.get());
  }

}

//
// These functions instantly modify the matrix data.
//

template<typename GkoMatType>
void Zero(GinkgoMatrix<GkoMatType>& A)
{
  const auto zero = gko::initialize<GkoVecType>({0.0}, A.gkoexec());
  // Reset because we are going to zero out the linear operator / matrix.
  A.Reset();
  A.gkomtx()->scale(zero.get());
}

template<typename GkoMatType>
void Fill(GinkgoMatrix<GkoMatType>& A, const realtype c)
{
  // Reset because we are going to fill out the linear operator / matrix.
  A.Reset();
  A.gkomtx()->fill(c);
}

template<typename GkoMatType>
void Copy(GinkgoMatrix<GkoMatType>& A, GinkgoMatrix<GkoMatType>& B)
{
  // We need to make sure the linear operator form of A is evaluated before
  // copying the matrix data, otherwise the matrix and linear operator might
  // not be equivalent.
  A.Evaluate();
  B.gkomtx()->copy_from(A.gkomtx().get());
}

template<typename GkoMatType>
void Print(GinkgoMatrix<GkoMatType>& A)
{
  gko::write(std::cout, A.gkomtx().get());
}

//
// Class methods
//

template<typename GkoMatType>
inline std::shared_ptr<GkoMatType> GinkgoMatrix<GkoMatType>::CreateIdentity(std::shared_ptr<const gko::Executor> exec, const gko::dim<2>& size)
{
  auto I = gko::share(GkoMatType::create(exec, size));
  I->fill(1.0);
  auto Idia = I->extract_diagonal();
  auto Icsr = gko::matrix::Csr<realtype>::create(exec, size);
  Idia->move_to(Icsr.get());
  Icsr->move_to(I.get());
  return I;
}

template<typename GkoMatType>
inline GinkgoMatrix<GkoMatType>::~GinkgoMatrix()
{
  if (sunmtx_) free(sunmtx_);
}

template<typename GkoMatType>
inline GinkgoMatrix<GkoMatType>::GinkgoMatrix(std::shared_ptr<gko::Executor> gko_exec, sunindextype M, sunindextype N, SUNContext sunctx)
  : sunmtx_(SUNMatNewEmpty(sunctx)),
    gkomtx_(gko::share(GkoMatType::create(gko_exec, gko::dim<2>(M, N)))),
    gkolinop_(gkomtx_),
    linop_was_updated_(false)
{
  sunmtx_->content = this;

  sunmtx_->ops->clone       = NULL;
  sunmtx_->ops->zero        = SUNMatZero_Ginkgo;
  sunmtx_->ops->copy        = SUNMatCopy_Ginkgo;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_Ginkgo;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_Ginkgo;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_Ginkgo;
  sunmtx_->ops->matvec      = SUNMatMatvec_Ginkgo;
  sunmtx_->ops->destroy     = SUNMatDestroy_Ginkgo;
  sunmtx_->ops->space       = SUNMatSpace_Ginkgo;
}

template<typename GkoMatType>
inline const gko::dim<2>& GinkgoMatrix<GkoMatType>::Size() const { return gkomtx_->get_size(); }

template<typename GkoMatType>
inline void GinkgoMatrix<GkoMatType>::Reset()
{
  // Resets the linear operator to be the matrix.
  gkolinop_ = gkomtx_;
}

template<typename GkoMatType>
inline void GinkgoMatrix<GkoMatType>::Evaluate()
{
  // If the linear operator form of the matrix was updated, then we are going
  // to evaluate it by applying it to the identity, and then reset the
  // linear operator to be the matrix.
  if (linop_was_updated_)
  {
    // TODO: should we keep I around?
    const auto I = CreateIdentity(gkoexec(), Size());
    gkolinop_->apply(I.get(), gkomtx_.get());
    gkolinop_ = gkomtx_;
    linop_was_updated_ = false;
  }
}

template<typename GkoMatType>
inline void GinkgoMatrix<GkoMatType>::ScaleAdd(const realtype c, GinkgoMatrix& B)
{
  // B = c*A + 1.0*B = c*A + B
  // B is the linear operator c*A+B, it is not computed upon creation.
  B.gkolinop_ =
    gko::share(
      gko::Combination<realtype>::create(
        gko::initialize<GkoVecType>({1.0}, gkoexec()),
        B.gkolinop_,
        gko::initialize<GkoVecType>({c}, gkoexec()),
        gkolinop_
      )
    );
  B.linop_was_updated_ = true;
}

template<typename GkoMatType>
inline void GinkgoMatrix<GkoMatType>::ScaleAddI(const realtype c)
{
  // A = c*A + 1.0*I = c*A + I
  // A is the linear operator c*A+I, it is not computed upon creation.
  gkolinop_ =
    gko::share(
      gko::Combination<realtype>::create(
        gko::initialize<GkoVecType>({c}, gkoexec()),
        gkolinop_,
        gko::initialize<GkoVecType>({1.0}, gkoexec()),
        gko::matrix::Identity<realtype>::create(gkoexec(), Size()[0])
      )
    );
  linop_was_updated_ = true;
}

template<typename GkoMatType>
inline long int GinkgoMatrix<GkoMatType>::WorkspaceSize()
{
  // TODO: if we keep I around, update this calculation
  return Size()[0] * Size()[1];
}


}//ginkgo
}//sundials
