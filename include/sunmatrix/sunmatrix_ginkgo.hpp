
#include <ginkgo/ginkgo.hpp>
#include <memory>
#include <sundials/sundials_matrix.h>

#ifndef _SUNMATRIX_GINKGO_HPP
#define _SUNMATRIX_GINKGO_HPP

namespace sundials {
namespace ginkgo {

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;
using GkoCsrMat   = gko::matrix::Csr<sunrealtype, sunindextype>;
using GkoVecType  = GkoDenseMat;

//
// Common matrix class
//

template<typename GkoMatType> class BaseMatrix {
public:
  BaseMatrix(std::shared_ptr<GkoMatType> gko_mat) : gkomtx_(gko_mat) {}

  std::shared_ptr<GkoMatType> gkomtx() const { return gkomtx_; }

  virtual bool isBlockDiagonal() const   = 0;
  virtual long int workspaceSize() const = 0;
  virtual gko::LinOp* gkolinop() { return nullptr; }

protected:
  std::shared_ptr<GkoMatType> gkomtx_;
};

//
// Standard matrix class
//

template<typename GkoMatType> class Matrix : public BaseMatrix<GkoMatType> {
public:
  Matrix(sunindextype M, sunindextype N,
         std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  {
    std::runtime_error(
        "Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  Matrix(sunindextype M, sunindextype N, sunindextype NNZ,
         std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
  {
    std::runtime_error(
        "Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  Matrix(std::shared_ptr<GkoMatType> gko_mat, SUNContext sunctx)
  {
    std::runtime_error(
        "Constructor is not implemented for the Ginkgo matrix type provided\n");
  }

  Matrix(const Matrix<GkoMatType>& A)
      : BaseMatrix<GkoMatType>(
            gko::share(GkoMatType::create(A.gkoexec(), A.gkodim())))
  {
    const SUNMatrix Asun = A;
    sunmtx_.reset(SUNMatNewEmpty(Asun->sunctx));
    sunmtx_->content = this;
    SUNMatCopyOps(Asun, sunmtx_.get());
  }

  Matrix(Matrix<GkoMatType>&& A)
      : BaseMatrix<GkoMatType>(A.gkomtx_), sunmtx_(std::move(A.sunmtx_))
  {}

  Matrix& operator=(const Matrix& A) { return *this = Matrix<GkoMatType>(A); }

  Matrix& operator=(Matrix&& A)
  {
    sunmtx_ = std::move(A.sunmtx_);
    return *this;
  }

  ~Matrix() {}

  operator SUNMatrix() { return sunmtx_.get(); }
  operator SUNMatrix() const { return sunmtx_.get(); }

  using BaseMatrix<GkoMatType>::gkomtx;
  std::shared_ptr<const gko::Executor> gkoexec() const
  {
    return gkomtx()->get_executor();
  }
  const gko::dim<2>& gkodim() const { return gkomtx()->get_size(); }

  bool isBlockDiagonal() const override { return false; }
  long int workspaceSize() const override { return gkodim()[0] * gkodim()[1]; }

  gko::LinOp* gkolinop() override
  {
    return static_cast<gko::LinOp*>(gkomtx().get());
  }

private:
  std::unique_ptr<struct _generic_SUNMatrix> sunmtx_;
};

//
// Methods that operate on SUNMatrix
//

template<typename GkoMatType> SUNMatrix_ID SUNMatGetID_Ginkgo(SUNMatrix A)
{
  return SUNMATRIX_GINKGODENSE;
}

template<typename GkoMatType> SUNMatrix SUNMatClone_Ginkgo(SUNMatrix A)
{
  auto Amat    = static_cast<Matrix<GkoMatType>*>(A->content);
  auto new_mat = new ginkgo::Matrix<GkoMatType>(*Amat);
  return ((SUNMatrix)*new_mat);
}

template<typename GkoMatType> void SUNMatDestroy_Ginkgo(SUNMatrix A)
{
  auto Amat = static_cast<Matrix<GkoMatType>*>(A->content);
  delete Amat;
}

template<typename GkoMatType> int SUNMatZero_Ginkgo(SUNMatrix A)
{
  auto Amat = static_cast<Matrix<GkoMatType>*>(A->content);
  Zero(*Amat);
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType> int SUNMatCopy_Ginkgo(SUNMatrix A, SUNMatrix B)
{
  auto Amat = static_cast<Matrix<GkoMatType>*>(A->content);
  Copy(*Amat, *static_cast<ginkgo::Matrix<GkoMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatScaleAdd_Ginkgo(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto Amat = static_cast<Matrix<GkoMatType>*>(A->content);
  ScaleAdd(c, *Amat, *static_cast<ginkgo::Matrix<GkoMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatScaleAddI_Ginkgo(sunrealtype c, SUNMatrix A)
{
  auto Amat = static_cast<Matrix<GkoMatType>*>(A->content);
  ScaleAddI(c, *Amat);
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType> int SUNMatMatvecSetup_Ginkgo(SUNMatrix A)
{
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatMatvec_Ginkgo(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto Amat = static_cast<Matrix<GkoMatType>*>(A->content);
  Matvec(*Amat, x, y);
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatSpace_Ginkgo(SUNMatrix A, long int* lenrw, long int* leniw)
{
  auto Amat = static_cast<ginkgo::Matrix<GkoMatType>*>(A->content);
  *lenrw    = Amat->workspaceSize();
  *leniw    = 0;
  return SUNMAT_SUCCESS;
}

//
// Class method specializations.
//

template<>
inline Matrix<GkoDenseMat>::Matrix(sunindextype M, sunindextype N,
                                   std::shared_ptr<gko::Executor> gko_exec,
                                   SUNContext sunctx)
    : BaseMatrix<GkoDenseMat>(
          gko::share(GkoDenseMat::create(gko_exec, gko::dim<2>(M, N)))),
      sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->clone       = SUNMatClone_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->zero        = SUNMatZero_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->copy        = SUNMatCopy_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->matvec      = SUNMatMatvec_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->destroy     = SUNMatDestroy_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->space       = SUNMatSpace_Ginkgo<GkoDenseMat>;
}

template<>
inline Matrix<GkoCsrMat>::Matrix(sunindextype M, sunindextype N,
                                 std::shared_ptr<gko::Executor> gko_exec,
                                 SUNContext sunctx)
    : BaseMatrix<GkoCsrMat>(
          gko::share(GkoCsrMat::create(gko_exec, gko::dim<2>(M, N)))),
      sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->clone       = SUNMatClone_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->zero        = SUNMatZero_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->copy        = SUNMatCopy_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->matvec      = SUNMatMatvec_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->destroy     = SUNMatDestroy_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->space       = SUNMatSpace_Ginkgo<GkoCsrMat>;
}

template<>
inline Matrix<GkoCsrMat>::Matrix(sunindextype M, sunindextype N, sunindextype NNZ,
                                 std::shared_ptr<gko::Executor> gko_exec,
                                 SUNContext sunctx)
    : BaseMatrix<GkoCsrMat>(
          gko::share(GkoCsrMat::create(gko_exec, gko::dim<2>(M, N), NNZ))),
      sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->clone       = SUNMatClone_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->zero        = SUNMatZero_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->copy        = SUNMatCopy_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->matvec      = SUNMatMatvec_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->destroy     = SUNMatDestroy_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->space       = SUNMatSpace_Ginkgo<GkoCsrMat>;
}

template<>
inline Matrix<GkoDenseMat>::Matrix(std::shared_ptr<GkoDenseMat> gko_mat,
                                   SUNContext sunctx)
    : BaseMatrix<GkoDenseMat>(gko_mat), sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->clone       = SUNMatClone_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->zero        = SUNMatZero_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->copy        = SUNMatCopy_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->matvec      = SUNMatMatvec_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->destroy     = SUNMatDestroy_Ginkgo<GkoDenseMat>;
  sunmtx_->ops->space       = SUNMatSpace_Ginkgo<GkoDenseMat>;
}

template<>
inline Matrix<GkoCsrMat>::Matrix(std::shared_ptr<GkoCsrMat> gko_mat,
                                 SUNContext sunctx)
    : BaseMatrix<GkoCsrMat>(gko_mat), sunmtx_(SUNMatNewEmpty(sunctx))
{
  sunmtx_->content = this;

  sunmtx_->ops->getid       = SUNMatGetID_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->clone       = SUNMatClone_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->zero        = SUNMatZero_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->copy        = SUNMatCopy_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->scaleadd    = SUNMatScaleAdd_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->scaleaddi   = SUNMatScaleAddI_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->matvecsetup = SUNMatMatvecSetup_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->matvec      = SUNMatMatvec_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->destroy     = SUNMatDestroy_Ginkgo<GkoCsrMat>;
  sunmtx_->ops->space       = SUNMatSpace_Ginkgo<GkoCsrMat>;
}

//
// Non-class methods
//

inline std::unique_ptr<GkoVecType> WrapVector(
    std::shared_ptr<const gko::Executor> gko_exec, N_Vector x)
{
  sunrealtype* x_arr       = (x->ops->nvgetdevicearraypointer)
                                 ? N_VGetDeviceArrayPointer(x)
                                 : N_VGetArrayPointer(x);
  const sunindextype x_len = N_VGetLength(x);
  return GkoVecType::create(gko_exec, gko::dim<2>(x_len, 1),
                            gko::Array<sunrealtype>::view(gko_exec, x_len, x_arr),
                            1);
}

inline std::unique_ptr<const GkoVecType> WrapConstVector(
    std::shared_ptr<const gko::Executor> gko_exec, N_Vector x)
{
  sunrealtype* x_arr       = (x->ops->nvgetdevicearraypointer)
                                 ? N_VGetDeviceArrayPointer(x)
                                 : N_VGetArrayPointer(x);
  const sunindextype x_len = N_VGetLength(x);
  return GkoVecType::create_const(gko_exec, gko::dim<2>(x_len, 1),
                                  gko::Array<sunrealtype>::const_view(gko_exec,
                                                                      x_len,
                                                                      x_arr),
                                  1);
}

template<typename GkoMatType>
void Print(Matrix<GkoMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.gkomtx().get());
}

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, GkoVecType* x, GkoVecType* y)
{
  A.gkomtx()->apply(x, y);
}

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y)
  {
    auto x_vec = WrapConstVector(A.gkoexec(), x);
    auto y_vec = WrapVector(A.gkoexec(), y);

    // y = Ax
    A.gkomtx()->apply(gko::lend(x_vec), gko::lend(y_vec));
  }
  else
  {
    auto x_vec = WrapVector(A.gkoexec(), x);

    // x = Ax
    A.gkomtx()->apply(gko::lend(x_vec), gko::lend(x_vec));
  }
}

template<typename GkoMatType>
std::shared_ptr<GkoMatType> CreateIdentity(
    std::shared_ptr<const gko::Executor> gko_exec, const gko::dim<2>& gko_dim)
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
inline void ScaleAdd(const sunrealtype c, Matrix<GkoDenseMat>& A,
                     Matrix<GkoDenseMat>& B)
{
  const auto I    = CreateIdentity<GkoDenseMat>(A.gkoexec(), A.gkodim());
  const auto one  = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  const auto cmat = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  // A = B + cA
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<>
inline void ScaleAdd(const sunrealtype c,
                     Matrix<gko::matrix::Csr<sunrealtype>>& A,
                     Matrix<gko::matrix::Csr<sunrealtype>>& B)
{
  const auto I =
      gko::matrix::Identity<sunrealtype>::create(A.gkoexec(), A.gkodim());
  const auto one  = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  const auto cmat = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  // A = B + cA
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<typename GkoMatType>
void ScaleAddI(const sunrealtype c, Matrix<GkoMatType>& A)
{
  const auto one  = gko::initialize<GkoDenseMat>({1.0}, A.gkoexec());
  const auto cmat = gko::initialize<GkoDenseMat>({c}, A.gkoexec());
  // A = 1*I + c*A = cA + I
  A.gkomtx()->add_scaled_identity(one.get(), cmat.get());
}

template<typename GkoMatType> void Zero(Matrix<GkoMatType>& A)
{
  A.gkomtx()->scale(gko::initialize<GkoDenseMat>({0.0}, A.gkoexec()).get());
}

template<> inline void Zero(Matrix<GkoDenseMat>& A) { A.gkomtx()->fill(0.0); }

template<typename GkoMatType>
void Copy(Matrix<GkoMatType>& A, Matrix<GkoMatType>& B)
{
  B.gkomtx()->copy_from(A.gkomtx().get());
}

} // namespace ginkgo
} // namespace sundials

#endif