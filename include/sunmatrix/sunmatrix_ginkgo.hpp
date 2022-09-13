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

#ifndef _SUNMATRIX_GINKGO_HPP
#define _SUNMATRIX_GINKGO_HPP

#include <memory>
#include <utility>
#include <ginkgo/ginkgo.hpp>
#include <sundials/core/sundials_base.hpp>
#include <sundials/sundials_matrix.hpp>

namespace sundials {
namespace ginkgo {

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;
using GkoCsrMat   = gko::matrix::Csr<sunrealtype, sunindextype>;
using GkoVecType  = GkoDenseMat;

// Forward decalaration of regular Matrix class
template<typename GkoMatType>
class Matrix;

//
// Methods that operate on SUNMatrix
//
template<typename GkoMatType>
SUNMatrix_ID SUNMatGetID_Ginkgo(SUNMatrix A)
{
  return SUNMATRIX_GINKGO;
}

template<typename GkoMatType>
SUNMatrix SUNMatClone_Ginkgo(SUNMatrix A)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  auto new_mat{new Matrix<GkoMatType>(*Amat)}; // NOLINT
  return new_mat->get();
}

template<typename GkoMatType>
void SUNMatDestroy_Ginkgo(SUNMatrix A)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  delete Amat; // NOLINT
  return;
}

template<typename GkoMatType>
int SUNMatZero_Ginkgo(SUNMatrix A)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  Zero(*Amat);
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatCopy_Ginkgo(SUNMatrix A, SUNMatrix B)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  Copy(*Amat, *static_cast<ginkgo::Matrix<GkoMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatScaleAdd_Ginkgo(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  ScaleAdd(c, *Amat, *static_cast<ginkgo::Matrix<GkoMatType>*>(B->content));
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatScaleAddI_Ginkgo(sunrealtype c, SUNMatrix A)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  ScaleAddI(c, *Amat);
  return SUNMAT_SUCCESS;
}

template<typename GkoMatType>
int SUNMatMatvec_Ginkgo(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto Amat{static_cast<Matrix<GkoMatType>*>(A->content)};
  Matvec(*Amat, x, y);
  return SUNMAT_SUCCESS;
}

namespace impl {

// //
// // Common base matrix class which makes RAII work.
// //
// template<class GkoMatType>
// class BaseMatrix : public BaseMatrix, public ConvertibleTo<SUNMatrix>
// {
// public:
//   BaseMatrix() : gkomtx_(nullptr), BaseMatrix() {}

//   explicit BaseMatrix(SUNContext sunctx) : gkomtx_(nullptr), BaseMatrix() {}

//   BaseMatrix(std::shared_ptr<GkoMatType> gko_mat, SUNContext sunctx)
//       : gkomtx_(std::move(gko_mat)), BaseMatrix(sunctx)
//   {}

//   // Move constructor
//   BaseMatrix(BaseMatrix&& that_matrix) noexcept : gkomtx_(std::move(that_matrix.gkomtx_)),
//   BaseMatrix(std::forward<BaseMatrix>(that_matrix)) {}

//   // Copy constructor clones the gko::matrix and SUNMatrix
//   BaseMatrix(const BaseMatrix& that_matrix) : gkomtx_(gko::clone(that_matrix.gkomtx_)), BaseMatrix(that_matrix) {}

//   // Move assignment
//   BaseMatrix& operator=(BaseMatrix&& rhs) noexcept
//   {
//     gkomtx_ = std::move(rhs.gkomtx_);
//     BaseMatrix::operator=(std::forward<BaseMatrix>(rhs));
//     return *this;
//   }

//   // Copy assignment clones the gko::matrix and SUNMatrix
//   BaseMatrix& operator=(const BaseMatrix& rhs)
//   {
//     gkomtx_ = gko::clone(rhs.gkomtx_);
//     BaseMatrix::operator=(rhs);
//     return *this;
//   }

//   // We have a pure virtual destructor to make this an asbtract class
//   virtual ~BaseMatrix() = 0;

//   // Override the ConvertibleTo methods
//   operator SUNMatrix() override { return object_.get(); }
//   operator SUNMatrix() const override { return object_.get(); }
//   SUNMatrix get() override { return object_.get(); }
//   SUNMatrix get() const override { return object_.get(); }

//   // Getters
//   std::shared_ptr<GkoMatType> gkomtx() const { return gkomtx_; }
//   std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx()->get_executor(); }
//   using BaseMatrix::sunctx;

// protected:
//   // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
//   std::shared_ptr<GkoMatType> gkomtx_;
// };

// // Pure virtual destructor requires implementation
// template<class GkoMatType>
// BaseMatrix<GkoMatType>::~BaseMatrix() = default;

} // namespace impl

//
// Standard matrix class
//
template<typename GkoMatType>
class Matrix : public sundials::impl::BaseMatrix, public sundials::ConvertibleTo<SUNMatrix>
{
public:
  // Default constructor means the matrix must be copied or moved to
  Matrix() = default;

  // We do not have implementations of these two constructors for general GkoMatType
  Matrix(sunindextype num_rows, sunindextype num_cols, std::shared_ptr<const gko::Executor> gko_exec,
         SUNContext sunctx);
  Matrix(sunindextype num_rows, sunindextype num_cols, sunindextype num_nonzeros,
         std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx);

  Matrix(std::shared_ptr<GkoMatType> gko_mat, SUNContext sunctx)
      : gkomtx_(gko_mat), sundials::impl::BaseMatrix(sunctx)
  {
    initSUNMatrix();
  }

  // Move constructor
  Matrix(Matrix&& that_matrix) noexcept
      : gkomtx_(std::move(that_matrix.gkomtx_)), sundials::impl::BaseMatrix(std::forward<Matrix>(that_matrix))
  {}

  // Copy constructor clones the gko::matrix and SUNMatrix
  Matrix(const Matrix& that_matrix)
      : gkomtx_(gko::clone(that_matrix.gkomtx_)), sundials::impl::BaseMatrix(that_matrix)
  {}

  // Move assignment
  Matrix& operator=(Matrix&& rhs) noexcept
  {
    gkomtx_ = std::move(rhs.gkomtx_);
    sundials::impl::BaseMatrix::operator=(std::forward<Matrix>(rhs));
    return *this;
  }

  // Copy assignment clones the gko::matrix and SUNMatrix
  Matrix& operator=(const Matrix& rhs)
  {
    gkomtx_ = gko::clone(rhs.gkomtx_);
    sundials::impl::BaseMatrix::operator=(rhs);
    return *this;
  }

  // Default destructor is fine since all members are RAII
  virtual ~Matrix() = default;

  // Getters
  std::shared_ptr<GkoMatType> gkomtx() const { return gkomtx_; }
  std::shared_ptr<const gko::Executor> gkoexec() const { return gkomtx()->get_executor(); }
  const gko::dim<2>& gkoSize() const { return gkomtx()->get_size(); }
  sunindextype gkodim(sunindextype dim) const { return gkomtx()->get_size()[dim]; }
  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods
  operator SUNMatrix() override { return object_.get(); }
  operator SUNMatrix() const override { return object_.get(); }
  SUNMatrix get() override { return object_.get(); }
  SUNMatrix get() const override { return object_.get(); }

private:
  std::shared_ptr<GkoMatType> gkomtx_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid     = SUNMatGetID_Ginkgo<GkoMatType>;
    this->object_->ops->clone     = SUNMatClone_Ginkgo<GkoMatType>;
    this->object_->ops->zero      = SUNMatZero_Ginkgo<GkoMatType>;
    this->object_->ops->copy      = SUNMatCopy_Ginkgo<GkoMatType>;
    this->object_->ops->scaleadd  = SUNMatScaleAdd_Ginkgo<GkoMatType>;
    this->object_->ops->scaleaddi = SUNMatScaleAddI_Ginkgo<GkoMatType>;
    this->object_->ops->matvec    = SUNMatMatvec_Ginkgo<GkoMatType>;
    this->object_->ops->destroy   = SUNMatDestroy_Ginkgo<GkoMatType>;
  }
};

//
// Specialized constructors
//

template<>
inline Matrix<GkoDenseMat>::Matrix(sunindextype num_rows, sunindextype num_cols,
                                   std::shared_ptr<const gko::Executor> gko_exec,
                                   SUNContext sunctx)
    : gkomtx_(GkoDenseMat::create(gko_exec, gko::dim<2>(num_rows, num_cols))), sundials::impl::BaseMatrix(sunctx)
{
  initSUNMatrix();
}

template<>
inline Matrix<GkoCsrMat>::Matrix(sunindextype num_rows, sunindextype num_cols, sunindextype num_nonzeros,
                                 std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
    : gkomtx_(GkoCsrMat::create(gko_exec, gko::dim<2>(num_rows, num_cols), num_nonzeros)),
      sundials::impl::BaseMatrix(sunctx)
{
  initSUNMatrix();
}

//
// Non-class methods
//

inline std::unique_ptr<GkoVecType> WrapVector(std::shared_ptr<const gko::Executor> gko_exec, N_Vector x)
{
  sunrealtype* x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x)};
  const sunindextype x_len{N_VGetLength(x)};
  return GkoVecType::create(gko_exec, gko::dim<2>(x_len, 1), gko::Array<sunrealtype>::view(gko_exec, x_len, x_arr), 1);
}

inline std::unique_ptr<const GkoVecType> WrapConstVector(std::shared_ptr<const gko::Executor> gko_exec, N_Vector x)
{
  sunrealtype* x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x)};
  const sunindextype x_len{N_VGetLength(x)};
  return GkoVecType::create_const(gko_exec, gko::dim<2>(x_len, 1),
                                  gko::Array<sunrealtype>::const_view(gko_exec, x_len, x_arr), 1);
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
  if (x != y) {
    auto x_vec{WrapConstVector(A.gkoexec(), x)};
    auto y_vec{WrapVector(A.gkoexec(), y)};

    // y = Ax
    A.gkomtx()->apply(x_vec.get(), y_vec.get());
  }
  else {
    auto x_vec{WrapVector(A.gkoexec(), x)};

    // x = Ax
    A.gkomtx()->apply(x_vec.get(), x_vec.get());
  }
}

template<typename GkoMatType>
void ScaleAdd(const sunrealtype c, Matrix<GkoMatType>& A, Matrix<GkoMatType>& B)
{
  const auto I{gko::matrix::Identity<sunrealtype>::create(A.gkoexec(), A.gkoSize())};
  const auto one{gko::initialize<GkoDenseMat>({1.0}, A.gkoexec())};
  const auto cmat{gko::initialize<GkoDenseMat>({c}, A.gkoexec())};
  // A = B + cA
  B.gkomtx()->apply(one.get(), I.get(), cmat.get(), A.gkomtx().get());
}

template<typename GkoMatType>
void ScaleAddI(const sunrealtype c, Matrix<GkoMatType>& A)
{
  const auto one{gko::initialize<GkoDenseMat>({1.0}, A.gkoexec())};
  const auto cmat{gko::initialize<GkoDenseMat>({c}, A.gkoexec())};
  // A = 1*I + c*A = cA + I
  A.gkomtx()->add_scaled_identity(one.get(), cmat.get());
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

} // namespace ginkgo
} // namespace sundials

#endif