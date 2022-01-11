
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>

using namespace sundials;
using GkoMatType = gko::matrix::Dense<sunrealtype>;
using GkoBatchMatType = gko::matrix::BatchDense<sunrealtype>;

#define GET_CONTENT(A) ((ginkgo::BaseMatrix*) A->content)

//
// Implementation specific methods
//

SUNMatrix SUNMatrix_GinkgoDense(std::shared_ptr<gko::Executor> gko_exec,
                                sunindextype M, sunindextype N, SUNContext sunctx)
{
  auto mat = new ginkgo::Matrix<GkoMatType>(gko_exec, M, N, sunctx);
  return ((SUNMatrix) *mat);
}

SUNMatrix SUNMatrix_GinkgoDenseBlock(std::shared_ptr<gko::Executor> gko_exec,
                                     sunindextype nblocks, sunindextype M,
                                     sunindextype N, SUNContext sunctx)
{
  auto mat = new ginkgo::BlockMatrix<GkoBatchMatType>(gko_exec, nblocks, M, N, sunctx);
  return ((SUNMatrix) *mat);
}

int SUNMatFill_GinkgoDense(sunrealtype c, SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    ginkgo::Fill(*Amat, c);
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    ginkgo::Fill(*Ablock, c);
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
}

int SUNMatPrint_GinkgoDense(SUNMatrix A)
{
  // if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(A->content))
  // {
  //   ginkgo::Print(Amat);
  // }
  // else
  // {
  //   auto Amat = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(A->content);
  //   ginkgo::Print(Amat);
  // }
  return SUNMAT_SUCCESS;
}

//
// SUNMatrix overrides
//

SUNMatrix_ID SUNMatGetID_GinkgoDense(SUNMatrix A)
{
  return SUNMATRIX_GINKGO;
}

SUNMatrix SUNMatClone_GinkgoDense(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    auto new_mat = new ginkgo::Matrix<GkoMatType>(*Amat);
    return ((SUNMatrix) *new_mat);
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    auto new_mat = new ginkgo::BlockMatrix<GkoBatchMatType>(*Ablock);
    return ((SUNMatrix) *new_mat);
  }
  else
  {
    return NULL;
  }
}

void SUNMatDestroy_GinkgoDense(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    delete Amat;
    return;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    delete Ablock;
    return;
  }
  else
  {
    return;
  }
}

int SUNMatZero_GinkgoDense(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    Amat->Zero();
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    Ablock->Zero();
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
}

int SUNMatCopy_GinkgoDense(SUNMatrix A, SUNMatrix B)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    Copy(*Amat, *static_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(B)));
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    Copy(*Ablock, *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(B)));
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_GinkgoDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    ScaleAdd(c, *Amat, *static_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(B)));
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    ScaleAdd(c, *Ablock, *static_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(B)));
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_GinkgoDense(sunrealtype c, SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    ScaleAddI(c, *Amat);
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    ScaleAddI(c, *Ablock);
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
  return SUNMAT_SUCCESS;
}

int SUNMatMatvecSetup_GinkgoDense(SUNMatrix A)
{
  return SUNMAT_SUCCESS;
}

int SUNMatMatvec_GinkgoDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    Matvec(*Amat, x, y);
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    Matvec(*Ablock, x, y);
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
  return SUNMAT_SUCCESS;
}

int SUNMatSpace_GinkgoDense(SUNMatrix A, long int *lenrw, long int *leniw)
{
  ginkgo::BaseMatrix* Amat = GET_CONTENT(A);
  *lenrw = Amat->workspaceSize();
  *leniw = 0;
  return SUNMAT_SUCCESS;
}
