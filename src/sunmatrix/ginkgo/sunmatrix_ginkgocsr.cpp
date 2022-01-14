
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>

using namespace sundials;
using GkoMatType = gko::matrix::Csr<sunrealtype>;
using GkoBatchMatType = gko::matrix::BatchCsr<sunrealtype>;

#define GET_CONTENT(A) ((ginkgo::BaseMatrix*) A->content)

//
// Implementation specific methods
//

SUNMatrix SUNMatrix_GinkgoCsr(sunindextype M, sunindextype N, sunindextype NNZ, std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
{
  auto mat = new ginkgo::Matrix<GkoMatType>(M, N, NNZ, gko_exec, sunctx);
  return ((SUNMatrix) *mat);
}

SUNMatrix SUNMatrix_GinkgoCsrBlock(sunindextype nblocks, sunindextype M, sunindextype N,
                                   sunindextype NNZ, std::shared_ptr<gko::Executor> gko_exec,
                                   SUNContext sunctx)
{
  auto mat = new ginkgo::BlockMatrix<GkoBatchMatType>(nblocks, M, N, NNZ, gko_exec, sunctx);
  return ((SUNMatrix) *mat);
}

int SUNMatPrint_GinkgoCsr(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    ginkgo::Print(*Amat);
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    // ginkgo::Print(*Ablock);
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
}

//
// SUNMatrix overrides
//

SUNMatrix_ID SUNMatGetID_GinkgoCsr(SUNMatrix A)
{
  return SUNMATRIX_GINKGO;
}

SUNMatrix SUNMatClone_GinkgoCsr(SUNMatrix A)
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

void SUNMatDestroy_GinkgoCsr(SUNMatrix A)
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

int SUNMatZero_GinkgoCsr(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    Zero(*Amat);
    return SUNMAT_SUCCESS;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    Zero(*Ablock);
    return SUNMAT_SUCCESS;
  }
  else
  {
    return SUNMAT_OPERATION_FAIL;
  }
}

int SUNMatCopy_GinkgoCsr(SUNMatrix A, SUNMatrix B)
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

int SUNMatScaleAdd_GinkgoCsr(sunrealtype c, SUNMatrix A, SUNMatrix B)
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

int SUNMatScaleAddI_GinkgoCsr(sunrealtype c, SUNMatrix A)
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

int SUNMatMatvecSetup_GinkgoCsr(SUNMatrix A)
{
  return SUNMAT_SUCCESS;
}

int SUNMatMatvec_GinkgoCsr(SUNMatrix A, N_Vector x, N_Vector y)
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

int SUNMatSpace_GinkgoCsr(SUNMatrix A, long int *lenrw, long int *leniw)
{
  ginkgo::BaseMatrix* Amat = GET_CONTENT(A);
  *lenrw = Amat->workspaceSize();
  *leniw = 0;
  return SUNMAT_SUCCESS;
}
