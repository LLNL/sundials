
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

SUNMatrix SUNMatrix_GinkgoDense(sunindextype M, sunindextype N,
                                std::shared_ptr<gko::Executor> gko_exec, SUNContext sunctx)
{
  auto mat = new ginkgo::Matrix<GkoMatType>(M, N, gko_exec, sunctx);
  return ((SUNMatrix) *mat);
}

SUNMatrix SUNMatrix_GinkgoDenseBlock(sunindextype nblocks, sunindextype M,
                                     sunindextype N, std::shared_ptr<gko::Executor> gko_exec,
                                     SUNContext sunctx)
{
  auto mat = new ginkgo::BlockMatrix<GkoBatchMatType>(nblocks, M, N, gko_exec, sunctx);
  return ((SUNMatrix) *mat);
}

int SUNMatPrint_GinkgoDense(SUNMatrix A)
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

GkoMatType* SUNMatrix_GinkgoDense_GetGkoMat(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    return Amat->gkomtx().get();
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    return NULL;
  }
  else
  {
    return NULL;
  }
}

GkoBatchMatType* SUNMatrix_GinkgoDense_GetGkoBatchMat(SUNMatrix A)
{
  if (auto Amat = dynamic_cast<ginkgo::Matrix<GkoMatType>*>(GET_CONTENT(A)))
  {
    return NULL;
  }
  else if (auto Ablock = dynamic_cast<ginkgo::BlockMatrix<GkoBatchMatType>*>(GET_CONTENT(A)))
  {
    return Ablock->gkomtx().get();
  }
  else
  {
    return NULL;
  }
}

// sunindextype SUNMatrix_GinkgoDense_LData(SUNMatrix A)
// {

// }

// sunindextype SUNMatrix_GinkgoDense_Rows(SUNMatrix A)
// {

// }

// sunindextype SUNMatrix_GinkgoDense_Columns(SUNMatrix A)
// {

// }

// sunindextype SUNMatrix_GinkgoDense_BlockRows(SUNMatrix A)
// {

// }

// sunindextype SUNMatrix_GinkgoDense_BlockColumns(SUNMatrix A)
// {

// }

// sunindextype SUNMatrix_GinkgoDense_BlockLData(SUNMatrix A)
// {

// }

// sunindextype SUNMatrix_GinkgoDense_NumBlocks(SUNMatrix A)
// {

// }

//
// SUNMatrix overrides
//

SUNMatrix_ID SUNMatGetID_GinkgoDense(SUNMatrix A)
{
  return SUNMATRIX_GINKGODENSE;
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
