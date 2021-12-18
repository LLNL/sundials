
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

using namespace sundials::ginkgo::matrix;
using GkoMatType = gko::matrix::Dense<sunrealtype>;

#define GET_CONTENT(A) ((GinkgoMatrix<GkoMatType>*) A->content)


//
// Implementation specific methods
//

SUNMatrix SUNMatrix_GinkgoDense(std::shared_ptr<gko::Executor> gko_exec,
                                sunindextype M, sunindextype N, SUNContext sunctx)
{
  auto mat = new GinkgoMatrix<GkoMatType>(gko_exec, M, N, sunctx);
  return ((SUNMatrix) *mat);
}

int SUNMatFill_GinkgoDense(sunrealtype c, SUNMatrix A)
{
  Fill(*GET_CONTENT(A), c);
  return SUNMAT_SUCCESS;
}

int SUNMatPrint_GinkgoDense(SUNMatrix A)
{
  Print(*GET_CONTENT(A));
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
  GinkgoMatrix<GkoMatType>* mat = GET_CONTENT(A);
  auto new_mat = new GinkgoMatrix<GkoMatType>(*mat);
  return ((SUNMatrix) *new_mat);
}

void SUNMatDestroy_GinkgoDense(SUNMatrix A)
{
  // We just delete mat because A is actually just stored in mat.
  GinkgoMatrix<GkoMatType>* mat = GET_CONTENT(A);
  delete mat;
}

int SUNMatZero_GinkgoDense(SUNMatrix A)
{
  Zero(*GET_CONTENT(A));
  return SUNMAT_SUCCESS;
}

int SUNMatCopy_GinkgoDense(SUNMatrix A, SUNMatrix B)
{
  Copy(*GET_CONTENT(A), *GET_CONTENT(B));
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_GinkgoDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  GinkgoMatrix<GkoMatType>* Amat = GET_CONTENT(A);
  GinkgoMatrix<GkoMatType>* Bmat = GET_CONTENT(B);
  ScaleAdd(c, *Amat, *Bmat);
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_GinkgoDense(sunrealtype c, SUNMatrix A)
{
  GinkgoMatrix<GkoMatType>* Amat = GET_CONTENT(A);
  ScaleAddI(c, *Amat);
  return SUNMAT_SUCCESS;
}

int SUNMatMatvecSetup_GinkgoDense(SUNMatrix A)
{
  return SUNMAT_SUCCESS;
}

int SUNMatMatvec_GinkgoDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  Matvec(*GET_CONTENT(A), x, y);
  return SUNMAT_SUCCESS;
}

int SUNMatSpace_GinkgoDense(SUNMatrix A, long int *lenrw, long int *leniw)
{
  GinkgoMatrix<GkoMatType>* Amat = GET_CONTENT(A);
  *lenrw = Amat->WorkspaceSize();
  *leniw = 0;
  return SUNMAT_SUCCESS;
}
