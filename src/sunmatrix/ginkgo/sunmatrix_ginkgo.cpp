
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

using namespace sundials::ginkgo;
using GkoMatType = gko::matrix::Dense<realtype>;

#define GET_CONTENT(A) ((GinkgoMatrix<GkoMatType>*) A->content)


//
// Implementation specific methods
//

SUNMatrix SUNMatrix_Ginkgo(std::shared_ptr<gko::Executor> gko_exec, sunindextype M, sunindextype N, SUNContext sunctx)
{
  auto mat = new GinkgoMatrix<GkoMatType>(gko_exec, M, N, sunctx);
  return ((SUNMatrix) *mat);
}

int SUNMatFill_Ginkgo(realtype c, SUNMatrix A)
{
  Fill(*GET_CONTENT(A), c);
  return SUNMAT_SUCCESS;
}

int SUNMatPrint_Ginkgo(SUNMatrix A)
{
  Print(*GET_CONTENT(A));
  return SUNMAT_SUCCESS;
}

//
// SUNMatrix overrides
//

SUNMatrix_ID SUNMatGetID_Ginkgo(SUNMatrix A)
{
  return SUNMATRIX_GINKGO;
}

SUNMatrix SUNMatClone_Ginkgo(SUNMatrix A)
{
  return NULL;
}

void SUNMatDestroy_Ginkgo(SUNMatrix A)
{
  // We just delete mat because A is actually just stored in mat.
  GinkgoMatrix<GkoMatType>* mat = GET_CONTENT(A);
  delete mat;
}

int SUNMatZero_Ginkgo(SUNMatrix A)
{
  Zero(*GET_CONTENT(A));
  return SUNMAT_SUCCESS;
}

int SUNMatCopy_Ginkgo(SUNMatrix A, SUNMatrix B)
{
  Copy(*GET_CONTENT(A), *GET_CONTENT(B));
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_Ginkgo(realtype c, SUNMatrix A, SUNMatrix B)
{
  GinkgoMatrix<GkoMatType>* Amat = GET_CONTENT(A);
  GinkgoMatrix<GkoMatType>* Bmat = GET_CONTENT(B);
  Amat->ScaleAdd(c, *Bmat);
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_Ginkgo(realtype c, SUNMatrix A)
{
  GinkgoMatrix<GkoMatType>* Amat = GET_CONTENT(A);
  Amat->ScaleAddI(c);
  return SUNMAT_SUCCESS;
}

int SUNMatMatvecSetup_Ginkgo(SUNMatrix A)
{
  return SUNMAT_SUCCESS;
}

int SUNMatMatvec_Ginkgo(SUNMatrix A, N_Vector x, N_Vector y)
{
  Matvec(*GET_CONTENT(A), x, y);
  return SUNMAT_SUCCESS;
}

int SUNMatSpace_Ginkgo(SUNMatrix A, long int *lenrw, long int *leniw)
{
  GinkgoMatrix<GkoMatType>* Amat = GET_CONTENT(A);
  *lenrw = Amat->WorkspaceSize();
  *leniw = 0;
  return SUNMAT_SUCCESS;
}
