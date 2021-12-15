
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>
#include "test_sunmatrix.h"

using namespace sundials::ginkgo;
using MtxType = GinkgoMatrix<gko::matrix::Dense<realtype>>;
using VecType = gko::matrix::Dense<realtype>;

int main(int argc, char* argv[])
{
  sundials::Context sunctx;

  auto gko_exec = gko::ReferenceExecutor::create();

  // Test the C++ interface
  {
    MtxType A(gko_exec, 2, 2, sunctx);
    MtxType B(gko_exec, 2, 2, sunctx);
    auto x = GkoVecType::create(gko_exec, gko::dim<2>(2, 1));
    x->fill(1.0);
    auto b = GkoVecType::create(gko_exec, gko::dim<2>(2, 1));
    b->fill(0.0);

    Zero(A);
    A.ScaleAddI(1.0);
    A.ScaleAddI(1.0);
    A.gkolinop()->apply(x.get(), b.get());

    std::cout << "x = \n";
    gko::write(std::cout, x.get());
    std::cout << "b = \n";
    gko::write(std::cout, b.get());

    Fill(B, 1.0);
    A.ScaleAdd(1.0, B);
    B.gkolinop()->apply(x.get(), b.get());

    std::cout << "x = \n";
    gko::write(std::cout, x.get());
    std::cout << "b = \n";
    gko::write(std::cout, b.get());

    Copy(A, B);

    std::cout << "A, B\n";
    gko::write(std::cout, A.gkomtx().get());
    gko::write(std::cout, B.gkomtx().get());
    std::cout << "\n";

    Matvec(A, x.get(), b.get());

    std::cout << "x = \n";
    gko::write(std::cout, x.get());
    std::cout << "b = \n";
    gko::write(std::cout, b.get());

    long int lenrw = A.WorkspaceSize();

    std::cout << "lenrw = " << lenrw << "\n";
  }

  std::cout << "----------------------------------------------------------\n";

  // Test the C interface
  {
    SUNMatrix A = SUNMatrix_Ginkgo(gko_exec, 2, 2, sunctx);
    SUNMatrix B = SUNMatrix_Ginkgo(gko_exec, 2, 2, sunctx);

    N_Vector x = N_VNew_Serial(2, sunctx);
    N_Vector b = N_VClone(x);

    N_VConst(1.0, x);
    N_VConst(0.0, b);

    SUNMatZero(A);
    SUNMatScaleAddI(1.0, A);
    SUNMatScaleAddI(1.0, A);
    SUNMatMatvec(A, x, b);

    std::cout << "x = \n";
    N_VPrint(x);
    std::cout << "b = \n";
    N_VPrint(b);

    SUNMatFill_Ginkgo(1.0, B);
    SUNMatScaleAdd(1.0, A, B);
    SUNMatMatvec(B, x, b);

    std::cout << "x = \n";
    N_VPrint(x);
    std::cout << "b = \n";
    N_VPrint(b);

    SUNMatCopy(A, B);

    std::cout << "A, B\n";
    SUNMatPrint_Ginkgo(A);
    SUNMatPrint_Ginkgo(B);
    std::cout << "\n";

    SUNMatMatvec(A, x, b);

    std::cout << "x = \n";
    N_VPrint(x);
    std::cout << "b = \n";
    N_VPrint(b);

    long int lenrw, leniw;
    SUNMatSpace(A, &lenrw, &leniw);

    std::cout << "lenrw = " << lenrw << "\n";

    SUNMatDestroy(A);
    SUNMatDestroy(B);
  }

  return 0;
}

int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  // int failure = 0;
  // sunindextype i = 0;
  // sunindextype Aldata = SUNMatrix_MagmaDense_LData(A);
  // sunindextype Bldata = SUNMatrix_MagmaDense_LData(B);
  // realtype *Adata = (realtype*) malloc(sizeof(realtype)*Aldata);
  // realtype *Bdata = (realtype*) malloc(sizeof(realtype)*Bldata);

  // /* copy data to host */
  // SUNMatrix_MagmaDense_CopyFromDevice(A, Adata);
  // SUNMatrix_MagmaDense_CopyFromDevice(B, Bdata);

  // /* check lengths */
  // if (Aldata != Bldata) {
  //   printf(">>> ERROR: check_matrix: Different data array lengths \n");
  //   return(1);
  // }

  // /* compare data */
  // for(i=0; i < Aldata; i++) {
  //   failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
  // }

  // free(Adata);
  // free(Bdata);

  // if (failure > ZERO)
  //   return(1);
  // else
  //   return(0);
  return 0;
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  // int failure = 0;
  // sunindextype i = 0;
  // sunindextype Aldata = SUNMatrix_MagmaDense_LData(A);
  // realtype *Adata = (realtype*) malloc(sizeof(realtype)*Aldata);

  // /* copy data to host */
  // SUNMatrix_MagmaDense_CopyFromDevice(A, Adata);

  // /* compare data */
  // for(i=0; i < Aldata; i++) {
  //   int check = SUNRCompareTol(Adata[i], val, tol);
  //   if (check) {
  //     printf("failed at %d\n", i);
  //     failure += check;
  //   }
  // }

  // free(Adata);

  // if (failure > ZERO)
  //   return(1);
  // else
  //   return(0);
  return 0;
}

int check_vector(N_Vector actual, N_Vector expected, realtype tol)
{
  // int failure = 0;
  // realtype *xdata, *ydata;
  // sunindextype xldata, yldata;
  // sunindextype i;

  // /* copy vectors to host */
  // HIP_OR_CUDA( N_VCopyFromDevice_Hip(actual);,
  //              N_VCopyFromDevice_Cuda(actual); )
  // HIP_OR_CUDA( N_VCopyFromDevice_Hip(expected);,
  //              N_VCopyFromDevice_Cuda(expected); )

  // /* get vector data */
  // xdata = N_VGetArrayPointer(actual);
  // ydata = N_VGetArrayPointer(expected);

  // /* check data lengths */
  // xldata = N_VGetLength(actual);
  // yldata = N_VGetLength(expected);


  // if (xldata != yldata) {
  //   printf(">>> ERROR: check_vector: Different data array lengths \n");
  //   return(1);
  // }

  // /* check vector data */
  // for(i=0; i < xldata; i++)
  //   failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  // if (failure > ZERO) {
  //   printf("Check_vector failures:\n");
  //   for(i=0; i < xldata; i++)
  //     if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
  //       printf("  actual[%ld] = %g != %e (err = %g)\n", (long int) i,
  //              xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  // }

  // if (failure > ZERO)
  //   return(1);
  // else
  //   return(0);
  return 0;
}

booleantype has_data(SUNMatrix A)
{
  // realtype *Adata = SUNMatrix_MagmaDense_Data(A);
  // if (Adata == NULL)
  //   return SUNFALSE;
  // else
  //   return SUNTRUE;
  return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  // if (SUNMatrix_MagmaDense_Rows(A) == SUNMatrix_MagmaDense_Columns(A))
  //   return SUNTRUE;
  // else
  //   return SUNFALSE;
  return SUNTRUE;
}

void sync_device(SUNMatrix A)
{
}