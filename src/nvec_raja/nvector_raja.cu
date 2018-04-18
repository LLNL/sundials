/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/raja/Vector.hpp>
#include <RAJA/RAJA.hpp>
#include <sundials/sundials_mpi_types.h>


#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

using namespace sunrajavec;

static constexpr sunindextype zeroIdx = 0;
static realtype VAllReduce_Raja(realtype d, int op, SUNDIALS_Comm comm);

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Raja(N_Vector v)
{
  return SUNDIALS_NVEC_RAJA;
}

N_Vector N_VNewEmpty_Raja(sunindextype length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Raja content;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Raja;
  ops->nvclone           = N_VClone_Raja;
  ops->nvcloneempty      = N_VCloneEmpty_Raja;
  ops->nvdestroy         = N_VDestroy_Raja;
  ops->nvspace           = N_VSpace_Raja;
  ops->nvgetarraypointer = NULL; //N_VGetArrayPointer_Raja;
  ops->nvsetarraypointer = NULL; //N_VSetArrayPointer_Raja;
  ops->nvlinearsum       = N_VLinearSum_Raja;
  ops->nvconst           = N_VConst_Raja;
  ops->nvprod            = N_VProd_Raja;
  ops->nvdiv             = N_VDiv_Raja;
  ops->nvscale           = N_VScale_Raja;
  ops->nvabs             = N_VAbs_Raja;
  ops->nvinv             = N_VInv_Raja;
  ops->nvaddconst        = N_VAddConst_Raja;
  ops->nvdotprod         = N_VDotProd_Raja;
  ops->nvmaxnorm         = N_VMaxNorm_Raja;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Raja;
  ops->nvwrmsnorm        = N_VWrmsNorm_Raja;
  ops->nvmin             = N_VMin_Raja;
  ops->nvwl2norm         = N_VWL2Norm_Raja;
  ops->nvl1norm          = N_VL1Norm_Raja;
  ops->nvcompare         = N_VCompare_Raja;
  ops->nvinvtest         = N_VInvTest_Raja;
  ops->nvconstrmask      = N_VConstrMask_Raja;
  ops->nvminquotient     = N_VMinQuotient_Raja;

  /* Attach ops and set content to NULL */
  v->content = NULL;
  v->ops     = ops;

  return(v);
}


N_Vector N_VNew_Raja(SUNDIALS_Comm comm,
                     sunindextype local_length,
                     sunindextype global_length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja(local_length);
  if (v == NULL) return(NULL);

  v->content = new Vector<realtype, sunindextype>(comm, local_length, global_length);

  return(v);
}


N_Vector N_VMake_Raja(N_VectorContent_Raja c)
{
  N_Vector v;
  Vector<realtype, sunindextype>* x = static_cast<Vector<realtype, sunindextype>*>(c);
  sunindextype length = x->size();

  v = NULL;
  v = N_VNewEmpty_Raja(length);
  if (v == NULL) return(NULL);

  v->content = c;

  return(v);
}


/* -----------------------------------------------------------------
 * Function to return the length of the vector.
 */
sunindextype N_VGetLength_Raja(N_Vector v)
{
  Vector<realtype, sunindextype>* xd = static_cast<Vector<realtype, sunindextype>*>(v->content);
  return xd->size();
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Raja(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  return (xv->host());
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Raja(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  return (xv->device());
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Raja(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  xv->copyToDev();
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Raja(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  xv->copyFromDev();
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Raja(N_Vector X)
{
  N_VPrintFile_Raja(X, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Raja(N_Vector X, FILE *outfile)
{
  const realtype *xd = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  sunindextype i;

  for (i = 0; i < N; ++i) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", xd[i]);
#else
    fprintf(outfile, "%11.8g\n", xd[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Raja(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;
  ops->nvprod            = w->ops->nvprod;
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale;
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */
  v->content = NULL;
  v->ops  = ops;

  return(v);
}

N_Vector N_VClone_Raja(N_Vector w)
{
  N_Vector v;
  Vector<realtype, sunindextype>* wdat = static_cast<Vector<realtype, sunindextype>*>(w->content);
  Vector<realtype, sunindextype>* vdat = new Vector<realtype, sunindextype>(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Raja(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Raja(N_Vector v)
{
  Vector<realtype, sunindextype>* x = static_cast<Vector<realtype, sunindextype>*>(v->content);
  if (x != NULL) {
    delete x;
    v->content = NULL;
  }

  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Raja(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
    *lrw = getSize<realtype, sunindextype>(X);
    *liw = 1;
}

void N_VConst_Raja(realtype c, N_Vector Z)
{
  const sunindextype N = getSize<realtype, sunindextype>(Z);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = a*xdata[i] + b*ydata[i];
  });
}

void N_VProd_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = xdata[i] * ydata[i];
  });
}

void N_VDiv_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = xdata[i] / ydata[i];
  });
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = c * xdata[i];
  });
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = abs(xdata[i]);
  });
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = RCONST(1.0) / xdata[i];
  });
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = xdata[i] + b;
  });
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result += xdata[i] * ydata[i] ;
  });

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  realtype gsum = VAllReduce_Raja(sum, 1, comm);
  return gsum;
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceMax<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result.max(abs(xdata[i]));
  });

  /* Reduce across MPI processes */
  realtype maximum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return VAllReduce_Raja(maximum, 2, comm);
}

realtype N_VWrmsNorm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *wdata = getDevData<realtype, sunindextype>(W);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  const sunindextype Nglobal = getGlobalSize<realtype, sunindextype>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
  });

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return std::sqrt(VAllReduce_Raja(sum, 1, comm)/Nglobal);
}

realtype N_VWrmsNormMask_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *wdata = getDevData<realtype, sunindextype>(W);
  const realtype *iddata = getDevData<realtype, sunindextype>(ID);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  const sunindextype Nglobal = getGlobalSize<realtype, sunindextype>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i] * iddata[i]);
  });

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return std::sqrt(VAllReduce_Raja(sum, 1, comm)/Nglobal);
}

realtype N_VMin_Raja(N_Vector X)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceMin<RAJA::cuda_reduce<128>, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result.min(xdata[i]);
  });

  /* Reduce across MPI processes */
  realtype minumum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return VAllReduce_Raja(minumum, 3, comm);
}

realtype N_VWL2Norm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *wdata = getDevData<realtype, sunindextype>(W);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
  });

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return std::sqrt(VAllReduce_Raja(sum, 1, comm));
}

realtype N_VL1Norm_Raja(N_Vector X)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    gpu_result += (abs(xdata[i]));
  });

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return VAllReduce_Raja(sum, 1, comm);
}

void N_VCompare_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(zeroIdx, N, [=] __device__(sunindextype i) {
     zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
  });
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(x);
  const sunindextype N = getSize<realtype, sunindextype>(x);
  realtype *zdata = getDevData<realtype, sunindextype>(z);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(ZERO);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    if (xdata[i] == ZERO) {
      gpu_result += ONE;
    } else {
      zdata[i] = ONE/xdata[i];
    }
  });

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(x);
  realtype global_minimum = VAllReduce_Raja(minimum, 3, comm);

  return (global_minimum < HALF);
}

booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m)
{
  const realtype *cdata = getDevData<realtype, sunindextype>(c);
  const realtype *xdata = getDevData<realtype, sunindextype>(x);
  const sunindextype N = getSize<realtype, sunindextype>(x);
  realtype *mdata = getDevData<realtype, sunindextype>(m);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(ZERO);
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    bool test = (abs(cdata[i]) > ONEPT5 && cdata[i]*xdata[i] <= ZERO) ||
                (abs(cdata[i]) > HALF   && cdata[i]*xdata[i] <  ZERO);
    mdata[i] = test ? ONE : ZERO;
    gpu_result += mdata[i];
  });

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(x);
  realtype global_minimum = VAllReduce_Raja(minimum, 3, comm);

  return (global_minimum < HALF);
}

realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom)
{
  const realtype *ndata = getDevData<realtype, sunindextype>(num);
  const realtype *ddata = getDevData<realtype, sunindextype>(denom);
  const sunindextype N = getSize<realtype, sunindextype>(num);

  RAJA::ReduceMin<RAJA::cuda_reduce<128>, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall<RAJA::cuda_exec<128> >(zeroIdx, N, [=] __device__(sunindextype i) {
    if (ddata[i] != ZERO)
      gpu_result.min(ndata[i]/ddata[i]);
  });

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNDIALS_Comm comm = getMPIComm<realtype, sunindextype>(num);
  return VAllReduce_Raja(minimum, 3, comm);
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static realtype VAllReduce_Raja(realtype d, int op, SUNDIALS_Comm comm)
{
  /*
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator
   */

#ifdef SUNDIALS_MPI_ENABLED

  realtype out;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);

#else

  /* If MPI is not enabled don't do reduction */
  return d;

#endif // ifdef SUNDIALS_MPI_ENABLED
}


} // extern "C"
