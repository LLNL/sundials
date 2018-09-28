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
 * -----------------------------------------------------------------
 * This is the implementation file for a MPI+RAJA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/raja/Vector.hpp>
#include <sundials/sundials_mpi.h>
#include <RAJA/RAJA.hpp>


#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

// RAJA defines
#define CUDA_BLOCK_SIZE 256
#define RAJA_NODE_TYPE RAJA::cuda_exec< CUDA_BLOCK_SIZE >
#define RAJA_REDUCE_TYPE RAJA::cuda_reduce< CUDA_BLOCK_SIZE >
#define RAJA_LAMBDA [=] __device__

extern "C" {

using namespace sunrajavec;

// Type defines
typedef sunrajavec::Vector<realtype, sunindextype> vector_type;

// Static constants
static constexpr sunindextype zeroIdx = 0;

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

  /* standard vector operations */
  ops->nvlinearsum    = N_VLinearSum_Raja;
  ops->nvconst        = N_VConst_Raja;
  ops->nvprod         = N_VProd_Raja;
  ops->nvdiv          = N_VDiv_Raja;
  ops->nvscale        = N_VScale_Raja;
  ops->nvabs          = N_VAbs_Raja;
  ops->nvinv          = N_VInv_Raja;
  ops->nvaddconst     = N_VAddConst_Raja;
  ops->nvdotprod      = N_VDotProd_Raja;
  ops->nvmaxnorm      = N_VMaxNorm_Raja;
  ops->nvwrmsnormmask = N_VWrmsNormMask_Raja;
  ops->nvwrmsnorm     = N_VWrmsNorm_Raja;
  ops->nvmin          = N_VMin_Raja;
  ops->nvwl2norm      = N_VWL2Norm_Raja;
  ops->nvl1norm       = N_VL1Norm_Raja;
  ops->nvcompare      = N_VCompare_Raja;
  ops->nvinvtest      = N_VInvTest_Raja;
  ops->nvconstrmask   = N_VConstrMask_Raja;
  ops->nvminquotient  = N_VMinQuotient_Raja;

  /* fused vector operations */
  ops->nvlinearcombination = N_VLinearCombination_Raja;
  ops->nvscaleaddmulti     = N_VScaleAddMulti_Raja;
  ops->nvdotprodmulti      = NULL;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Raja;
  ops->nvscalevectorarray             = N_VScaleVectorArray_Raja;
  ops->nvconstvectorarray             = N_VConstVectorArray_Raja;
  ops->nvwrmsnormvectorarray          = NULL;
  ops->nvwrmsnormmaskvectorarray      = NULL;
  ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Raja;
  ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;

  /* Attach ops and set content to NULL */
  v->content = NULL;
  v->ops     = ops;

  return(v);
}


#if SUNDIALS_MPI_ENABLED
N_Vector N_VNew_Raja(MPI_Comm comm,
                     sunindextype local_length,
                     sunindextype global_length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja(local_length);
  if (v == NULL) return(NULL);

  v->content = new vector_type(comm, local_length, global_length);

  return(v);
}
#else
N_Vector N_VNew_Raja(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja(length);
  if (v == NULL) return(NULL);

  v->content = new vector_type(SUNMPI_COMM_WORLD, length, length);

  return(v);
}
#endif


N_Vector N_VMake_Raja(N_VectorContent_Raja c)
{
  N_Vector v;
  vector_type* x = static_cast<vector_type*>(c);
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
  vector_type* xd = static_cast<vector_type*>(v->content);
  return xd->size();
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->host());
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->device());
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyToDev();
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
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

  /* standard vector operations */
  ops->nvlinearsum    = w->ops->nvlinearsum;
  ops->nvconst        = w->ops->nvconst;
  ops->nvprod         = w->ops->nvprod;
  ops->nvdiv          = w->ops->nvdiv;
  ops->nvscale        = w->ops->nvscale;
  ops->nvabs          = w->ops->nvabs;
  ops->nvinv          = w->ops->nvinv;
  ops->nvaddconst     = w->ops->nvaddconst;
  ops->nvdotprod      = w->ops->nvdotprod;
  ops->nvmaxnorm      = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
  ops->nvmin          = w->ops->nvmin;
  ops->nvwl2norm      = w->ops->nvwl2norm;
  ops->nvl1norm       = w->ops->nvl1norm;
  ops->nvcompare      = w->ops->nvcompare;
  ops->nvinvtest      = w->ops->nvinvtest;
  ops->nvconstrmask   = w->ops->nvconstrmask;
  ops->nvminquotient  = w->ops->nvminquotient;

  /* fused vector operations */
  ops->nvlinearcombination = w->ops->nvlinearcombination;
  ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
  ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
  ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
  ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
  ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
  ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
  ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
  ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

  /* Create content */
  v->content = NULL;
  v->ops  = ops;

  return(v);
}

N_Vector N_VClone_Raja(N_Vector w)
{
  N_Vector v;
  vector_type* wdat = static_cast<vector_type*>(w->content);
  vector_type* vdat = new vector_type(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Raja(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Raja(N_Vector v)
{
  vector_type* x = static_cast<vector_type*>(v->content);
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
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  int npes;

  SUNMPI_Comm_size(comm, &npes);

  *lrw = getGlobalSize<realtype, sunindextype>(X);
  *liw = 2*npes;
}

void N_VConst_Raja(realtype c, N_Vector Z)
{
  const sunindextype N = getSize<realtype, sunindextype>(Z);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N), RAJA_LAMBDA(sunindextype i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = a*xdata[i] + b*ydata[i];
    }
  );
}

void N_VProd_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] * ydata[i];
    }
  );
}

void N_VDiv_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] / ydata[i];
    }
  );
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = c * xdata[i];
    }
  );
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]);
    }
  );
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = ONE / xdata[i];
    }
  );
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] + b;
    }
  );
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *ydata = getDevData<realtype, sunindextype>(Y);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += xdata[i] * ydata[i] ;
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return gsum;
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceMax< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.max(abs(xdata[i]));
    }
  );

  /* Reduce across MPI processes */
  realtype maximum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return SUNMPI_Allreduce_scalar(maximum, 2, comm);
}

realtype N_VWrmsNorm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *wdata = getDevData<realtype, sunindextype>(W);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  const sunindextype Nglobal = getGlobalSize<realtype, sunindextype>(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return std::sqrt(SUNMPI_Allreduce_scalar(sum, 1, comm)/Nglobal);
}

realtype N_VWrmsNormMask_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *wdata = getDevData<realtype, sunindextype>(W);
  const realtype *iddata = getDevData<realtype, sunindextype>(ID);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  const sunindextype Nglobal = getGlobalSize<realtype, sunindextype>(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (iddata[i] > ZERO)
        gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return std::sqrt(SUNMPI_Allreduce_scalar(sum, 1, comm)/Nglobal);
}

realtype N_VMin_Raja(N_Vector X)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.min(xdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype minumum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return SUNMPI_Allreduce_scalar(minumum, 3, comm);
}

realtype N_VWL2Norm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const realtype *wdata = getDevData<realtype, sunindextype>(W);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return std::sqrt(SUNMPI_Allreduce_scalar(sum, 1, comm));
}

realtype N_VL1Norm_Raja(N_Vector X)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (abs(xdata[i]));
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  return SUNMPI_Allreduce_scalar(sum, 1, comm);
}

void N_VCompare_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(X);
  const sunindextype N = getSize<realtype, sunindextype>(X);
  realtype *zdata = getDevData<realtype, sunindextype>(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
    }
  );
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = getDevData<realtype, sunindextype>(x);
  const sunindextype N = getSize<realtype, sunindextype>(x);
  realtype *zdata = getDevData<realtype, sunindextype>(z);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (xdata[i] == ZERO) {
        gpu_result += ONE;
      } else {
        zdata[i] = ONE/xdata[i];
      }
    }
  );

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(x);
  realtype global_minimum = SUNMPI_Allreduce_scalar(minimum, 3, comm);

  return (global_minimum < HALF);
}

booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m)
{
  const realtype *cdata = getDevData<realtype, sunindextype>(c);
  const realtype *xdata = getDevData<realtype, sunindextype>(x);
  const sunindextype N = getSize<realtype, sunindextype>(x);
  realtype *mdata = getDevData<realtype, sunindextype>(m);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      bool test = (abs(cdata[i]) > ONEPT5 && cdata[i]*xdata[i] <= ZERO) ||
                  (abs(cdata[i]) > HALF   && cdata[i]*xdata[i] <  ZERO);
      mdata[i] = test ? ONE : ZERO;
      gpu_result += mdata[i];
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(x);
  realtype global_sum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return (global_sum < HALF);
}

realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom)
{
  const realtype *ndata = getDevData<realtype, sunindextype>(num);
  const realtype *ddata = getDevData<realtype, sunindextype>(denom);
  const sunindextype N = getSize<realtype, sunindextype>(num);

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (ddata[i] != ZERO)
        gpu_result.min(ndata[i]/ddata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(num);
  return SUNMPI_Allreduce_scalar(minimum, 3, comm);
}


/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearCombination_Raja(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  cudaError_t  err;

  sunindextype N = getSize<realtype, sunindextype>(z);
  realtype* d_zd = getDevData<realtype, sunindextype>(z);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = getDevData<realtype, sunindextype>(X[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      d_zd[i] = d_c[0] * d_Xd[0][i];
      for (int j=1; j<nvec; j++)
        d_zd[i] += d_c[j] * d_Xd[j][i];
    }
  );

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VScaleAddMulti_Raja(int nvec, realtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getSize<realtype, sunindextype>(x);
  realtype* d_xd = getDevData<realtype, sunindextype>(x);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Yd[j] = getDevData<realtype, sunindextype>(Y[j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = getDevData<realtype, sunindextype>(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = d_c[j] * d_xd[i] + d_Yd[j][i];
    }
  );

  // Free host array
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Raja(int nvec,
                                 realtype a, N_Vector* X,
                                 realtype b, N_Vector* Y,
                                 N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getSize<realtype, sunindextype>(Z[0]);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = getDevData<realtype, sunindextype>(X[j]);

  realtype** h_Yd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Yd[j] = getDevData<realtype, sunindextype>(Y[j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = getDevData<realtype, sunindextype>(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = a * d_Xd[j][i] + b * d_Yd[j][i];
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VScaleVectorArray_Raja(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getSize<realtype, sunindextype>(Z[0]);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = getDevData<realtype, sunindextype>(X[j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = getDevData<realtype, sunindextype>(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = d_c[j] * d_Xd[j][i];
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VConstVectorArray_Raja(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getSize<realtype, sunindextype>(Z[0]);

  // Create array of device pointers on host
  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = getDevData<realtype, sunindextype>(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = c;
    }
  );

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VScaleAddMultiVectorArray_Raja(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;

  sunindextype N = getSize<realtype, sunindextype>(X[0]);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = getDevData<realtype, sunindextype>(X[j]);

  realtype** h_Yd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Yd[j*nsum+k] = getDevData<realtype, sunindextype>(Y[k][j]);

  realtype** h_Zd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Zd[j*nsum+k] = getDevData<realtype, sunindextype>(Z[k][j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nsum*nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nsum*nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        for (int k=0; k<nsum; k++)
          d_Zd[j*nsum+k][i] = d_c[k] * d_Xd[j][i] + d_Yd[j*nsum+k][i];
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VLinearCombinationVectorArray_Raja(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getSize<realtype, sunindextype>(Z[0]);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Xd[j*nsum+k] = getDevData<realtype, sunindextype>(X[k][j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = getDevData<realtype, sunindextype>(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++) {
        d_Zd[j][i] = d_c[0] * d_Xd[j*nsum][i];
        for (int k=1; k<nsum; k++) {
          d_Zd[j][i] += d_c[k] * d_Xd[j*nsum+k][i];
        }
      }
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}

} // extern "C"
