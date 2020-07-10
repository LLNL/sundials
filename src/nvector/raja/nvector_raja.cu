/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos, Daniel McGreer @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a RAJA+CUDA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_raja.h>
#include <RAJA/RAJA.hpp>

#include "sundials_debug.h"
#include "sundials_cuda.h"

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

// RAJA defines
#define CUDA_BLOCK_SIZE 256
#define RAJA_NODE_TYPE RAJA::cuda_exec< CUDA_BLOCK_SIZE >
#define RAJA_REDUCE_TYPE RAJA::cuda_reduce
#define RAJA_LAMBDA [=] __device__

extern "C" {

// Static constants
static constexpr sunindextype zeroIdx = 0;

// Helpful macros
#define NVEC_RAJA_CONTENT(x) ((N_VectorContent_Raja)(x->content))
#define NVEC_RAJA_MEMSIZE(x) (NVEC_RAJA_CONTENT(x)->length * sizeof(realtype))
#define NVEC_RAJA_PRIVATE(x) ((N_PrivateVectorContent_Raja)(NVEC_RAJA_CONTENT(x)->priv))

struct _N_PrivateVectorContent_Raja
{
  booleantype use_managed_mem; /* indicates if the data pointers and buffer pointers are managed memory */
};

typedef struct _N_PrivateVectorContent_Raja *N_PrivateVectorContent_Raja;

static int AllocateData(N_Vector v);

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Raja(N_Vector v)
{
  return SUNDIALS_NVEC_RAJA;
}

N_Vector N_VNewEmpty_Raja()
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Raja;
  v->ops->nvclone           = N_VClone_Raja;
  v->ops->nvcloneempty      = N_VCloneEmpty_Raja;
  v->ops->nvdestroy         = N_VDestroy_Raja;
  v->ops->nvspace           = N_VSpace_Raja;
  v->ops->nvgetlength       = N_VGetLength_Raja;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Raja;
  v->ops->nvconst        = N_VConst_Raja;
  v->ops->nvprod         = N_VProd_Raja;
  v->ops->nvdiv          = N_VDiv_Raja;
  v->ops->nvscale        = N_VScale_Raja;
  v->ops->nvabs          = N_VAbs_Raja;
  v->ops->nvinv          = N_VInv_Raja;
  v->ops->nvaddconst     = N_VAddConst_Raja;
  v->ops->nvdotprod      = N_VDotProd_Raja;
  v->ops->nvmaxnorm      = N_VMaxNorm_Raja;
  v->ops->nvmin          = N_VMin_Raja;
  v->ops->nvl1norm       = N_VL1Norm_Raja;
  v->ops->nvinvtest      = N_VInvTest_Raja;
  v->ops->nvconstrmask   = N_VConstrMask_Raja;
  v->ops->nvminquotient  = N_VMinQuotient_Raja;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Raja;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Raja;
  v->ops->nvwl2norm      = N_VWL2Norm_Raja;
  v->ops->nvcompare      = N_VCompare_Raja;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Raja;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Raja;
  v->ops->nvdotprodlocal     = N_VDotProd_Raja;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Raja;
  v->ops->nvminlocal         = N_VMin_Raja;
  v->ops->nvl1normlocal      = N_VL1Norm_Raja;
  v->ops->nvinvtestlocal     = N_VInvTest_Raja;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Raja;
  v->ops->nvminquotientlocal = N_VMinQuotient_Raja;

  v->content = (N_VectorContent_Raja) malloc(sizeof(_N_VectorContent_Raja));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  NVEC_RAJA_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Raja));
  if (NVEC_RAJA_CONTENT(v)->priv == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  NVEC_RAJA_CONTENT(v)->length          = 0;
  NVEC_RAJA_CONTENT(v)->own_data        = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  return(v);
}

N_Vector N_VNew_Raja(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->own_data        = SUNTRUE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return(v);
}

N_Vector N_VNewManaged_Raja(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  /* if using managed memory, we can attach an operation for
     nvgetarraypointer since the host and device pointers are the same */
  v->ops->nvgetarraypointer = N_VGetHostArrayPointer_Raja;

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->own_data        = SUNTRUE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return(v);
}


N_Vector N_VMake_Raja(sunindextype length, realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->own_data        = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data       = h_vdata;
  NVEC_RAJA_CONTENT(v)->device_data     = d_vdata;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  return(v);
}


N_Vector N_VMakeManaged_Raja(sunindextype length, realtype *vdata)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  /* if using managed memory, we can attach an operation for
     nvgetarraypointer since the host and device pointers are the same */
  v->ops->nvgetarraypointer = N_VGetHostArrayPointer_Raja;

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->own_data        = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data       = vdata;
  NVEC_RAJA_CONTENT(v)->device_data     = vdata;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNTRUE;

  return(v);
}

/* -----------------------------------------------------------------
 * Function to return the global length of the vector.
 */
sunindextype N_VGetLength_Raja(N_Vector v)
{
  return NVEC_RAJA_CONTENT(v)->length;
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Raja(N_Vector x)
{
  return NVEC_RAJA_CONTENT(x)->host_data;
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Raja(N_Vector x)
{
  return NVEC_RAJA_CONTENT(x)->device_data;
}

/* ----------------------------------------------------------------------------
 * Return a flag indicating if the memory for the vector data is managed
 */
booleantype N_VIsManagedMemory_Raja(N_Vector x)
{
  return NVEC_RAJA_PRIVATE(x)->use_managed_mem;
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Raja(N_Vector x)
{
  if(!NVEC_RAJA_PRIVATE(x)->use_managed_mem)
  {
    cudaError_t err = cudaMemcpyAsync(NVEC_RAJA_CONTENT(x)->device_data,
                                      NVEC_RAJA_CONTENT(x)->host_data,
                                      NVEC_RAJA_MEMSIZE(x),
                                      cudaMemcpyHostToDevice,
                                      0);
    SUNDIALS_CUDA_VERIFY(err);
  }
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Raja(N_Vector x)
{
  if(!NVEC_RAJA_PRIVATE(x)->use_managed_mem)
  {
    cudaError_t err = cudaMemcpyAsync(NVEC_RAJA_CONTENT(x)->host_data,
                                      NVEC_RAJA_CONTENT(x)->device_data,
                                      NVEC_RAJA_MEMSIZE(x),
                                      cudaMemcpyDeviceToHost,
                                      0);
   SUNDIALS_CUDA_VERIFY(err);
  }
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
  const realtype *xd = NVEC_RAJA_CONTENT(X)->host_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
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

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_RAJA_CONTENT(v)->length          = NVEC_RAJA_CONTENT(w)->length;
  NVEC_RAJA_CONTENT(v)->own_data        = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = NVEC_RAJA_PRIVATE(w)->use_managed_mem;


  return(v);
}

N_Vector N_VClone_Raja(N_Vector w)
{
  N_Vector v;
  v = NULL;
  v = N_VCloneEmpty_Raja(w);
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->own_data = SUNTRUE;
  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

return(v);

}


void N_VDestroy_Raja(N_Vector v)
{
  if (v == NULL) return;

  N_VectorContent_Raja vc = NVEC_RAJA_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

  if (vc->own_data)
  {
    if (vcp != NULL && !vcp->use_managed_mem) free(vc->host_data);
    cudaFree(vc->device_data);
    vc->device_data = NULL;
    vc->host_data = NULL;
  }

  /* free private content struct */
  if (vcp) free(vcp);
  vc->priv = NULL;

  /* free content struct */
  free(vc);
  v->content = NULL;

  /* free ops */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* free vector */
  free(v);
  v = NULL;

  return;
}

void N_VSpace_Raja(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_RAJA_CONTENT(X)->length;
  *liw = 2;
}

void N_VConst_Raja(realtype c, N_Vector Z)
{
  const sunindextype N = NVEC_RAJA_CONTENT(Z)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N), RAJA_LAMBDA(sunindextype i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const realtype *ydata = NVEC_RAJA_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = a*xdata[i] + b*ydata[i];
    }
  );
}

void N_VProd_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const realtype *ydata = NVEC_RAJA_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] * ydata[i];
    }
  );
}

void N_VDiv_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const realtype *ydata = NVEC_RAJA_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] / ydata[i];
    }
  );
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = c * xdata[i];
    }
  );
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]);
    }
  );
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = ONE / xdata[i];
    }
  );
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] + b;
    }
  );
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const realtype *ydata = NVEC_RAJA_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += xdata[i] * ydata[i] ;
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceMax< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.max(abs(xdata[i]));
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWSqrSumLocal_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const realtype *wdata = NVEC_RAJA_CONTENT(W)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWrmsNorm_Raja(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Raja(X, W);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  return std::sqrt(sum/N);
}

realtype N_VWSqrSumMaskLocal_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const realtype *wdata = NVEC_RAJA_CONTENT(W)->device_data;
  const realtype *iddata = NVEC_RAJA_CONTENT(ID)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (iddata[i] > ZERO)
        gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWrmsNormMask_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype sum = N_VWSqrSumMaskLocal_Raja(X, W, ID);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  return std::sqrt(sum/N);
}

realtype N_VMin_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.min(xdata[i]);
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWL2Norm_Raja(N_Vector X, N_Vector W)
{
  return std::sqrt(N_VWSqrSumLocal_Raja(X, W));
}

realtype N_VL1Norm_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (abs(xdata[i]));
    }
  );

  return (static_cast<realtype>(gpu_result));
}

void N_VCompare_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(X)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(Z)->device_data;

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
    }
  );
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = NVEC_RAJA_CONTENT(x)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype *zdata = NVEC_RAJA_CONTENT(z)->device_data;

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
  realtype minimum = static_cast<realtype>(gpu_result);
  return (minimum < HALF);
}

booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m)
{
  const realtype *cdata = NVEC_RAJA_CONTENT(c)->device_data;
  const realtype *xdata = NVEC_RAJA_CONTENT(x)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype *mdata = NVEC_RAJA_CONTENT(m)->device_data;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      bool test = (abs(cdata[i]) > ONEPT5 && cdata[i]*xdata[i] <= ZERO) ||
                  (abs(cdata[i]) > HALF   && cdata[i]*xdata[i] <  ZERO);
      mdata[i] = test ? ONE : ZERO;
      gpu_result += mdata[i];
    }
  );

  realtype sum = static_cast<realtype>(gpu_result);
  return(sum < HALF);
}

realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom)
{
  const realtype *ndata = NVEC_RAJA_CONTENT(num)->device_data;
  const realtype *ddata = NVEC_RAJA_CONTENT(denom)->device_data;
  const sunindextype N = NVEC_RAJA_CONTENT(num)->length;

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (ddata[i] != ZERO)
        gpu_result.min(ndata[i]/ddata[i]);
    }
  );
  return (static_cast<realtype>(gpu_result));
}


/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearCombination_Raja(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  cudaError_t  err;

  sunindextype N = NVEC_RAJA_CONTENT(z)->length;
  realtype* d_zd = NVEC_RAJA_CONTENT(z)->device_data;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = NVEC_RAJA_CONTENT(X[j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return(0);
}


int N_VScaleAddMulti_Raja(int nvec, realtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype* d_xd = NVEC_RAJA_CONTENT(x)->device_data;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Yd[j] = NVEC_RAJA_CONTENT(Y[j])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = NVEC_RAJA_CONTENT(Z[j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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

  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = NVEC_RAJA_CONTENT(X[j])->device_data;

  realtype** h_Yd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Yd[j] = NVEC_RAJA_CONTENT(Y[j])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = NVEC_RAJA_CONTENT(Z[j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return(0);
}


int N_VScaleVectorArray_Raja(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = NVEC_RAJA_CONTENT(X[j])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = NVEC_RAJA_CONTENT(Z[j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return(0);
}


int N_VConstVectorArray_Raja(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Create array of device pointers on host
  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = NVEC_RAJA_CONTENT(Z[j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return(0);
}


int N_VScaleAddMultiVectorArray_Raja(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;

  sunindextype N = NVEC_RAJA_CONTENT(X[0])->length;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = NVEC_RAJA_CONTENT(X[j])->device_data;

  realtype** h_Yd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Yd[j*nsum+k] = NVEC_RAJA_CONTENT(Y[k][j])->device_data;

  realtype** h_Zd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Zd[j*nsum+k] = NVEC_RAJA_CONTENT(Z[k][j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return(0);
}


int N_VLinearCombinationVectorArray_Raja(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Xd[j*nsum+k] = NVEC_RAJA_CONTENT(X[k][j])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = NVEC_RAJA_CONTENT(Z[j])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

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
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return(0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Raja;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Raja;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Raja;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Raja;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Raja;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Raja;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Raja;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Raja;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Raja;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Raja;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Raja;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Raja;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}


int AllocateData(N_Vector v)
{
  cudaError_t err;
  N_VectorContent_Raja vc = NVEC_RAJA_CONTENT(v);
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

  if (vcp->use_managed_mem)
  {
    err = cudaMallocManaged((void**) &vc->device_data, NVEC_RAJA_MEMSIZE(v));
    vc->host_data = vc->device_data;
    if (!SUNDIALS_CUDA_VERIFY(err))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: cudaMallocManaged failed\n");
      return -1;
    }
  }
  else
  {
    vc->host_data = (realtype*) malloc(NVEC_RAJA_MEMSIZE(v));
    if(vc->host_data == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: Malloc failed\n");
      return -1;
    }

    err = cudaMalloc((void**) &vc->device_data, NVEC_RAJA_MEMSIZE(v));
    if (!SUNDIALS_CUDA_VERIFY(err))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: cudaMalloc failed\n");
      return -1;
    }
  }

  return 0;

}

} // extern "C"
