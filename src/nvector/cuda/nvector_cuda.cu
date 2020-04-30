/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, and Cody J. Balos @ LLNL
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
 * This is the implementation file for a CUDA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <nvector/nvector_cuda.h>
#include "VectorKernels.cuh"
#include "VectorArrayKernels.cuh"

#include "sundials_cuda.h"
#include "sundials_debug.h"

#define HALF   RCONST(0.5)

extern "C" {

using namespace sundials;
using namespace sundials::nvector_cuda;

/*
 * Macro definitions
 */

#define NVEC_CUDA_CONTENT(x) ((N_VectorContent_Cuda)(x->content))
#define NVEC_CUDA_PRIVATE(x) ((N_PrivateVectorContent_Cuda)(NVEC_CUDA_CONTENT(x)->priv))
#define NVEC_CUDA_MEMSIZE(x) (NVEC_CUDA_CONTENT(x)->length * sizeof(realtype))
#define NVEC_CUDA_STREAM(x)  (NVEC_CUDA_CONTENT(x)->stream_exec_policy->stream())

/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Cuda
{
    booleantype use_managed_mem; /* indicates if the data pointers and buffer pointers are managed memory */
    size_t      reduce_buffer_allocated_bytes; /* current size of the reduction buffer */
    realtype*   reduce_buffer_dev;      /* device buffer used for reductions */
    realtype*   reduce_buffer_host;     /* host buffer used for reductions */
    void*       (*userallocfn)(size_t); /* a user provided allocator (assumes managed mem) */
    void        (*userfreefn)(void*);   /* a user provided free function */
    booleantype own_exec; /* indicates if the exec policy is owned by the vector */
};

typedef struct _N_PrivateVectorContent_Cuda *N_PrivateVectorContent_Cuda;

/*
 * Private function definitions
 */

static int AllocateData(N_Vector v);
static void AllocateReductionBuffer(N_Vector v);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n);
static void GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid, size_t& block,
                                size_t& shMemSize, cudaStream_t& stream, size_t n = 0);
static void PostKernelLaunch();                               

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Cuda(N_Vector v)
{
  return SUNDIALS_NVEC_CUDA;
}

N_Vector N_VNewEmpty_Cuda()
{
  N_Vector v;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Cuda;
  v->ops->nvclone           = N_VClone_Cuda;
  v->ops->nvcloneempty      = N_VCloneEmpty_Cuda;
  v->ops->nvdestroy         = N_VDestroy_Cuda;
  v->ops->nvspace           = N_VSpace_Cuda;
  v->ops->nvgetlength       = N_VGetLength_Cuda;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Cuda;
  v->ops->nvconst        = N_VConst_Cuda;
  v->ops->nvprod         = N_VProd_Cuda;
  v->ops->nvdiv          = N_VDiv_Cuda;
  v->ops->nvscale        = N_VScale_Cuda;
  v->ops->nvabs          = N_VAbs_Cuda;
  v->ops->nvinv          = N_VInv_Cuda;
  v->ops->nvaddconst     = N_VAddConst_Cuda;
  v->ops->nvdotprod      = N_VDotProd_Cuda;
  v->ops->nvmaxnorm      = N_VMaxNorm_Cuda;
  v->ops->nvmin          = N_VMin_Cuda;
  v->ops->nvl1norm       = N_VL1Norm_Cuda;
  v->ops->nvinvtest      = N_VInvTest_Cuda;
  v->ops->nvconstrmask   = N_VConstrMask_Cuda;
  v->ops->nvminquotient  = N_VMinQuotient_Cuda;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Cuda;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Cuda;
  v->ops->nvwl2norm      = N_VWL2Norm_Cuda;
  v->ops->nvcompare      = N_VCompare_Cuda;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProd_Cuda;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Cuda;
  v->ops->nvminlocal         = N_VMin_Cuda;
  v->ops->nvl1normlocal      = N_VL1Norm_Cuda;
  v->ops->nvinvtestlocal     = N_VInvTest_Cuda;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Cuda;
  v->ops->nvminquotientlocal = N_VMinQuotient_Cuda;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Cuda;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Cuda;

  /* Create content */

  v->content = (N_VectorContent_Cuda) malloc(sizeof(_N_VectorContent_Cuda));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  NVEC_CUDA_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Cuda));
  if (NVEC_CUDA_CONTENT(v)->priv == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  NVEC_CUDA_CONTENT(v)->length                        = 0;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaStreamingExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaReductionExecPolicy(256);
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = NULL;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;
  
  return(v);
}

N_Vector N_VNew_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaStreamingExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaReductionExecPolicy(256);
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = NULL;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  AllocateData(v);

  return(v);
}

N_Vector N_VNewManaged_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* if using managed memory, we can attach an operation for
     nvgetarraypointer since the host and device pointers are the same */
  v->ops->nvgetarraypointer = N_VGetHostArrayPointer_Cuda;

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaStreamingExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaReductionExecPolicy(256);
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = NULL;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  AllocateData(v);

  return(v);
}

N_Vector N_VMake_Cuda(sunindextype length, realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->host_data                     = h_vdata;
  NVEC_CUDA_CONTENT(v)->device_data                   = d_vdata;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaStreamingExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaReductionExecPolicy(256);
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = NULL;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  return(v);
}

N_Vector N_VMakeManaged_Cuda(sunindextype length, realtype *vdata)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* if using managed memory, we can attach an operation for
     nvgetarraypointer since the host and device pointers are the same */
  v->ops->nvgetarraypointer = N_VGetHostArrayPointer_Cuda;

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->host_data                     = vdata;
  NVEC_CUDA_CONTENT(v)->device_data                   = vdata;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaStreamingExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaReductionExecPolicy(256);
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = NULL;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  return(v);
}

N_Vector N_VMakeWithManagedAllocator_Cuda(sunindextype length,
                                          void* (*allocfn)(size_t),
                                          void (*freefn)(void*))
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* if using managed memory, we can attach an operation for
     nvgetarraypointer since the host and device pointers are the same */
  v->ops->nvgetarraypointer = N_VGetHostArrayPointer_Cuda;

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaStreamingExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaReductionExecPolicy(256);
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = allocfn;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = freefn;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  AllocateData(v);

  return(v);
}

/* -----------------------------------------------------------------
 * Function to return the global length of the vector.
 */
sunindextype N_VGetLength_Cuda(N_Vector x)
{
  return NVEC_CUDA_CONTENT(x)->length;
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Cuda(N_Vector x)
{
  return NVEC_CUDA_CONTENT(x)->host_data;
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Cuda(N_Vector x)
{
  return NVEC_CUDA_CONTENT(x)->device_data;
}

/* ----------------------------------------------------------------------------
 * Return a flag indicating if the memory for the vector data is managed
 */
booleantype N_VIsManagedMemory_Cuda(N_Vector x)
{
  return NVEC_CUDA_PRIVATE(x)->use_managed_mem;
}

int N_VSetKernelExecPolicy_Cuda(N_Vector x,
                                SUNCudaExecPolicy* stream_exec_policy,
                                SUNCudaExecPolicy* reduce_exec_policy)
{
  if (x == NULL || stream_exec_policy == NULL ||
      reduce_exec_policy == NULL)
    return -1;

  if (NVEC_CUDA_PRIVATE(x)->own_exec)
  {
    delete NVEC_CUDA_CONTENT(x)->stream_exec_policy;
    delete NVEC_CUDA_CONTENT(x)->reduce_exec_policy;
  }
  NVEC_CUDA_CONTENT(x)->stream_exec_policy = stream_exec_policy;
  NVEC_CUDA_CONTENT(x)->reduce_exec_policy = reduce_exec_policy;
  NVEC_CUDA_PRIVATE(x)->own_exec = SUNFALSE;
  return 0;
}

/*
 * ----------------------------------------------------------------------------
 * DEPRECATED: will be removed in SUNDIALS v6.
 * Sets the cudaStream_t to use for execution of the CUDA kernels.
 */
void N_VSetCudaStream_Cuda(N_Vector x, cudaStream_t *stream)
{
  CudaStreamingExecPolicy* s = 
    new CudaStreamingExecPolicy(NVEC_CUDA_CONTENT(x)->stream_exec_policy->blockSize(), *stream);
  CudaReductionExecPolicy* r =
    new CudaReductionExecPolicy(NVEC_CUDA_CONTENT(x)->reduce_exec_policy->blockSize(), *stream);
  N_VSetKernelExecPolicy_Cuda(x, s, r);
  NVEC_CUDA_PRIVATE(x)->own_exec = SUNTRUE;
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Cuda(N_Vector x)
{
  cudaError_t err;

  /* If the host and device pointers are the same, then we don't need
     to do a copy (this happens in the managed memory case), but we
     still need to synchronize the device to adhere to the unified
     memory access rules. */
  if (NVEC_CUDA_PRIVATE(x)->use_managed_mem)
  {
    err = cudaStreamSynchronize(NVEC_CUDA_STREAM(x));
    SUNDIALS_CUDA_VERIFY(err);
  }
  else
  {
    err = cudaMemcpyAsync(NVEC_CUDA_CONTENT(x)->device_data,
                          NVEC_CUDA_CONTENT(x)->host_data,
                          NVEC_CUDA_MEMSIZE(x),
                          cudaMemcpyHostToDevice,
                          NVEC_CUDA_STREAM(x));
    SUNDIALS_CUDA_VERIFY(err);
  }
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Cuda(N_Vector x)
{
  cudaError_t err;

  /* If the host and device pointers are the same, then we don't need
     to do a copy (this happens in the managed memory case), but we
     still need to synchronize the device to adhere to the unified
     memory access rules. */
  if (NVEC_CUDA_PRIVATE(x)->use_managed_mem)
  {
    err = cudaStreamSynchronize(NVEC_CUDA_STREAM(x));
    SUNDIALS_CUDA_VERIFY(err);
  }
  else
  {
    err = cudaMemcpyAsync(NVEC_CUDA_CONTENT(x)->host_data,
                          NVEC_CUDA_CONTENT(x)->device_data,
                          NVEC_CUDA_MEMSIZE(x),
                          cudaMemcpyDeviceToHost,
                          NVEC_CUDA_STREAM(x));
    SUNDIALS_CUDA_VERIFY(err);
  }
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to stdout
 */

void N_VPrint_Cuda(N_Vector x)
{
  N_VPrintFile_Cuda(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to outfile
 */

void N_VPrintFile_Cuda(N_Vector x, FILE *outfile)
{
  sunindextype i;

  for (i = 0; i < NVEC_CUDA_CONTENT(x)->length; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", NVEC_CUDA_CONTENT(x)->host_data[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", NVEC_CUDA_CONTENT(x)->host_data[i]);
#else
    fprintf(outfile, "%11.8g\n", NVEC_CUDA_CONTENT(x)->host_data[i]);
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

N_Vector N_VCloneEmpty_Cuda(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_CUDA_CONTENT(v)->length                        = NVEC_CUDA_CONTENT(w)->length;
  NVEC_CUDA_CONTENT(v)->own_data                      = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = NVEC_CUDA_CONTENT(w)->stream_exec_policy->clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = NVEC_CUDA_CONTENT(w)->reduce_exec_policy->clone();
  NVEC_CUDA_PRIVATE(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = NVEC_CUDA_PRIVATE(w)->use_managed_mem;
  NVEC_CUDA_PRIVATE(v)->userallocfn                   = NVEC_CUDA_PRIVATE(w)->userallocfn;
  NVEC_CUDA_PRIVATE(v)->userfreefn                    = NVEC_CUDA_PRIVATE(w)->userfreefn;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  return(v);
}

N_Vector N_VClone_Cuda(N_Vector w)
{
  N_Vector v;
  v = NULL;
  v = N_VCloneEmpty_Cuda(w);
  if (v == NULL) return(NULL);

  AllocateData(v);

  return(v);
}


void N_VDestroy_Cuda(N_Vector v)
{
  if (v == NULL) return;

  N_VectorContent_Cuda vc = NVEC_CUDA_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  N_PrivateVectorContent_Cuda priv = NVEC_CUDA_PRIVATE(v);

  /* free content */
  if (vc->own_data)
  {
    if (priv != NULL && priv->userfreefn)
    {
      priv->userfreefn(vc->device_data);
      vc->device_data = NULL;
      vc->host_data = NULL;
    }
    else
    {
      if (priv != NULL && !priv->use_managed_mem) free(vc->host_data);
      SUNDIALS_CUDA_VERIFY(cudaFree(vc->device_data));
      vc->device_data = NULL;
      vc->host_data = NULL;
    }
  }

  /* free reduction buffer */
  FreeReductionBuffer(v);

  /* free execution policies */
  if (priv != NULL && priv->own_exec)
  {
    if (vc->stream_exec_policy) delete vc->stream_exec_policy;
    if (vc->reduce_exec_policy) delete vc->reduce_exec_policy;
  }

  /* free private content */
  if (priv) free(priv);
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

void N_VSpace_Cuda(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_CUDA_CONTENT(X)->length;
  *liw = 2;
}

void N_VConst_Cuda(realtype a, N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  setConstKernel<<<grid, block, shMemSize, stream>>>
  (
    a,
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VLinearSum_Cuda(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  linearSumKernel<<<grid, block, shMemSize, stream>>>
  (
    a,
    NVEC_CUDA_CONTENT(X)->device_data,
    b,
    NVEC_CUDA_CONTENT(Y)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VProd_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  prodKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Y)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VDiv_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  divKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Y)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VScale_Cuda(realtype a, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  scaleKernel<<<grid, block, shMemSize, stream>>>
  (
    a,
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VAbs_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  absKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VInv_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  invKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VAddConst_Cuda(N_Vector X, realtype b, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  addConstKernel<<<grid, block, shMemSize, stream>>>
  (
    b,
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

realtype N_VDotProd_Cuda(N_Vector X, N_Vector Y)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  dotProdKernel<realtype,sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Y)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    sumReduceKernel<realtype,sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype gpu_result = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }

  return gpu_result;
}

realtype N_VMaxNorm_Cuda(N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  maxNormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    maxNormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype gpu_result = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    if (NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i] > gpu_result)
      gpu_result = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }
  return gpu_result;
}

realtype N_VWSqrSumLocal_Cuda(N_Vector X, N_Vector W)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  wL2NormSquareKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(W)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    sumReduceKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype gpu_result =  NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result +=  NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }
  return gpu_result;
}

realtype N_VWrmsNorm_Cuda(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Cuda(X, W);
  return std::sqrt(sum/NVEC_CUDA_CONTENT(X)->length);
}

realtype N_VWSqrSumMaskLocal_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  wL2NormSquareMaskKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(W)->device_data,
    NVEC_CUDA_CONTENT(Id)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    sumReduceKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype gpu_result = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result += NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }
  return gpu_result;
}

realtype N_VWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  const realtype sum = N_VWSqrSumMaskLocal_Cuda(X, W, Id);
  return std::sqrt(sum/NVEC_CUDA_CONTENT(X)->length);
}

realtype N_VMin_Cuda(N_Vector X)
{
  realtype maxVal = std::numeric_limits<realtype>::max();

  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  findMinKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    maxVal,
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    findMinKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      maxVal,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype gpu_result = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    if (NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i] < gpu_result)
      gpu_result = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }
  return gpu_result;
}

realtype N_VWL2Norm_Cuda(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Cuda(X, W);
  return std::sqrt(sum);
}

realtype N_VL1Norm_Cuda(N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  L1NormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    sumReduceKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype gpu_result =  NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    gpu_result +=  NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }
  return gpu_result;
}

void N_VCompare_Cuda(realtype c, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  compareKernel<<<grid, block, shMemSize, stream>>>
  (
    c,
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

booleantype N_VInvTest_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  invTestKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    sumReduceKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype locmin = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    locmin += NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }
  
  return (locmin < HALF);
}

booleantype N_VConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(X); /* an allocation occurs only if it is necessary */

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  constrMaskKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_CONTENT(C)->device_data,
    NVEC_CUDA_CONTENT(X)->device_data,
    NVEC_CUDA_CONTENT(M)->device_data,
    NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(X, true, grid, block, shMemSize, stream, n);
    sumReduceKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(X)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X, n);
  realtype loc_sum = NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    loc_sum += NVEC_CUDA_PRIVATE(X)->reduce_buffer_host[i];
  }

  return (loc_sum < HALF);
}

realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  // Starting value for min reduction
  const realtype maxVal = std::numeric_limits<realtype>::max();

  size_t grid, block, shMemSize;
  cudaStream_t stream;

  AllocateReductionBuffer(num); /* an allocation occurs only if it is necessary */

  GetKernelParameters(num, true, grid, block, shMemSize, stream);
  minQuotientKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    maxVal,
    NVEC_CUDA_CONTENT(num)->device_data,
    NVEC_CUDA_CONTENT(denom)->device_data,
    NVEC_CUDA_PRIVATE(num)->reduce_buffer_dev,
    NVEC_CUDA_CONTENT(num)->length
  );
  PostKernelLaunch();

  // All quotients are computed by now. Find the minimum.
  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // (Re)run reduction kernel
    GetKernelParameters(num, true, grid, block, shMemSize, stream, n);
    findMinKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
    (
      maxVal,
      NVEC_CUDA_PRIVATE(num)->reduce_buffer_dev,
      NVEC_CUDA_PRIVATE(num)->reduce_buffer_dev,
      n
    );
    PostKernelLaunch();
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(num, n);
  realtype gpu_result = NVEC_CUDA_PRIVATE(num)->reduce_buffer_host[0];
  for (unsigned int i=1; i<n; i++)
  {
    if (NVEC_CUDA_PRIVATE(num)->reduce_buffer_host[i] < gpu_result)
      gpu_result = NVEC_CUDA_PRIVATE(num)->reduce_buffer_host[i];
  }
  return gpu_result;
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_CONTENT(X[i])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Set kernel parameters and launch
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X[0], false, grid, block, shMemSize, stream);
  linearCombinationKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_c,
    d_Xd,
    NVEC_CUDA_CONTENT(Z)->device_data,
    NVEC_CUDA_CONTENT(Z)->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = cudaFree(d_c);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return 0;
}

int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector X, N_Vector* Y,
                          N_Vector* Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_CUDA_CONTENT(Y[i])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_CONTENT(Z[i])->device_data;

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

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  scaleAddMultiKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_c,
    NVEC_CUDA_CONTENT(X)->device_data,
    d_Yd,
    d_Zd,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

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

  return 0;
}


int N_VDotProdMulti_Cuda(int nvec, N_Vector X, N_Vector* Y, realtype* dots)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_CUDA_CONTENT(Y[i])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;
  
  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  // shMemSize has to be computed differently for this case
  shMemSize = nvec*block*sizeof(realtype);

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  dotProdMultiKernel<realtype,sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    NVEC_CUDA_CONTENT(X)->device_data,
    d_Yd,
    d_buff,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax) 
  {
    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    sumReduceVectorKernel<realtype,sunindextype><<<grid, block, shMemSize, stream>>>
    (
      nvec,
      d_buff,
      d_buff,
      n
    );
    PostKernelLaunch();

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  realtype* h_buff = new realtype[nvec*n*sizeof(realtype)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(realtype), cudaMemcpyDeviceToHost);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  for (int k=0; k<nvec; k++) {
    dots[k] = h_buff[k*n];
    for (unsigned i=1; i<n; i++){
      dots[k] += h_buff[i + k*n];
    }
  }

  // Free host array
  delete[] h_Yd;
  delete[] h_buff;

  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_buff);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return 0;
}



/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Cuda(int nvec, realtype a, N_Vector* X, realtype b,
                                 N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_CONTENT(X[i])->device_data;

  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_CUDA_CONTENT(Y[i])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_CONTENT(Z[i])->device_data;

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

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(Z[0], false, grid, block, shMemSize, stream);
  linearSumVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    a,
    d_Xd,
    b,
    d_Yd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

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

  return 0;
}


int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_CONTENT(X[i])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_CONTENT(Z[i])->device_data;

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

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(Z[0], false, grid, block, shMemSize, stream);
  scaleVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_c,
    d_Xd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return 0;
}


int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_CONTENT(Z[i])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(Z[0], false, grid, block, shMemSize, stream);
  constVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    c,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return 0;
}

int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                realtype* norms)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_CONTENT(X[i])->device_data;

  realtype** h_Wd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = NVEC_CUDA_CONTENT(W[i])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;
  
  GetKernelParameters(X[0], true, grid, block, shMemSize, stream);
  // shMemSize has to be computed differently for this case
  shMemSize = nvec*block*sizeof(realtype);

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  wL2NormSquareVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_Xd,
    d_Wd,
    d_buff,
    NVEC_CUDA_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    sumReduceVectorKernel<realtype,sunindextype><<<grid, block, shMemSize, stream>>>
    (
      nvec,
      d_buff,
      d_buff,
      n
    );
    PostKernelLaunch();

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  realtype* h_buff = new realtype[nvec*n*sizeof(realtype)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(realtype), cudaMemcpyDeviceToHost);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  for (int k=0; k<nvec; k++) {
    norms[k] = h_buff[k*n];
    for (unsigned i=1; i<n; i++){
      norms[k] += h_buff[i + k*n];
    }
  }

  for (int k=0; k<nvec; ++k)
    norms[k] = std::sqrt(norms[k]/NVEC_CUDA_CONTENT(X[0])->length);

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;
  delete[] h_buff;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Wd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_buff);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return 0;
}

int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                    N_Vector id, realtype* norms)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_CONTENT(X[i])->device_data;

  realtype** h_Wd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = NVEC_CUDA_CONTENT(W[i])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  realtype** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Get kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;
  
  GetKernelParameters(X[0], true, grid, block, shMemSize, stream);
  // shMemSize has to be computed differently for this case
  shMemSize = nvec*block*sizeof(realtype);

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = cudaMalloc((void**) &d_buff, nvec*grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  wL2NormSquareMaskVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_Xd,
    d_Wd,
    NVEC_CUDA_CONTENT(id)->device_data,
    d_buff,
    NVEC_CUDA_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  unsigned n = grid;
  unsigned nmax = 2*block;
  while (n > nmax)
  {
    // Recompute partitioning
    grid = (n + block - 1)/block;

    // Rerun reduction kernel
    sumReduceVectorKernel<realtype,sunindextype><<<grid, block, shMemSize, stream>>>
    (
      nvec,
      d_buff,
      d_buff,
      n
    );
    PostKernelLaunch();

    // update buffer array working length
    n = grid;
  }

  // Finish reduction on CPU if there are less than two blocks of data left.
  realtype* h_buff = new realtype[nvec*n*sizeof(realtype)];
  err = cudaMemcpy(h_buff, d_buff, nvec*n*sizeof(realtype), cudaMemcpyDeviceToHost);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  for (int k=0; k<nvec; k++) {
    norms[k] = h_buff[k*n];
    for (unsigned i=1; i<n; i++){
      norms[k] += h_buff[i + k*n];
    }
  }

  for (int k=0; k<nvec; ++k)
    norms[k] = std::sqrt(norms[k]/NVEC_CUDA_CONTENT(X[0])->length);

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;
  delete[] h_buff;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_Wd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaFree(d_buff);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  return 0;
}


int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return -1;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_CONTENT(X[i])->device_data;

  realtype** h_Yd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Yd[j*nsum+i] = NVEC_CUDA_CONTENT(Y[i][j])->device_data;

  realtype** h_Zd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Zd[j*nsum+i] = NVEC_CUDA_CONTENT(Z[i][j])->device_data;

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

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(Z[0][0], false, grid, block, shMemSize, stream);
  scaleAddMultiVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    nsum,
    d_c,
    d_Xd,
    d_Yd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0][0])->length
  );
  PostKernelLaunch();

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

  return 0;
}


int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Xd[j*nsum+i] = NVEC_CUDA_CONTENT(X[i][j])->device_data;

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_CONTENT(Z[i])->device_data;

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(Z[0], false, grid, block, shMemSize, stream);
  linearCombinationVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    nsum,
    d_c,
    d_Xd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return err;

  return cudaGetLastError();
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Cuda;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Cuda;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Cuda;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Cuda;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Cuda;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Cuda;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Cuda;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Cuda;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Cuda;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;
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


int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Cuda;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Cuda;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_Cuda;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Cuda;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Cuda;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Cuda;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Cuda;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Cuda;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Cuda;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}

/*
 * Private helper functions.
 */

int AllocateData(N_Vector v)
{
  int err;

  if (NVEC_CUDA_PRIVATE(v)->userallocfn)
  {
    /* We assume managed memory when a custom allocator is provided */
    NVEC_CUDA_CONTENT(v)->device_data = 
      (realtype *) NVEC_CUDA_PRIVATE(v)->userallocfn(NVEC_CUDA_MEMSIZE(v));
    NVEC_CUDA_CONTENT(v)->host_data = NVEC_CUDA_CONTENT(v)->device_data;
    err = NVEC_CUDA_CONTENT(v)->device_data != NULL;
  } 
  else if (NVEC_CUDA_PRIVATE(v)->use_managed_mem)
  {
    err = !SUNDIALS_CUDA_VERIFY(
      cudaMallocManaged((void**) &NVEC_CUDA_CONTENT(v)->device_data, NVEC_CUDA_MEMSIZE(v))
    );
    NVEC_CUDA_CONTENT(v)->host_data = NVEC_CUDA_CONTENT(v)->device_data;
  }
  else
  {
    NVEC_CUDA_CONTENT(v)->host_data = (realtype*) malloc(NVEC_CUDA_MEMSIZE(v));
    err = !SUNDIALS_CUDA_VERIFY(
      cudaMalloc((void**) &NVEC_CUDA_CONTENT(v)->device_data, NVEC_CUDA_MEMSIZE(v))
    );
  }

  return err ? -1 : 0;
}

/* 
 * Computes the number of bytes needed for the reduction buffer
 * based on the current length.
 */
size_t ReduceBufferBytesNeeded(N_Vector v)
{
  N_VectorContent_Cuda vc = NVEC_CUDA_CONTENT(v);
  return vc->reduce_exec_policy->gridSize(vc->length) * sizeof(realtype);
}

/* 
 * Allocates the internal buffer used for reductions.
 * If the buffer is already allocated, it will only be reallocated
 * if it is no longer large enough. This may occur if the length
 * of the vector is increased.
 */
void AllocateReductionBuffer(N_Vector v)
{
  cudaError_t err;
  size_t bytes = ReduceBufferBytesNeeded(v);
  N_PrivateVectorContent_Cuda priv = NVEC_CUDA_PRIVATE(v);

  /* we allocate if the existing reduction buffer is not large enough */
  if (priv->reduce_buffer_allocated_bytes < bytes)
  {
    if (priv->reduce_buffer_allocated_bytes) FreeReductionBuffer(v);
  }
  else
  {
    return;
  }

  if (priv->userallocfn != nullptr) {

    priv->reduce_buffer_dev = (realtype*) priv->userallocfn(bytes);
    if(priv->reduce_buffer_dev == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateReductionBuffer: could not allocate device data with user allocator\n");
    }
    priv->reduce_buffer_host = priv->reduce_buffer_dev;

  } else if (priv->use_managed_mem) {

    err = cudaMallocManaged((void**) &priv->reduce_buffer_dev, bytes);
    SUNDIALS_CUDA_VERIFY(err);
    priv->reduce_buffer_host = priv->reduce_buffer_dev;

  } else {

    priv->reduce_buffer_host = (realtype *) malloc(bytes);
    if(priv->reduce_buffer_host == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateReductionBuffer: could not allocate host data\n");
    }
    err = cudaMalloc((void**) &priv->reduce_buffer_dev, bytes);
    SUNDIALS_CUDA_VERIFY(err);

  }

  priv->reduce_buffer_allocated_bytes = bytes;
}

/* Free the reduction buffer
 */
void FreeReductionBuffer(N_Vector v)
{
  cudaError_t err;

  N_PrivateVectorContent_Cuda priv = NVEC_CUDA_PRIVATE(v);
  if (priv == NULL) return;

  if (priv->use_managed_mem) {
    /* managed memory */
    if (priv->userfreefn) {
      if (priv->reduce_buffer_dev != NULL) 
        priv->userfreefn(priv->reduce_buffer_dev);
    } else {
      if (priv->reduce_buffer_dev != NULL) {
        err = cudaFree(priv->reduce_buffer_dev);
        SUNDIALS_CUDA_VERIFY(err);
      }
    }
    priv->reduce_buffer_dev = priv->reduce_buffer_host = NULL;
  } else {
    /* unmanaged memory */
    if (priv->reduce_buffer_dev != NULL) {
      err = cudaFree(priv->reduce_buffer_dev);
      SUNDIALS_CUDA_VERIFY(err);
    }
    if (priv->reduce_buffer_host != NULL) free(priv->reduce_buffer_host);
    priv->reduce_buffer_dev = NULL;
    priv->reduce_buffer_host = NULL;
  }
}

/* Copy the reduction buffer from the device to the host.
 */
int CopyReductionBufferFromDevice(N_Vector v, size_t n)
{
  cudaError_t err;

  /* If using managed memory, then we don't need to do a copy, but we
      still need to synchronize the device to adhere to the unified
      memory access rules. */
  if (NVEC_CUDA_PRIVATE(v)->use_managed_mem)
  { 
    err = cudaStreamSynchronize(NVEC_CUDA_STREAM(v));
  }
  else
  {
    err = cudaMemcpyAsync(NVEC_CUDA_PRIVATE(v)->reduce_buffer_host,
                          NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev,
                          n*sizeof(realtype),
                          cudaMemcpyDeviceToHost,
                          NVEC_CUDA_STREAM(v));
  }
  return (!SUNDIALS_CUDA_VERIFY(err)) ? -1 : 0;
}

/* Get the kernel launch parameters based on the kernel type (reduction or not),
 * using the appropriate kernel execution policy.
 */
void GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid, size_t& block,
                         size_t& shMemSize, cudaStream_t& stream, size_t n)
{
  n = (n == 0) ? NVEC_CUDA_CONTENT(v)->length : n;
  if (reduction)
  {
    SUNCudaExecPolicy* reduce_exec_policy = NVEC_CUDA_CONTENT(v)->reduce_exec_policy;
    grid      = reduce_exec_policy->gridSize(n);
    block     = reduce_exec_policy->blockSize();
    shMemSize = reduce_exec_policy->blockSize() * sizeof(realtype);
    stream    = reduce_exec_policy->stream();
  }
  else
  {
    SUNCudaExecPolicy* stream_exec_policy = NVEC_CUDA_CONTENT(v)->stream_exec_policy;
    grid      = stream_exec_policy->gridSize(n);
    block     = stream_exec_policy->blockSize();
    shMemSize = 0;
    stream    = stream_exec_policy->stream();
  }

#ifdef SUNDIALS_DEBUG
  if (grid == 0)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in NVECTOR_CUDA: kernel grid size must be > 0\n");
  }

  if (block == 0)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in NVECTOR_CUDA: kernel block size must be > 0\n");
  }
#endif
}

/* Should be called after a kernel launch.
 * If SUNDIALS_DEBUG_CUDA_LASTERROR is not defined, then the function does nothing.
 * If it is defined, the function will synchronize and check the last CUDA error.
 */
void PostKernelLaunch()
{
#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  SUNDIALS_CUDA_VERIFY(cudaGetLastError());
#endif
}


} // extern "C"
