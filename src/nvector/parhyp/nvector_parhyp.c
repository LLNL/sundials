/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL and Jean M. Sexton @ SMU
 * -----------------------------------------------------------------
 * Based on work by Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                  and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a HYPRE ParVector wrapper
 * for the NVECTOR package.
 * -----------------------------------------------------------------*/
#include <nvector/nvector_parhyp.h>

/* --- Backend-specific headers --- */

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
#pragma message "hypre backend SERIAL confirmed from nvector_parhyp.c"
#include <stdio.h>
#include <stdlib.h>
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "sundials_debug.h"           /* located in src/nvector/sundials */
#endif

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
#pragma message "hypre backend CUDA confirmed from nvector_parhyp.c"
#include "sundials_cuda.h"            /* located in src/nvector/sundials */
#include "VectorKernels.cuh"          /* located in src/nvector/cuda     */
#include "VectorArrayKernels.cuh"     /* located in src/nvector/cuda     */

#elif defined(SUNDIALS_HYPRE_BACKENDS_HIP)
#pragma message "hypre backend HIP confirmed from nvector_parhyp.c"
#include "sundials_hip.h"             /* located in src/nvector/sundials */
#include "VectorKernels.hip.hpp"      /* located in src/nvector/hip      */
#include "VectorArrayKernels.hip.hpp" /* located in src/nvector/hip      */
#endif

/* --- Backend-specific namespaces --- */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
using namespace sundials;
using namespace sundials::cuda;
using namespace sundials::cuda::impl;

#elif defined(SUNDIALS_HYPRE_BACKENDS_HIP)
using namespace sundials;
using namespace sundials::hip;
using namespace sundials::hip::impl;
#endif

/* --- Defined constants --- */

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/*
 * -----------------------------------------------------------------
 * Simplifying macros NV_PH_CONTENT, NV_PH_DATA, NV_PH_LOCLENGTH,
 *                    NV_PH_GLOBLENGTH, and NV_PH_COMM
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * sunindextype v_len, s_len, i;
 *
 * (1) NV_PH_CONTENT
 *
 *     This routines gives access to the contents of the HYPRE
 *     vector wrapper (the N_Vector).
 *
 *     The assignment v_cont = NV_PH_CONTENT(v) sets v_cont to be
 *     a pointer to the N_Vector content structure.
 *
 * (2) NV_PH_DATA, NV_PH_LOCLENGTH, NV_PH_GLOBLENGTH, and NV_PH_COMM
 *
 *     These routines give access to the individual parts of
 *     the content structure of a parhyp N_Vector.
 *
 *     The assignment v_llen = NV_PH_LOCLENGTH(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_PH_LOCLENGTH(v) = llen_v generally should NOT be used! It
 *     will change locally stored value with the HYPRE local vector
 *     length, but it will NOT change the length of the actual HYPRE
 *     local vector.
 *
 *     The assignment v_glen = NV_PH_GLOBLENGTH(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_PH_GLOBLENGTH(v) = glen_v generally should NOT be used! It
 *     will change locally stored value with the HYPRE parallel vector
 *     length, but it will NOT change the length of the actual HYPRE
 *     parallel vector.
 *
 *     The assignment v_comm = NV_PH_COMM(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_PH_COMM_C(v) = comm_v sets the MPI communicator of v to be
 *     NV_PH_COMM(v) = comm_v generally should NOT be used! It
 *     will change locally stored value with the HYPRE parallel vector
 *     communicator, but it will NOT change the communicator of the
 *     actual HYPRE parallel vector.
 *
 * (3) NV_PH_DATA, NV_PH_HYPRE_PARVEC
 *
 *     The assignment v_data = NV_PH_DATA(v) sets v_data to be
 *     a pointer to the first component of the data inside the
 *     local vector of the HYPRE_parhyp vector for the vector v.
 *     The assignment NV_PH_DATA(v) = data_v should NOT be used.
 *     Instead, use NV_PH_HYPRE_PARVEC to obtain pointer to HYPRE
 *     vector and then use HYPRE functions to manipulate vector data.
 *
 *     The assignment v_parhyp = NV_PH_HYPRE_PARVEC(v) sets v_parhyp
 *     to be a pointer to HYPRE_ParVector of vector v. The assignment
 *     NV_PH_HYPRE_PARVEC(v) = parhyp_v sets pointer to
 *     HYPRE_ParVector of vector v to be parhyp_v.
 *
 * -----------------------------------------------------------------
 */

/* --- (TEMPORARY) Reference accesor macros defined by hypre --- */

// For hypre ParVector
// #define hypre_ParVectorComm(vector)             ((vector) -> comm)
// #define hypre_ParVectorGlobalSize(vector)       ((vector) -> global_size)
// #define hypre_ParVectorFirstIndex(vector)       ((vector) -> first_index)
// #define hypre_ParVectorLastIndex(vector)        ((vector) -> last_index)
// #define hypre_ParVectorPartitioning(vector)     ((vector) -> partitioning)
// #define hypre_ParVectorActualLocalSize(vector)  ((vector) -> actual_local_size)
// #define hypre_ParVectorLocalVector(vector)      ((vector) -> local_vector)
// #define hypre_ParVectorOwnsData(vector)         ((vector) -> owns_data)
// #define hypre_ParVectorAllZeros(vector)         ((vector) -> all_zeros)
// #define hypre_ParVectorNumVectors(vector)       (hypre_VectorNumVectors(hypre_ParVectorLocalVector(vector)))
// #define hypre_ParVectorAssumedPartition(vector) ((vector) -> assumed_partition)
// static inline HYPRE_MemoryLocation
// hypre_ParVectorMemoryLocation(hypre_ParVector *vector)
// {
//    return hypre_VectorMemoryLocation(hypre_ParVectorLocalVector(vector));
// }

// For hypre hypreVector (local_vector inside a ParVector)
// #define hypre_VectorData(vector)                  ((vector) -> data)
// #define hypre_VectorSize(vector)                  ((vector) -> size)
// #define hypre_VectorComponent(vector)             ((vector) -> component)
// #define hypre_VectorOwnsData(vector)              ((vector) -> owns_data)
// #define hypre_VectorMemoryLocation(vector)        ((vector) -> memory_location)
// #define hypre_VectorNumVectors(vector)            ((vector) -> num_vectors)
// #define hypre_VectorMultiVecStorageMethod(vector) ((vector) -> multivec_storage_method)
// #define hypre_VectorVectorStride(vector)          ((vector) -> vecstride)
// #define hypre_VectorIndexStride(vector)           ((vector) -> idxstride)

/* --- Common accessor macros --- */

#define NV_PH_CONTENT(v)      ( (N_VectorContent_ParHyp)(v->content) )
#define NV_PH_LOCLENGTH(v)    ( NV_PH_CONTENT(v)->local_length )
#define NV_PH_GLOBLENGTH(v)   ( NV_PH_CONTENT(v)->global_length )
#define NV_PH_OWN_PARVEC(v)   ( NV_PH_CONTENT(v)->own_parvector )
#define NV_PH_COMM(v)         ( NV_PH_CONTENT(v)->comm )
// hypre ParVector accessor macros
#define NV_PH_HYPRE_PARVEC(v) ( NV_PH_CONTENT(v)->x )
#define NV_PH_HYPRE_MEMLOC(v) ( (HYPRE_MemoryLocation) hypre_ParVectorMemoryLocation(NV_PH_HYPRE_PARVEC(v)) )

/* --- Backend-dependent accessor macros --- */

#define NV_PH_DATA(v)         ( NV_PH_HYPRE_PARVEC(v) == NULL ? NULL : hypre_VectorData(hypre_ParVectorLocalVector(NV_PH_HYPRE_PARVEC(v))) )

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
#define NV_PH_MEMHELP(v)      (NV_PH_CONTENT(v)->mem_helper)
#define NV_PH_MEMSIZE(v)      (NV_PH_CONTENT(v)->length * sizeof(realtype))
#define NV_PH_STREAM(v)       (NV_PH_CONTENT(v)->stream_exec_policy->stream())
// Private content accessor macros
#define NV_PH_PRIVATE(v)      ((N_PrivateVectorContent_ParHyp)(NV_PH_CONTENT(v)->priv))
#define NV_PH_HBUFFERp(v)     ((realtype*) NV_PH_PRIVATE(v)->reduce_buffer_host->ptr)
#define NV_PH_DBUFFERp(v)     ((realtype*) NV_PH_PRIVATE(v)->reduce_buffer_dev->ptr)
#define NV_PH_DCOUNTERp(v)    ((unsigned int*) NV_PH_PRIVATE(v)->device_counter->ptr)
#endif

/* --- Private structure definition --- */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
struct _N_PrivateVectorContent_ParHyp
{
  booleantype use_managed_mem; /* do data pointers use managed memory */

  // Reduction workspace
  SUNMemory device_counter;      // device memory for a counter (used in LDS reductions)
  SUNMemory reduce_buffer_dev;   // device memory for reductions
  SUNMemory reduce_buffer_host;  // host memory for reductions
  size_t    reduce_buffer_bytes; // current size of reduction buffers

  // // Fused vector operations workspace
  // SUNMemory fused_buffer_dev;    // device memory for fused ops
  // SUNMemory fused_buffer_host;   // host memory for fused ops
  // size_t    fused_buffer_bytes;  // current size of the buffers
  // size_t    fused_buffer_offset; // current offset into the buffer
};
typedef struct _N_PrivateVectorContent_ParHyp *N_PrivateVectorContent_ParHyp;
#endif

/* --- Default execution policies --- */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
ThreadDirectExecPolicy DEFAULT_STREAMING_EXECPOLICY(256);
BlockReduceAtomicExecPolicy DEFAULT_REDUCTION_EXECPOLICY(256);
#elif defined(SUNDIALS_HYPRE_BACKENDS_HIP)
ThreadDirectExecPolicy DEFAULT_STREAMING_EXECPOLICY(512);
BlockReduceExecPolicy DEFAULT_REDUCTION_EXECPOLICY(512);
#endif

/* --- Private function prototypes --- */

/* z=x+y */
static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=ax+y */
static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_ParHyp(N_Vector v)
{
  return SUNDIALS_NVEC_PARHYP;
}


/* ----------------------------------------------------------------
 * Function to create a new parhyp vector without underlying
 * HYPRE vector.
 */
N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm,
                            sunindextype local_length,
                            sunindextype global_length,
                            SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_ParHyp content;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_ParHyp;
  v->ops->nvclone           = N_VClone_ParHyp;
  v->ops->nvcloneempty      = N_VCloneEmpty_ParHyp;
  v->ops->nvdestroy         = N_VDestroy_ParHyp;
  v->ops->nvspace           = N_VSpace_ParHyp;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_ParHyp;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_ParHyp;
  v->ops->nvgetcommunicator = N_VGetCommunicator_ParHyp;
  v->ops->nvgetlength       = N_VGetLength_ParHyp;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_ParHyp;
  v->ops->nvconst        = N_VConst_ParHyp;
  v->ops->nvprod         = N_VProd_ParHyp;
  v->ops->nvdiv          = N_VDiv_ParHyp;
  v->ops->nvscale        = N_VScale_ParHyp;
  v->ops->nvabs          = N_VAbs_ParHyp;
  v->ops->nvinv          = N_VInv_ParHyp;
  v->ops->nvaddconst     = N_VAddConst_ParHyp;
  v->ops->nvdotprod      = N_VDotProd_ParHyp;
  v->ops->nvmaxnorm      = N_VMaxNorm_ParHyp;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_ParHyp;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_ParHyp;
  v->ops->nvmin          = N_VMin_ParHyp;
  v->ops->nvwl2norm      = N_VWL2Norm_ParHyp;
  v->ops->nvl1norm       = N_VL1Norm_ParHyp;
  v->ops->nvcompare      = N_VCompare_ParHyp;
  v->ops->nvinvtest      = N_VInvTest_ParHyp;
  v->ops->nvconstrmask   = N_VConstrMask_ParHyp;
  v->ops->nvminquotient  = N_VMinQuotient_ParHyp;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProdLocal_ParHyp;
  v->ops->nvmaxnormlocal     = N_VMaxNormLocal_ParHyp;
  v->ops->nvminlocal         = N_VMinLocal_ParHyp;
  v->ops->nvl1normlocal      = N_VL1NormLocal_ParHyp;
  v->ops->nvinvtestlocal     = N_VInvTestLocal_ParHyp;
  v->ops->nvconstrmasklocal  = N_VConstrMaskLocal_ParHyp;
  v->ops->nvminquotientlocal = N_VMinQuotientLocal_ParHyp;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_ParHyp;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_ParHyp;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = N_VDotProdMultiLocal_ParHyp;
  v->ops->nvdotprodmultiallreduce = N_VDotProdMultiAllReduce_ParHyp;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_ParHyp;
  v->ops->nvbufpack   = N_VBufPack_ParHyp;
  v->ops->nvbufunpack = N_VBufUnpack_ParHyp;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ParHyp) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Attach lengths and communicator */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_parvector = SUNFALSE;
  content->x             = NULL;

  return(v);
}


/* ----------------------------------------------------------------
 * Function to create a parhyp N_Vector wrapper around user
 * supplie HYPRE vector.
 */

N_Vector N_VMake_ParHyp(HYPRE_ParVector x, SUNContext sunctx)
{
  N_Vector v;
  MPI_Comm comm = hypre_ParVectorComm(x);
  sunindextype global_length = (sunindextype) hypre_ParVectorGlobalSize(x);
  sunindextype local_begin   = (sunindextype) hypre_ParVectorFirstIndex(x);
  sunindextype local_end     = (sunindextype) hypre_ParVectorLastIndex(x);
  sunindextype local_length  = local_end - local_begin + 1;

  v = NULL;
  v = N_VNewEmpty_ParHyp(comm, local_length, global_length, sunctx);
  if (v == NULL)
    return(NULL);

  NV_PH_OWN_PARVEC(v)   = SUNFALSE;
  NV_PH_HYPRE_PARVEC(v) = x;

  return(v);
}


/* ----------------------------------------------------------------
 * Function to create an array of new parhyp vectors.
 */

N_Vector *N_VCloneVectorArray_ParHyp(int count, N_Vector w)
{
  return(N_VCloneVectorArray(count, w));
}

/* ----------------------------------------------------------------
 * Function to create an array of new parhyp vector wrappers
 * without uderlying HYPRE vectors.
 */

N_Vector *N_VCloneVectorArrayEmpty_ParHyp(int count, N_Vector w)
{
  return(N_VCloneEmptyVectorArray(count, w));
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_ParHyp
 */

void N_VDestroyVectorArray_ParHyp(N_Vector *vs, int count)
{
  N_VDestroyVectorArray(vs, count);
  return;
}


/* ----------------------------------------------------------------
 * Extract HYPRE vector
 */

HYPRE_ParVector N_VGetVector_ParHyp(N_Vector v)
{
  return NV_PH_HYPRE_PARVEC(v);
}

/* ----------------------------------------------------------------
 * Function to print a parhyp vector.
 * TODO: Consider using a HYPRE function for this.
 */

void N_VPrint_ParHyp(N_Vector x)
{
  N_VPrintFile_ParHyp(x, stdout);
}

/* ----------------------------------------------------------------
 * Function to print a parhyp vector.
 * TODO: Consider using a HYPRE function for this.
 */

void N_VPrintFile_ParHyp(N_Vector x, FILE *outfile)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%g\n", xd[i]);
#else
    fprintf(outfile, "%g\n", xd[i]);
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

N_Vector N_VCloneEmpty_ParHyp(N_Vector w)
{
  N_Vector v;
  N_VectorContent_ParHyp content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ParHyp) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = NV_PH_LOCLENGTH(w);
  content->global_length = NV_PH_GLOBLENGTH(w);
  content->comm          = NV_PH_COMM(w);
  content->own_parvector = SUNFALSE;
  content->x             = NULL;

  return(v);
}

/*
 * Clone HYPRE vector wrapper.
 *
 */
N_Vector N_VClone_ParHyp(N_Vector w)
{
  N_Vector v;
  HYPRE_ParVector vx;
  const HYPRE_ParVector wx = NV_PH_HYPRE_PARVEC(w);

  v = NULL;
  v = N_VCloneEmpty_ParHyp(w);
  if (v==NULL) return(NULL);

  vx = hypre_ParVectorCreate(wx->comm, wx->global_size, wx->partitioning);
  hypre_ParVectorInitialize(vx);

#ifdef hypre_ParVectorOwnsPartitioning
  hypre_ParVectorSetPartitioningOwner(vx, 0);
#endif
  hypre_ParVectorSetDataOwner(vx, 1);
  hypre_SeqVectorSetDataOwner(hypre_ParVectorLocalVector(vx), 1);

  NV_PH_HYPRE_PARVEC(v) = vx;
  NV_PH_OWN_PARVEC(v)   = SUNTRUE;

  return(v);
}

void N_VDestroy_ParHyp(N_Vector v)
{
  if (v == NULL) return;

  /* free content */
  if (v->content != NULL) {
    /* free the hypre parvector if it's owned by the vector wrapper */
    if (NV_PH_OWN_PARVEC(v) && NV_PH_HYPRE_PARVEC(v) != NULL) {
      hypre_ParVectorDestroy(NV_PH_HYPRE_PARVEC(v));
      NV_PH_HYPRE_PARVEC(v) = NULL;
    }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}


void N_VSpace_ParHyp(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_PH_COMM(v);
  MPI_Comm_size(comm, &npes);

  *lrw = NV_PH_GLOBLENGTH(v);
  *liw = 2*npes;

  return;
}


/*
 * This function is disabled in ParHyp implementation and returns NULL.
 * The user should extract HYPRE vector using N_VGetVector_ParHyp and
 * then use HYPRE functions to get pointer to raw data of the local HYPRE
 * vector.
 */
realtype *N_VGetArrayPointer_ParHyp(N_Vector v)
{
  return NULL; /* ((realtype *) NV_PH_DATA(v)); */
}


/*
 * This method is not implemented for HYPRE vector wrapper.
 * TODO: Put error handler in the function body.
 */
void N_VSetArrayPointer_ParHyp(realtype *v_data, N_Vector v)
{
  /* Not implemented for Hypre vector */
}


void *N_VGetCommunicator_ParHyp(N_Vector v)
{
  return((void *) &(NV_PH_COMM(v)));
}

sunindextype N_VGetLength_ParHyp(N_Vector v)
{
  return(NV_PH_GLOBLENGTH(v));
}

/*
 * Computes z[i] = a*x[i] + b*y[i]
 *
 */
void N_VLinearSum_ParHyp(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    HYPRE_Complex   alpha=a;
    HYPRE_ParVectorAxpy( alpha, (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(x),
                                (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(y));
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    HYPRE_Complex   beta=b;
    HYPRE_ParVectorAxpy( beta, (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(y),
                               (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(x));
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_ParHyp(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_ParHyp(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_ParHyp(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_ParHyp(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_ParHyp(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_ParHyp(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}

void N_VConst_ParHyp(realtype c, N_Vector z)
{
  HYPRE_Complex value = c;
  HYPRE_ParVectorSetConstantValues( (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(z), value);
  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  // if Serial
  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];
  // if CUDA kernal launch... divKernel<<<...>>>

  return;
}


void N_VScale_ParHyp(realtype c, N_Vector x, N_Vector z)
{
  HYPRE_Complex value = c;

  if (x != z) {
     HYPRE_ParVectorCopy((HYPRE_ParVector) NV_PH_HYPRE_PARVEC(x), (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(z));
  }
  HYPRE_ParVectorScale(value, (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(z));

  return;
}


void N_VAbs_ParHyp(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = SUNRabs(xd[i]);

  return;
}

void N_VInv_ParHyp(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}

void N_VAddConst_ParHyp(N_Vector x, realtype b, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
     zd[i] = xd[i] + b;

  return;
}

realtype N_VDotProdLocal_ParHyp(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype sum, *xd, *yd;
  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);

  sum = ZERO;
  for (i = 0; i < N; i++)
    sum += xd[i]*yd[i];
  return(sum);
}

realtype N_VDotProd_ParHyp(N_Vector x, N_Vector y)
{
  HYPRE_Real gsum;
  HYPRE_ParVectorInnerProd( (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(x),
                            (HYPRE_ParVector) NV_PH_HYPRE_PARVEC(y), &gsum);

  return(gsum);
}

realtype N_VMaxNormLocal_ParHyp(N_Vector x)
{
  sunindextype i, N;
  realtype max, *xd;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);

  max = ZERO;
  for (i = 0; i < N; i++)
    if (SUNRabs(xd[i]) > max) max = SUNRabs(xd[i]);
  return(max);
}

realtype N_VMaxNorm_ParHyp(N_Vector x)
{
  realtype lmax, gmax;
  lmax = N_VMaxNormLocal_ParHyp(x);
  MPI_Allreduce(&lmax, &gmax, 1, MPI_SUNREALTYPE, MPI_MAX, NV_PH_COMM(x));
  return(gmax);
}

realtype N_VWSqrSumLocal_ParHyp(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  wd = NV_PH_DATA(w);

  sum = ZERO;
  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }
  return(sum);
}

realtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w)
{
  realtype lsum, gsum;
  lsum = N_VWSqrSumLocal_ParHyp(x, w);
  MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_PH_COMM(x));
  return(SUNRsqrt(gsum/(NV_PH_GLOBLENGTH(x))));
}

realtype N_VWSqrSumMaskLocal_ParHyp(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd, *idd;

  N   = NV_PH_LOCLENGTH(x);
  xd  = NV_PH_DATA(x);
  wd  = NV_PH_DATA(w);
  idd = NV_PH_DATA(id);

  sum = ZERO;
  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i]*wd[i];
      sum += SUNSQR(prodi);
    }
  }
  return(sum);
}

realtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w, N_Vector id)
{
  realtype lsum, gsum;
  lsum = N_VWSqrSumMaskLocal_ParHyp(x, w, id);
  MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_PH_COMM(x));
  return(SUNRsqrt(gsum/(NV_PH_GLOBLENGTH(x))));
}

realtype N_VMinLocal_ParHyp(N_Vector x)
{
  sunindextype i, N;
  realtype min, *xd;

  xd  = NULL;
  N   = NV_PH_LOCLENGTH(x);
  min = BIG_REAL;

  if (N > 0) {
    xd = NV_PH_DATA(x);
    min = xd[0];
    for (i = 1; i < N; i++)
      if (xd[i] < min) min = xd[i];
  }
  return(min);
}

realtype N_VMin_ParHyp(N_Vector x)
{
  realtype lmin, gmin;
  lmin = N_VMinLocal_ParHyp(x);
  MPI_Allreduce(&lmin, &gmin, 1, MPI_SUNREALTYPE, MPI_MIN, NV_PH_COMM(x));
  return(gmin);
}

realtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w)
{
  realtype lsum, gsum;
  lsum = N_VWSqrSumLocal_ParHyp(x, w);
  MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_PH_COMM(x));
  return(SUNRsqrt(gsum));
}

realtype N_VL1NormLocal_ParHyp(N_Vector x)
{
  sunindextype i, N;
  realtype sum, *xd;

  N   = NV_PH_LOCLENGTH(x);
  xd  = NV_PH_DATA(x);
  sum = ZERO;

  for (i = 0; i<N; i++)  sum += SUNRabs(xd[i]);
  return(sum);
}

realtype N_VL1Norm_ParHyp(N_Vector x)
{
  realtype lsum, gsum;
  lsum = N_VL1NormLocal_ParHyp(x);
  MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_PH_COMM(x));
  return(gsum);
}

void N_VCompare_ParHyp(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++) {
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}

booleantype N_VInvTestLocal_ParHyp(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd, val;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  zd = NV_PH_DATA(z);

  val = ONE;
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO)
      val = ZERO;
    else
      zd[i] = ONE/xd[i];
  }
  if (val == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z)
{
  realtype val, gval;
  val = (N_VInvTestLocal_ParHyp(x, z)) ? ONE : ZERO;
  MPI_Allreduce(&val, &gval, 1, MPI_SUNREALTYPE, MPI_MIN, NV_PH_COMM(x));
  if (gval == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VConstrMaskLocal_ParHyp(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  booleantype test;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  cd = NV_PH_DATA(c);
  md = NV_PH_DATA(m);

  temp = ZERO;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i]*cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF   && xd[i]*cd[i] <  ZERO);
    if (test) {
      temp = md[i] = ONE;
    }
  }

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

booleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x, N_Vector m)
{
  realtype temp, temp2;
  temp = (N_VConstrMaskLocal_ParHyp(c, x, m)) ? ZERO : ONE;
  MPI_Allreduce(&temp, &temp2, 1, MPI_SUNREALTYPE, MPI_MAX, NV_PH_COMM(x));
  return (temp2 == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotientLocal_ParHyp(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  sunindextype i, N;
  realtype *nd, *dd, min;

  nd = dd = NULL;

  N  = NV_PH_LOCLENGTH(num);
  nd = NV_PH_DATA(num);
  dd = NV_PH_DATA(denom);

  notEvenOnce = SUNTRUE;
  min = BIG_REAL;

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
      else {
        min = nd[i]/dd[i];
        notEvenOnce = SUNFALSE;
      }
    }
  }
  return(min);
}

realtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom)
{
  realtype lmin, gmin;
  lmin = N_VMinQuotientLocal_ParHyp(num, denom);
  MPI_Allreduce(&lmin, &gmin, 1, MPI_SUNREALTYPE, MPI_MIN, NV_PH_COMM(num));
  return(gmin);
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */


int N_VLinearCombination_ParHyp(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_ParHyp(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_ParHyp(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_PH_LOCLENGTH(z);
  zd = NV_PH_DATA(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
    for (i=1; i<nvec; i++) {
      xd = NV_PH_DATA(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z) {
    for (j=0; j<N; j++) {
      zd[j] *= c[0];
    }
    for (i=1; i<nvec; i++) {
      xd = NV_PH_DATA(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_PH_DATA(X[0]);
  for (j=0; j<N; j++) {
    zd[j] = c[0] * xd[j];
  }
  for (i=1; i<nvec; i++) {
    xd = NV_PH_DATA(X[i]);
    for (j=0; j<N; j++) {
      zd[j] += c[i] * xd[j];
    }
  }
  return(0);
}


int N_VScaleAddMulti_ParHyp(int nvec, realtype* a, N_Vector x, N_Vector* Y,
                             N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_ParHyp(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      yd = NV_PH_DATA(Y[i]);
      for (j=0; j<N; j++) {
        yd[j] += a[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    yd = NV_PH_DATA(Y[i]);
    zd = NV_PH_DATA(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a[i] * xd[j] + yd[j];
    }
  }
  return(0);
}


int N_VDotProdMulti_ParHyp(int nvec, N_Vector x, N_Vector* Y,
                           realtype* dotprods)
{
  int          i, retval;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_ParHyp(x, Y[0]);
    return(0);
  }

  /* get vector length, data array, and communicator */
  N    = NV_PH_LOCLENGTH(x);
  xd   = NV_PH_DATA(x);
  comm = NV_PH_COMM(x);

  /* compute multiple dot products */
  for (i=0; i<nvec; i++) {
    yd = NV_PH_DATA(Y[i]);
    dotprods[i] = ZERO;
    for (j=0; j<N; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }
  retval = MPI_Allreduce(MPI_IN_PLACE, dotprods, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  return retval == MPI_SUCCESS ? 0 : -1;
}


/*
 * -----------------------------------------------------------------
 * single buffer reduction operations
 * -----------------------------------------------------------------
 */


int N_VDotProdMultiLocal_ParHyp(int nvec, N_Vector x, N_Vector* Y,
                                realtype* dotprods)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* get vector length, data array, and communicator */
  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);

  /* compute multiple dot products */
  for (i=0; i<nvec; i++) {
    yd = NV_PH_DATA(Y[i]);
    dotprods[i] = ZERO;
    for (j=0; j<N; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }

  return 0;
}


int N_VDotProdMultiAllReduce_ParHyp(int nvec, N_Vector x, realtype* sum)
{
  int retval;
  retval = MPI_Allreduce(MPI_IN_PLACE, sum, nvec, MPI_SUNREALTYPE, MPI_SUM,
                         NV_PH_COMM(x));
  return retval == MPI_SUCCESS ? 0 : -1;
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */


int N_VLinearSumVectorArray_ParHyp(int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_ParHyp(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_PH_LOCLENGTH(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i=0; i<nvec; i++) {
    xd = NV_PH_DATA(X[i]);
    yd = NV_PH_DATA(Y[i]);
    zd = NV_PH_DATA(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a * xd[j] + b * yd[j];
    }
  }

  return(0);
}


int N_VScaleVectorArray_ParHyp(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_ParHyp(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_PH_LOCLENGTH(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_PH_DATA(X[i]);
      for (j=0; j<N; j++) {
        xd[j] *= c[i];
      }
    }
    return(0);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_PH_DATA(X[i]);
    zd = NV_PH_DATA(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c[i] * xd[j];
    }
  }

  return(0);
}


int N_VConstVectorArray_ParHyp(int nvec, realtype c, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_ParHyp(c, Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_PH_LOCLENGTH(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i=0; i<nvec; i++) {
    zd = NV_PH_DATA(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c;
    }
  }

  return(0);
}


int N_VWrmsNormVectorArray_ParHyp(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  int          i, retval;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_ParHyp(X[0], W[0]);
    return(0);
  }

  /* get vector lengths and communicator */
  Nl   = NV_PH_LOCLENGTH(X[0]);
  Ng   = NV_PH_GLOBLENGTH(X[0]);
  comm = NV_PH_COMM(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_PH_DATA(X[i]);
    wd = NV_PH_DATA(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
  }
  retval = MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return retval == MPI_SUCCESS ? 0 : -1;
}


int N_VWrmsNormMaskVectorArray_ParHyp(int nvec, N_Vector* X, N_Vector* W,
                                        N_Vector id, realtype* nrm)
{
  int          i, retval;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  realtype*    idd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_ParHyp(X[0], W[0], id);
    return(0);
  }

  /* get vector lengths, communicator, and mask data */
  Nl   = NV_PH_LOCLENGTH(X[0]);
  Ng   = NV_PH_GLOBLENGTH(X[0]);
  comm = NV_PH_COMM(X[0]);
  idd  = NV_PH_DATA(id);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_PH_DATA(X[i]);
    wd = NV_PH_DATA(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      if (idd[j] > ZERO)
        nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
  }
  retval = MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return retval == MPI_SUCCESS ? 0 : -1;
}


int N_VScaleAddMultiVectorArray_ParHyp(int nvec, int nsum, realtype* a,
                                          N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int          i, j;
  sunindextype k, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  int          retval;
  N_Vector*    YY;
  N_Vector*    ZZ;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VLinearSum */
    if (nsum == 1) {
      N_VLinearSum_ParHyp(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector *) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_ParHyp(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_ParHyp(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N  = NV_PH_LOCLENGTH(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_PH_DATA(X[i]);
      for (j=0; j<nsum; j++){
        yd = NV_PH_DATA(Y[j][i]);
        for (k=0; k<N; k++) {
          yd[k] += a[j] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_PH_DATA(X[i]);
    for (j=0; j<nsum; j++) {
      yd = NV_PH_DATA(Y[j][i]);
      zd = NV_PH_DATA(Z[j][i]);
      for (k=0; k<N; k++) {
        zd[k] = a[j] * xd[k] + yd[k];
      }
    }
  }
  return(0);
}


int N_VLinearCombinationVectorArray_ParHyp(int nvec, int nsum,
                                             realtype* c,
                                             N_Vector** X,
                                             N_Vector* Z)
{
  int          i; /* vector arrays index in summation [0,nsum) */
  int          j; /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  realtype*    ctmp;
  N_Vector*    Y;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VScale */
    if (nsum == 1) {
      N_VScale_ParHyp(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_ParHyp(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    N_VLinearCombination_ParHyp(nsum, c, Y, Z[0]);

    free(Y);
    return(0);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1) {

    ctmp = (realtype*) malloc(nvec * sizeof(realtype));

    for (j=0; j<nvec; j++) {
      ctmp[j] = c[0];
    }

    N_VScaleVectorArray_ParHyp(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(0);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    N_VLinearSumVectorArray_ParHyp(nvec, c[0], X[0], c[1], X[1], Z);
    return(0);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_PH_LOCLENGTH(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE)) {
    for (j=0; j<nvec; j++) {
      zd = NV_PH_DATA(Z[j]);
      for (i=1; i<nsum; i++) {
        xd = NV_PH_DATA(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z) {
    for (j=0; j<nvec; j++) {
      zd = NV_PH_DATA(Z[j]);
      for (k=0; k<N; k++) {
        zd[k] *= c[0];
      }
      for (i=1; i<nsum; i++) {
        xd = NV_PH_DATA(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j=0; j<nvec; j++) {
    xd = NV_PH_DATA(X[0][j]);
    zd = NV_PH_DATA(Z[j]);
    for (k=0; k<N; k++) {
      zd[k] = c[0] * xd[k];
    }
    for (i=1; i<nsum; i++) {
      xd = NV_PH_DATA(X[i][j]);
      for (k=0; k<N; k++) {
        zd[k] += c[i] * xd[k];
      }
    }
  }
  return(0);
}


/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */


int N_VBufSize_ParHyp(N_Vector x, sunindextype *size)
{
  if (x == NULL) return(-1);
  *size = NV_PH_LOCLENGTH(x) * ((sunindextype)sizeof(realtype));
  return(0);
}


int N_VBufPack_ParHyp(N_Vector x, void *buf)
{
  sunindextype i, N;
  realtype     *xd = NULL;
  realtype     *bd = NULL;

  if (x == NULL || buf == NULL) return(-1);

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  bd = (realtype*) buf;

  for (i = 0; i < N; i++)
    bd[i] = xd[i];

  return(0);
}


int N_VBufUnpack_ParHyp(N_Vector x, void *buf)
{
  sunindextype i, N;
  realtype     *xd = NULL;
  realtype     *bd = NULL;

  if (x == NULL || buf == NULL) return(-1);

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  bd = (realtype*) buf;

  for (i = 0; i < N; i++)
    xd[i] = bd[i];

  return(0);
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}


static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]+yd[i]);

  return;
}

static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+yd[i];

  return;
}

static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_PH_LOCLENGTH(x);
  xd = NV_PH_DATA(x);
  yd = NV_PH_DATA(y);
  zd = NV_PH_DATA(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])-yd[i];

  return;
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_ParHyp;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_ParHyp;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_ParHyp;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_ParHyp;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_ParHyp;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_ParHyp;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_ParHyp;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_ParHyp;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_ParHyp;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_ParHyp;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_ParHyp;
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
    /* disable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_ParHyp;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_ParHyp;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_ParHyp;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_ParHyp;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_ParHyp;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_ParHyp;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_ParHyp;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_ParHyp;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_ParHyp;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_ParHyp;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMultiLocal_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_ParHyp;
  else
    v->ops->nvdotprodmultilocal = NULL;

  /* return success */
  return(0);
}
