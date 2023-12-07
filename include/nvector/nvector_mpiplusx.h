/* -----------------------------------------------------------------
 * Programmer(s): Cody Balos @ LLNL
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
 * This is the header file for the MPI+X implementation of the
 * NVECTOR module. The MPIPlusX NVECTOR is really just an extension
 * of the ManyVector.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_MPIPLUSX_H
#define _NVECTOR_MPIPLUSX_H

#include <mpi.h>
#include <sundials/sundials_core.h>
#include <nvector/nvector_mpimanyvector.h>
#include <sundials/priv/sundials_errors_impl.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef N_VectorContent_MPIManyVector N_VectorContent_MPIPlusX;

SUNDIALS_EXPORT
N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X, SUNContext sunctx);

SUNDIALS_STATIC_INLINE
N_Vector_ID N_VGetVectorID_MPIPlusX(N_Vector v)
{
  return SUNDIALS_NVEC_MPIPLUSX;
}

SUNDIALS_STATIC_INLINE
sunrealtype* N_VGetArrayPointer_MPIPlusX(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  sunrealtype* arr = N_VGetSubvectorArrayPointer_MPIManyVector(v, 0); SUNCheckLastErrNoRet();
  return arr;
}

SUNDIALS_STATIC_INLINE
void N_VSetArrayPointer_MPIPlusX(sunrealtype* vdata, N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  N_VSetSubvectorArrayPointer_MPIManyVector(vdata, v, 0); SUNCheckLastErrNoRet();
}

SUNDIALS_STATIC_INLINE
N_Vector N_VGetLocalVector_MPIPlusX(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  N_Vector result = N_VGetSubvector_MPIManyVector(v, 0); SUNCheckLastErrNoRet();
  return result;
}

SUNDIALS_STATIC_INLINE
sunindextype N_VGetLocalLength_MPIPlusX(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  sunindextype len = N_VGetLength(N_VGetLocalVector_MPIPlusX(v)); SUNCheckLastErrNoRet();
  return len;
}

SUNDIALS_STATIC_INLINE
SUNErrCode N_VEnableFusedOps_MPIPlusX(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);
  SUNCheckCall(N_VEnableFusedOps_MPIManyVector(v, tf));
  return SUN_SUCCESS;
}

SUNDIALS_EXPORT
void N_VPrint_MPIPlusX(N_Vector x);

SUNDIALS_EXPORT 
void N_VPrintFile_MPIPlusX(N_Vector x, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
