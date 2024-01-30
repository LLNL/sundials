/* -----------------------------------------------------------------
 * Programmer(s): Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
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
#include <nvector/nvector_mpimanyvector.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef N_VectorContent_MPIManyVector N_VectorContent_MPIPlusX;

SUNDIALS_EXPORT
N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector_ID N_VGetVectorID_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT
sunrealtype* N_VGetArrayPointer_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT
void N_VSetArrayPointer_MPIPlusX(sunrealtype* vdata, N_Vector v);

SUNDIALS_EXPORT
N_Vector N_VGetLocalVector_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT
sunindextype N_VGetLocalLength_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_MPIPlusX(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
void N_VPrint_MPIPlusX(N_Vector x);

SUNDIALS_EXPORT
void N_VPrintFile_MPIPlusX(N_Vector x, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
