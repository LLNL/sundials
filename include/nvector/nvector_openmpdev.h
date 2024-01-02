/* -------------------------------------------------------------------
 * Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
 * -------------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
 * -------------------------------------------------------------------
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
 * This is the header file for the OpenMP 4.5+ implementation of the
 * NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definition of the type 'sunrealtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'sunbooleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_OpenMPDEV(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_OPENMPDEV_H
#define _NVECTOR_OPENMPDEV_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * OpenMPDEV implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_OpenMPDEV
{
  sunindextype length;     /* vector length       */
  sunbooleantype own_data; /* data ownership flag */
  sunrealtype* host_data;  /* host data array     */
  sunrealtype* dev_data;   /* device data array   */
};

typedef struct _N_VectorContent_OpenMPDEV* N_VectorContent_OpenMPDEV;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_OMPDEV, NV_DATA_HOST_OMPDEV, NV_OWN_DATA_OMPDEV,
 *        NV_LENGTH_OMPDEV, and NV_Ith_OMPDEV
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_OMPDEV(v) ((N_VectorContent_OpenMPDEV)(v->content))

#define NV_LENGTH_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->length)

#define NV_OWN_DATA_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->own_data)

#define NV_DATA_HOST_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->host_data)

#define NV_DATA_DEV_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->dev_data)

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_openmpdev
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_OpenMPDEV(sunindextype vec_length,
                                          SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VNewEmpty_OpenMPDEV(sunindextype vec_length,
                                               SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_OpenMPDEV(sunindextype vec_length,
                                           sunrealtype* h_data,
                                           sunrealtype* v_data,
                                           SUNContext sunctx);

SUNDIALS_EXPORT sunindextype N_VGetLength_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT sunrealtype* N_VGetHostArrayPointer_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT sunrealtype* N_VGetDeviceArrayPointer_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_OpenMPDEV(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT void N_VCopyToDevice_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT void N_VCopyFromDevice_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_OpenMPDEV(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_OpenMPDEV(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_OpenMPDEV(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_OpenMPDEV(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_OpenMPDEV(N_Vector v, sunindextype* lrw,
                                        sunindextype* liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_OpenMPDEV(sunrealtype a, N_Vector x,
                                            sunrealtype b, N_Vector y,
                                            N_Vector z);
SUNDIALS_EXPORT void N_VConst_OpenMPDEV(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_OpenMPDEV(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_OpenMPDEV(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_OpenMPDEV(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_OpenMPDEV(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_OpenMPDEV(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_OpenMPDEV(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_OpenMPDEV(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_OpenMPDEV(N_Vector x, N_Vector w,
                                                      N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_OpenMPDEV(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_OpenMPDEV(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_OpenMPDEV(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_OpenMPDEV(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_OpenMPDEV(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_OpenMPDEV(N_Vector c, N_Vector x,
                                                       N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_OpenMPDEV(N_Vector num,
                                                     N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearCombination_OpenMPDEV(int nvec,
                                                          sunrealtype* c,
                                                          N_Vector* V,
                                                          N_Vector z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMulti_OpenMPDEV(int nvec, sunrealtype* a,
                                                      N_Vector x, N_Vector* Y,
                                                      N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMulti_OpenMPDEV(int nvec, N_Vector x,
                                                     N_Vector* Y,
                                                     sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearSumVectorArray_OpenMPDEV(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleVectorArray_OpenMPDEV(int nvec, sunrealtype* c,
                                                         N_Vector* X,
                                                         N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VConstVectorArray_OpenMPDEV(int nvecs, sunrealtype c,
                                                         N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormVectorArray_OpenMPDEV(int nvecs,
                                                            N_Vector* X,
                                                            N_Vector* W,
                                                            sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormMaskVectorArray_OpenMPDEV(
  int nvecs, N_Vector* X, N_Vector* W, N_Vector id, sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMultiVectorArray_OpenMPDEV(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT SUNErrCode N_VLinearCombinationVectorArray_OpenMPDEV(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal_OpenMPDEV(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal_OpenMPDEV(N_Vector x, N_Vector w,
                                                          N_Vector id);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_OpenMPDEV(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_OpenMPDEV(N_Vector v,
                                                                sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_OpenMPDEV(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_OpenMPDEV(N_Vector v,
                                                           sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_OpenMPDEV(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_OpenMPDEV(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_OpenMPDEV(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_OpenMPDEV(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_OpenMPDEV(N_Vector v,
                                                      sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_OpenMPDEV(N_Vector v,
                                                       sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_OpenMPDEV(N_Vector v,
                                                           sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
