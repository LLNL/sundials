/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file contains the prototypes for functions to
 * test an NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#ifndef _TEST_NVECTOR_COMPLEX_H
#define _TEST_NVECTOR_COMPLEX_H

#include <math.h>
#include <sundials/sundials_types.h>

/* define constants */
#define NEG_TWO  SUN_RCONST(-2.0)
#define NEG_ONE  SUN_RCONST(-1.0)
#define NEG_HALF SUN_RCONST(-0.5)
#define ZERO     SUN_RCONST(0.0)
#define HALF     SUN_RCONST(0.5)
#define ONE      SUN_RCONST(1.0)
#define TWO      SUN_RCONST(2.0)

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

extern SUNContext sunctx;

/* Forward declarations for implementation specific utility functions */
int check_ans_Z(sunscalartype ans, N_Vector X, sunindextype local_length);
void set_element_Z(N_Vector X, sunindextype i, sunscalartype val);
void set_element_range_Z(N_Vector X, sunindextype is, sunindextype ie,
                         sunscalartype val);
sunscalartype get_element_Z(N_Vector X, sunindextype i);

/* Vector ID test */
int Test_N_VGetVectorID_Z(N_Vector X, N_Vector_ID ID, int myid);

/* Vector constructor tests */
int Test_N_VMake_Z(N_Vector X, sunindextype local_length, int myid);

/* Vector clone tests */
int Test_N_VCloneVectorArray_Z(int count, N_Vector W, sunindextype local_length,
                             int myid);
int Test_N_VCloneEmptyVectorArray_Z(int count, N_Vector W, int myid);
int Test_N_VCloneEmpty_Z(N_Vector W, int myid);
int Test_N_VClone_Z(N_Vector W, sunindextype local_length, int myid);

/* Vector get/set function tests */
int Test_N_VGetArrayPointer_Z(N_Vector W, sunindextype local_length, int myid);
int Test_N_VSetArrayPointer_Z(N_Vector W, sunindextype local_length, int myid);
int Test_N_VGetLength_Z(N_Vector W, int myid);
int Test_N_VGetCommunicator_Z(N_Vector W, SUNComm comm, int myid);
int Test_N_VGetCommunicatorMPI_Z(N_Vector W, SUNComm comm, int myid);

/* Standard vector operation tests */
int Test_N_VLinearSum_Z(N_Vector X, N_Vector Y, N_Vector Z,
                      sunindextype local_length, int myid);
int Test_N_VConst_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VProd_Z(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length,
                 int myid);
int Test_N_VDiv_Z(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length,
                int myid);
int Test_N_VScale_Z(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
int Test_N_VAbs_Z(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
int Test_N_VInv_Z(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
int Test_N_VAddConst_Z(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
int Test_N_VDotProd_Z(N_Vector X, N_Vector Y, sunindextype local_length, int myid);
int Test_N_VMaxNorm_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VWrmsNorm_Z(N_Vector X, N_Vector W, sunindextype local_length, int myid);
int Test_N_VWrmsNormMask_Z(N_Vector X, N_Vector W, N_Vector ID,
                         sunindextype local_length, int myid);
int Test_N_VMin_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VWL2Norm_Z(N_Vector X, N_Vector W, sunindextype local_length, int myid);
int Test_N_VL1Norm_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VCompare_Z(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
int Test_N_VInvTest_Z(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
int Test_N_VConstrMask_Z(N_Vector C, N_Vector X, N_Vector M,
                       sunindextype local_length, int myid);
int Test_N_VMinQuotient_Z(N_Vector NUM, N_Vector DENOM, sunindextype local_length,
                        int myid);

/* Fused vector operation tests */
int Test_N_VLinearCombination_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VScaleAddMulti_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VDotProdMulti_Z(N_Vector X, sunindextype local_length, int myid);

/* Vector array operation tests */
int Test_N_VLinearSumVectorArray_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VScaleVectorArray_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VConstVectorArray_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VWrmsNormVectorArray_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VWrmsNormMaskVectorArray_Z(N_Vector X, sunindextype local_length,
                                    int myid);
int Test_N_VScaleAddMultiVectorArray_Z(N_Vector X, sunindextype local_length,
                                     int myid);
int Test_N_VLinearCombinationVectorArray_Z(N_Vector X, sunindextype local_length,
                                         int myid);

/* Local reduction operation tests */
int Test_N_VDotProdLocal_Z(N_Vector X, N_Vector Y, sunindextype local_length,
                         int myid);
int Test_N_VMaxNormLocal_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VMinLocal_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VL1NormLocal_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VWSqrSumLocal_Z(N_Vector X, N_Vector W, sunindextype local_length,
                         int myid);
int Test_N_VWSqrSumMaskLocal_Z(N_Vector X, N_Vector W, N_Vector ID,
                             sunindextype local_length, int myid);
int Test_N_VInvTestLocal_Z(N_Vector X, N_Vector Z, sunindextype local_length,
                         int myid);
int Test_N_VConstrMaskLocal_Z(N_Vector C, N_Vector X, N_Vector M,
                            sunindextype local_length, int myid);
int Test_N_VMinQuotientLocal_Z(N_Vector NUM, N_Vector DENOM,
                             sunindextype local_length, int myid);

/* Single buffer reduction tests */
int Test_N_VDotProdMultiLocal_Z(N_Vector X, sunindextype local_length, int myid);
int Test_N_VDotProdMultiAllReduce_Z(N_Vector X, sunindextype local_length,
                                  int myid);

/* XBraid interface operations */
int Test_N_VBufSize_Z(N_Vector x, sunindextype local_length, int myid);
int Test_N_VBufPack_Z(N_Vector x, sunindextype local_length, int myid);
int Test_N_VBufUnpack_Z(N_Vector x, sunindextype local_length, int myid);


#ifdef __cplusplus
}
#endif

#endif
