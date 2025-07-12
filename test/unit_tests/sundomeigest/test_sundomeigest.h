/*
 * -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * test SUNDomEigEstimator module implementations.
 * -----------------------------------------------------------------
 */

#include <sundials/sundials_core.h>

/* define constants */
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Test function declarations */
int Test_SUNDomEigEst_SetATimes(SUNDomEigEstimator DEE, void* ATdata,
                               SUNATimesFn ATimes, int myid);
int Test_SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE,
                                 sunindextype max_iters, int myid);
int Test_SUNDomEigEst_SetNumPreProcess(SUNDomEigEstimator DEE, int numwarmups,
                                       int myid);
int Test_SUNDomEigEst_SetTol(SUNDomEigEstimator DEE, sunrealtype tol, int myid);
int Test_SUNDomEigEst_Initialize(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEigEst_PreProcess(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEigEst_ComputeHess(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                           sunrealtype* lambdaI, int myid);
int Test_SUNDomEigEst_GetCurRes(SUNDomEigEstimator DEE, sunrealtype* curres, int myid);
int Test_SUNDomEigEst_GetCurNumIters(SUNDomEigEstimator DEE, int* curniter, int myid);
int Test_SUNDomEigEst_GetMaxNumIters(SUNDomEigEstimator DEE, int* maxniter, int myid);
int Test_SUNDomEigEst_GetMinNumIters(SUNDomEigEstimator DEE, int* minniter, int myid);
int Test_SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE, long int* nATimes, int myid);
int Test_SUNDomEigEst_PrintStats(SUNDomEigEstimator DEE, int myid);

/* Timing function */
void SetTiming(int onoff);

#ifdef __cplusplus
}
#endif
