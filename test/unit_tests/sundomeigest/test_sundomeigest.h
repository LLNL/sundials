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

#include <math.h>
#include <sundials/sundials_core.h>

/* define constants */
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif
/* Forward declarations for implementation specific utility functions */
// void sync_device(void);

/* Test function declarations */
int Test_SUNDomEigEstGetID(SUNDomEigEstimator DEE,
                           SUNDomEigEstimator_ID suntype, int myid);
int Test_SUNDomEigEstSetATimes(SUNDomEigEstimator DEE, void* ATdata,
                               SUNATimesFn ATimes, int myid);
int Test_SUNDomEigEstSetMaxPowerIter(SUNDomEigEstimator DEE,
                                     sunindextype max_powiter, int myid);
int Test_SUNDomEigEstSetNumPreProcess(SUNDomEigEstimator DEE, int numwarmups, int myid);
int Test_SUNDomEigEstInitialize(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEigEstPreProcess(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEigEstComputeHess(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEigEstimate(SUNDomEigEstimator DEE, suncomplextype* dom_eig, int myid);
int Test_SUNDomEigEstNumIters(SUNDomEigEstimator DEE, int* niter, int myid);
int Test_SUNDomEigEstRes(SUNDomEigEstimator DEE, sunrealtype* res, int myid);

/* Timing function */
void SetTiming(int onoff);

#ifdef __cplusplus
}
#endif
