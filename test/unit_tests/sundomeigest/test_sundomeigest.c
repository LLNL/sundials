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
 * These test functions are designed to check a SUNDomEigEstimator
 * module implementation.
 * -----------------------------------------------------------------
 */

#include "test_sundomeigest.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_domeigestimator.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

/* private functions */
static double get_time(void);

int print_time = 0;

#define PRINT_TIME(format, time) \
  if (print_time) printf(format, time)

/* ----------------------------------------------------------------------
 * SUNDomEigEstGetID Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstGetID(SUNDomEigEstimator DEE,
                             SUNDomEigEstimator_ID suntype, int myid)
{
  double start_time, stop_time;
  SUNDomEigEstimator_ID myestid;

  start_time = get_time();
  myestid  = SUNDomEigEstGetID(DEE);
  // sync_device();
  stop_time = get_time();

  if (suntype != myestid)
  {
    printf(">>> FAILED test -- SUNDomEigEstGetID, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstGetID Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstGetID \n");
    PRINT_TIME("    SUNDomEigEstGetID Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstSetATimes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstSetATimes(SUNDomEigEstimator DEE, void* ATdata,
                               SUNATimesFn ATimes, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetATimes routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstSetATimes(DEE, ATdata, ATimes);
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstSetATimes returned %d on Proc %d \n",
           failure, myid);
    PRINT_TIME("    SUNDomEigEstSetATimes Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstSetATimes \n");
    PRINT_TIME("    SUNDomEigEstSetATimes Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}


/* ----------------------------------------------------------------------
 * SUNDomEigEstSetMaxPowerIter Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstSetMaxPowerIter(SUNDomEigEstimator DEE, sunindextype max_powiter, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetMaxPowerIter routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstSetMaxPowerIter(DEE, max_powiter);
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstSetMaxPowerIter returned %d on Proc %d \n",
           failure, myid);
    PRINT_TIME("    SUNDomEigEstSetMaxPowerIter Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstSetMaxPowerIter \n");
    PRINT_TIME("    SUNDomEigEstSetMaxPowerIter Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

int Test_SUNDomEigEstSetNumPreProcess(SUNDomEigEstimator DEE, int numwarmups, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEstSetNumPreProcess routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstSetNumPreProcess(DEE, numwarmups);
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstSetNumPreProcess check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstSetNumPreProcess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstSetNumPreProcess \n");
    PRINT_TIME("    SUNDomEigEstSetNumPreProcess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstInitialize Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstInitialize(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstInitialize(DEE);
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstInitialize check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstInitialize Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstInitialize \n");
    PRINT_TIME("    SUNDomEigEstInitialize Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstPreProcess Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstPreProcess(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstPreProcess(DEE);
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstPreProcess check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstPreProcess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstPreProcess \n");
    PRINT_TIME("    SUNDomEigEstPreProcess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstComputeHess Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstComputeHess(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstComputeHess(DEE);
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstComputeHess check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstComputeHess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstComputeHess \n");
    PRINT_TIME("    SUNDomEigEstComputeHess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimate Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimate(SUNDomEigEstimator DEE, suncomplextype* dom_eig, int myid)
{
  int failure;
  suncomplextype estimated_dom_eig;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstimate(DEE, &estimated_dom_eig);
  *dom_eig = estimated_dom_eig;
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimate check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstimate Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimate \n");
    PRINT_TIME("    SUNDomEigEstimate Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstNumIters(SUNDomEigEstimator DEE, int* niter, int myid)
{
  int failure;
  int num_iters;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstNumIters(DEE, &num_iters);
  *niter     = num_iters;
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstNumIters check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstNumIters Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstNumIters \n");
    PRINT_TIME("    SUNDomEigEstNumIters Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstRes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstRes(SUNDomEigEstimator DEE, sunrealtype* res, int myid)
{
  int failure;
  sunrealtype residual;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstRes(DEE, &residual);
  *res       = residual;
  // sync_device();
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstRes check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstRes Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstRes \n");
    PRINT_TIME("    SUNDomEigEstRes Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ======================================================================
 * Private functions
 * ====================================================================*/

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#endif

void SetTiming(int onoff)
{
  print_time = onoff;

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  base_time_tv_sec = spec.tv_sec;
#endif
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time(void)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  double time = (double)(spec.tv_sec - base_time_tv_sec) +
                ((double)(spec.tv_nsec) / 1E9);
#else
  double time = 0;
#endif
  return time;
}
