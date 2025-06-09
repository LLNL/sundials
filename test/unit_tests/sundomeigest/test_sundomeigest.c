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
 * SUNDomEigEstGetType Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstGetType(SUNDomEigEstimator DEE,
                             SUNLinearSolver_Type suntype, int myid)
{
  double start_time, stop_time;
  SUNLinearSolver_Type mysuntype;

  start_time = get_time();
  mysuntype  = SUNDomEigEstGetType(DEE);
  // sync_device();
  stop_time = get_time();

  if (suntype != mysuntype)
  {
    printf(">>> FAILED test -- SUNDomEigEstGetType, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstGetType Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstGetType \n");
    PRINT_TIME("    SUNDomEigEstGetType Time: %22.15e \n \n",
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
  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstComputeHess Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstComputeHess(SUNDomEigEstimator DEE, int myid)
{
  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimate Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimate(SUNDomEigEstimator DEE, int myid) { return (0); }

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
