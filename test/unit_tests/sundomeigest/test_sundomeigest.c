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
 * SUNDomEigEst_SetATimes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_SetATimes(SUNDomEigEstimator DEE, void* ATdata,
                               SUNATimesFn ATimes, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetATimes routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetATimes(DEE, ATdata, ATimes);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_SetATimes returned %d on Proc %d \n",
           failure, myid);
    PRINT_TIME("    SUNDomEigEst_SetATimes Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_SetATimes \n");
    PRINT_TIME("    SUNDomEigEst_SetATimes Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_SetMaxIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE,
                                 int max_iters, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetMaxIters routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetMaxIters(DEE, max_iters);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_SetMaxIters returned %d on Proc "
           "%d \n",
           failure, myid);
    PRINT_TIME("    SUNDomEigEst_SetMaxIters Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_SetMaxIters \n");
    PRINT_TIME("    SUNDomEigEst_SetMaxIters Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

int Test_SUNDomEigEst_SetNumPreProcess(SUNDomEigEstimator DEE, int numwarmups,
                                      int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEst_SetNumPreProcess routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetNumPreProcess(DEE, numwarmups);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_SetNumPreProcess check, Proc %d \n",
           myid);
    PRINT_TIME("    SUNDomEigEst_SetNumPreProcess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_SetNumPreProcess \n");
    PRINT_TIME("    SUNDomEigEst_SetNumPreProcess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_SetTol Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_SetTol(SUNDomEigEstimator DEE, sunrealtype tol, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEst_SetTol routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetTol(DEE, tol);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_SetTol check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_SetTol Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_SetTol \n");
    PRINT_TIME("    SUNDomEigEst_SetTol Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_Initialize Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_Initialize(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_Initialize(DEE);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_Initialize check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_Initialize Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_Initialize \n");
    PRINT_TIME("    SUNDomEigEst_Initialize Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_PreProcess Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_PreProcess(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_PreProcess(DEE);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_PreProcess check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_PreProcess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_PreProcess \n");
    PRINT_TIME("    SUNDomEigEst_PreProcess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_ComputeHess Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_ComputeHess(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_ComputeHess(DEE);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_ComputeHess check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_ComputeHess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_ComputeHess \n");
    PRINT_TIME("    SUNDomEigEst_ComputeHess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEig_Estimate Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                           sunrealtype* lambdaI, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEig_Estimate(DEE, lambdaR, lambdaI);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEig_Estimate check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEig_Estimate Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEig_Estimate \n");
    PRINT_TIME("    SUNDomEig_Estimate Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetCurRes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetCurRes(SUNDomEigEstimator DEE, sunrealtype* curres, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetCurRes(DEE, curres);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetCurRes check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetCurRes Time: %22.15e \n \n", stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetCurRes \n");
    PRINT_TIME("    SUNDomEigEst_GetCurRes Time: %22.15e \n \n", stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetCurNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetCurNumIters(SUNDomEigEstimator DEE, int* curniter, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetCurNumIters(DEE, curniter);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetCurNumIters check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetCurNumIters Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetCurNumIters \n");
    PRINT_TIME("    SUNDomEigEst_GetCurNumIters Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetMaxNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetMaxNumIters(SUNDomEigEstimator DEE, int* maxniter, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetMaxNumIters(DEE, maxniter);
  stop_time = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetMaxNumIters check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetMaxNumIters Time: %22.15e \n \n", stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetMaxNumIters \n");
    PRINT_TIME("    SUNDomEigEst_GetMaxNumIters Time: %22.15e \n \n", stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetMinNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetMinNumIters(SUNDomEigEstimator DEE, int* minniter, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetMinNumIters(DEE, minniter);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetMinNumIters check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetMinNumIters Time: %22.15e \n \n", stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetMinNumIters \n");
    PRINT_TIME("    SUNDomEigEst_GetMinNumIters Time: %22.15e \n \n", stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetNumATimesCalls Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE, long int* nATimes, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetNumATimesCalls(DEE, nATimes);
  stop_time  = get_time();
  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetNumATimesCalls check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetNumATimesCalls Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetNumATimesCalls \n");
    PRINT_TIME("    SUNDomEigEst_GetNumATimesCalls Time: %22.15e \n \n",
               stop_time - start_time);
  }
  
  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_PrintStats Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_PrintStats(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();

  printf("    SUNDomEigEst_PrintStats: \n");
  printf("    --------------------------------------------\n");
  failure    = SUNDomEigEst_PrintStats(DEE, stdout);
  printf("    --------------------------------------------\n");
  stop_time = get_time();
  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_PrintStats check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_PrintStats Time: %22.15e \n \n", stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_PrintStats \n");
    PRINT_TIME("    SUNDomEigEst_PrintStats Time: %22.15e \n \n", stop_time - start_time);
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
