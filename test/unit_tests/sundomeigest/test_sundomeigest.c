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
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_SetATimes returned %d on Proc %d "
           "\n",
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
int Test_SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE, long int max_iters,
                                  int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetMaxIters routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetMaxIters(DEE, max_iters);
  stop_time  = get_time();

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

int Test_SUNDomEigEst_SetNumPreProcess(SUNDomEigEstimator DEE, int num_warmups,
                                       int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEst_SetNumPreProcess routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetNumPreProcess(DEE, num_warmups);
  stop_time  = get_time();

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
 * SUNDomEigEst_SetRelTol Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_SetRelTol(SUNDomEigEstimator DEE, sunrealtype tol, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEst_SetRelTol routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEst_SetRelTol(DEE, tol);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_SetRelTol check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_SetRelTol Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_SetRelTol \n");
    PRINT_TIME("    SUNDomEigEst_SetRelTol Time: %22.15e \n \n",
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
  stop_time  = get_time();

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
 * SUNDomEig_Estimate Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                            sunrealtype* lambdaI, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEig_Estimate(DEE, lambdaR, lambdaI);
  stop_time  = get_time();

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
 * SUNDomEigEst_GetRes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetRes(SUNDomEigEstimator DEE, sunrealtype* cur_res,
                             int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetRes(DEE, cur_res);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetRes check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetRes Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetRes \n");
    PRINT_TIME("    SUNDomEigEst_GetRes Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetNumIters(SUNDomEigEstimator DEE, long int* curniter,
                                  int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetNumIters(DEE, curniter);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetNumIters check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_GetNumIters Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_GetNumIters \n");
    PRINT_TIME("    SUNDomEigEst_GetNumIters Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEst_GetNumATimesCalls Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE,
                                        long int* num_ATimes, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEst_GetNumATimesCalls(DEE, num_ATimes);
  stop_time  = get_time();
  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_GetNumATimesCalls check, Proc %d "
           "\n",
           myid);
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
 * SUNDomEigEst_Write Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEst_Write(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();

  printf("    SUNDomEigEst_Write: \n");
  printf("    --------------------------------------------\n");
  failure = SUNDomEigEst_Write(DEE, stdout);
  printf("    --------------------------------------------\n");
  stop_time = get_time();
  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEst_Write check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEst_Write Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEst_Write \n");
    PRINT_TIME("    SUNDomEigEst_Write Time: %22.15e \n \n",
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
