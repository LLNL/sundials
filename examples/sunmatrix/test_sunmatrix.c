/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * These test functions are designed to check a SUNMatrix module
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <math.h> /* include isnan */
#include <stdio.h>
#include <stdlib.h>

#include "test_sunmatrix.h"

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif


/* private functions */
static double get_time();

int print_time = 0;

#define PRINT_TIME(format, time) if(print_time) printf(format, time)


/* ----------------------------------------------------------------------
 * SUNMatGetID Test
 * --------------------------------------------------------------------*/
int Test_SUNMatGetID(SUNMatrix A, SUNMatrix_ID sunid, int myid)
{
  double       start_time, stop_time;
  SUNMatrix_ID mysunid;

  start_time = get_time();   
  sunid = SUNMatGetID(A);
  stop_time = get_time();   

  if (sunid != mysunid) {
    printf(">>> FAILED test -- SUNMatGetID, Proc %d \n", myid);
    PRINT_TIME("    SUNMatGetID Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  } else if (myid == 0) {
    printf("    PASSED test -- SUNMatGetID \n");
    PRINT_TIME("    SUNMatGetID Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VClone Test
 *
 * NOTE: This routine depends on SUNMatCopy to check matrix data.
 * --------------------------------------------------------------------*/
int Test_SUNMatClone(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B;

  /* clone vector */
  start_time = get_time();   
  B = SUNMatClone(A);
  stop_time = get_time();   

  /* check cloned vector */
  if (B == NULL) {
    printf(">>> FAILED test -- SUNMatClone, Proc %d \n", myid);
    printf("    After SUNMatClone, B == NULL \n \n");
    return(1);
  } 

  /* check cloned vector data */
  if (!has_data(B)) {
    printf(">>> FAILED test -- SUNMatClone, Proc %d \n", myid);
    printf("    Matrix data == NULL \n \n");
    SUNMatDestroy(B);
    return(1);
  }    

  SUNMatCopy(B, A);
  failure = check_matrix(B, A);
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy, Proc %d \n", myid);
    printf("    Failed SUNMatCopy check \n \n");
    SUNMatDestroy(B);
    return(1);
  }    

  SUNMatDestroy(B); 

  if (myid == 0) {
    printf("    PASSED test -- N_VClone \n");
    PRINT_TIME("    SUNMatClone Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}



/* ----------------------------------------------------------------------
 * SUNMatZero Test
 * --------------------------------------------------------------------*/
int Test_SUNMatZero(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B;

  /* protect A */
  B = SUNMatClone(A);

  /* set matrix data to zero */
  start_time = get_time();
  failure = SUNMatZero(B);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNMatZero returned %d on Proc %d \n", 
           failure, myid);
    PRINT_TIME("    SUNMatZero Time: %22.15e \n \n", stop_time - start_time);
  }

  /* A data should be a vector of zeros */
  failure = check_matrix_entry(B, ZERO);

  if (failure) {
    printf(">>> FAILED test -- SUNMatZero, Proc %d \n", myid);
    PRINT_TIME("    SUNMatZero Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNMatZero \n");
    PRINT_TIME("    SUNMatZero Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(failure);
}


/* ----------------------------------------------------------------------
 * SUNMatCopy Test
 * --------------------------------------------------------------------*/
int Test_SUNMatCopy(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B; 

  B = SUNMatClone(A);

  /* copy matrix data */
  start_time = get_time();
  failure = SUNMatCopy(B, A);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d on Proc %d \n", 
           failure, myid);
    PRINT_TIME("    SUNMatCopy Time: %22.15e \n \n", stop_time - start_time);
  }

  /* check matrix entries */
  failure = check_matrix(B, A);

  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy, Proc %d \n", myid);
    PRINT_TIME("    SUNMatCopy Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VConst \n");
    PRINT_TIME("    SUNMatCopy Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(failure);
}



/* ----------------------------------------------------------------------
 * N_VLinearSum Tests: A = c * A + B 
 *
 * NOTE: Sparse matrices will need additional testing for possibly
 * different sparsity patterns
 * --------------------------------------------------------------------*/
int Test_SUNMatScaleAdd(SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B;

  /* protect A */
  B = SUNMatClone(A);
  SUNMatCopy(B,A);
  
  /* fill vector data */
  start_time = get_time(); 
  failure = SUNMatrixScaleAdd(NEG_ONE, B, B);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd returned %d on Proc %d \n", 
           failure, myid);
    PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);
  }
  
  /* check matrix entries */
  failure = check_matrix_entry(B, ZERO);

  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd, Proc %d \n", myid);
    PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNMatScaleAdd Case 1a \n");
    PRINT_TIME("    SUNMatScaleAdd Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(failure);
}


/* ----------------------------------------------------------------------
 * SUNMatScaleAddI Tests
 *
 * NOTE: Sparse matrices will need additional testing for possibly
 * different sparsity patterns
 * --------------------------------------------------------------------*/
int Test_SUNMatScaleAddI(SUNMatrix A, SUNMatrix I, int myid)
{
  int       failure;
  double    start_time, stop_time;
  SUNMatrix B;

  /* protect A */
  B = SUNMatClone(A);
  SUNMatCopy(B,I);
  
  /* fill vector data */
  start_time = get_time(); 
  failure = SUNMatScaleAddI(NEG_ONE, B);
  stop_time = get_time(); 
  
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI returned %d on Proc %d \n",
           failure, myid);
    PRINT_TIME("    SUNMatScaleAddI Time: %22.15e \n \n", stop_time - start_time);
  }

  /* check matrix */
  failure = check_matrix_entry(B, ZERO);

  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI, Proc %d \n", myid);
    PRINT_TIME("    SUNMatScaleAddI Time: %22.15e \n \n", stop_time - start_time);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNMatScaleAddI \n");
    PRINT_TIME("    SUNMatScaleAddI Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(failure);
}



/* ----------------------------------------------------------------------
 * SUNMatMatvec Test
 * --------------------------------------------------------------------*/
int Test_SUNMatMatvec(SUNMatix A, N_Vector x, N_Vector y, int myid)
{
  int      failure;
  double   start_time, stop_time;
  SUNMatrix B;
  N_Vector  z, w;

  /* protect A */
  B = SUNMatClone(A);
  SUNMatCopy(B,A);

  /* compute matrix vector product */
  SUNMatScaleAddI(THREE,B);

  z = N_VClone(y);
  w = N_VClone(y);

  start_time = get_time();
  SUNMatMatvec(B,x,z)
  stop_time = get_time(); 

  w = N_VLinearSum(THREE,y,ONE,x,w);
  
  failure = check_vector(w,z);

  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec, Proc %d \n", myid);
    PRINT_TIME("    SUNMatMatvec Time: %22.15e \n \n", stop_time - start_time);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNMatMatvec \n");
    PRINT_TIME("    SUNMatMatvec Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(failure);
}






/* ======================================================================
 * Private functions
 * ====================================================================*/

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#endif

void SetTiming(int onoff)
{
   print_time = onoff;

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;  
  clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
  base_time_tv_sec = spec.tv_sec;
#endif
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time()
{
#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;  
  clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
  double time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
  double time = 0;
#endif
  return time;
}


