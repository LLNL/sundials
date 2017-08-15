/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNMatrix Diagonal module 
 * implementation in serial.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_diagonal.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"


/* private functions */
/*    checks function return values  */
static int check_flag(void *flagvalue, const char *funcname, int opt);
/*    uniform random number generator in [0,1] */
static realtype urand();

/* global copy of diagonal length (for check_vector routines) */
sunindextype diag_length;

/* ----------------------------------------------------------------------
 * Diagonal SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int          fails=0;          /* counter for test failures  */
  N_Vector     x, y;             /* test vectors               */
  realtype     *xdata, *ydata;   /* pointers to vector data    */
  realtype     *Adata;           /* pointer to matrix data     */
  SUNMatrix    A, I;             /* test matrices              */
  int          print_timing;
  sunindextype i, N;

  /* check input and set vector length */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: matrix rows, print timing \n");
    return(-1);
  }
  N = atol(argv[1]);
  diag_length = N;
  if (N <= 0) {
    printf("ERROR: local problem size must be a positive integer\n");
    return 1; 
  }
  
  print_timing = atoi(argv[2]);
  SetTiming(print_timing);
  
  printf("\nDiagonal matrix test:\n");
  printf("  problem size = %ld\n", (long int) N);
  printf("  timing output flag = %i\n\n", print_timing);
  
  /* Create vectors and matrices */
  x = N_VNew_Serial(N);
  if (check_flag(x, "N_VNew_Serial", 0)) return 1;
  y = N_VNew_Serial(N);
  if (check_flag(y, "N_VNew_Serial", 0)) return 1;

  A = SUNDiagonalMatrix(x);
  I = SUNDiagonalMatrix(x);
  
  /* Fill matrices and vectors */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);
  Adata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(A));
  for (i=0; i<N; i++) {
    xdata[i] = ONE + urand();
    Adata[i] = ONE + urand();
    ydata[i] = xdata[i] * Adata[i];
  }
  N_VConst(ONE, SUNDiagonalMatrix_Diag(I));
    
  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_DIAGONAL, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  fails += Test_SUNMatScaleAdd(A, I, 0);
  fails += Test_SUNMatScaleAddI(A, I, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: Diagonal SUNMatrix module failed %i tests \n \n", fails);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  /* Free vectors and matrices */
  N_VDestroy(x);
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNMatDestroy(I);

  return(0);
}

/* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*/


/* uniform random number generator */
static realtype urand()
{
  return (random() / (pow(RCONST(2.0),RCONST(31.0)) - ONE));
}

/* Check function return value based on "opt" input:
     0:  function allocates memory so check for NULL pointer
     1:  function returns a flag so check for flag != 0 */
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if function returned NULL pointer - no memory allocated */
  if (opt==0 && flagvalue==NULL) {
    fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return 1; }

  /* Check if flag != 0 */
  if (opt==1) {
    errflag = (int *) flagvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return 1; }}

  return 0;
}


/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  realtype *Adata, *Bdata;
  sunindextype i;
  
  /* get data pointers */
  Adata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(A));
  Bdata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(B));

  /* compare data */
  for(i=0; i < diag_length; i++){
    failure += FNEQ(Adata[i], Bdata[i], tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  realtype *Adata;
  sunindextype i;
  
  /* get data pointer */
  Adata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(A));

  /* compare data */
  for(i=0; i < diag_length; i++)
    failure += FNEQ(Adata[i], val, tol);

  if (failure > ZERO) {
    printf("Check_matrix_entry failures:\n");
    for(i=0; i < diag_length; i++)
      if (FNEQ(Adata[i], val, tol) != 0)
        printf("  Adata[%ld] = %g != %g (err = %g)\n", (long int) i,
               Adata[i], val, SUNRabs(Adata[i]-val));
  }
  
  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector x, N_Vector y, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype i;

  /* get vector data */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);

  /* check vector data */
  for(i=0; i < diag_length; i++)
    failure += FNEQ(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for(i=0; i < diag_length; i++)
      if (FNEQ(xdata[i], ydata[i], tol) != 0)
        printf("  xdata[%ld] = %g != %g (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  }
  
  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  N_Vector Adiag  = SUNDiagonalMatrix_Diag(A);
  realtype *Adata = N_VGetArrayPointer(Adiag);
  if ((Adiag == NULL) || (Adata == NULL))
    return FALSE;
  else
    return TRUE;
}

booleantype is_square(SUNMatrix A)
{
  return TRUE;
}
