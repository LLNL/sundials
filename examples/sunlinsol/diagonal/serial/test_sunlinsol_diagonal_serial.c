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
 * This is the testing routine to check the SUNLinSol Diagonal
 * module implementation in serial.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_diagonal.h>
#include <sunmatrix/sunmatrix_diagonal.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"

/* constants */
#define TWO       RCONST(2.0)
#define FIVE      RCONST(5.0)
#define THOUSAND  RCONST(1000.0)

/* private functions */
/*    checks function return values  */
static int check_flag(void *flagvalue, const char *funcname, int opt);
/*    uniform random number generator in [0,1] */
static realtype urand();

/* global copy of diagonal length (for check_vector routine) */
sunindextype diag_length;

/* ----------------------------------------------------------------------
 * Diagonal Linear Solver Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int             fails=0;          /* counter for test failures */
  SUNLinearSolver LS;               /* linear solver object      */
  SUNMatrix       A, I;             /* test matrices             */
  N_Vector        x, b;             /* test vectors              */
  int             print_timing;
  sunindextype    i, N;
  realtype        *xdata, *Adata;

  /* check inputs: local problem size, timing flag */
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

  printf("\nDiagonal linear solver test:\n");
  printf("  problem size = %ld\n", (long int) N);
  printf("  timing output flag = %i\n\n", print_timing);
  
  /* Create vectors */
  x = N_VNew_Serial(N);
  if (check_flag(x, "N_VNew_Serial", 0)) return 1;
  b = N_VNew_Serial(N);
  if (check_flag(b, "N_VNew_Serial", 0)) return 1;

  A = SUNDiagonalMatrix(x);
  I = SUNDiagonalMatrix(x);

  /* Fill matrices and vectors */
  xdata = N_VGetArrayPointer(x);
  Adata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(A));
  for (i=0; i<N; i++) {
    xdata[i] = ONE + urand();
    Adata[i] = TWO + urand();
  }
  N_VConst(ONE, SUNDiagonalMatrix_Diag(I));

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  
  /* Create diagonal linear solver */
  LS = SUNDiagonalLinearSolver(x, A);

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, RCONST(1.0e-15), 0);
 
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);
  fails += Test_SUNLinSolNumIters(LS, 0);
  fails += Test_SUNLinSolResNorm(LS, 0);
  fails += Test_SUNLinSolSetATimes(LS, NULL, NULL, NULL, FALSE, 0);
  fails += Test_SUNLinSolSetPreconditioner(LS, NULL, NULL, NULL, FALSE, 0);
  fails += Test_SUNLinSolSetScalingVectors(LS, x, b, FALSE, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNSPGMR module, problem 6, failed %i tests\n\n", fails);
  } else {
    printf("SUCCESS: Diagonal SUNLinSol module passed all tests \n \n");
  }

  /* Free solver and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(x);
  N_VDestroy(b);

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
