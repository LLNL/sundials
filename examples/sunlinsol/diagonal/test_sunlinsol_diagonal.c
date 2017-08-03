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
 * module implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_diagonal.h>
#include <sunmatrix/sunmatrix_diagonal.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"
#include "mpi.h"

/* constants */
#define TWO       RCONST(2.0)
#define FIVE      RCONST(5.0)
#define THOUSAND  RCONST(1000.0)

/* private functions */
/*    checks function return values  */
static int check_flag(void *flagvalue, const char *funcname, int opt);
/*    uniform random number generator in [0,1] */
static realtype urand();

/* global copy of Nloc (for check_vector routine) */
sunindextype local_length;

/* ----------------------------------------------------------------------
 * Diagonal Linear Solver Testing Routine
 * 
 * If this is run with 1 MPI task, our tests use the serial N_Vector 
 * module; otherwise we use the MPI-parallel N_Vector module.
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int             fails=0;          /* counter for test failures */
  SUNLinearSolver LS;               /* linear solver object      */
  SUNMatrix       A, I;             /* test matrices             */
  N_Vector        x, b;             /* test vectors              */
  int             print_timing;
  sunindextype    i, Nloc;
  MPI_Comm        comm;
  int             myid;
  int             nprocs;
  realtype        *xdata, *Adata;

  /* Set up MPI environment */
  fails = MPI_Init(&argc, &argv);
  if (check_flag(&fails, "MPI_Init", 1)) return 1;
  comm = MPI_COMM_WORLD;
  fails = MPI_Comm_size(comm, &nprocs);
  if (check_flag(&fails, "MPI_Comm_size", 1)) return 1;
  fails = MPI_Comm_rank(comm, &myid);
  if (check_flag(&fails, "MPI_Comm_rank", 1)) return 1;

  /* check inputs: local problem size, timing flag */
  if (argc < 3){
    printf("ERROR: THREE (2) Input required: local matrix rows, print timing \n");
    return(-1);
  }
  Nloc = atol(argv[1]);
  local_length = Nloc;
  if (Nloc <= 0) {
    printf("ERROR: local problem size must be a positive integer\n");
    return 1; 
  }
  
  print_timing = atoi(argv[2]);
  SetTiming(print_timing);

  if (myid == 0) {
    printf("\nDiagonal linear solver test:\n");
    printf("  nprocs = %i\n", nprocs);
    printf("  local/global problem sizes = %ld/%ld\n", (long int) Nloc,
           (long int) (nprocs * Nloc));
    printf("  timing output flag = %i\n\n", print_timing);
  }
  
  /* Create vectors */
  if (nprocs == 1) {
    x = N_VNew_Serial(Nloc);
    if (check_flag(x, "N_VNew_Serial", 0)) return 1;
    b = N_VNew_Serial(Nloc);
    if (check_flag(b, "N_VNew_Serial", 0)) return 1;
  } else {
    x = N_VNew_Parallel(comm, Nloc, nprocs*Nloc);
    if (check_flag(x, "N_VNew_Parallel", 0)) return 1;
    b = N_VNew_Parallel(comm, Nloc, nprocs*Nloc);
    if (check_flag(b, "N_VNew_Parallel", 0)) return 1;
  }
  A = SUNDiagonalMatrix(x);
  I = SUNDiagonalMatrix(x);

  /* Fill matrices and vectors */
  xdata = N_VGetArrayPointer(x);
  Adata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(A));
  for (i=0; i<Nloc; i++) {
    xdata[i] = ONE + urand();
    Adata[i] = TWO + urand();
  }
  N_VConst(ONE, SUNDiagonalMatrix_Diag(I));

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  
  /* Create diagonal linear solver */
  LS = SUNDiagonalLinearSolver(x, A);

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, myid);
  fails += Test_SUNLinSolSetup(LS, A, myid);
  fails += Test_SUNLinSolSolve(LS, A, x, b, RCONST(1.0e-15), myid);
 
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, myid);
  fails += Test_SUNLinSolLastFlag(LS, myid);
  fails += Test_SUNLinSolSpace(LS, myid);
  fails += Test_SUNLinSolNumIters(LS, myid);
  fails += Test_SUNLinSolResNorm(LS, myid);
  fails += Test_SUNLinSolSetATimes(LS, NULL, NULL, NULL, FALSE, myid);
  fails += Test_SUNLinSolSetPreconditioner(LS, NULL, NULL, NULL, FALSE, myid);
  fails += Test_SUNLinSolSetScalingVectors(LS, x, b, FALSE, myid);

  /* Print result */
  if (fails) 
    printf("FAIL: SUNSPGMR module, problem 6, failed %i tests\n\n", fails);
  else if (myid == 0)
    printf("SUCCESS: Diagonal SUNLinSol module passed all tests \n \n");

  /* Free solver and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(x);
  N_VDestroy(b);

  MPI_Finalize();
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
  for(i=0; i < local_length; i++)
    failure += FNEQ(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for(i=0; i < local_length; i++)
      if (FNEQ(xdata[i], ydata[i], tol) != 0)
        printf("  xdata[%ld] = %g != %g (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  }
  
  if (failure > ZERO)
    return(1);
  else
    return(0);
}
