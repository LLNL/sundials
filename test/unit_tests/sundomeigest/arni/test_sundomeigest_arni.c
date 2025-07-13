/* -----------------------------------------------------------------
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
 * These test functions check some components of Arnoldi Iteration
 * module implementation.
 * -----------------------------------------------------------------
 */

#include <nvector/nvector_serial.h>
#include <sundomeigest/sundomeigest_arni.h>
#include "../test_sundomeigest.h"

/* constants */
#define ZERO SUN_RCONST(0.0)

#define factor   SUN_RCONST(-100.0)
#define realpart SUN_RCONST(-30000.0)
#define imagpart SUN_RCONST(+40000.0)

/* user data structure */
typedef struct
{
  sunindextype N; /* problem size */
  N_Vector diag;  /* matrix diagonal */

  /* nondiagonal entries of the matrix that lead to the complex conjugate eigenvalues */
  sunrealtype real_part;
  sunrealtype imag_part;
} UserData;

/* private functions */
/*    matrix-vector product  */
int ATimes(void* ProbData, N_Vector v, N_Vector z);
/*    checks function return values  */
int check_flag(void* flagvalue, const char* funcname, int opt);

/* ----------------------------------------------------------------------
 * DomEig Module Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails              = 0;     /* counter for test failures  */
  int passfail           = 0;     /* overall pass/fail flag     */
  SUNDomEigEstimator DEE = NULL;  /* domeig estimator object    */
  UserData ProbData;              /* problem data structure     */
  int numwarmups;                 /* number of the preprocessing warmups */
  int max_iters;                  /* max power iteration        */
  int krydim;                     /* Krylov subspace dimension  */
  int curniter;                   /* cur. number of iterations  */
  int maxniter;                   /* max. number of iterations  */
  int minniter;                   /* min. number of iterations  */
  long int nATimes;               /* number of ATimes calls     */
  int print_timing;               /* timing output flag         */
  sunrealtype curres;             /* current residual           */
  sunrealtype lambdaR, lambdaI;   /* computed domeig parts      */
  sunrealtype tlambdaR, tlambdaI; /* true domeig parts          */
  SUNContext sunctx;
  sunrealtype rel_tol = SUN_RCONST(1.0e-2); /* relative tol for pass/fail */
  sunrealtype rel_error;
  N_Vector q; /* random initial eigenvector */

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }

  /* check inputs: local problem size, Krylov dimension, preprocessing items, timing flag */
  if (argc < 5)
  {
    printf("ERROR: FOUR (4) Inputs required:\n");
    printf("  Problem size should be >= 2\n");
    printf("  Krylov subspace dimension should be >0\n");
    printf("  Number of preprocessing should be >= 0\n");
    printf("  Include timers for calculation (0=off, 1=on)\n");
    return 1;
  }
  ProbData.N = (sunindextype)atol(argv[1]);
  if (ProbData.N <= 0)
  {
    printf("ERROR: Problem size must be a positive integer\n");
    return 1;
  }
  krydim = atoi(argv[2]);
  if (krydim <= 0)
  {
    printf("ERROR: Krylov subspace dimension must be a positive integer\n");
    return 1;
  }
  numwarmups = atoi(argv[3]);
  if (numwarmups < 0)
  {
    printf("ERROR: Number of preprocessing must be a nonnegative integer\n");
    return 1;
  }
  print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  printf("\nDomEig module test:\n");
  printf("  Problem size = %ld\n", (long int)ProbData.N);
  printf("  Krylov subspace dimension = %i\n", krydim);
  printf("  Number of preprocessing = %i\n", numwarmups);
  printf("  Timing output flag = %i\n\n", print_timing);

  /* Create vectors */
  ProbData.diag = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(ProbData.diag, "N_VNew_Serial", 0)) { return 1; }

  q = N_VClone(ProbData.diag);
  if (check_flag(q, "N_VClone", 0)) { return 1; }

  sunrealtype* qd = N_VGetArrayPointer(q);
  for (int i = 0; i < ProbData.N; i++)
  {
    qd[i] = (sunrealtype)rand() / (sunrealtype)RAND_MAX;
  }

  /* Fill matrix diagonal and problem data */
  // real diag is [3 4 5 ... N 0 0]*factor
  // 2x2 block matrix attached to the last two diagonals is
  // [ realpart   imagpart;
  // [-imagpart   realpart]
  // This setup allows two types of dominant eigenvalues (real and complex)
  // based on the "factor" and the problem dimension N.
  sunrealtype* v = N_VGetArrayPointer(ProbData.diag);
  for (int i = 0; i < ProbData.N - 2; i++) { v[i] = factor * (i + 3); }

  // Set the problem data corresponding to 2x2 block matrix
  ProbData.real_part = realpart;
  ProbData.imag_part = imagpart;

  /* Create Arnoldi Iteration Dominant Eigvalue Estimator (DEE)*/
  DEE = SUNDomEigEst_ArnI(q, krydim, sunctx);
  if (check_flag(DEE, "SUNDomEigEst_ArnI", 0)) { return 1; }

  fails += Test_SUNDomEigEst_SetATimes(DEE, &ProbData, ATimes, 0);
  // SUNDomEigEst_SetMaxIters is not an option for Arnoldi iteration.
  // It should return with SUN_SUCCESS
  max_iters = krydim;
  fails += Test_SUNDomEigEst_SetMaxIters(DEE, max_iters, 0);
  fails += Test_SUNDomEigEst_SetNumPreProcess(DEE, numwarmups, 0);
  fails += Test_SUNDomEigEst_SetTol(DEE, rel_tol, 0);
  fails += Test_SUNDomEigEst_Initialize(DEE, 0);
  fails += Test_SUNDomEigEst_PreProcess(DEE, 0);
  fails += Test_SUNDomEigEst_ComputeHess(DEE, 0);
  fails += Test_SUNDomEig_Estimate(DEE, &lambdaR, &lambdaI, 0);
  // SUNDomEigEst_GetCurRes, SUNDomEigEst_GetCurNumIters, SUNDomEigEst_GetMaxNumIters
  // and SUNDomEigEst_GetMinNumIters are not options for Arnoldi iteration.
  // They should return with 0.
  fails += Test_SUNDomEigEst_GetCurRes(DEE, &curres, 0);
  if (curres > SUN_SMALL_REAL)
  {
    printf("    >>> FAILED test -- SUNDomEigEst_GetCurRes return value\n");
    fails++;
  }
  fails += Test_SUNDomEigEst_GetCurNumIters(DEE, &curniter, 0);
  if (curniter != 0)
  {
    printf("    >>> FAILED test -- SUNDomEigEst_GetCurNumIters return value\n");
    fails++;
  }
  fails += Test_SUNDomEigEst_GetMaxNumIters(DEE, &maxniter, 0);
  if (maxniter != 0)
  {
    printf(
      "    >>> FAILED test -- SUNDomEigEst_GetMaxNumIters return  value\n");
    fails++;
  }
  fails += Test_SUNDomEigEst_GetMinNumIters(DEE, &minniter, 0);
  if (minniter != 0)
  {
    printf("    >>> FAILED test -- SUNDomEigEst_GetMinNumIters return value\n");
    fails++;
  }
  fails += Test_SUNDomEigEst_GetNumATimesCalls(DEE, &nATimes, 0);
  fails += Test_SUNDomEigEst_PrintStats(DEE, 0);

  if (fails)
  {
    printf("FAIL: SUNDomEigEst_ArnI module failed %i initialization tests\n\n",
           fails);
    return 1;
  }
  else
  {
    printf(
      "SUCCESS: SUNDomEigEst_ArnI module passed all initialization tests\n\n");
  }

  /* First check if the computed eigenvalue has a nonzero magnitute */
  sunrealtype norm_of_dom_eig = SUNRsqrt(lambdaR * lambdaR + lambdaI * lambdaI);
  if (norm_of_dom_eig < SUN_SMALL_REAL)
  {
    printf("FAIL: Dominant Eigenvalue Test Failed\n\n");
    return 1;
  }

  /* Identify the tlambdaR and tlambdaI based on given parameters*/
  if (SUNRsqrt(realpart * realpart + imagpart * imagpart) > -factor * ProbData.N)
  {
    /* Dominant eigenvalue corresponds to the 2x2 block matrix */
    tlambdaR = realpart;
    tlambdaI = imagpart;
  }
  else
  {
    /* Dominant eigenvalue corresponds to the maximum real value at the diagonal */
    tlambdaR = factor * ProbData.N;
    tlambdaI = ZERO;
  }

  printf("\ncomputed dominant eigenvalue = " SUN_FORMAT_G " + " SUN_FORMAT_G
         " i\n",
         lambdaR, lambdaI);
  printf("    true dominant eigenvalue = " SUN_FORMAT_G " + " SUN_FORMAT_G
         " i\n",
         tlambdaR, tlambdaI);

  /* Compare the estimated dom_eig with the tlambdaR and tlambdaI*/
  rel_error = SUNRsqrt((lambdaR - tlambdaR) * (lambdaR - tlambdaR) +
                       (lambdaI - tlambdaI) * (lambdaI - tlambdaI));

  rel_error /= norm_of_dom_eig;

  if (rel_error < rel_tol)
  {
    printf("\n\nPASS:   relative error = " SUN_FORMAT_G " \n\n", rel_error);
  }
  else
  {
    printf("\n\nFAIL:   relative error = " SUN_FORMAT_G " \n\n", rel_error);
    passfail += 1;
  }

  /* Free solver and vectors */
  N_VDestroy(ProbData.diag);
  SUNContext_Free(&sunctx);
  N_VDestroy(q);
  DEE->ops->free(DEE);

  return (passfail);
}

/* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*/

/* matrix-vector product  */
int ATimes(void* Data, N_Vector v_vec, N_Vector z_vec)
{
  /* local variables */
  sunrealtype *v, *z, *diag, real_part, imag_part;
  sunindextype i, N;
  UserData* ProbData;

  /* access user data structure and vector data */
  ProbData = (UserData*)Data;
  v        = N_VGetArrayPointer(v_vec);
  if (check_flag(v, "N_VGetArrayPointer", 0)) { return 1; }
  z = N_VGetArrayPointer(z_vec);
  if (check_flag(z, "N_VGetArrayPointer", 0)) { return 1; }
  N         = ProbData->N;
  real_part = ProbData->real_part;
  imag_part = ProbData->imag_part;
  diag      = N_VGetArrayPointer(ProbData->diag);
  if (check_flag(diag, "N_VGetArrayPointer", 0)) { return 1; }

  /* perform product on the diagonal part of the matrix */
  for (i = 0; i < N - 2; i++) { z[i] = diag[i] * v[i]; }

  /* perform product at the non-diagonal last two rows */
  z[N - 2] = v[N - 2] * real_part + v[N - 1] * imag_part;
  z[N - 1] = v[N - 1] * real_part - v[N - 2] * imag_part;
  /* return with success */
  return 0;
}

/* Check function return value based on "opt" input:
     0:  function allocates memory so check for NULL pointer
     1:  function returns a flag so check for flag != 0 */
int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return 1;
  }

  /* Check if flag != 0 */
  if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag != 0)
    {
      fprintf(stderr, "\nERROR: %s() failed with flag = %d\n\n", funcname,
              *errflag);
      return 1;
    }
  }

  return 0;
}
