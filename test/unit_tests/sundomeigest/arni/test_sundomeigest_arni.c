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
#include <stdlib.h>
#include <sundials/sundials_math.h>

#include <sundomeigest/sundomeigest_arni.h>

#include "../test_sundomeigest.h"

/* constants */
#define ZERO SUN_RCONST(0.0)

#define factor   (-100.0)
#define realpart (-30000.0)
#define imagpart (+40000.0)

/* user data structure */
typedef struct
{
  sunindextype N;        /* problem size */
  N_Vector diag;         /* matrix diagonal */

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
  int fails    = 0;              /* counter for test failures */
  int passfail = 0;              /* overall pass/fail flag     */
  SUNDomEigEstimator DEE = NULL; /* domeig estimator object    */
  N_Vector q;                    /* test vectors               */
  UserData ProbData;             /* problem data structure     */
  int numwarmups;                /* Number of the preprocessing warmups */
  int max_powiter;               /* max power iteration        */
  int krydim;                    /* Krylov subspace dimension  */
  int niter;                     /* number of iterations       */
  int print_timing;              /* timing output flag         */
  sunrealtype res;               /* current residual           */
  suncomplextype dom_eig;        /* computed domeig value      */
  suncomplextype true_dom_eig;   /* true domeig value          */
  SUNContext sunctx;
  sunrealtype rel_tol = 1.0e-2;  /* relative tol for pass/fail */
  sunrealtype rel_error;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }

  /* check inputs: local problem size, timing flag */
  if (argc < 5)
  {
    printf("ERROR: TWO (4) Inputs required:\n");
    printf("  Problem size should be >0\n");
    printf("  Krylov subspace dimension should be >0\n");
    printf("  Number of preprocessing should be >= 0\n");
    return 1;
  }
  ProbData.N   = (sunindextype)atol(argv[1]);
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
    printf("ERROR: Number of preprocessing warmups must be a nonnegative integer\n");
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
  q = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(q, "N_VNew_Serial", 0)) { return 1; }
  ProbData.diag = N_VClone(q);
  if (check_flag(ProbData.diag, "N_VClone", 0)) { return 1; }

  /* Fill matrix diagonal and problem data */
  // real diag is [3 4 5 ... N 0 0]*factor
  // 2x2 block marix attached to the last two diagonals is
  // [ realpart   imagpart;
  // [-imagpart   realpart]
  // This setup allows two types of dominant eigenvalues (real and complex)
  // based on the "factor" and the problem dimension N.
  sunrealtype* v = N_VGetArrayPointer(ProbData.diag);
  int i;
  for (i = 0; i < ProbData.N - 2; i++)
  {
    v[i] = factor * (i + 3);
  }
  // Set the probem data corresponding to 2x2 block marix
  ProbData.real_part = realpart;
  ProbData.imag_part = imagpart;

  /* Create Arnoldi DomEig estimator*/
  DEE = SUNDomEigEst_ArnI(q, krydim, sunctx);
  if (check_flag(DEE, "SUNDomEigEst_ArnI", 0)) { return 1; }

  fails += Test_SUNDomEigEstGetID(DEE, SUNDSOMEIGESTIMATOR_ARNOLDI, 0);
  fails += Test_SUNDomEigEstSetATimes(DEE, &ProbData, ATimes, 0);
  fails += Test_SUNDomEigEstSetNumPreProcess(DEE, numwarmups, 0);
  // SUNDomEigEstSetMaxPowerIter is not an option for Arnoldi iteration.
  // It should return with SUN_SUCCESS
  max_powiter = krydim;
  fails += Test_SUNDomEigEstSetMaxPowerIter(DEE, max_powiter, 0);
  fails += Test_SUNDomEigEstInitialize(DEE, 0);
  fails += Test_SUNDomEigEstPreProcess(DEE, 0);
  fails += Test_SUNDomEigEstComputeHess(DEE, 0);
  fails += Test_SUNDomEigEstimate(DEE, &dom_eig, 0);
  // SUNDomEigEstNumIters and SUNDomEigEstRes are not options for
  // Arnoldi iteration. They should return with 0
  fails += Test_SUNDomEigEstNumIters(DEE, &niter, 0);
  if(niter != 0)
  {
    printf("    >>> FAILED test -- SUNDomEigEstNumIters return value\n");
    fails++;
  }
  fails += Test_SUNDomEigEstRes(DEE, &res, 0);
  if(res > SUN_SMALL_REAL)
  {
    printf("    >>> FAILED test -- Test_SUNDomEigEstRes return value\n");
    fails++;
  }

  if (fails)
  {
    printf("FAIL: SUNDomEigEst_ArnI module failed %i initialization tests\n\n",
           fails);
    return 1;
  }
  else
  {
    printf("SUCCESS: SUNDomEigEst_ArnI module passed all initialization tests\n\n");
  }

  /* First check if the computed eigenvalue has a nonzero magnitute */
  sunrealtype norm_of_dom_eig =
    SUNRsqrt(dom_eig.real * dom_eig.real + dom_eig.imag * dom_eig.imag);
  if (norm_of_dom_eig < SUN_SMALL_REAL)
  {
    printf("FAIL: Dominant Eigenvalue Test Failed\n\n");
    return 1;
  }

  /* Identify the true_dom_eig based on given parameters*/
  if (SUNRsqrt(realpart * realpart + imagpart * imagpart) > -factor * ProbData.N)
  {
    /* Dominant eigenvalue corresponds to the 2x2 block marix */
    true_dom_eig.real = realpart;
    true_dom_eig.imag = imagpart;
  }
  else
  {
    /* Dominant eigenvalue corresponds to the maximum real value at the diagonal */
    true_dom_eig.real = factor * ProbData.N;
    true_dom_eig.imag = ZERO;
  }

  printf("\ncomputed dominant eigenvalue = %20.4lf + %20.4lfi\n", dom_eig.real,
         dom_eig.imag);
  printf("    true dominant eigenvalue = %20.4lf + %20.4lfi\n",
         true_dom_eig.real, true_dom_eig.imag);

  /* Compare the estimated dom_eig with the true_dom_eig*/
  rel_error = SUNRsqrt(
    (dom_eig.real - true_dom_eig.real) * (dom_eig.real - true_dom_eig.real) +
    (dom_eig.imag - true_dom_eig.imag) * (dom_eig.imag - true_dom_eig.imag));

  rel_error /= norm_of_dom_eig;

  if (rel_error < rel_tol && passfail == 0)
  {
    printf("\n\nPASS:   relative error = %lf\n\n", rel_error);
  }
  else
  {
    printf("\n\nFAIL:   relative error = %lf\n\n", rel_error);
    passfail += 1;
  }

  /* Free solver and vectors */
  N_VDestroy(q);
  N_VDestroy(ProbData.diag);
  SUNContext_Free(&sunctx);
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
