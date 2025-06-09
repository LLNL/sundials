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
 * These test functions check some components of an DomEig iteration
 * module implementation.
 * -----------------------------------------------------------------
 */

#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <sundials/sundials_domeigestimator.h>
#include <sundomeigest/sundomeigest_pi.h>

#include "../test_sundomeigest.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* constants */
#define ZERO     SUN_RCONST(0.0)

#define factor    (-100.0)
#define realpart  (-30000.0)
#define imagpart  (0.0)

/* user data structure */
typedef struct
{
  sunindextype N;         /* problem size */
  N_Vector d;             /* matrix diagonal */
  sunrealtype real_part;  /* nondiagonal entries of the matrix */
  sunrealtype imag_part;  /* nondiagonal entries of the matrix */
} UserData;

/* private functions */
/*    matrix-vector product  */
int ATimes(void* ProbData, N_Vector v, N_Vector z);
/*    checks function return values  */
int check_flag(void* flagvalue, const char* funcname, int opt);
/*    uniform random number generator in [0,1] */
int check_vector(N_Vector X, N_Vector Y, sunrealtype tol);

/* global copy of the problem size (for check_vector routine) */
sunindextype problem_size;


/* ----------------------------------------------------------------------
 * DomEig Module Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails    = 0;    /* counter for test failures */
  int passfail = 0;    /* overall pass/fail flag    */
  N_Vector q;          /* test vectors              */
  UserData ProbData;   /* problem data structure    */
  int maxl, failure, power_of_A;
  SUNContext sunctx;
  sunrealtype tolerans = 1.0e-2;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }

  /* check inputs: local problem size */
  if (argc < 4)
  {
    printf("ERROR: TWO (2) Inputs required:\n");
    printf("  Problem size should be >0\n");
    printf("  Krylov subspace dimension should be >0\n");
    printf("  Number of preprocessing should be >= 0\n");
    return 1;
  }
  ProbData.N   = (sunindextype)atol(argv[1]);
  problem_size = ProbData.N;
  if (ProbData.N <= 0)
  {
    printf("ERROR: Problem size must be a positive integer\n");
    return 1;
  }
  maxl = atoi(argv[2]);
  if (maxl <= 0)
  {
    printf(
      "ERROR: Krylov subspace dimension must be a positive integer\n");
    return 1;
  }
  power_of_A = atoi(argv[3]);
  if (power_of_A < 0)
  {
    printf(
      "ERROR: Number of preprocessing must be a nonnegative integer\n");
    return 1;
  }
  printf("\nDomEig module test:\n");
  printf("  Problem size = %ld\n", (long int)ProbData.N);
  printf("  Krylov subspace dimension = %i\n", maxl);
  printf("  Number of preprocessing = %i\n", power_of_A);

  /* Create vectors */
  q = N_VNew_Serial(ProbData.N, sunctx);
  if (check_flag(q, "N_VNew_Serial", 0)) { return 1; }
  ProbData.d = N_VClone(q);
  if (check_flag(ProbData.d, "N_VClone", 0)) { return 1; }

  /* Fill q vector with 1s */
  N_VConst(SUN_RCONST(1.0), q);

  /* Fill matrix diagonal and problem data */
  // real diag is factor*[3 4 5 ... N 0 0]
  // suncomplextype diag is [ realpart   imagpart;
  //                         -imagpart   realpart]
  sunrealtype *v = N_VGetArrayPointer(ProbData.d);
  int i, j;
  for(i = 0; i < ProbData.N-2; i++)
  {
    v[i] = factor*(i + 3);
  }

  ProbData.real_part = realpart;
  ProbData.imag_part = imagpart;

  /* Create DomEig estimator*/
  SUNDomEigEstimator DEE = NULL;

  DEE = SUNDomEigEst_PI(q, maxl, sunctx);
  if (check_flag(DEE, "SUNDomEigEst_PI", 0)) { return 1; }

  /* Set Atimes*/
  passfail = DEE->ops->setatimes(DEE, &ProbData, ATimes);
  if (check_flag(&passfail, "setatimes", 1)) { return 1; }

  /* Set the number of preprocessings */
  if(DEE->ops->setnumofperprocess != NULL)
  {
    passfail = DEE->ops->setnumofperprocess(DEE, power_of_A);
    if (check_flag(&passfail, "setnumofperprocess", 1)) { return 1; }
  }

  /* Initialize the estimator */
  passfail = DEE->ops->initialize(DEE);
  if (check_flag(&passfail, "initialize", 1)) { return 1; }

  /* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
  passfail = DEE->ops->preprocess(DEE);
  if (check_flag(&passfail, "preprocess", 1)) { return 1; }

  /* Estimate the dominant eigenvalue */
  suncomplextype dom_eig;
  passfail = DEE->ops->estimate(DEE, &dom_eig);
  if (check_flag(&passfail, "estimate", 1)) { return 1; }

  sunrealtype norm_of_dom_eig = SUNRsqrt(dom_eig.real * dom_eig.real + dom_eig.imag * dom_eig.imag);
  if(norm_of_dom_eig < SUN_SMALL_REAL) {
    printf("FAIL:   Dominant Eigenvalue Test Failed");
    return 1;
  }

  suncomplextype true_dom_eig;
  /* Identify true_dom_eig based on given parameters*/
  if(SUNRsqrt(realpart * realpart + imagpart * imagpart) > -1.0 * factor * ProbData.N)
  {
    true_dom_eig.real = realpart;
    true_dom_eig.imag = imagpart;
  }
  else
  {
    true_dom_eig.real = factor*ProbData.N;
    true_dom_eig.imag = ZERO;
  }

  printf("\ncomputed dominant eigenvalue = %20.4lf + %20.4lfi\n", dom_eig.real, dom_eig.imag);
  printf("    true dominant eigenvalue = %20.4lf + %20.4lfi\n", true_dom_eig.real, true_dom_eig.imag);

  /* Compare the estimated dom_eig with true_dom_eig*/
  sunrealtype error = SUNRsqrt((dom_eig.real - true_dom_eig.real) * (dom_eig.real - true_dom_eig.real)
                             + (dom_eig.imag - true_dom_eig.imag) * (dom_eig.imag - true_dom_eig.imag));

  error /= norm_of_dom_eig;

  if(error < tolerans) {
    printf("\n\nPASS:   relative error = %lf\n", error);
    return 0;
  }
  else {
    printf("\n\nFAIL:   relative error = %lf\n", error);
    return 0;
  }

  /* Free solver and vectors */
  N_VDestroy(q);
  N_VDestroy(ProbData.d);
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
  N = ProbData->N;
  real_part = ProbData->real_part;
  imag_part = ProbData->imag_part;
  diag = N_VGetArrayPointer(ProbData->d);
  if (check_flag(diag, "N_VGetArrayPointer", 0)) { return 1; }

  /* perform product on the diagonal part of the matrix */
  for (i = 0; i < N - 2; i++)
  {
    z[i] = diag[i] * v[i];
  }

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
