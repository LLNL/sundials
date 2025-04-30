/* -----------------------------------------------------------------
 * Programmer(s): Sylvia Amihere @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * These test functions check some components of a complex-valued
 * SUNLINEARSOLVER module implementation (for more thorough tests,
 * see the main SUNDIALS repository, inside examples/sunlinsol/).

 * The solvers tested are the SUNQRFACT (QR factorization which uses Givens
 * rotations and QRsol (solves the system Rx = Q^{T}b or the Ax=b, for the solve x).
 * ----------------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sundials/sundials_math.h>
#include <sundials/sundials_iterative.h>


#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

int SUNQRfact(int n, sunscalartype** h, sunscalartype* q, int job);
int SUNQRsol(int n, sunscalartype** h, sunscalartype* q, sunscalartype* b);

int main(int argc, char* argv[])
{
  sunscalartype* givens;
  sunscalartype* yg;
  sunscalartype* exactSol;
  SUNContext sunctx;
  sunscalartype** H;
  int ndim = 5; //number of columns of the matrix
  int k,l, job, n, krydim;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    std::cerr << "ERROR: SUNContext_Create failed\n";
    return (-1);
  }

  /* Create vectors */

  H = (sunscalartype**)malloc((ndim+1) * sizeof(sunscalartype*));
  for (k = 0; k <= ndim; k++)
  {
    H[k] = NULL;
    H[k] = (sunscalartype*)malloc(ndim * sizeof(sunscalartype));
  }

  givens = (sunscalartype*)malloc(2 * (ndim) * sizeof(sunscalartype));
  yg = (sunscalartype*)malloc((ndim +1) * sizeof(sunscalartype));
  exactSol = (sunscalartype*)malloc((ndim) * sizeof(sunscalartype));


  /* set up matrix */
  #if defined(SUNDIALS_SCALAR_TYPE_REAL)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~ 6-by-5 real-valued matrix, ndim = 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  H[0][0] = SUN_RCONST(4.0);  H[0][1] = SUN_RCONST(3.0);  H[0][2] = SUN_RCONST(-1.0); H[0][3] = SUN_RCONST(2.0);  H[0][4] = SUN_RCONST(6.0);
  H[1][0] = SUN_RCONST(-1.0); H[1][1] = SUN_RCONST(9.0);  H[1][2] = SUN_RCONST(5.0);  H[1][3] = SUN_RCONST(2.0);  H[1][4] = SUN_RCONST(3.0);
  H[2][0] = SUN_RCONST(0.0);  H[2][1] = SUN_RCONST(-4.0); H[2][2] = SUN_RCONST(3.0);  H[2][3] = SUN_RCONST(1.0);  H[2][4] = SUN_RCONST(-7.0);
  H[3][0] = SUN_RCONST(0.0);  H[3][1] = SUN_RCONST(0.0);  H[3][2] = SUN_RCONST(1.0);  H[3][3] = SUN_RCONST(10.0); H[3][4] = SUN_RCONST(4.0);
  H[4][0] = SUN_RCONST(0.0);  H[4][1] = SUN_RCONST(0.0);  H[4][2] = SUN_RCONST(0.0);  H[4][3] = SUN_RCONST(3.0);  H[4][4] = SUN_RCONST(-2.0);
  H[5][0] = SUN_RCONST(0.0);  H[5][1] = SUN_RCONST(0.0);  H[5][2] = SUN_RCONST(0.0);  H[5][3] = SUN_RCONST(0.0);  H[5][4] = SUN_RCONST(8.0);

    /* right-hand side vector */
    yg[0] = SUN_RCONST(1.0);
    yg[1] = SUN_RCONST(-4.0);
    yg[2] = SUN_RCONST(3.0);
    yg[3] = SUN_RCONST(2.0);
    yg[4] = SUN_RCONST(8.0);
    yg[5] = SUN_RCONST(5.0);

    /* exact solution, x, from the system Ax = b */
    exactSol[0] = SUN_RCONST(0.367281678716080);
    exactSol[1] = SUN_RCONST(-0.715472518804680);
    exactSol[2] = SUN_RCONST(0.268391562142706);
    exactSol[3] = SUN_RCONST(0.322709872729248);
    exactSol[4] = SUN_RCONST(0.208661611483312);

  #elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~ 6-by-5 complex-valued matrix, ndim = 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    H[0][0] = SUN_CCONST(4.0, -1.0); H[0][1] = SUN_CCONST(3.0, 2.0);  H[0][2] = SUN_CCONST(-1.0, -1.0); H[0][3] = SUN_CCONST(2.0, 5.0);  H[0][4] = SUN_CCONST(6.0, -8.0);
    H[1][0] = SUN_CCONST(-1.0, 3.0); H[1][1] = SUN_CCONST(9.0, 8.0);  H[1][2] = SUN_CCONST(5.0, 2.0);   H[1][3] = SUN_CCONST(2.0, 4.0);  H[1][4] = SUN_CCONST(3.0, 7.0);
    H[2][0] = SUN_RCONST(0.0);       H[2][1] = SUN_CCONST(-4.0, 1.0); H[2][2] = SUN_CCONST(3.0, -2.0);  H[2][3] = SUN_CCONST(1.0, 1.0);  H[2][4] = SUN_CCONST(-7.0, -1.0);
    H[3][0] = SUN_RCONST(0.0);       H[3][1] = SUN_RCONST(0.0);       H[3][2] = SUN_CCONST(1.0, -2.0);  H[3][3] = SUN_CCONST(10.0, 2.0); H[3][4] = SUN_CCONST(4.0, 3.0);
    H[4][0] = SUN_RCONST(0.0);       H[4][1] = SUN_RCONST(0.0);       H[4][2] = SUN_RCONST(0.0);        H[4][3] = SUN_CCONST(3.0, -4.0); H[4][4] = SUN_CCONST(-2.0, 3.0);
    H[5][0] = SUN_RCONST(0.0);       H[5][1] = SUN_RCONST(0.0);       H[5][2] = SUN_RCONST(0.0);        H[5][3] = SUN_RCONST(0.0);       H[5][4] = SUN_CCONST(8.0, + 2.0);

    /* right-hand side vector */
    yg[0] = SUN_CCONST(1.0, 4.0);
    yg[1] = SUN_CCONST(-4.0, -1.0);
    yg[2] = SUN_CCONST(3.0, 7.0);
    yg[3] = SUN_CCONST(2.0, 4.0);
    yg[4] = SUN_CCONST(8.0, 3.0);
    yg[5] = SUN_CCONST(5.0, -6.0);

    /* exact solution, x, from the system Ax = b */
    exactSol[0] = SUN_CCONST(2.443458511228827, 3.513895963512003);
    exactSol[1] = SUN_CCONST(-0.404998816794420, -0.513911027696170);
    exactSol[2] = SUN_CCONST(1.118249944274408, -0.001711121348630);
    exactSol[3] = SUN_CCONST(0.079261343295451, 0.964160077719036);
    exactSol[4] = SUN_CCONST(0.150470525675297, -0.886295308913532);

  #else
  #error                                                                  \
    "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
  #endif


  /* perform QR factorization using Givens rotation */
  for (k=0; k<ndim; k++)
  {
    krydim = k+1;
    SUNQRfact(krydim, H, givens, k);
  }

  /* Solve the system Ax=b, for x, using Rx = Q^{T}b */
  SUNQRsol(ndim, H, givens, yg);

  /* Compare the approximated solution with the actual solution */
  sunrealtype tolerance = SUN_RCONST(1e-14);
  int checkSolReal = 0;
  int checkSolImag = 0;
  for (k=0; k<ndim; k++){
    if ((SUNabs(SUN_REAL(yg[k])-SUN_REAL(exactSol[k])))>tolerance){checkSolReal=1;}
    if ((SUNabs(SUN_IMAG(yg[k])-SUN_IMAG(exactSol[k])))>tolerance){checkSolImag=1;}
  }
  if (checkSolReal==0 && checkSolImag==0){
    std::cout << "Test Passed!" << "\n";
  }
  else {
    std::cout << "Test Failed!" << "\n";
    return 1;
  }

  /* return with success */
  return 0;
}
