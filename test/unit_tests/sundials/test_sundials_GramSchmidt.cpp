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
 * ------------------------------------------------------------------------
 * These test functions check some components of a complex-valued
 * SUNLINEARSOLVER module implementation (for more thorough tests,
 * see the main SUNDIALS repository, inside examples/sunlinsol/).

 * The solvers tested are SUNClassicalGS (Classical Gram-Schmidt) and 
 * SUNModifiedGS (Modified Gram-Schmidt).
 * ------------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sundials_iterative.h"
#include "sundials_iterative_impl.h" 


#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

SUNErrCode SUNClassicalGS(N_Vector* v, sunscalartype** h, int k, int p,
                          sunrealtype* new_vk_norm, sunscalartype* stemp,
                          N_Vector* vtemp);

SUNErrCode SUNModifiedGS(N_Vector* v, sunscalartype** h, int k, int p,
                         sunrealtype* new_vk_norm);

int main(int argc, char* argv[])
{
  N_Vector* V;
  N_Vector x;
  N_Vector* vtemp;
  sunscalartype* stemp;
  sunscalartype* vdata;
  SUNContext sunctx;
  sunscalartype** H;
  int k, l;
  sunrealtype vnorm;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }
  #if defined(SUNDIALS_SCALAR_TYPE_REAL)
    /* Create vectors */
    x = N_VNew(3, sunctx);
    V = N_VCloneVectorArray(3, x);
    N_Vector CheckV =  N_VClone(x);

    H = (sunscalartype**)malloc((3+1) * sizeof(sunscalartype*)); 
     for (k = 0; k <= 3; k++)
     {
       H[k] = NULL;
       H[k] = (sunscalartype*)malloc((3) * sizeof(sunscalartype)); 
     }

    vtemp = (N_Vector*)malloc((3) * sizeof(N_Vector)); 
    stemp = (sunscalartype*)malloc((3) * sizeof(sunscalartype)); 

    /* set up matrix */
    vdata = N_VGetArrayPointer(V[0]);
    vdata[0] = SUN_RCONST(12.0); 
    vdata[1] = 6.0;
    vdata[2] = -4.0;
    vdata = N_VGetArrayPointer(V[1]);
    vdata[0] = -51.0;
    vdata[1] = 167.0;
    vdata[2] = 24.0;
    vdata = N_VGetArrayPointer(V[2]);
    vdata[0] = 4.0;
    vdata[1] = -68.0;
    vdata[2] = -41.0;

    /* perform Gram-Schmidt process for all vectors in V */
    char functionName[100] = "ClassicalGS"; // default function name

    if (argc > 1) {
      strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
    }

    if (strcmp(functionName, "ClassicalGS") == 0) {
      // printf("Using Classical Gram Schmidt \n");
      cout << "Using Classical Gram Schmidt!" << "\n";
      for (k=0; k<3; k++){
        SUNClassicalGS(V, H, k, 3, &vnorm, stemp, vtemp);
        N_VScale(1.0/vnorm, V[k], V[k]);
      }
    }
    else if (strcmp(functionName, "ModifiedGS") == 0) {
      // printf("Using Modified Gram Schmidt \n");
      cout << "Using Modified Gram Schmidt!" << "\n";
      for (k=0; k<3; k++){
        SUNModifiedGS(V, H, k, 3, &vnorm);
        N_VScale(1.0/vnorm, V[k], V[k]);
      }
    }
    else {
      // printf("Incorrect function name, use: ClassicalGS or ModifiedGS. \nUsing default: ClassicalGS\n");
      cout << "Incorrect function name, use: ClassicalGS or ModifiedGS" << endl;
      cout << "Using default: ClassicalGS!" << endl;
      for (k=0; k<3; k++){
        SUNClassicalGS(V, H, k, 3, &vnorm, stemp, vtemp);// Default function
        N_VScale(1.0/vnorm, V[k], V[k]);
      }
    }

    /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
    sunrealtype tolerance = 1e-14;

    /* check dot product results */
    int unit_vectorsReal = 0;
    int unit_vectorsImag = 0;
    int orthogonalReal = 0;
    int orthogonalImag = 0;

    /* check dot product results */
    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        float vnorm = N_VDotProd(V[k],V[l]);
        if ((k==l) && (SUNabs(SUNabs(creal(vnorm))-SUN_RCONST(1.0)))>tolerance){unit_vectorsReal = 1;} //unit vectors
        if ((k==l) && (SUNabs(cimag(vnorm))>tolerance)){unit_vectorsImag = 1;}
        if ((k!=l) && (SUNabs(creal(vnorm))>tolerance)) {orthogonalReal = 1;}//orthogonal vectors
        if ((k!=l) && (SUNabs(cimag(vnorm))>tolerance)) {orthogonalImag = 1;}
      }
    }

  #elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
    /* Create vectors */
    x = N_VNew(5, sunctx);
    V = N_VCloneVectorArray(5, x);
    N_Vector CheckV =  N_VClone(x);

    H = (sunscalartype**)malloc((5+1) * sizeof(sunscalartype*)); 
    for (k = 0; k <= 5; k++) 
    {
      H[k] = NULL;
      H[k] = (sunscalartype*)malloc((5) * sizeof(sunscalartype)); 
    }

    vtemp = (N_Vector*)malloc((5) * sizeof(N_Vector)); 
    stemp = (sunscalartype*)malloc((5) * sizeof(sunscalartype)); 

    /* set up matrix */
    vdata = N_VGetArrayPointer(V[0]);
    vdata[0] = SUN_RCONST(1.0) + 1.0*I; 
    vdata[1] = 2.0 - 2.0*I;
    vdata[2] = 3.0 + 3.0*I;
    vdata[3] = 2.0 - 3.0*I;
    vdata[4] = 6.0 + 1.0*I;
    vdata = N_VGetArrayPointer(V[1]);
    vdata[0] = 2.0 - 1.0*I;
    vdata[1] = 1.0 + 1.0*I;
    vdata[2] = 4.0 - 3.0*I;
    vdata[3] = 1.0 - 2.0*I;
    vdata[4] = 3.0 + 5.0*I; 
    vdata = N_VGetArrayPointer(V[2]);
    vdata[0] = 3.0 + 2.0*I;
    vdata[1] = 1.0 - 1.0*I;
    vdata[2] = 5.0 + 4.0*I;
    vdata[3] = 4.0 + 3.0*I;
    vdata[4] = 5.0 - 1.0*I;
    vdata = N_VGetArrayPointer(V[3]);
    vdata[0] = 1.0 - 4.0*I;
    vdata[1] = 5.0 + 2.0*I;
    vdata[2] = 3.0 - 4.0*I;
    vdata[3] = 1.0 - 3.0*I;
    vdata[4] = 1.0 + 1.0*I;
    vdata = N_VGetArrayPointer(V[4]);
    vdata[0] = 4.0 + 5.0*I;
    vdata[1] = 6.0 + 3.0*I;
    vdata[2] = 2.0 - 1.0*I;
    vdata[3] = 2.0 + 5.0*I;
    vdata[4] = 3.0 - 5.0*I;

    /* perform Gram-Schmidt process for all vectors in V */
    char functionName[100] = "ClassicalGS"; // default function name

    if (argc > 1) {
      strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
    }

    if (strcmp(functionName, "ClassicalGS") == 0) {
      // printf("Using Classical Gram Schmidt \n");
      cout << "Using Classical Gram Schmidt!" << "\n";
      for (k=0; k<5; k++){
        SUNClassicalGS(V, H, k, 5, &vnorm, stemp, vtemp);
        N_VScale(1.0/vnorm, V[k], V[k]);
      }
    }
    else if (strcmp(functionName, "ModifiedGS") == 0) {
      // printf("Using Modified Gram Schmidt \n");
      cout << "Using Modified Gram Schmidt!" << "\n";
      for (k=0; k<5; k++){
        SUNModifiedGS(V, H, k, 5, &vnorm);
        N_VScale(1.0/vnorm, V[k], V[k]);
      }
    }
    else {
      // printf("Incorrect function name, use: ClassicalGS or ModifiedGS. \nUsing default: ClassicalGS\n");
      cout << "Incorrect function name, use: ClassicalGS or ModifiedGS." << endl;
      cout << "Using default: ClassicalGS!" << endl;
      for (k=0; k<5; k++){
        SUNClassicalGS(V, H, k, 5, &vnorm, stemp, vtemp);// Default function
        N_VScale(1.0/vnorm, V[k], V[k]);
      }
    }

    /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
    sunrealtype tolerance = 1e-14;

    /* check dot product results */
    int unit_vectorsReal = 0;
    int unit_vectorsImag = 0;
    int orthogonalReal = 0;
    int orthogonalImag = 0;

    /* check dot product results */
    for (k=0; k<5; k++) {
      for (l=0; l<5; l++) {
        float vnorm = N_VDotProd(V[k],V[l]);
        if ((k==l) && (SUNabs(SUNabs(creal(vnorm))-SUN_RCONST(1.0)))>tolerance){unit_vectorsReal = 1;} //unit vectors
        if ((k==l) && (SUNabs(cimag(vnorm))>tolerance)){unit_vectorsImag = 1;}
        if ((k!=l) && (SUNabs(creal(vnorm))>tolerance)) {orthogonalReal = 1;}//orthogonal vectors
        if ((k!=l) && (SUNabs(cimag(vnorm))>tolerance)) {orthogonalImag = 1;}
      }
    }

  #else
  #error                                                                  \
    "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
  #endif

  /* Check if the columns of Q are orthonormal. */
  if ((orthogonalReal==0) && (orthogonalImag==0) && (unit_vectorsReal==0) && (unit_vectorsImag==0)){
    // printf("The columns of Q are orthonormal!\nTest Passed!\n");
    cout << "Test Passed!" << "\n";
  } 
  else {
    // printf("The columns of Q are not orthonormal!\nTest failed!\n");
    cout << "Test Failed!" << "\n";
    return 1;
  } 

  /* return with success */
  return 0;
}
