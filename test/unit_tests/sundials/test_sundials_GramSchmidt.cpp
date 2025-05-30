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

#include <iostream>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>

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
    std::cerr << "ERROR: SUNContext_Create failed\n";
    return (-1);
  }
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
  /* Create vectors */
  x               = N_VNew_Serial(3, sunctx);
  V               = N_VCloneVectorArray(3, x);
  N_Vector CheckV = N_VClone(x);

  H = (sunscalartype**)malloc((3 + 1) * sizeof(sunscalartype*));
  for (k = 0; k <= 3; k++)
  {
    H[k] = NULL;
    H[k] = (sunscalartype*)malloc((3) * sizeof(sunscalartype));
  }

  vtemp = (N_Vector*)malloc((3) * sizeof(N_Vector));
  stemp = (sunscalartype*)malloc((3) * sizeof(sunscalartype));

  /* set up matrix */
  vdata    = N_VGetArrayPointer(V[0]);
  vdata[0] = SUN_RCONST(12.0);
  vdata[1] = SUN_RCONST(6.0);
  vdata[2] = SUN_RCONST(-4.0);
  vdata    = N_VGetArrayPointer(V[1]);
  vdata[0] = SUN_RCONST(-51.0);
  vdata[1] = SUN_RCONST(167.0);
  vdata[2] = SUN_RCONST(24.0);
  vdata    = N_VGetArrayPointer(V[2]);
  vdata[0] = SUN_RCONST(4.0);
  vdata[1] = SUN_RCONST(-68.0);
  vdata[2] = SUN_RCONST(-41.0);

  /* perform Gram-Schmidt process for all vectors in V */
  char functionName[100] = "ClassicalGS"; // default function name

  if (argc > 1)
  {
    strcpy(functionName,
           argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "ClassicalGS") == 0)
  {
    std::cout << "Using Classical Gram Schmidt!"
              << "\n";
    for (k = 0; k < 3; k++)
    {
      SUNClassicalGS(V, H, k, 3, &vnorm, stemp, vtemp);
      N_VScale(SUN_RCONST(1.0) / vnorm, V[k], V[k]);
    }
  }
  else if (strcmp(functionName, "ModifiedGS") == 0)
  {
    std::cout << "Using Modified Gram Schmidt!"
              << "\n";
    for (k = 0; k < 3; k++)
    {
      SUNModifiedGS(V, H, k, 3, &vnorm);
      N_VScale(SUN_RCONST(1.0) / vnorm, V[k], V[k]);
    }
  }
  else
  {
    std::cout << "Incorrect function name, use: ClassicalGS or ModifiedGS"
              << std::endl;
    std::cout << "Using default: ClassicalGS!" << std::endl;
    for (k = 0; k < 3; k++)
    {
      SUNClassicalGS(V, H, k, 3, &vnorm, stemp, vtemp); // Default function
      N_VScale(SUN_RCONST(1.0) / vnorm, V[k], V[k]);
    }
  }

  /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
  sunrealtype tolerance = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* check dot product results */
  int unit_vectorsReal = 0;
  int unit_vectorsImag = 0;
  int orthogonalReal   = 0;
  int orthogonalImag   = 0;

  /* check dot product results */
  for (k = 0; k < 3; k++)
  {
    for (l = 0; l < 3; l++)
    {
      float vnorm = N_VDotProd(V[k], V[l]);
      if ((k == l) &&
          (SUNabs(SUNabs(SUN_REAL(vnorm)) - SUN_RCONST(1.0))) > tolerance)
      {
        unit_vectorsReal = 1;
      } //unit vectors
      if ((k == l) && (SUNabs(SUN_IMAG(vnorm)) > tolerance))
      {
        unit_vectorsImag = 1;
      }
      if ((k != l) && (SUNabs(SUN_REAL(vnorm)) > tolerance))
      {
        orthogonalReal = 1;
      } //orthogonal vectors
      if ((k != l) && (SUNabs(SUN_IMAG(vnorm)) > tolerance))
      {
        orthogonalImag = 1;
      }
    }
  }

#elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
  /* Create vectors */
  x               = N_VNew_Serial(5, sunctx);
  V               = N_VCloneVectorArray(5, x);
  N_Vector CheckV = N_VClone(x);

  H = (sunscalartype**)malloc((5 + 1) * sizeof(sunscalartype*));
  for (k = 0; k <= 5; k++)
  {
    H[k] = NULL;
    H[k] = (sunscalartype*)malloc((5) * sizeof(sunscalartype));
  }

  vtemp = (N_Vector*)malloc((5) * sizeof(N_Vector));
  stemp = (sunscalartype*)malloc((5) * sizeof(sunscalartype));

  /* set up matrix */
  vdata    = N_VGetArrayPointer(V[0]);
  vdata[0] = SUN_CCONST(1.0, 1.0);
  vdata[1] = SUN_CCONST(2.0, -2.0);
  vdata[2] = SUN_CCONST(3.0, 3.0);
  vdata[3] = SUN_CCONST(2.0, -3.0);
  vdata[4] = SUN_CCONST(6.0, 1.0);
  vdata    = N_VGetArrayPointer(V[1]);
  vdata[0] = SUN_CCONST(2.0, -1.0);
  vdata[1] = SUN_CCONST(1.0, 1.0);
  vdata[2] = SUN_CCONST(4.0, -3.0);
  vdata[3] = SUN_CCONST(1.0, -2.0);
  vdata[4] = SUN_CCONST(3.0, 5.0);
  vdata    = N_VGetArrayPointer(V[2]);
  vdata[0] = SUN_CCONST(3.0, 2.0);
  vdata[1] = SUN_CCONST(1.0, -1.0);
  vdata[2] = SUN_CCONST(5.0, 4.0);
  vdata[3] = SUN_CCONST(4.0, 3.0);
  vdata[4] = SUN_CCONST(5.0, -1.0);
  vdata    = N_VGetArrayPointer(V[3]);
  vdata[0] = SUN_CCONST(1.0, -4.0);
  vdata[1] = SUN_CCONST(5.0, 2.0);
  vdata[2] = SUN_CCONST(3.0, -4.0);
  vdata[3] = SUN_CCONST(1.0, -3.0);
  vdata[4] = SUN_CCONST(1.0, 1.0);
  vdata    = N_VGetArrayPointer(V[4]);
  vdata[0] = SUN_CCONST(4.0, 5.0);
  vdata[1] = SUN_CCONST(6.0, 3.0);
  vdata[2] = SUN_CCONST(2.0, -1.0);
  vdata[3] = SUN_CCONST(2.0, 5.0);
  vdata[4] = SUN_CCONST(3.0, -5.0);

  /* perform Gram-Schmidt process for all vectors in V */
  char functionName[100] = "ClassicalGS"; // default function name

  if (argc > 1)
  {
    strcpy(functionName,
           argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "ClassicalGS") == 0)
  {
    std::cout << "Using Classical Gram Schmidt!"
              << "\n";
    for (k = 0; k < 5; k++)
    {
      SUNClassicalGS(V, H, k, 5, &vnorm, stemp, vtemp);
      N_VScale(SUN_RCONST(1.0) / vnorm, V[k], V[k]);
    }
  }
  else if (strcmp(functionName, "ModifiedGS") == 0)
  {
    std::cout << "Using Modified Gram Schmidt!"
              << "\n";
    for (k = 0; k < 5; k++)
    {
      SUNModifiedGS(V, H, k, 5, &vnorm);
      N_VScale(SUN_RCONST(1.0) / vnorm, V[k], V[k]);
    }
  }
  else
  {
    std::cout << "Incorrect function name, use: ClassicalGS or ModifiedGS."
              << std::endl;
    std::cout << "Using default: ClassicalGS!" << std::endl;
    for (k = 0; k < 5; k++)
    {
      SUNClassicalGS(V, H, k, 5, &vnorm, stemp, vtemp); // Default function
      N_VScale(SUN_RCONST(1.0) / vnorm, V[k], V[k]);
    }
  }

  /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
  sunrealtype tolerance = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* check dot product results */
  int unit_vectorsReal = 0;
  int unit_vectorsImag = 0;
  int orthogonalReal   = 0;
  int orthogonalImag   = 0;

  /* check dot product results */
  for (k = 0; k < 5; k++)
  {
    for (l = 0; l < 5; l++)
    {
      float vnorm = N_VDotProd(V[k], V[l]);
      if ((k == l) &&
          (SUNabs(SUNabs(SUN_REAL(vnorm)) - SUN_RCONST(1.0))) > tolerance)
      {
        unit_vectorsReal = 1;
      } //unit vectors
      if ((k == l) && (SUNabs(SUN_IMAG(vnorm)) > tolerance))
      {
        unit_vectorsImag = 1;
      }
      if ((k != l) && (SUNabs(SUN_REAL(vnorm)) > tolerance))
      {
        orthogonalReal = 1;
      } //orthogonal vectors
      if ((k != l) && (SUNabs(SUN_IMAG(vnorm)) > tolerance))
      {
        orthogonalImag = 1;
      }
    }
  }

#else
#error \
  "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
#endif

  /* Check if the columns of Q are orthonormal. */
  if ((orthogonalReal == 0) && (orthogonalImag == 0) &&
      (unit_vectorsReal == 0) && (unit_vectorsImag == 0))
  {
    std::cout << "Test Passed!"
              << "\n";
  }
  else
  {
    std::cout << "Test Failed!"
              << "\n";
    return 1;
  }

  /* return with success */
  return 0;
}
