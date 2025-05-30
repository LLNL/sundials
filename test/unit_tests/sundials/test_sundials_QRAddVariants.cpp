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

 * The solvers tested are the different QR decomposition variants
 * used in Anderson Acceleration.
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

/* Declaration of SUNQRData structure (copied from src/sundials/sundials_iterative_impl.h) */
typedef struct _SUNQRData* SUNQRData;

struct _SUNQRData
{
  N_Vector vtemp;
  N_Vector vtemp2;
  sunscalartype* temp_array;
};

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

SUNErrCode SUNQRAdd_MGS(N_Vector* Q, sunscalartype* R, N_Vector df, int m,
                        int mMax, void* QRdata);

SUNErrCode SUNQRAdd_ICWY(N_Vector* Q, sunscalartype* R, N_Vector df, int m,
                         int mMax, void* QRdata);

SUNErrCode SUNQRAdd_CGS2(N_Vector* Q, sunscalartype* R, N_Vector df, int m,
                         int mMax, void* QRdata);

SUNErrCode SUNQRAdd_DCGS2(N_Vector* Q, sunscalartype* R, N_Vector df, int m,
                          int mMax, void* QRdata);

int main(int argc, char* argv[])
{
  N_Vector* Q;
  N_Vector x;
  SUNContext sunctx;
  int k, l, m, mMax;
  sunrealtype vnorm;
  sunscalartype* temp_array;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    std::cerr << "ERROR: SUNContext_Create failed\n";
    return (-1);
  }
  SUNQRData qrdata;
  qrdata = (SUNQRData)malloc(sizeof(*qrdata));

#if defined(SUNDIALS_SCALAR_TYPE_REAL)
  /* Create vectors */
  x                  = N_VNew_Serial(3, sunctx);
  Q                  = N_VCloneVectorArray(3, x);
  N_Vector vtemp     = N_VClone(x);
  N_Vector vtemp2    = N_VClone(x);
  N_Vector df_True   = N_VClone(x);
  N_Vector R_Approx  = N_VClone(x);
  N_Vector df_Approx = N_VClone(x);

  qrdata->vtemp = vtemp;
  qrdata->vtemp2 =
    vtemp2; //together with temp_array, they are used in all the other QRAdd variants except QRAdd_MGS

  /* this stores the elements of the correction matrix (square matrix) as a one column vector by stacking the columns together starting with the first column */
  qrdata->temp_array = (sunscalartype*)malloc((9) * sizeof(sunscalartype));

  m    = 2; //number of vectors already orthogonalised (and are othornormal)
  mMax = 3; //number of rows = number of columns of the matrix Q

  /* the vector to orthogonalise */
  sunscalartype* dfdata = N_VGetArrayPointer(df_True);
  dfdata[0]             = SUN_RCONST(4.0);
  dfdata[1]             = SUN_RCONST(-68.0);
  dfdata[2]             = SUN_RCONST(-41.0);

  /* set up matrix, last column is obtained from any of the QRAdd functions */
  sunscalartype* vdata = N_VGetArrayPointer(Q[0]);
  vdata[0]             = SUN_RCONST(6.0) / SUN_RCONST(7.0);
  vdata[1]             = SUN_RCONST(3.0) / SUN_RCONST(7.0);
  vdata[2]             = SUN_RCONST(-2.0) / SUN_RCONST(7.0);

  vdata    = N_VGetArrayPointer(Q[1]);
  vdata[0] = SUN_RCONST(-69.0) / SUN_RCONST(175.0);
  vdata[1] = SUN_RCONST(158.0) / SUN_RCONST(175.0);
  vdata[2] = SUN_RCONST(6.0) / SUN_RCONST(35.0);

  vdata    = N_VGetArrayPointer(Q[2]);
  vdata[0] = SUN_RCONST(0.0);
  vdata[1] = SUN_RCONST(0.0);
  vdata[2] = SUN_RCONST(0.0);

  /* upper triangular matrix R (stored column-wise), the last column is obtained from any of the QRAdd functions*/
  sunscalartype R[9] = {SUN_RCONST(14.0), SUN_RCONST(0.0),   SUN_RCONST(0.0),
                        SUN_RCONST(21.0), SUN_RCONST(175.0), SUN_RCONST(0.0),
                        SUN_RCONST(0.0),  SUN_RCONST(0.0),   SUN_RCONST(0.0)};

  char functionName[100] = "sunqradd_mgs"; // default function name

  if (argc > 1)
  {
    strcpy(functionName,
           argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "sunqradd_mgs") == 0)
  {
    SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata);
    std::cout << "Using SUNQRAdd_MGS!"
              << "\n";
  }
  else if (strcmp(functionName, "sunqradd_icwy") == 0)
  {
    SUNQRAdd_ICWY(Q, R, df_True, m, mMax, qrdata);
    std::cout << "Using SUNQRAdd_ICWY!"
              << "\n";
  }
  else if (strcmp(functionName, "sunqradd_cgs2") == 0)
  {
    SUNQRAdd_CGS2(Q, R, df_True, m, mMax, qrdata);
    std::cout << "Using SUNQRAdd_CGS2!"
              << "\n";
  }
  else if (strcmp(functionName, "sunqradd_dcgs2") == 0)
  {
    SUNQRAdd_DCGS2(Q, R, df_True, m, mMax, qrdata);
    std::cout << "Using SUNQRAdd_DCGS2!"
              << "\n";
  }
  else
  {
    std::cout << "Incorrect function name, use: sunqradd_mgs or sunqradd_icwy "
                 "or sunqradd_cgs2 or sunqradd_dcgs2"
              << std::endl;
    std::cout << "Using default: sunqradd_mgs!" << std::endl;
    SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata); // Default function
  }

  /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
  sunrealtype tolerance = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* check dot product results */
  int unit_vectorsReal = 0;
  int unit_vectorsImag = 0;
  int orthogonalReal   = 0;
  int orthogonalImag   = 0;
  int solnCheckReal    = 0;
  int solnCheckImag    = 0;

  for (k = 0; k < 3; k++)
  {
    for (l = 0; l < 3; l++)
    {
      float vnorm = N_VDotProd(Q[k], Q[l]);
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

  /* the last column in R */
  sunscalartype* Rdata = N_VGetArrayPointer(R_Approx);
  Rdata[0]             = R[6];
  Rdata[1]             = R[7];
  Rdata[2]             = R[8];

  /* use the last column in R to check if the product of the last column of Q and R gives df_True */
  sunscalartype* finalR = N_VGetArrayPointer(df_Approx);
  finalR[0]             = SUN_RCONST(0.0);
  finalR[1]             = SUN_RCONST(0.0);
  finalR[2]             = SUN_RCONST(0.0);

  /* multiply Q by the last column of R (the result) and the final answer should be df */
  N_VLinearCombination(3, Rdata, Q, df_Approx);
  for (l = 0; l < 3; l++)
  {
    if (SUNabs(SUN_REAL(dfdata[l]) - SUN_REAL(finalR[l])) > tolerance)
    {
      solnCheckReal = 1;
    }
    if (SUNabs(SUN_IMAG(dfdata[l]) - SUN_IMAG(finalR[l])) > tolerance)
    {
      solnCheckImag = 1;
    }
  }

#elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
  /* Create vectors */
  x                  = N_VNew_Serial(5, sunctx);
  Q                  = N_VCloneVectorArray(5, x);
  N_Vector vtemp     = N_VClone(x);
  N_Vector vtemp2    = N_VClone(x);
  N_Vector df_True   = N_VClone(x);
  N_Vector R_Approx  = N_VClone(x);
  N_Vector df_Approx = N_VClone(x);

  qrdata->vtemp = vtemp;
  qrdata->vtemp2 =
    vtemp2; //together with temp_array, they are used in all the other QRAdd variants except QRAdd_MGS

  /* this stores the elements of the correction matrix (square matrix) as a one column vector by stacking the columns together starting with the first column */
  qrdata->temp_array = (sunscalartype*)malloc((16) * sizeof(sunscalartype));

  m    = 4; //number of vectors already orthogonalised (and are othornormal)
  mMax = 5; //number of rows = number of columns of the matrix Q

  /* the vector to orthogonalise */
  sunscalartype* dfdata = N_VGetArrayPointer(df_True);
  dfdata[0]             = SUN_CCONST(4.0, 5.0);
  dfdata[1]             = SUN_CCONST(6.0, 3.0);
  dfdata[2]             = SUN_CCONST(2.0, -1.0);
  dfdata[3]             = SUN_CCONST(2.0, 5.0);
  dfdata[4]             = SUN_CCONST(3.0, -5.0);

  /* set up matrix, last column is obtained from any of the QRAdd functions */
  sunscalartype* vdata = N_VGetArrayPointer(Q[0]);
  vdata[0]             = SUN_CCONST(0.113227703414460, 0.226455406828919);
  vdata[1]             = SUN_CCONST(0.226455406828919, 0.339683110243379);
  vdata[2]             = SUN_CCONST(0.339683110243379, 0.113227703414460);
  vdata[3]             = SUN_CCONST(0.226455406828919, -0.339683110243379);
  vdata[4]             = SUN_CCONST(0.679366220486758, 0.113227703414460);

  vdata    = N_VGetArrayPointer(Q[1]);
  vdata[0] = SUN_CCONST(0.358047898868247, -0.209079065032553);
  vdata[1] = SUN_CCONST(0.449519989819989, 0.164649763713135);
  vdata[2] = SUN_CCONST(0.352820922242433, 0.044429301319417);
  vdata[3] = SUN_CCONST(-0.334526504052085, -0.007840464938720);
  vdata[4] = SUN_CCONST(-0.376342317058596, 0.467814408010337);

  vdata    = N_VGetArrayPointer(Q[2]);
  vdata[0] = SUN_CCONST(0.368417696619559, 0.108720463349240);
  vdata[1] = SUN_CCONST(0.382885110056019, -0.076920802132466);
  vdata[2] = SUN_CCONST(-0.108648842490643, 0.349438169091529);
  vdata[3] = SUN_CCONST(0.326877598633683, 0.632698664840043);
  vdata[4] = SUN_CCONST(0.056007511422335, -0.236062349933527);

  vdata    = N_VGetArrayPointer(Q[3]);
  vdata[0] = SUN_CCONST(-0.173120531596438, -0.317326783017719);
  vdata[1] = SUN_CCONST(0.305340355271806, 0.559154947299423);
  vdata[2] = SUN_CCONST(-0.270428892207880, -0.452932177178935);
  vdata[3] = SUN_CCONST(0.395829946721830, 0.018686126433033);
  vdata[4] = SUN_CCONST(-0.144400733195965, -0.085349715976197);

  vdata    = N_VGetArrayPointer(Q[4]);
  vdata[0] = SUN_CCONST(0.0, 0.0);
  vdata[1] = SUN_CCONST(0.0, 0.0);
  vdata[2] = SUN_CCONST(0.0, 0.0);
  vdata[3] = SUN_CCONST(0.0, 0.0);
  vdata[4] = SUN_CCONST(0.0, 0.0);

  /* upper triangular matrix R (stored column-wise), the last column is obtained from any of the QRAdd functions */
  sunscalartype R[25] = {SUN_CCONST(8.831760866327848, 0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_CCONST(7.586256128768794, 2.717464881947030),
                         SUN_RCONST(4.905517563326268),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_CCONST(5.887840577551898, 0.113227703414459),
                         SUN_CCONST(-0.606329288594409, -0.786659982184980),
                         SUN_CCONST(7.438685615532769, 0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_CCONST(3.623286509262707, -3.396831102433787),
                         SUN_CCONST(4.432476178690121, -2.524629710268074),
                         SUN_CCONST(-1.780852648997921, -2.366997755750347),
                         SUN_CCONST(4.851661421774982, 0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0),
                         SUN_RCONST(0.0)};

  /* perform QR decomposition using Gram-Schmidt process for df, enter any QRAdd function to use */
  char functionName[100] = "sunqradd_mgs"; // default function name

  if (argc > 1)
  {
    strcpy(functionName,
           argv[1]); //if a function name (2nd argument) is provided after executable name
  }

  if (strcmp(functionName, "sunqradd_mgs") == 0)
  {
    std::cout << "Using SUNQRAdd_MGS!"
              << "\n";
    SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata);
  }
  else if (strcmp(functionName, "sunqradd_icwy") == 0)
  {
    std::cout << "Using SUNQRAdd_ICWY!"
              << "\n";
    SUNQRAdd_ICWY(Q, R, df_True, m, mMax, qrdata);
  }
  else if (strcmp(functionName, "sunqradd_cgs2") == 0)
  {
    std::cout << "Using SUNQRAdd_CGS2!"
              << "\n";
    SUNQRAdd_CGS2(Q, R, df_True, m, mMax, qrdata);
  }
  else if (strcmp(functionName, "sunqradd_dcgs2") == 0)
  {
    std::cout << "Using SUNQRAdd_DCGS2!"
              << "\n";
    SUNQRAdd_DCGS2(Q, R, df_True, m, mMax, qrdata);
  }
  else
  {
    std::cout << "Incorrect function name, use: sunqradd_mgs or sunqradd_icwy "
                 "or sunqradd_cgs2 or sunqradd_dcgs2"
              << std::endl;
    std::cout << "Using default: sunqradd_mgs!" << std::endl;
    SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata); // Default function
  }

  /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
  sunrealtype tolerance = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* check dot product results */
  int unit_vectorsReal = 0;
  int unit_vectorsImag = 0;
  int orthogonalReal   = 0;
  int orthogonalImag   = 0;
  int solnCheckReal    = 0;
  int solnCheckImag    = 0;

  for (k = 0; k < 5; k++)
  {
    for (l = 0; l < 5; l++)
    {
      float vnorm = N_VDotProd(Q[k], Q[l]);
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

  /* the last column in R */
  sunscalartype* Rdata = N_VGetArrayPointer(R_Approx);
  Rdata[0]             = R[20];
  Rdata[1]             = R[21];
  Rdata[2]             = R[22];
  Rdata[3]             = R[23];
  Rdata[4]             = R[24];

  /* use the last column in R to check if the product of the last column of Q and R gives df_True */
  sunscalartype* finalR = N_VGetArrayPointer(df_Approx);
  finalR[0]             = SUN_RCONST(0.0);
  finalR[1]             = SUN_RCONST(0.0);
  finalR[2]             = SUN_RCONST(0.0);
  finalR[3]             = SUN_RCONST(0.0);
  finalR[4]             = SUN_RCONST(0.0);

  /* multiply Q by the last column of R (the result) and the final answer should be df */
  N_VLinearCombination(5, Rdata, Q, df_Approx);
  for (l = 0; l < 5; l++)
  {
    if (SUNabs(SUN_REAL(dfdata[l]) - SUN_REAL(finalR[l])) > tolerance)
    {
      solnCheckReal = 1;
    }
    if (SUNabs(SUN_IMAG(dfdata[l]) - SUN_IMAG(finalR[l])) > tolerance)
    {
      solnCheckImag = 1;
    }
  }

#else
#error \
  "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
#endif

  /* Check if the computed last columns of Q and R are correct. */
  if ((solnCheckReal == 0) && (solnCheckImag == 0) && (orthogonalReal == 0) &&
      (orthogonalImag == 0) && (unit_vectorsReal == 0) && (unit_vectorsImag == 0))
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
