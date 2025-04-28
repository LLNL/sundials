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
 * The matrix and vectors used here are complex-valued.
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
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }
  SUNQRData qrdata;
  qrdata = (SUNQRData)malloc(sizeof(*qrdata));

  #if defined(SUNDIALS_SCALAR_TYPE_REAL)
    /* Create vectors */
    x = N_VNew(3, sunctx);
    Q = N_VCloneVectorArray(3, x);
    N_Vector vtemp = N_VClone(x);
    N_Vector vtemp2 = N_VClone(x);
    N_Vector df_True = N_VClone(x);
    N_Vector R_Approx =  N_VClone(x);
    N_Vector df_Approx =  N_VClone(x);

    qrdata->vtemp = vtemp;
    qrdata->vtemp2 = vtemp2;//together with temp_array, they are used in all the other QRAdd variants except QRAdd_MGS

    /* this stores the elements of the correction matrix (square matrix) as a one column vector by stacking the columns together starting with the first column */
    qrdata->temp_array = (sunscalartype*)malloc((9) * sizeof(sunscalartype));

    m = 2;  //number of vectors already orthogonalised (and are othornormal)
    mMax = 3; //number of rows = number of columns of the matrix Q

    /* the vector to orthogonalise */
    sunscalartype* dfdata = N_VGetArrayPointer(df_True);
    dfdata[0] = 4.0;
    dfdata[1] = -68.0;
    dfdata[2] = -41.0;

    /* set up matrix, last column is obtained from any of the QRAdd functions */
    sunscalartype* vdata = N_VGetArrayPointer(Q[0]);
    vdata[0] = SUN_RCONST(6.0 / 7.0);
    vdata[1] = 3.0 / 7.0;
    vdata[2] = -2.0 / 7.0;

    vdata = N_VGetArrayPointer(Q[1]);
    vdata[0] = -69.0 / 175.0;
    vdata[1] = 158.0 / 175.0;
    vdata[2] = 6.0 / 35.0;

    vdata = N_VGetArrayPointer(Q[2]);
    vdata[0] = 0.0;
    vdata[1] = 0.0;
    vdata[2] = 0.0;

    /* upper trinagular matrix R, the last column is obtained from any of the QRAdd functions*/
    sunscalartype R[9] = {14.0, 0.0, 0.0,
                         21.0, 175.0, 0.0,
                         0.0, 0.0, 0.0 };//R matrix stored as a column vector by stacking columns together starting with the first column 
 

    char functionName[100] = "sunqradd_mgs"; // default function name

    if (argc > 1) {
      strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
    }

    if (strcmp(functionName, "sunqradd_mgs") == 0) {
      SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata);
      // printf("Using SUNQRAdd_MGS \n");
      cout << "Using SUNQRAdd_MGS!" << "\n";
    }
    else if (strcmp(functionName, "sunqradd_icwy") == 0) {
      SUNQRAdd_ICWY(Q, R, df_True, m, mMax, qrdata);
      // printf("Using SUNQRAdd_ICWY \n");
      cout << "Using SUNQRAdd_ICWY!" << "\n";
    }
    else if (strcmp(functionName, "sunqradd_cgs2") == 0) {
      SUNQRAdd_CGS2(Q, R, df_True, m, mMax, qrdata);
      // printf("Using SUNQRAdd_CGS2 \n");
      cout << "Using SUNQRAdd_CGS2!" << "\n";
      
    }
    else if (strcmp(functionName, "sunqradd_dcgs2") == 0) {
      SUNQRAdd_DCGS2(Q, R, df_True, m, mMax, qrdata);
      // printf("Using SUNQRAdd_DCGS2 \n");
      cout << "Using SUNQRAdd_DCGS2!" << "\n";
    }
    else {
      // printf("Incorrect function name, use: sunqradd_mgs or sunqradd_icwy or sunqradd_cgs2 or sunqradd_dcgs2. \nUsing default: sunqradd_mgs\n");
      cout << "Incorrect function name, use: sunqradd_mgs or sunqradd_icwy or sunqradd_cgs2 or sunqradd_dcgs2" << endl;
      cout << "Using default: sunqradd_mgs!" << endl;
      SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata); // Default function
    }

    /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
    sunrealtype tolerance = 1e-14;

    /* check dot product results */
    int unit_vectorsReal = 0;
    int unit_vectorsImag = 0;
    int orthogonalReal = 0;
    int orthogonalImag = 0;
    int solnCheckReal = 0;
    int solnCheckImag = 0;

    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        float vnorm = N_VDotProd(Q[k],Q[l]);
        if ((k==l) && (SUNabs(SUNabs(creal(vnorm))-SUN_RCONST(1.0)))>tolerance){unit_vectorsReal = 1;} //unit vectors
        if ((k==l) && (SUNabs(cimag(vnorm))>tolerance)){unit_vectorsImag = 1;}
        if ((k!=l) && (SUNabs(creal(vnorm))>tolerance)) {orthogonalReal = 1;}//orthogonal vectors
        if ((k!=l) && (SUNabs(cimag(vnorm))>tolerance)) {orthogonalImag = 1;}
      }
    }

    /* the last column in R */
    sunscalartype* Rdata = N_VGetArrayPointer(R_Approx);
    Rdata[0] = R[6];
    Rdata[1] = R[7];
    Rdata[2] = R[8];

    /* use the last column in R to check if the product of the last column of Q and R gives df_True */
    sunscalartype* finalR = N_VGetArrayPointer(df_Approx);
    finalR[0] = 0.0;
    finalR[1] = 0.0;
    finalR[2] = 0.0;

    /* multiply Q by the last column of R (the result) and the final answer should be df */
    N_VLinearCombination(3, Rdata, Q, df_Approx);
    for (l=0;l<3;l++) {
      if (SUNabs(creal(dfdata[l]) - creal(finalR[l]))>tolerance ){solnCheckReal = 1;}
      if (SUNabs(cimag(dfdata[l]) - cimag(finalR[l]))>tolerance ){solnCheckImag = 1;}
    }

  #elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
    /* Create vectors */
    x = N_VNew(5, sunctx);
    Q = N_VCloneVectorArray(5, x);
    N_Vector vtemp = N_VClone(x);
    N_Vector vtemp2 = N_VClone(x);
    N_Vector df_True = N_VClone(x);
    N_Vector R_Approx =  N_VClone(x);
    N_Vector df_Approx =  N_VClone(x);

    qrdata->vtemp = vtemp;
    qrdata->vtemp2 = vtemp2;//together with temp_array, they are used in all the other QRAdd variants except QRAdd_MGS

    /* this stores the elements of the correction matrix (square matrix) as a one column vector by stacking the columns together starting with the first column */
    qrdata->temp_array = (sunscalartype*)malloc((16) * sizeof(sunscalartype));

    m = 4;  //number of vectors already orthogonalised (and are othornormal)
    mMax = 5; //number of rows = number of columns of the matrix Q

    /* the vector to orthogonalise */
    sunscalartype* dfdata = N_VGetArrayPointer(df_True);
    dfdata[0] = 4.0+5.0*I;
    dfdata[1] = 6.0+3.0*I;
    dfdata[2] = 2.0-1.0*I;
    dfdata[3] = 2.0+5.0*I;
    dfdata[4] = 3.0-5.0*I;

    /* set up matrix, last column is obtained from any of the QRAdd functions */
    sunscalartype* vdata = N_VGetArrayPointer(Q[0]);
    vdata[0] = SUN_RCONST(0.113227703414460) + 0.226455406828919*I;
    vdata[1] = 0.226455406828919 + 0.339683110243379*I;
    vdata[2] = 0.339683110243379 + 0.113227703414460*I;
    vdata[3] = 0.226455406828919 - 0.339683110243379*I;
    vdata[4] = 0.679366220486758 + 0.113227703414460*I;

    vdata = N_VGetArrayPointer(Q[1]);
    vdata[0] = 0.358047898868247 - 0.209079065032553*I;
    vdata[1] = 0.449519989819989 + 0.164649763713135*I;
    vdata[2] = 0.352820922242433 + 0.044429301319417*I;
    vdata[3] = -0.334526504052085 - 0.007840464938720*I;
    vdata[4] = -0.376342317058596 + 0.467814408010337*I;

    vdata = N_VGetArrayPointer(Q[2]);
    vdata[0] = 0.368417696619559 + 0.108720463349240*I;
    vdata[1] = 0.382885110056019 - 0.076920802132466*I;
    vdata[2] = -0.108648842490643 + 0.349438169091529*I;
    vdata[3] = 0.326877598633683 + 0.632698664840043*I;
    vdata[4] = 0.056007511422335 - 0.236062349933527*I;

    vdata = N_VGetArrayPointer(Q[3]);
    vdata[0] = -0.173120531596438 - 0.317326783017719*I;
    vdata[1] = 0.305340355271806 + 0.559154947299423*I;
    vdata[2] = -0.270428892207880 - 0.452932177178935*I;
    vdata[3] = 0.395829946721830 + 0.018686126433033*I;
    vdata[4] = -0.144400733195965 - 0.085349715976197*I;

    vdata = N_VGetArrayPointer(Q[4]);
    vdata[0] = 0.0 + 0.0*I;
    vdata[1] = 0.0 + 0.0*I;
    vdata[2] = 0.0 + 0.0*I;
    vdata[3] = 0.0 + 0.0*I;
    vdata[4] = 0.0 + 0.0*I;

    /* upper trinagular matrix R, the last column is obtained from any of the QRAdd functions*/
    sunscalartype R[25] = {8.831760866327848+0.0*I, 0.0, 0.0, 0.0, 0.0,
                         7.586256128768794+2.717464881947030*I, 4.905517563326268-0.000000000000002*I, 0.0, 0.0, 0.0,
                         5.887840577551898+0.113227703414459*I, -0.606329288594409-0.786659982184980*I, 7.438685615532769+0.000000000000000*I, 0.0, 0.0,
                         3.623286509262707 - 3.396831102433787*I, 4.432476178690121 - 2.524629710268074*I, -1.780852648997921 - 2.366997755750347*I,  4.851661421774982 + 0.000000000000001*I, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0 };//R matrix stored as a column vector by stacking columns together starting with the first column


    /* perform QR decomposition using Gram-Schmidt process for df, enter any QRAdd function to use */                       
    char functionName[100] = "sunqradd_mgs"; // default function name

    if (argc > 1) {
      strcpy(functionName, argv[1]); //if a function name (2nd argument) is provided after executable name
    }

    if (strcmp(functionName, "sunqradd_mgs") == 0) {
      // printf("Using SUNQRAdd_MGS\n");
      cout << "Using SUNQRAdd_MGS!" << "\n";
      SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata);
    }
    else if (strcmp(functionName, "sunqradd_icwy") == 0) {
      // printf("Using SUNQRAdd_ICWY \n");
      cout << "Using SUNQRAdd_ICWY!" << "\n";
      SUNQRAdd_ICWY(Q, R, df_True, m, mMax, qrdata);
    }
    else if (strcmp(functionName, "sunqradd_cgs2") == 0) {
      // printf("Using SUNQRAdd_CGS2 \n");
      cout << "Using SUNQRAdd_CGS2!" << "\n";
      SUNQRAdd_CGS2(Q, R, df_True, m, mMax, qrdata);
    }
    else if (strcmp(functionName, "sunqradd_dcgs2") == 0) {
      // printf("Using SUNQRAdd_DCGS2 \n");
      cout << "Using SUNQRAdd_DCGS2!" << "\n";
      SUNQRAdd_DCGS2(Q, R, df_True, m, mMax, qrdata);
    }
    else {
      // printf("Incorrect function name, use: sunqradd_mgs or sunqradd_icwy or sunqradd_cgs2 or sunqradd_dcgs2. \nUsing default: sunqradd_mgs\n");
      cout << "Incorrect function name, use: sunqradd_mgs or sunqradd_icwy or sunqradd_cgs2 or sunqradd_dcgs2" << endl;
      cout << "Using default: sunqradd_mgs!" << endl;
      SUNQRAdd_MGS(Q, R, df_True, m, mMax, qrdata); // Default function
    }

    /* Threshold for orthogonal vectors in matrix Q and imaginary component for the norm of a column vector in Q */
    sunrealtype tolerance = 1e-14;

    /* check dot product results */
    int unit_vectorsReal = 0;
    int unit_vectorsImag = 0;
    int orthogonalReal = 0;
    int orthogonalImag = 0;
    int solnCheckReal = 0;
    int solnCheckImag = 0;

    for (k=0; k<5; k++) {
      for (l=0; l<5; l++) {
        float vnorm = N_VDotProd(Q[k],Q[l]);
        if ((k==l) && (SUNabs(SUNabs(creal(vnorm))-SUN_RCONST(1.0)))>tolerance){unit_vectorsReal = 1;} //unit vectors
        if ((k==l) && (SUNabs(cimag(vnorm))>tolerance)){unit_vectorsImag = 1;}
        if ((k!=l) && (SUNabs(creal(vnorm))>tolerance)) {orthogonalReal = 1;}//orthogonal vectors
        if ((k!=l) && (SUNabs(cimag(vnorm))>tolerance)) {orthogonalImag = 1;}
      }
    }

    /* the last column in R */
    sunscalartype* Rdata = N_VGetArrayPointer(R_Approx);
    Rdata[0] = R[20];
    Rdata[1] = R[21];
    Rdata[2] = R[22];
    Rdata[3] = R[23];
    Rdata[4] = R[24];

    /* use the last column in R to check if the product of the last column of Q and R gives df_True */
    sunscalartype* finalR = N_VGetArrayPointer(df_Approx);
    finalR[0] = 0.0;
    finalR[1] = 0.0;
    finalR[2] = 0.0;
    finalR[3] = 0.0;
    finalR[4] = 0.0;

    /* multiply Q by the last column of R (the result) and the final answer should be df */
    N_VLinearCombination(5, Rdata, Q, df_Approx);
    for (l=0;l<5;l++){
      if (SUNabs(creal(dfdata[l]) - creal(finalR[l]))>tolerance ){solnCheckReal = 1;}
      if (SUNabs(cimag(dfdata[l]) - cimag(finalR[l]))>tolerance ){solnCheckImag = 1;}
    }

  #else
  #error                                                                  \
    "SUNDIALS scalar type not defined, report to github.com/LLNL/sundials/issues"
  #endif

  /* Check if the computed last columns of Q and R are correct. */
  if ((solnCheckReal==0) && (solnCheckImag==0) && (orthogonalReal==0) && (orthogonalImag==0) && (unit_vectorsReal==0) && (unit_vectorsImag==0)) {
    // printf("Test passed!\n");
    cout << "Test Passed!" << "\n";
  } 
  else {
    // printf("Test failed!\n");
    cout << "Test Failed!" << "\n";
    return 1;
  }

  /* return with success */
  return 0;
}
