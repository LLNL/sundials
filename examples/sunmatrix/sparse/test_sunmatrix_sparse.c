/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
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
 * This is the testing routine to check the SUNMatrix Sparse module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>


/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int       fails = 0;                /* counter for test failures  */
  int       mattype;                  /* matrix type (CSR or CSC)   */
  long int  matrows, matcols, matnnz; /* vector length              */
  N_Vector  x, y;                     /* test vectors               */
  realtype* xdata, ydata;             /* pointer to vector data     */
  SUNMatrix A, B, I;                  /* test matrices              */
  realtype* Adata, Bdata, Idata;      /* pointer to matrix data     */
  int       print_timing; 

  /* check input and set vector length */
  if (argc < 4){
    printf("ERROR: THREE (3) Input required: matrix rows, matrix cols, print timing \n");
    return(-1);
  }
  
  matrows = atol(argv[1]); 
  if (matrows < 5) {
    printf("ERROR: number of rows must be a positive integer greater than 4\n");
    return(-1); 
  }
  
  matcols = atol(argv[2]); 
  if (matcols < 5) {
    printf("ERROR: number of cols must be a positive integer greater than 4\n");
    return(-1); 
  }

  /* or add check around scale add I and skip if not square */
  if (matrows != matcols) {
    printf("ERROR: number of rows must equal the number of cols \n");
    return(-1); 
  }
  
  print_timing = atoi(argv[3]);
  SetTiming(print_timing);
  
  printf("\nRunning with matrix size, %ld by %ld \n \n", matrows, matcols);
    
  /* check creating sparse matrix from dense matrix */
  B = SUNDenseMatrix(5,5);
  
  Bdata = SUNDenseMatrix_Data(B);
  Bdata[2]  = ONE;
  Bdata[5]  = ONE;
  Bdata[9]  = ONE;
  Bdata[11] = ONE;
  Bdata[13] = ONE;
  Bdata[18] = ONE;
  Bdata[20] = ONE;
  Bdata[21] = ONE;

  /* Check CSR */
  A = SUNSparseFromDenseMatrix(B, 1e-15, CSR_MAT);
  fails += check_dense_conversion(A, B, 1e-15);

  if (fails) {
    printf("FAIL: SUNMatrix Sparse CSR conversion failed\n");
  }

  /* Check CSC */
  A = SUNSparseFromDenseMatrix(B, 1e-15, CSC_MAT);
  fails += check_dense_conversion(A, B, 1e-15);

  if (fails) {
    printf("FAIL: SUNMatrix Sparse CSC conversion failed\n");
  }

  /* check creating sparse matrix from banded matrix */

  /* CSR */

  /* CSC */


  /* Create vectors and matrices */
  x = N_VNew_Serial(matcols);
  y = N_VNew_Serial(matrows);
  A = SUNSparseMatrix(matrows, matcols, matnnz, mattype);
  I = SUNSparseMatrix(matrows, matcols, matcols, mattype);

  /* Fill matrices and vectors */
  Idata      = SUNSparseMatrix_Data(I);
  Iindexvals = SUNSparseMatrix_IndexValues(I);
  Iindexptrs = SUNSparseMatrix_IndexPointers(I);

  for(i=0, i < matrows; i++) {
    Idata[i] = ONE;
    Iindexvals[i] = i;
    Iindexptrs[i] = i;
  }

  xdata = N_VGetArrayPointer(x);
  for(i=0; i < matcols; i++) {
    xdata[i] = ONE / (i+1);
  }

  ydata = N_VGetArrayPointer(y);
  for(i=0; i < matrows; i++) {
    m = i;
    n = m + matcols - 1;
    ydata[i] = HALF*(n+1-m)*(n+m);
  }
    
  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_DENSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  fails += Test_SUNScaleAdd(A, 0);
  fails += Test_SUNScaleAddI(A, I, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);

  /* Free vectors and matrices */
  N_VDestroy_Serial(x);
  N_VDestroy_Serial(y);
  SUNMatDestroy(A);
  SUNMatDestroy(I);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  return(0);
}

/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_dense_conversion(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  realtype *Adata, *Bdata;

  Adata   = SUNSparseMatrix_Data(A);
  Aidxval = SUNSparseMatrix_IndexValues(A);
  Aidxptr = SUNSparseMatrix_IndexPointers(A);

  Bdata = SUNDenseMatrix_Data(B);
  
  for(j=0; j < matcols; j++) {
    for(i=0; i < matrows; i++) {
      if (Bdata[j*matrows + i] > ZERO) {
        for(k=Aidxptr[i]; k < Aidxptr[i+1]; i++) {
          if (Aidxval[k] == j) {
            failure += FNEQ(Bdata[j*matrows + i], Adata[Aidxptr[i]+k], tol);
            break;
          }
        }
      }
    }
  }

  return(failure);
}


int check_band_conversion(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  realtype *Adata, *Bdata;


  return(failure);
}


int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int      failure = 0;
  realtype *Adata, *Bdata;
  long int Aldata, Bldata;
  long int i;
  
  /* get data pointers */
  Adata = SUNSparseMatrix_Data(A);
  Bdata = SUNSparseMatrix_Data(B);

  /* get and check data lengths */
  Aldata = SUNSparseMatrix_LData(A);
  Bldata = SUNSparseMatrix_LData(B);

  if (Aldata != Bldata) {
    printf(">>> ERROR: check_matrix: Different data array lengths \n");
    return(1);
  }
  
  /* compare data */
  for(i=0; i < Aldata; i++){
    failure += FNEQ(Adata[i], Bdata[i], tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int      failure = 0;
  realtype *Adata;
  long int Aldata;
  long int i;
  
  /* get data pointer */
  Adata = SUNSparseMatrix_Data(A);

  /* compare data */
  Aldata = SUNSparseMatrix_LData(A);
  for(i=0; i < Aldata; i++){
    failure += FNEQ(Adata[i], val, tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector x, N_Vector y, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  long int xldata, yldata;
  long int i;

  /* get vector data */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);

  /* check data lengths */
  xldata = N_VGetLength_Serial(x);
  yldata = N_VGetLength_Serial(y);

  if (xldata != yldata) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++){
    failure += FNEQ(xdata[i], ydata[i], tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  realtype *Adata = SUNSparseMatrix_Data(A);
  if (Adata == NULL)
    return FALSE;
  else
    return TRUE;
}
