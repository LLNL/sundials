/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * This is the header file for the dense implementation of the 
 * SUNMATRIX module.
 * 
 * Part I contains declarations specific to the dense implementation
 * of the supplied SUNMATRIX module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNMatrix type without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor 
 * SUNMatrixNew_Dense as well as implementation-specific prototypes 
 * for various useful matrix operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_DENSE_H
#define _SUNMATRIX_DENSE_H

#include <sundials/sundials_matrix.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Dense implementation of SUNMatrix
 *
 * The dense implementation of the SUNMatrix 'content' structure
 * contains:
 *   M     - number of rows
 *   N     - number of columns
 *   ldim  - leading dimension (ldim >= M)
 *   data  - pointer to a contiguous block of realtype variables
 *   ldata - length of the data array = ldim*N
 *   cols  - array of pointers. cols[j] points to the first element 
 *           of the j-th column of the matrix in the array data.
 * The elements of a dense matrix are stored columnwise (i.e. columns 
 * are stored one on top of the other in memory); i.e. if A is a 
 * SUNMatrix_Dense object, then the (i,j)th element of A (with 
 * 0 <= i < M and 0 <= j < N) is given by (A->data)[j*ldim+i].
 *
 * -----------------------------------------------------------------
 */
  
struct _SUNMatrixContent_Dense {
  long int M;
  long int N;
  long int ldim;
  realtype *data;
  long int ldata;
  realtype **cols;
};

typedef struct _SUNMatrixContent_Dense *SUNMatrixContent_Dense;

/*
 * -----------------------------------------------------------------
 * PART II: macros SM_CONTENT_D, SM_DATA_D, SM_ROWS_D, SM_COLUMNS_D, 
 *          SM_COLUMN_D, and SM_ELEMENT_D
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNMatrix A;
 * SUNMatrixContent_Dense A_cont;
 * realtype *A_col_j, *A_data, A_ij;
 * long int i, j, A_rows, A_cols;
 *
 * (1) SM_CONTENT_D
 *
 *     This macro gives access to the contents of the dense
 *     SUNMatrix
 *
 *     The assignment A_cont = SM_CONTENT_D(A) sets A_cont to be
 *     a pointer to the dense SUNMatrix content structure.
 *
 * (2) SM_DATA_D, SM_ROWS_D, SM_COLUMNS_D
 *
 *     These macros give access to the individual parts of
 *     the content structure of a dense SUNMatrix.
 *
 *     The assignment A_data = SM_DATA_D(A) sets A_data to be
 *     a pointer to the first component of A. 
 *
 *     The assignment A_rows = SM_ROWS_D(A) sets A_rows to be
 *     the number of rows in A.
 *
 *     The assignment A_cols = SM_COLUMNS_D(A) sets A_cols to be
 *     the number of columns in A.
 *
 * (3) SM_COLUMN_D and SM_ELEMENT_D
 *
 *     These macros give access to the individual columns and 
 *     elements of a dense SUNMatrix, respectively.  In the 
 *     following, the entries of a SUNMatrix are indexed (i,j) 
 *     where i=0,...,M-1 and j=0,...,N-1.
 *     
 *     The assignment A_col_j = SM_COLUMN_D(A,j) sets A_col_j to 
 *     be a pointer to the jth column of the M-by-N dense
 *     matrix A, 0 <= j < N.  After the assignment, A_col_j may 
 *     be treated as an array indexed from 0 to M-1.
 *     The (i,j)-th element of A is thus referenced by col_j[i].
 *
 *     The assignment A_ij = SM_ELEMENT_D(A,i,j) sets A_ij to 
 *     the value of the (i,j)th element of the dense M-by-N matrix
 *     A, 0 <= i < M ; 0 <= j < N.  Similarly, the assignment 
 *     SM_ELEMENT_D(A,i,j) = A_ij sets the value of A_ij into the 
 *     (i,j) location of the matrix A.
 *
 * -----------------------------------------------------------------
 */

#define SM_CONTENT_D(A)     ( (SUNMatrixContent_Dense)(A->content) )

#define SM_ROWS_D(A)        ( SM_CONTENT_D(A)->M )

#define SM_COLUMNS_D(A)     ( SM_CONTENT_D(A)->N )

#define SM_DATA_D(A)        ( SM_CONTENT_D(A)->data )

#define SM_COLUMN_D(A,j)    ( (SM_CONTENT_D(A)->cols)[j] )

#define SM_ELEMENT_D(A,i,j) ( (SM_CONTENT_D(A)->cols)[j][i] )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunmatrix_dense
 * 
 * CONSTRUCTORS:
 *    SUNMatrixNew_Dense
 * DESTRUCTORS:
 *    SUNMatrixDestroy_Dense
 * OTHER:
 *    SUNMatrixPrint_Dense
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function: SUNMatrixNew_Dense
 * -----------------------------------------------------------------
 * SUNMatrixNew_Dense creates and allocates memory for an M-by-N 
 * dense SUNMatrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNMatrix SUNMatrixNew_Dense(long int M, long int N);

/*
 * -----------------------------------------------------------------
 * Functions: SUNMatrixPrint_Dense
 * -----------------------------------------------------------------
 * This function prints the content of a M-by-N dense matrix A to
 * standard output as it would normally appear on paper.
 * It is intended as debugging tools with small values of M and N.
 * The elements are printed using the %g/%lg/%Lg option. 
 * A blank line is printed before and after the matrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void SUNMatrixPrint_Dense(SUNMatrix A);


/*
 * -----------------------------------------------------------------
 * dense implementations of various useful matrix operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNMatrix_ID SUNMatrixGetID_Dense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatrixClone_Dense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatrixDestroy_Dense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatrixZero_Dense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatrixScale_Dense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatrixCopy_Dense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatrixAddIdentity_Dense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatrixAdd_Dense(SUNMatrix A, SUNMatrix B);

  
#ifdef __cplusplus
}
#endif

#endif
