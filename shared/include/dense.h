/*******************************************************************
 *                                                                 *
 * File          : dense.h                                         *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for a generic DENSE linear solver       *
 * package.  The routines listed in this file all use type         *
 * DenseMat, defined below, for matrices.  These routines in turn  *
 * call routines in the smalldense.h/smalldense.c module, which    *
 * use the type realtype** for matrices.  This separation allows   *
 * for possible modifications in which matrices of type DenseMat   *
 * may not be stored contiguously, while small matrices can still  *
 * be treated with the routines in smalldense.                     *
 *                                                                 * 
 * Routines that work with the type DenseMat begin with "Dense".   *
 * The DenseAllocMat function allocates a dense matrix for use in  *
 * the other DenseMat routines listed in this file. Matrix         *
 * storage details are given in the documentation for the type     *
 * DenseMat. The DenseAllocPiv function allocates memory for       *
 * pivot information. The storage allocated by DenseAllocMat and   *
 * DenseAllocPiv is deallocated by the routines DenseFreeMat and   *
 * DenseFreePiv, respectively. The DenseFactor and DenseBacksolve  *
 * routines perform the actual solution of a dense linear system.  *
 *                                                                 *
 * Routines that work with realtype** begin with "den" (except for *
 * the factor and solve routines which are called gefa and gesl,   *
 * respectively). The underlying matrix storage is described in    *
 * the documentation for denalloc in smalldense.h                  *
 *                                                                 *
 *******************************************************************/
 
#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _dense_h
#define _dense_h


#include "sundialstypes.h"
#include "smalldense.h"

 
/******************************************************************
 *                                                                *
 * Type: DenseMat                                                 *
 *----------------------------------------------------------------*
 * The type DenseMat is defined to be a pointer to a structure    *
 * with a size and a data field. The size field indicates the     *
 * number of columns (== number of rows) of a dense matrix, while *
 * the data field is a two dimensional array used for component   *
 * storage. The elements of a dense matrix are stored columnwise  *
 * (i.e columns are stored one on top of the other in memory). If *
 * A is of type DenseMat, then the (i,j)th element of A (with     *
 * 0 <= i,j <= size-1) is given by the expression (A->data)[j][i] *
 * or by the expression (A->data)[0][j*n+i]. The macros below     *
 * allow a user to access efficiently individual matrix           *
 * elements without writing out explicit data structure           *
 * references and without knowing too much about the underlying   *
 * element storage. The only storage assumption needed is that    *
 * elements are stored columnwise and that a pointer to the jth   *
 * column of elements can be obtained via the DENSE_COL macro.    *
 * Users should use these macros whenever possible.               *
 *                                                                *
 ******************************************************************/

typedef struct _DenseMat {
  integertype size;
  realtype  **data;
} *DenseMat;
 

/* DenseMat accessor macros */

 
/******************************************************************
 *                                                                *
 * Macro : DENSE_ELEM                                             *
 * Usage : DENSE_ELEM(A,i,j) = a_ij;  OR                          *
 *         a_ij = DENSE_ELEM(A,i,j);                              *
 *----------------------------------------------------------------*
 * DENSE_ELEM(A,i,j) references the (i,j)th element of the N by N *
 * DenseMat A, 0 <= i,j <= N-1.                                   *
 *                                                                *
 ******************************************************************/

#define DENSE_ELEM(A,i,j) ((A->data)[j][i])


/******************************************************************
 *                                                                *
 * Macro : DENSE_COL                                              *
 * Usage : col_j = DENSE_COL(A,j);                                *
 *----------------------------------------------------------------*
 * DENSE_COL(A,j) references the jth column of the N by N         *
 * DenseMat A, 0 <= j <= N-1. The type of the expression          *
 * DENSE_COL(A,j) is realtype *. After the assignment in the usage*
 * above, col_j may be treated as an array indexed from 0 to N-1. *
 * The (i,j)th element of A is referenced by col_j[i].            *
 *                                                                *
 ******************************************************************/

#define DENSE_COL(A,j) ((A->data)[j])
 

/* Functions that use the DenseMat representation for a dense matrix */

 
/******************************************************************
 *                                                                *
 * Function : DenseAllocMat                                       *
 * Usage    : A = DenseAllocMat(N);                               *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * DenseAllocMat allocates memory for an N by N dense matrix and  *
 * returns the storage allocated (type DenseMat). DenseAllocMat   *
 * returns NULL if the request for matrix storage cannot be       *
 * satisfied. See the above documentation for the type DenseMat   *
 * for matrix storage details.                                    * 
 *                                                                * 
 ******************************************************************/

DenseMat DenseAllocMat(integertype N);


/******************************************************************
 *                                                                *
 * Function : DenseAllocPiv                                       *
 * Usage    : p = DenseAllocPiv(N);                               *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * DenseAllocPiv allocates memory for pivot information to be     *
 * filled in by the DenseFactor routine during the factorization  *
 * of an N by N dense matrix. The underlying type for pivot       *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, DenseAllocPiv returns NULL. *
 *                                                                * 
 ******************************************************************/

integertype *DenseAllocPiv(integertype N);


/******************************************************************
 *                                                                *
 * Function : DenseFactor                                         *
 * Usage    : ier = DenseFactor(A, p);                            *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * DenseFactor performs the LU factorization of the N by N dense  *
 * matrix A. This is done using standard Gaussian elimination     *
 * with partial pivoting.                                         *
 *                                                                *
 * A successful LU factorization leaves the matrix A and the      *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.  *
 *                                                                *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * DenseFactor returns 0 if successful. Otherwise it encountered  *
 * a zero diagonal element during the factorization. In this case *
 * it returns the column index (numbered from one) at which       *
 * it encountered the zero.                                       *
 *                                                                *
 ******************************************************************/

integertype DenseFactor(DenseMat A, integertype *p);


/******************************************************************
 *                                                                *
 * Function : DenseBacksolve                                      *
 * Usage    : DenseBacksolve(A, p, b);                            *
 *----------------------------------------------------------------*
 * DenseBacksolve solves the N-dimensional system A x = b using   *
 * the LU factorization in A and the pivot information in p       *
 * computed in DenseFactor. The solution x is returned in b. This *
 * routine cannot fail if the corresponding call to DenseFactor   *
 * did not fail.                                                  *
 *                                                                *
 ******************************************************************/

void DenseBacksolve(DenseMat A, integertype *p, realtype *b);


/******************************************************************
 *                                                                *
 * Function : DenseZero                                           *
 * Usage    : DenseZero(A);                                       *
 *----------------------------------------------------------------*
 * DenseZero sets all the elements of the N by N matrix A to 0.0. *
 *                                                                *
 ******************************************************************/

void DenseZero(DenseMat A);


/******************************************************************
 *                                                                *
 * Function : DenseCopy                                           *
 * Usage    : DenseCopy(A, B);                                    *
 *----------------------------------------------------------------*
 * DenseCopy copies the contents of the N by N matrix A into the  *
 * N by N matrix B.                                               *
 *                                                                *
 ******************************************************************/

void DenseCopy(DenseMat A, DenseMat B);


/******************************************************************
 *                                                                *
 * Function: DenseScale                                           *
 * Usage   : DenseScale(c, A);                                    *
 *----------------------------------------------------------------*
 * DenseScale scales the elements of the N by N matrix A by the   *
 * constant c and stores the result back in A.                    *
 *                                                                *
 ******************************************************************/

void DenseScale(realtype c, DenseMat A);


/******************************************************************
 *                                                                *
 * Function : DenseAddI                                           *
 * Usage    : DenseAddI(A);                                       *
 *----------------------------------------------------------------*
 * DenseAddI adds the identity matrix to A and stores the result  *
 * back in A.                                                     *
 *                                                                *
 ******************************************************************/

void DenseAddI(DenseMat A);


/******************************************************************
 *                                                                *
 * Function : DenseFreeMat                                        *
 * Usage    : DenseFreeMat(A);                                    *
 *----------------------------------------------------------------*
 * DenseFreeMat frees the memory allocated by DenseAllocMat for   *
 * the N by N matrix A.                                           *
 *                                                                *
 ******************************************************************/

void DenseFreeMat(DenseMat A);


/******************************************************************
 *                                                                *
 * Function : DenseFreePiv                                        *
 * Usage    : DenseFreePiv(p);                                    *
 *----------------------------------------------------------------*
 * DenseFreePiv frees the memory allocated by DenseAllocPiv for   *
 * the pivot information array p.                                 *
 *                                                                *
 ******************************************************************/

void DenseFreePiv(integertype *p);


/******************************************************************
 *                                                                *
 * Function : DensePrint                                          *
 * Usage    : DensePrint(A);                                      *
 *----------------------------------------------------------------*
 * This routine prints the N by N dense matrix A to standard      *
 * output as it would normally appear on paper. It is intended    *
 * as a debugging tool with small values of N. The elements are   *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/

void DensePrint(DenseMat A);
 

#endif
#ifdef __cplusplus
}
#endif
