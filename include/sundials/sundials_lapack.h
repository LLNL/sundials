/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:17 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for a generic package of dense matrix
 * operations for use with BLAS/LAPACK.  The routines listed in this 
 * file all use type LapackMat, defined below, for M-by-N matrices.
 *
 * The LapackAllocDenseMat function allocates a dense matrix for use
 * in the other LapackMat routines listed in this file. Matrix storage 
 * details are given in the documentation for the type LapackMat.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_LAPACK_H
#define _SUNDIALS_LAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_types.h>

  /*
   * =================================================================
   *              L A P A C K     C O N S T A N T S
   * =================================================================
   */

  /*
   *  LAPACK_DENSE: use dense matrices
   *  LAPACK_BAND:  use banded matrices
   */

#define LAPACK_DENSE 1
#define LAPACK_BAND  2

  /*
   * ==================================================================
   * Type definitions
   * ==================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Type : LapackMat
   * -----------------------------------------------------------------
   * The type LapackMat is defined to be a pointer to a structure
   * with various sizes, a data field, and an array of pointers to
   * the columns. The M and N fields indicates the number of rows and 
   * columns, respectively. The data field is a one dimensional array 
   * used for component storage. The cols field stores the pointers in 
   * data for the beginning of each column.
   * 
   * The LapackMat type can be used to store either a dense or a 
   * banded matrix.
   *
   * Dense Matrices
   * --------------
   *
   * The fields in LapackMat relevant for a dense matrix are:
   *    type  = LAPACK_DENSE
   *    M     - number of rows
   *    N     - number of columns
   *    ldim  - leading dimension (ldim >= M)
   *    data  - pointer to a contiguous block of realtype variables
   *    ldata - length of the data array =ldim*N
   *    cols  - array of pointers. cols[j] points to the first element 
   *            of the j-th column of the matrix in the array data.
   *
   * The elements of a dense matrix are stored columnwise (i.e columns 
   * are stored one on top of the other in memory). 
   * If A is of type LapackMat, then the (i,j)th element of A (with 
   * 0 <= i < M and 0 <= j < N) is given by (A->data)[j*n+i]. 
   *
   * The LAPACK_DENSE_COL and LAPACK_DENSE_ELEM macros below allow a 
   * user to access efficiently individual matrix elements without 
   * writing out explicit data structure references and without knowing 
   * too much about the underlying element storage. The only storage 
   * assumption needed is that elements are stored columnwise and that 
   * a pointer to the jth column of elements can be obtained via the 
   * LAPACK_DENSE_COL macro.
   *
   * Banded Matrices
   * ---------------
   *
   * The fields in LapackMat relevant for a banded matrix are:
   *    type  = LAPACK_BAND
   *    M     - number of rows
   *    N     - number of columns
   *    mu    - upper bandwidth, 0 <= mu <= min(M,N)
   *    ml    - lower bandwidth, 0 <= ml <= min(M,N)
   *    storage_mu - storage upper bandwidth, mu <= storage_mu <= N-1.
   *            The dgbtrf routine writes the LU factors into the storage 
   *            for A. The upper triangular factor U, however, may have 
   *            an upper bandwidth as big as MIN(N-1,mu+ml) because of 
   *            partial pivoting. The storage_mu field holds the upper 
   *            bandwidth allocated for A.
   *    ldim  - leading dimension (ldim >= storage_mu)
   *    data  - pointer to a contiguous block of realtype variables
   *    ldata - length of the data array =ldim*(storage_mu+ml+1)
   *    cols  - array of pointers. cols[j] points to the first element 
   *            of the j-th column of the matrix in the array data.
   *
   * The LAPACK_BAND_COL, LAPACK_BAND_COL_ELEM, and LAPACK_BAND_ELEM
   * macros below allow a user to access individual matrix
   * elements without writing out explicit data structure
   * references and without knowing too much about the underlying
   * element storage. The only storage assumption needed is that
   * elements are stored columnwise and that a pointer into the jth
   * column of elements can be obtained via the LAPACK_BAND_COL macro.
   * The LAPACK_BAND_COL_ELEM macro selects an element from a column
   * which has already been isolated via LAPACK_BAND_COL. The macro
   * LAPACK_BAND_COL_ELEM allows the user to avoid the translation 
   * from the matrix location (i,j) to the index in the array returned 
   * by LAPACK_BAND_COL at which the (i,j)th element is stored. 
   * -----------------------------------------------------------------
   */

  typedef struct _LapackMat {
    int type;
    int M;
    int N;
    int ldim;
    int mu;
    int ml;
    int storage_mu;
    realtype *data;
    int ldata;
    realtype **cols;
  } *LapackMat;

  /*
   * ==================================================================
   * Data accessor macros
   * ==================================================================
   */

  /*
   * -----------------------------------------------------------------
   * LAPACK_DENSE_COL and LAPACK_DENSE_ELEM
   * -----------------------------------------------------------------
   *
   * LAPACK_DENSE_COL(A,j) references the jth column of the M-by-N
   * dense LapackMat A, 0 <= j < N. The type of the expression 
   * LAPACK_DENSE_COL(A,j) is (realtype *). After the assignment in the
   * usage above, col_j may be treated as an array indexed from 0 to M-1. 
   * The (i,j)-th element of A is thus referenced by col_j[i].
   *
   * LAPACK_DENSE_ELEM(A,i,j) references the (i,j)th element of the 
   * dense M-by-N LapackMat A, 0 <= i < M ; 0 <= j < N.
   *
   * -----------------------------------------------------------------
   */

#define LAPACK_DENSE_COL(A,j) ((A->cols)[j])
#define LAPACK_DENSE_ELEM(A,i,j) ((A->cols)[j][i])

   /*
   * -----------------------------------------------------------------
   * LAPACK_BAND_COL, LAPACK_BAND_COL_ELEM, and LAPACK_BAND_ELEM
   * -----------------------------------------------------------------
   *  
   * LAPACK_BAND_COL(A,j) references the diagonal element of the jth
   * column of the N by N band matrix A, 0 <= j <= N-1. The type of
   * the expression LAPACK_BAND_COL(A,j) is realtype *. The pointer
   * returned by the call LAPACK_BAND_COL(A,j) can be treated as an
   * array which is indexed from -(A->mu) to (A->ml).
   * 
   * LAPACK_BAND_COL_ELEM references the (i,j)th entry of the band
   * matrix A when used in conjunction with LAPACK_BAND_COL. The 
   * index (i,j) should satisfy j-(A->mu) <= i <= j+(A->ml).
   *
   * LAPACK_BAND_ELEM(A,i,j) references the (i,j)th element of the
   * N by N band matrix A, where 0 <= i,j <= N-1. The location
   * (i,j) should further satisfy j-(A->mu) <= i <= j+(A->ml). 
   *
   * -----------------------------------------------------------------
   */
 
#define LAPACK_BAND_COL(A,j) (((A->cols)[j])+(A->storage_mu))
#define LAPACK_BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
#define LAPACK_BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->storage_mu)])

  /*
   * ==================================================================
   * Exported function prototypes
   * ==================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : LapackAllocDenseMat
   * -----------------------------------------------------------------
   * LapackAllocDenseMat allocates memory for an M-by-N dense matrix and
   * returns the storage allocated (type LapackMat). LapackAllocDenseMat
   * returns NULL if the request for matrix storage cannot be
   * satisfied. See the above documentation for the type LapackMat
   * for matrix storage details.
   * -----------------------------------------------------------------
   */

  LapackMat LapackAllocDenseMat(int M, int N);

  /*
   * -----------------------------------------------------------------
   * Function : LapackAllocBandMat
   * -----------------------------------------------------------------
   * LapackAllocBandMat allocates memory for an N by N band matrix 
   * with upper bandwidth mu, lower bandwidth ml, and storage upper
   * bandwidth storage_mu. Pass storage_mu as follows depending on 
   * whether A will be factored by BandGBTRF:
   *
   * (1) Pass storage_mu = mu if A will not be factored.
   *
   * (2) Pass storage_mu = MIN(N-1,mu+ml) if A will be factored.
   *
   * LapackAllocBandMat returns the storage allocated (type BandMat) or
   * NULL if the request for matrix storage cannot be satisfied.
   * See the documentation for the type BandMat for matrix storage
   * details.
   * -----------------------------------------------------------------
   */

  LapackMat LapackAllocBandMat(int N, int mu, int ml, int storage_mu);

  /*
   * -----------------------------------------------------------------
   * Function : LapackDenseMat
   * -----------------------------------------------------------------
   * LapackFreeMat frees the memory allocated by LapackAllocDenseMat
   * or LapackAllocBandMat.
   * -----------------------------------------------------------------
   */

  void LapackFreeMat(LapackMat A);

  /*
   * -----------------------------------------------------------------
   * Function : LapackAllocIntArray
   * -----------------------------------------------------------------
   * Usage : p = LapackAllocIntArray(N);
   *         if (p == NULL) ... memory request failed
   * -----------------------------------------------------------------
   * LapackAllocIntArray allocates memory an array of N integers and
   * returns the pointer to the memory it allocates. If the request for
   * memory storage cannot be satisfied, it returns NULL.
   * -----------------------------------------------------------------
   */

  int *LapackAllocIntArray(int N);

  /*
   * -----------------------------------------------------------------
   * Function : LapackAllocRealArray
   * -----------------------------------------------------------------
   * LapackAllocRealArray allocates memory an array of N realtype and
   * returns the pointer to the memory it allocates. If the request for
   * memory storage cannot be satisfied, it returns NULL.
   * -----------------------------------------------------------------
   */

  realtype *LapackAllocRealArray(int N);

  /*
   * -----------------------------------------------------------------
   * Function : LapackFreeArray
   * -----------------------------------------------------------------
   * LapackFreeArray frees memory allocated by LapackAllocRealArray or
   * LapackAllocIntArray.
   * -----------------------------------------------------------------
   */

  void LapackFreeArray(void *p);

  /*
   * -----------------------------------------------------------------
   * Function : LapackPrintMat
   * -----------------------------------------------------------------
   * This routine prints the M-by-N dense matrix A to standard output
   * as it would normally appear on paper. It is intended as a 
   * debugging tool with small values of M and N. The elements are
   * printed using the %g/%lg/%Lg option. A blank line is printed
   * before and after the matrix.
   * -----------------------------------------------------------------
   */

  void LapackPrintMat(LapackMat A);

  /*
   * ==================================================================
   * Blas and Lapack functions
   * ==================================================================
   */

#if defined(F77_FUNC)

#define dcopy_f77       F77_FUNC(dcopy, DCOPY)
#define dscal_f77       F77_FUNC(dscal, DSCAL)
#define dgemv_f77       F77_FUNC(dgemv, DGEMV)
#define dtrsv_f77       F77_FUNC(dtrsv, DTRSV)
#define dsyrk_f77       F77_FUNC(dsyrk, DSKYR)

#define dgbtrf_f77      F77_FUNC(dgbtrf, DGBTRF)
#define dgbtrs_f77      F77_FUNC(dgbtrs, DGBTRS)
#define dgetrf_f77      F77_FUNC(dgetrf, DGETRF)
#define dgetrs_f77      F77_FUNC(dgetrs, DGETRS)
#define dgeqp3_f77      F77_FUNC(dgeqp3, DGEQP3)
#define dgeqrf_f77      F77_FUNC(dgeqrf, DGEQRF)
#define dormqr_f77      F77_FUNC(dormqr, DORMQR)
#define dpotrf_f77      F77_FUNC(dpotrf, DPOTRF)
#define dpotrs_f77      F77_FUNC(dpotrs, DPOTRS)

#else

#define dcopy_f77       dcopy_
#define dscal_f77       dscal_
#define dgemv_f77       dgemv_
#define dtrsv_f77       dtrsv_
#define dsyrk_f77       dsyrk_

#define dgbtrf_f77      dgbtrf_
#define dgbtrs_f77      dgbtrs_
#define dgeqp3_f77      dgeqp3_
#define dgeqrf_f77      dgeqrf_
#define dgetrf_f77      dgetrf_
#define dgetrs_f77      dgetrs_
#define dormqr_f77      dormqr_
#define dpotrf_f77      dpotrf_
#define dpotrs_f77      dpotrs_

#endif

  /* Level-1 BLAS */
  
  extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);
  extern void dscal_f77(int *n, const double *alpha, double *x, const int *inc_x);

  /* Level-2 BLAS */

  extern void dgemv_f77(const char *trans, int *m, int *n, const double *alpha, const double *a, 
                        int *lda, const double *x, int *inc_x, const double *beta, double *y, int *inc_y, 
                        int len_trans);

  extern void dtrsv_f77(const char *uplo, const char *trans, const char *diag, const int *n, 
                        const double *a, const int *lda, double *x, const int *inc_x, 
                        int len_uplo, int len_trans, int len_diag);

  /* Level-3 BLAS */

  extern void dsyrk_f77(const char *uplo, const char *trans, const int *n, const int *k, 
                        const double *alpha, const double *a, const int *lda, const double *beta, 
                        const double *c, const int *ldc, int len_uplo, int len_trans);
  
  /* LAPACK */

  extern void dgbtrf_f77(const int *m, const int *n, const int *kl, const int *ku, 
                         double *ab, int *ldab, int *ipiv, int *info);

  extern void dgbtrs_f77(const char *trans, const int *n, const int *kl, const int *ku, const int *nrhs, 
                         double *ab, const int *ldab, int *ipiv, double *b, const int *ldb, 
                         int *info, int len_trans);


  extern void dgeqp3_f77(const int *m, const int *n, double *a, const int *lda, int *jpvt, double *tau, 
                         double *work, const int *lwork, int *info);

  extern void dgeqrf_f77(const int *m, const int *n, double *a, const int *lda, double *tau, double *work, 
                         const int *lwork, int *info);

  extern void dgetrf_f77(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);

  extern void dgetrs_f77(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, 
                         int *ipiv, double *b, const int *ldb, int *info, int len_trans);


  extern void dormqr_f77(const char *side, const char *trans, const int *m, const int *n, const int *k, 
                         double *a, const int *lda, double *tau, double *c, const int *ldc, 
                         double *work, const int *lwork, int *info, int len_side, int len_trans);

  extern void dpotrf_f77(const char *uplo, const int *n, double *a, int *lda, int *info, int len_uplo);

  extern void dpotrs_f77(const char *uplo, const int *n, const int *nrhs, double *a, const int *lda, 
                         double *b, const int *ldb, int * info, int len_uplo);


#ifdef __cplusplus
}
#endif

#endif
