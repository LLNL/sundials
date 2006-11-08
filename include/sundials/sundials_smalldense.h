/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-11-08 01:01:17 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for a generic DENSE linear solver
 * package, intended for small m-by-n dense matrices, with m<=n.
 * These routines use the type realtype** for dense matrix arguments.
 * -----------------------------------------------------------------
 */

#ifndef _SMALLDENSE_H
#define _SMALLDENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_types.h>

  /*
   * -----------------------------------------------------------------
   * Functions: denalloc and denfree
   * -----------------------------------------------------------------
   * Usage : realtype **a;
   *         a = denalloc(m, n);
   *         if (a == NULL) ... memory request failed
   *
   *         denfree(a);
   * -----------------------------------------------------------------
   * denalloc(m, n) allocates storage for an m-by-n dense matrix.
   * It returns a pointer to the newly allocated storage if successful. 
   * If the memory request cannot be satisfied, then denalloc returns 
   * NULL. The underlying type of the dense matrix returned is a double
   * pointer to realtype. If we allocate a dense matrix a by 
   * a = denalloc(m, n), then a[j][i] references the (i,j)-th element 
   * of the matrix a, 0 <= i < m, 0 <= j < n,  and a[j] is a pointer 
   * to the first element in the jth column of a. The location a[0] 
   * contains a pointer to m*n contiguous locations which contain the 
   * elements of a.
   * -----------------------------------------------------------------
   * denfree(a) frees the dense matrix a allocated by denalloc.
   * -----------------------------------------------------------------
   */

  realtype **denalloc(long int m, long int n);
  void denfree(realtype **a);

  /*
   * -----------------------------------------------------------------
   * Functions: denallocpiv and denfreepiv
   * -----------------------------------------------------------------
   * Usage : long int *pivot;
   *         pivot = denallocpiv(n);
   *         if (pivot == NULL) ... memory request failed
   *
   *         denfreepiv(p);
   * -----------------------------------------------------------------
   * denallocpiv(n) allocates an array of n long int. It returns
   * a pointer to the first element in the array if successful.
   * It returns NULL if the memory request could not be satisfied.
   * -----------------------------------------------------------------
   * denfreepiv(p) frees the pivot array p allocated by denallocpiv.
   * -----------------------------------------------------------------
   */

  long int *denallocpiv(long int n);
  void denfreepiv(long int *p);

  /*
   * -----------------------------------------------------------------
   * Functions: denseallocbeta and densefreebeta
   * -----------------------------------------------------------------
   * denseallocBeta allocates memory for the beta array to be
   * filled in by the denseGEQRF routine during the factorization
   * of an m-by-n dense matrix. The underlying type for beta
   * information is an array of n realtype and this routine returns
   * the pointer to the memory it allocates. If the request for
   * beta storage cannot be satisfied, denseallocBeta returns NULL.
   * -----------------------------------------------------------------
   * densefreeBeta frees the memory allocated by denseallocBeta for
   * the array beta.
   * -----------------------------------------------------------------
   */

  realtype *denallocbeta(long int m);
  void denfreebeta(realtype *beta);

  /*
   * -----------------------------------------------------------------
   * Functions: denGETRF and denGETRS
   * -----------------------------------------------------------------
   * Usage : long int ier;
   *         ier = denGETRF(a,m,n,p);
   *         if (ier > 0) ... zero element encountered during
   *                          the factorization
   *
   *         realtype *b;
   *         ier = denGETRF(a,n,n,p);
   *         if (ier == 0) denGETRS(a,n,p,b);
   * -----------------------------------------------------------------
   * denGETRF(a,m,n,p) factors the m-by-n dense matrix a, m>=n.
   * It overwrites the elements of a with its LU factors and keeps 
   * track of the pivot rows chosen in the pivot array p.
   *
   * A successful LU factorization leaves the matrix a and the
   * pivot array p with the following information:
   *
   * (1) p[k] contains the row number of the pivot element chosen
   *     at the beginning of elimination step k, k=0, 1, ..., n-1.
   *
   * (2) If the unique LU factorization of a is given by Pa = LU,
   *     where P is a permutation matrix, L is a lower trapezoidal
   *     matrix with all 1.0 on the diagonal, and U is an upper
   *     triangular matrix, then the upper triangular part of a
   *     (including its diagonal) contains U and the strictly lower
   *     trapezoidal part of a contains the multipliers, I-L.
   *
   * Note that for square matrices (m=n), L is unit lower triangular.
   *
   * denGETRF returns 0 if successful. Otherwise it encountered a
   * zero diagonal element during the factorization. In this case
   * it returns the column index (numbered from one) at which it
   * encountered the zero.
   * -----------------------------------------------------------------
   * denGETRS(a,n,p,b) solves the n by n linear system a*x = b.
   * It assumes that a has been LU factored and the pivot array p has
   * been set by a successful call to denGETRF(a,n,n,p). 
   * denGETRS does not check whether a is square!
   * The solution x is written into the b array.
   * -----------------------------------------------------------------
   */

  long int denGETRF(realtype **a, long int m, long int n, long int *p);
  void denGETRS(realtype **a, long int n, long int *p, realtype *b);

  /*
   * -----------------------------------------------------------------
   * Functions : densePOTRF and densePOTRS
   * -----------------------------------------------------------------
   * densePOTRF computes the Cholesky factorization of a real symmetric
   * positive definite matrix a.
   * -----------------------------------------------------------------
   * densePOTRS solves a system of linear equations a*x = b with a 
   * symmetric positive definite matrix a using the Cholesky factorization
   * a = L*L**T computed by densePOTRF.
   * -----------------------------------------------------------------
   */

  long int denPOTRF(realtype **a, long int m);
  void denPOTRS(realtype **a, long int m, realtype *b);

  /*
   * -----------------------------------------------------------------
   * Functions : denseGEQRF and denseORMQR
   * -----------------------------------------------------------------
   * denseGEQRF computes a QR factorization of a real m-by-n matrix a:
   * a = Q * R (with m>=n)
   *
   * denseGEQRF reuires a temporary work vector wrk of length m.
   * -----------------------------------------------------------------
   * denseORMQR computes the product w = Q * v where Q is a real 
   * orthogonal matrix defined as the product of k elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   * as returned by denseGEQRF. Q is an m-by-n matrix, v is a vector
   * of length n and w is a vector of length m (with m>=n).
   *
   * denseORMQR reuires a temporary work vector wrk of length m.
   * -----------------------------------------------------------------
   */

  int denGEQRF(realtype **a, long int m, long int n, realtype *beta, realtype *v);
  int denORMQR(realtype **a, long int m, long int n, realtype *beta,
               realtype *v, realtype *w, realtype *wrk);

  /*
   * -----------------------------------------------------------------
   * Function : denzero
   * -----------------------------------------------------------------
   * denzero(a,m,n) sets all the elements of the m-by-n dense matrix
   * a to be 0.0.
   * -----------------------------------------------------------------
   */

  void denzero(realtype **a, long int m, long int n);

  /*
   * -----------------------------------------------------------------
   * Function : dencopy
   * -----------------------------------------------------------------
   * dencopy(a,b,m,n) copies the m-by-n dense matrix a into the
   * m-by-n dense matrix b.
   * -----------------------------------------------------------------
   */

  void dencopy(realtype **a, realtype **b, long int m, long int n);

  /*
   * -----------------------------------------------------------------
   * Function : denscale
   * -----------------------------------------------------------------
   * denscale(c,a,m,n) scales every element in the m-by-n dense
   * matrix a by c.
   * -----------------------------------------------------------------
   */

  void denscale(realtype c, realtype **a, long int m, long int n);

  /*
   * -----------------------------------------------------------------
   * Function : denaddI
   * -----------------------------------------------------------------
   * denaddI(a,n) increments the diagonal elements of the dense 
   * m-by-n matrix a by 1.0. (a_ii <= a_ii + 1, i=1,2,..n-1.)
   * denaddI is typically used with square matrices.
   * denaddI does NOT check for m >= n! Therefore, a segmentation 
   * fault will occur if m<n!
   * -----------------------------------------------------------------
   */

  void denaddI(realtype **a, long int n);

  /*
   * -----------------------------------------------------------------
   * Function : denprint
   * -----------------------------------------------------------------
   * denprint(a,m,n) prints the m-by-n dense matrix a to standard
   * output as it would normally appear on paper. It is intended as
   * a debugging tool with small values of m and n. The elements are
   * printed using the %g/%lg/%Lg option. A blank line is printed
   * before and after the matrix.
   * -----------------------------------------------------------------
   */

  void denprint(realtype **a, long int m, long int n);

#ifdef __cplusplus
}
#endif

#endif
