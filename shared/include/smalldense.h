/*******************************************************************
 *                                                                 *
 * File          : smalldense.h                                    *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for a generic DENSE linear solver       *
 * package, intended for small dense matrices.  These routines     *
 * use the type realtype** for dense matrix arguments.             *
 *                                                                 *
 * These routines begin with "den" (except for the factor and      *
 * solve routines which are called gefa and gesl, respectively).   *
 * The underlying matrix storage is described in the               *
 * documentation for denalloc.                                     *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _smalldense_h
#define _smalldense_h

#include "sundialstypes.h"

 
/******************************************************************
 *                                                                *
 * Function : denalloc                                            *
 * Usage    : realtype **a;                                       *
 *            a = denalloc(n);                                    *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * denalloc(n) allocates storage for an n by n dense matrix. It   *
 * returns a pointer to the newly allocated storage if            *
 * successful. If the memory request cannot be satisfied, then    *
 * denalloc returns NULL. The underlying type of the dense matrix *
 * returned is realtype **. If we allocate a dense matrix         *
 * realtype **a by a = denalloc(n), then a[j][i] references the   *
 * (i,j)th element of the matrix a, 0 <= i,j <= n-1, and a[j] is  *
 * a pointer to the first element in the jth column of a.         *
 * The location a[0] contains a pointer to n^2 contiguous         *
 * locations which contain the elements of a.                     *
 *                                                                *
 ******************************************************************/

realtype **denalloc(integertype n);


/******************************************************************
 *                                                                *
 * Function : denallocpiv                                         *
 * Usage    : integertype *pivot;                                 *
 *            pivot = denallocpiv(n);                             *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * denallocpiv(n) allocates an array of n integertype. It returns *
 * a pointer to the first element in the array if successful.     *
 * It returns NULL if the memory request could not be satisfied.  *
 *                                                                *
 ******************************************************************/

integertype *denallocpiv(integertype n);


/******************************************************************
 *                                                                *
 * Function : gefa                                                *
 * Usage    : integertype ier;                                    *
 *            ier = gefa(a,n,p);                                  *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * gefa(a,n,p) factors the n by n dense matrix a. It overwrites   *
 * the elements of a with its LU factors and keeps track of the   *
 * pivot rows chosen in the pivot array p.                        *
 *                                                                *
 * A successful LU factorization leaves the matrix a and the      *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., n-1.  *
 *                                                                *
 * (2) If the unique LU factorization of a is given by Pa = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of a     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of a contains the multipliers, I-L.        *
 *                                                                * 
 * gefa returns 0 if successful. Otherwise it encountered a zero  *
 * diagonal element during the factorization. In this case it     *
 * returns the column index (numbered from one) at which it       *
 * encountered the zero.                                          *
 *                                                                *
 ******************************************************************/

integertype gefa(realtype **a, integertype n, integertype *p);


/******************************************************************
 *                                                                *
 * Function : gesl                                                *
 * Usage    : realtype *b;                                        *
 *            ier = gefa(a,n,p);                                  *
 *            if (ier == 0) gesl(a,n,p,b);                        *
 *----------------------------------------------------------------*
 * gesl(a,n,p,b) solves the n by n linear system ax = b. It       *
 * assumes that a has been LU factored and the pivot array p has  *
 * been set by a successful call to gefa(a,n,p). The solution x   *
 * is written into the b array.                                   *
 *                                                                *
 ******************************************************************/

void gesl(realtype **a, integertype n, integertype *p, realtype *b);


/******************************************************************
 *                                                                *
 * Function : denzero                                             *
 * Usage    : denzero(a,n);                                       *
 *----------------------------------------------------------------*
 * denzero(a,n) sets all the elements of the n by n dense matrix  *
 * a to be 0.0.                                                   *
 *                                                                *
 ******************************************************************/

void denzero(realtype **a, integertype n);


/******************************************************************
 *                                                                *
 * Function : dencopy                                             *
 * Usage    : dencopy(a,b,n);                                     *
 *----------------------------------------------------------------*
 * dencopy(a,b,n) copies the n by n dense matrix a into the       *
 * n by n dense matrix b.                                         *
 *                                                                *
 ******************************************************************/

void dencopy(realtype **a, realtype **b, integertype n);


/******************************************************************
 *                                                                *
 * Function : denscale                                            *
 * Usage    : denscale(c,a,n);                                    *
 *----------------------------------------------------------------*
 * denscale(c,a,n) scales every element in the n by n dense       *
 * matrix a by c.                                                 *
 *                                                                *
 ******************************************************************/

void denscale(realtype c, realtype **a, integertype n);


/******************************************************************
 *                                                                *
 * Function : denaddI                                             *
 * Usage    : denaddI(a,n);                                       *
 *----------------------------------------------------------------*
 * denaddI(a,n) increments the n by n dense matrix a by the       *
 * identity matrix.                                               *
 *                                                                *
 ******************************************************************/

void denaddI(realtype **a, integertype n);


/******************************************************************
 *                                                                *
 * Function : denfreepiv                                          *
 * Usage    : denfreepiv(p);                                      *
 *----------------------------------------------------------------*
 * denfreepiv(p) frees the pivot array p allocated by             *
 * denallocpiv.                                                   *
 *                                                                *
 ******************************************************************/

void denfreepiv(integertype *p);


/******************************************************************
 *                                                                *
 * Function : denfree                                             *
 * Usage    : denfree(a);                                         *
 *----------------------------------------------------------------*
 * denfree(a) frees the dense matrix a allocated by denalloc.     *
 *                                                                *
 ******************************************************************/

void denfree(realtype **a);


/******************************************************************
 *                                                                *
 * Function : denprint                                            *
 * Usage    : denprint(a,n);                                      *
 *----------------------------------------------------------------*
 * denprint(a,n) prints the n by n dense matrix a to standard     *
 * output as it would normally appear on paper. It is intended as *
 * a debugging tool with small values of n. The elements are      *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/

void denprint(realtype **a, integertype n);
 

#endif

#ifdef __cplusplus
}
#endif
