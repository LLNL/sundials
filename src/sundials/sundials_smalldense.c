/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-10-11 16:34:20 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic DENSE linear
 * solver package, intended for small dense matrices.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_smalldense.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Implementation */


realtype **denalloc(long int m, long int n)
{
  long int j;
  realtype **a;

  if ( (n <= 0) || (m <= 0) ) return(NULL);

  a = NULL;
  a = (realtype **) malloc(n * sizeof(realtype *));
  if (a == NULL) return(NULL);

  a[0] = NULL;
  a[0] = (realtype *) malloc(m * n * sizeof(realtype));
  if (a[0] == NULL) {
    free(a); a = NULL;
    return(NULL);
  }

  for (j=1; j < n; j++) a[j] = a[0] + j * m;

  return(a);
}

long int *denallocpiv(long int m)
{
  long int *piv;

  if (m <= 0) return(NULL);

  piv = NULL;
  piv = (long int *) malloc(m * sizeof(long int));

  return(piv);
}

long int denGETRF(realtype **a, long int m, long int n, long int *p)
{
  long int i, j, k, l;
  realtype *col_j, *col_k, *diag_k;
  realtype temp, mult, a_kj;
  booleantype swap;

  /* k = elimination step number */

  for (k=0; k < m-1; k++, p++) {

    col_k     = a[k];
    diag_k    = col_k + k;

    /* find l = pivot row number */

    l=k;
    for (i=k+1; i < m; i++)
      if (ABS(col_k[i]) > ABS(col_k[l])) l=i;
    *p = l;

    /* check for zero pivot element */

    if (col_k[l] == ZERO) return(k+1);
    
    /* swap a(l,k) and a(k,k) if necessary */
    
    if ( (swap = (l != k) )) {
      temp = col_k[l];
      col_k[l] = *diag_k;
      *diag_k = temp;
    }

    /* Scale the elements below the diagonal in         */
    /* column k by -1.0 / a(k,k). After the above swap, */
    /* a(k,k) holds the pivot element. This scaling     */
    /* stores the pivot row multipliers -a(i,k)/a(k,k)  */
    /* in a(i,k), i=k+1, ..., m-1.                      */

    mult = -ONE / (*diag_k);
    for(i=k+1; i < m; i++)
      col_k[i] *= mult;

    /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1 */
    /* row k is the pivot row after swapping with row l.      */
    /* The computation is done one column at a time,          */
    /* column j=k+1, ..., n-1.                                */

    for (j=k+1; j < n; j++) {

      col_j = a[j];
      a_kj = col_j[l];

      /* Swap the elements a(k,j) and a(k,l) if l!=k. */

      if (swap) {
	col_j[l] = col_j[k];
	col_j[k] = a_kj;
      }

      /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
      /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */

      if (a_kj != ZERO) {
	for (i=k+1; i < m; i++)
	  col_j[i] += a_kj * col_k[i];
      }
    }
  }

  /* set the last pivot row to be m-1 and check for a zero pivot */

  *p = m-1;
  if (a[m-1][m-1] == ZERO) return(m);

  /* return 0 to indicate success */

  return(0);
}

void denGETRS(realtype **a, long int m, long int *p, realtype *b)
{
  long int k, l, i;
  realtype mult, *col_k;

  /* Solve Ly = Pb, store solution y in b */

  for (k=0; k < m-1; k++) {
    l = p[k];
    mult = b[l];
    if (l != k) {
      b[l] = b[k];
      b[k] = mult;
    }
    col_k = a[k];
    for (i=k+1; i < m; i++)
      b[i] += mult*col_k[i];
  }
  
  /* Solve Ux = y, store solution x in b */
  
  for (k=m-1; k >= 0; k--) {
    col_k = a[k];
    b[k] /= col_k[k];
    mult = -b[k];
    for (i=0; i < k; i++)
      b[i] += mult*col_k[i];
  }
}

void denzero(realtype **a, long int m, long int n)
{
  long int i, j;
  realtype *col_j;

  for (j=0; j < n; j++) {
    col_j = a[j];
    for (i=0; i < m; i++)
      col_j[i] =  ZERO;
  }
}

void dencopy(realtype **a, realtype **b, long int m, long int n)
{
  long int i, j;
  realtype *a_col_j, *b_col_j;

  for (j=0; j < n; j++) {
    a_col_j = a[j];
    b_col_j = b[j];
    for (i=0; i < m; i++)
      b_col_j[i] = a_col_j[i];
  }

}

void denscale(realtype c, realtype **a, long int m, long int n)
{
  long int i, j;
  realtype *col_j;

  for (j=0; j < n; j++) {
    col_j = a[j];
    for (i=0; i < m; i++)
      col_j[i] *= c;
  }
}

void denaddI(realtype **a, long int m)
{
  long int i;
  
  for (i=0; i < m; i++) a[i][i] += ONE;
}

void denfreepiv(long int *p)
{
  free(p); p = NULL;
}

void denfree(realtype **a)
{
  free(a[0]); a[0] = NULL;
  free(a); a = NULL;
}

void denprint(realtype **a, long int m, long int n)
{
  long int i, j;

  printf("\n");
  for (i=0; i < m; i++) {
    for (j=0; j < n; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%10Lg  ", a[j][i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%10lg  ", a[j][i]);
#else
      printf("%10g  ", a[j][i]);
#endif
    }
    printf("\n");
  }
  printf("\n");
}
