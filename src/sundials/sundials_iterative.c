/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 *                Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the iterative.h header
 * file. It contains the implementation of functions that may be
 * useful for many different iterative solvers of A x = b.
 * -----------------------------------------------------------------
 */

#include <stdio.h>

#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>

#define FACTOR RCONST(1000.0)
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Function : ModifiedGS
 * -----------------------------------------------------------------
 * This implementation of ModifiedGS is a slight modification of a
 * previous modified Gram-Schmidt routine (called mgs) written by
 * Milo Dorr.
 * -----------------------------------------------------------------
 */
 
int ModifiedGS(N_Vector *v, realtype **h, int k, int p, 
               realtype *new_vk_norm)
{
  int  i, k_minus_1, i0;
  realtype new_norm_2, new_product, vk_norm, temp;
  
  vk_norm = SUNRsqrt(N_VDotProd(v[k],v[k]));
  k_minus_1 = k - 1;
  i0 = SUNMAX(k-p, 0);
  
  /* Perform modified Gram-Schmidt */
  
  for (i=i0; i < k; i++) {
    h[i][k_minus_1] = N_VDotProd(v[i], v[k]);
    N_VLinearSum(ONE, v[k], -h[i][k_minus_1], v[i], v[k]);
  }

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = SUNRsqrt(N_VDotProd(v[k], v[k]));

  /* If the norm of the new vector at v[k] is less than
     FACTOR (== 1000) times unit roundoff times the norm of the
     input vector v[k], then the vector will be reorthogonalized
     in order to ensure that nonorthogonality is not being masked
     by a very small vector length. */

  temp = FACTOR * vk_norm;
  if ((temp + (*new_vk_norm)) != temp) return(0);
  
  new_norm_2 = ZERO;

  for (i=i0; i < k; i++) {
    new_product = N_VDotProd(v[i], v[k]);
    temp = FACTOR * h[i][k_minus_1];
    if ((temp + new_product) == temp) continue;
    h[i][k_minus_1] += new_product;
    N_VLinearSum(ONE, v[k],-new_product, v[i], v[k]);
    new_norm_2 += SUNSQR(new_product);
  }

  if (new_norm_2 != ZERO) {
    new_product = SUNSQR(*new_vk_norm) - new_norm_2;
    *new_vk_norm = (new_product > ZERO) ? SUNRsqrt(new_product) : ZERO;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : ClassicalGS
 * -----------------------------------------------------------------
 * This implementation of ClassicalGS was contributed by Homer Walker
 * and Peter Brown.
 * -----------------------------------------------------------------
 */

int ClassicalGS(N_Vector *v, realtype **h, int k, int p, realtype *new_vk_norm,
                realtype *stemp, N_Vector *vtemp)
{
  int  i, i0, k_minus_1, retval;
  realtype vk_norm;

  k_minus_1 = k - 1;
  i0 = SUNMAX(k-p,0);

  /* Perform Classical Gram-Schmidt */

  retval = N_VDotProdMulti(k-i0+1, v[k], v+i0, stemp);
  if (retval != 0) return(-1);

  vk_norm = SUNRsqrt(stemp[k-i0]);
  for (i=k-i0-1; i >= 0; i--) {
    h[i][k_minus_1] = stemp[i];
    stemp[i+1] = -stemp[i];
    vtemp[i+1] = v[i];
  }
  stemp[0] = ONE;
  vtemp[0] = v[k];

  retval = N_VLinearCombination(k-i0+1, stemp, vtemp, v[k]);
  if (retval != 0) return(-1);

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = SUNRsqrt(N_VDotProd(v[k], v[k]));

  /* Reorthogonalize if necessary */

  if ((FACTOR * (*new_vk_norm)) < vk_norm) {

    retval = N_VDotProdMulti(k-i0, v[k], v+i0, stemp+1);
    if (retval != 0) return(-1);

    stemp[0] = ONE;
    vtemp[0] = v[k];
    for (i=i0; i < k; i++) {
      h[i][k_minus_1] += stemp[i-i0+1];
      stemp[i-i0+1] = -stemp[i-i0+1];
      vtemp[i-i0+1] = v[i-i0];
    }

    retval = N_VLinearCombination(k+1, stemp, vtemp, v[k]);
    if (retval != 0) return(-1);

    *new_vk_norm = SUNRsqrt(N_VDotProd(v[k],v[k]));
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : QRfact
 * -----------------------------------------------------------------
 * This implementation of QRfact is a slight modification of a
 * previous routine (called qrfact) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

int QRfact(int n, realtype **h, realtype *q, int job)
{
  realtype c, s, temp1, temp2, temp3;
  int i, j, k, q_ptr, n_minus_1, code=0;

  switch (job) {
  case 0:

    /* Compute a new factorization of H */

    code = 0;
    for (k=0; k < n; k++) {
      
      /* Multiply column k by the previous k-1 Givens rotations */

      for (j=0; j < k-1; j++) {
	i = 2*j;
	temp1 = h[j][k];
	temp2 = h[j+1][k];
	c = q[i];
	s = q[i+1];
	h[j][k] = c*temp1 - s*temp2;
	h[j+1][k] = s*temp1 + c*temp2;
      }
      
      /* Compute the Givens rotation components c and s */

      q_ptr = 2*k;
      temp1 = h[k][k];
      temp2 = h[k+1][k];
      if( temp2 == ZERO) {
	c = ONE;
	s = ZERO;
      } else if (SUNRabs(temp2) >= SUNRabs(temp1)) {
	temp3 = temp1/temp2;
	s = -ONE/SUNRsqrt(ONE+SUNSQR(temp3));
	c = -s*temp3;
      } else {
	temp3 = temp2/temp1;
	c = ONE/SUNRsqrt(ONE+SUNSQR(temp3));
	s = -c*temp3;
      }
      q[q_ptr] = c;
      q[q_ptr+1] = s;
      if( (h[k][k] = c*temp1 - s*temp2) == ZERO) code = k+1;
    }
    break;

  default:

    /* Update the factored H to which a new column has been added */

    n_minus_1 = n - 1;
    code = 0;
    
    /* Multiply the new column by the previous n-1 Givens rotations */

    for (k=0; k < n_minus_1; k++) {
      i = 2*k;
      temp1 = h[k][n_minus_1];
      temp2 = h[k+1][n_minus_1];
      c = q[i];
      s = q[i+1];
      h[k][n_minus_1] = c*temp1 - s*temp2;
      h[k+1][n_minus_1] = s*temp1 + c*temp2;
    }
    
    /* Compute new Givens rotation and multiply it times the last two
       entries in the new column of H.  Note that the second entry of 
       this product will be 0, so it is not necessary to compute it. */

    temp1 = h[n_minus_1][n_minus_1];
    temp2 = h[n][n_minus_1];
    if (temp2 == ZERO) {
      c = ONE;
      s = ZERO;
    } else if (SUNRabs(temp2) >= SUNRabs(temp1)) {
      temp3 = temp1/temp2;
      s = -ONE/SUNRsqrt(ONE+SUNSQR(temp3));
      c = -s*temp3;
    } else {
      temp3 = temp2/temp1;
      c = ONE/SUNRsqrt(ONE+SUNSQR(temp3));
      s = -c*temp3;
    }
    q_ptr = 2*n_minus_1;
    q[q_ptr] = c;
    q[q_ptr+1] = s;
    if ((h[n_minus_1][n_minus_1] = c*temp1 - s*temp2) == ZERO)
      code = n;
  }
  
  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : QRsol
 * -----------------------------------------------------------------
 * This implementation of QRsol is a slight modification of a
 * previous routine (called qrsol) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

int QRsol(int n, realtype **h, realtype *q, realtype *b)
{
  realtype c, s, temp1, temp2;
  int i, k, q_ptr, code=0;

  /* Compute Q*b */
  
  for (k=0; k < n; k++) {
    q_ptr = 2*k;
    c = q[q_ptr];
    s = q[q_ptr+1];
    temp1 = b[k];
    temp2 = b[k+1];
    b[k] = c*temp1 - s*temp2;
    b[k+1] = s*temp1 + c*temp2;
  }

  /* Solve  R*x = Q*b */

  for (k=n-1; k >= 0; k--) {
    if (h[k][k] == ZERO) {
      code = k + 1;
      break;
    }
    b[k] /= h[k][k];
    for (i=0; i < k; i++) b[i] -= b[k]*h[i][k];
  }
  
  return (code);
}


/*
 * -----------------------------------------------------------------
 * Function : QRAdd_MGS
 * -----------------------------------------------------------------
 * Implementation of QRAdd to be called in Anderson Acceleration 
 * -----------------------------------------------------------------
 */

int QRAdd_MGS(N_Vector *Q, realtype *R, N_Vector df,
              int m, int mMax, long int index, void *QRdata)
{
    sunindextype j;
    QRData qrdata = (QRData) QRdata;    

    N_VScale(ONE, df, qrdata->vtemp);
    for (j=0; j < (m-1); j++) {
      R[(m-1)*mMax+j] = N_VDotProd(Q[j], qrdata->vtemp);
      N_VLinearSum(ONE, qrdata->vtemp, -R[(m-1)*mMax + j], Q[j], qrdata->vtemp);
    }
    R[(m-1)*mMax + m-1] = SUNRsqrt(N_VDotProd(qrdata->vtemp, qrdata->vtemp));
    N_VScale((1/R[(m-1)*mMax +m-1]), qrdata->vtemp, Q[index]);

    /* Return success */
    return 0;
}

/*
 * -----------------------------------------------------------------
 * Function : QRAdd_ICWY
 * -----------------------------------------------------------------
 * Low synchronous implementation of QRAdd to be called in 
 * Anderson Acceleration. 
 * -----------------------------------------------------------------
 */

int QRAdd_ICWY(N_Vector *Q, realtype *R, N_Vector df, 
               int m, int mMax, long int index, void *QRdata)
{
    sunindextype j, k;
    int t_offset, retval;
    QRData qrdata = (QRData) QRdata;    

    /* T matrix is stored in row-wise order with each row having an additional
     *   value because it is lower triangular */
    N_VScale(ONE, df, qrdata->vtemp); /* stores d_fi in temp */

    t_offset = (m-2)*(m-1)/2;
    if (m-2 > 0) {
      /* T(k-1,1:k-1) = Q(:,1:k-1)^T * Q(:,k-1) */
      N_VDotProdMulti(m-1, Q[m-2], Q, qrdata->temp_array + t_offset);
    }
    /* T(k-1,k-1) = 1.0 */
    qrdata->temp_array[t_offset+(m-2)] = 1.0; /* I think adding iter-2 is correct to get to diagonal */

    /* R(1:k-1,k) = Q_k-1^T * df */
    N_VDotProdMulti(m-1, qrdata->vtemp, Q, R + (m-1)*mMax );

    /* Solve T * R(1:k-1,k) = R(1:k-1,k) */
    for (k = 0; k < (m-1); k++) {
      t_offset = k*(k+1)/2;
      R[(m-1)*mMax + k] = R[(m-1)*mMax + k] / qrdata->temp_array[t_offset + k];
      for (j = 0; j < (m-1); j++) {
        t_offset = j*(j+1)/2;
        R[(m-1)*mMax + j] -= R[(m-1)*mMax + k] * qrdata->temp_array[t_offset + k];
      }
    }

    /* Q(:,k) = df - Q_k-1 R(1:k-1,k) */
    N_VLinearCombination(m-1, R + (m-1)*mMax, Q, qrdata->vtemp2); 
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, Q[m-1]);
    
    /* R(k,k) = \| df \| */
    R[(m-1)*mMax + m-1] = SUNRsqrt(N_VDotProd(qrdata->vtemp, qrdata->vtemp));
    /* Q(:,k) = df / \| df \| */
    N_VScale((1/R[(m-1)*mMax + m-1]), qrdata->vtemp, Q[index]);
    
    /* Return success */
    return 0;
}

/*
 * -----------------------------------------------------------------
 * Function : QRAdd_CGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration. 
 * -----------------------------------------------------------------
 */

int QRAdd_CGS2(N_Vector *Q, realtype *R, N_Vector df, 
               int m, int mMax, long int index, void *QRdata)
{
    sunindextype j;
    QRData qrdata = (QRData) QRdata;    

    N_VScale(ONE, df, qrdata->vtemp); /* temp = df */ 
    
    /* s_k = Q_k-1^T df_aa -- update with sdata as a realtype* array */
    N_VDotProdMulti(m-1, qrdata->vtemp, Q, qrdata->temp_array); 

    /* y = df - Q_k-1 s_k */
    N_VLinearCombination(m-1, qrdata->temp_array, Q, qrdata->vtemp2); 
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp2);

    /* z_k = Q_k-1^T y */
    N_VDotProdMulti(m-1, qrdata->vtemp2, Q, R + (m-1)*mMax);

    /* df = y - Q_k-1 z_k  -- update using N_VLinearCombination */
    N_VLinearCombination(m-1, R + (m-1)*mMax, Q, qrdata->vtemp2); 
    N_VLinearSum(ONE, qrdata->vtemp2, -ONE, qrdata->vtemp, qrdata->vtemp);

    /* R(1:k-1,k) = s_k + z_k */
    for (j = 0; j < (m-1); j++) {
      R[(m-1)*mMax + j] = R[(m-1)*mMax + j] + qrdata->temp_array[j];
    }
    /* R(k,k) = \| df \| */
    R[(m-1)*mMax + m-1] = SUNRsqrt(N_VDotProd(qrdata->vtemp, qrdata->vtemp));
    /* Q(:,k) = df / R(k,k) */
    N_VScale((1/R[(m-1)*mMax + m-1]), qrdata->vtemp, Q[index]);
    /* -- end kernel -- */

    /* Return success */
    return 0;
}

/*
 * -----------------------------------------------------------------
 * Function : QRAdd_DCGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration. 
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int QRAdd_DCGS2(N_Vector *Q, realtype *R, N_Vector df, 
                                int m, int mMax, long int index, void *QRdata)
{
    sunindextype j;
    QRData qrdata = (QRData) QRdata;    

    N_VScale(ONE, df, qrdata->vtemp); /* temp = df */ 
    
    /* R(1:k-1,k) = Q_k-1^T df_aa */
    N_VDotProdMulti(m-2, qrdata->vtemp, Q, R + (m-1)*mMax );
    /* Delayed reorthogonalization */
    if ((m-1) > 1) {
        /* s = Q_k-2^T Q(:,k-1) */
        N_VDotProdMulti(m-1, Q[m-2], Q, qrdata->temp_array);

        /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
        N_VLinearCombination(m-2, qrdata->temp_array, Q, qrdata->vtemp2); 
        N_VLinearSum(ONE, Q[m-2], -ONE, qrdata->vtemp2, Q[m-2]);
        
        /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
        for (j = 0; j < (m-2); j++) {
          R[(m-2)*mMax + j] = R[(m-2)*mMax + j] + qrdata->temp_array[j];
        }
    }

    /* df = df - Q(:,k-1) s */
    N_VLinearCombination(m-1, R + (m-1)*mMax, Q, qrdata->vtemp2); 
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);

    /* R(k,k) = \| df \| */
    R[(m-1)*mMax + m-1] = SUNRsqrt(N_VDotProd(qrdata->vtemp, qrdata->vtemp));
    /* Q(:,k) = df / R(k,k) */
    N_VScale((1/R[(m-1)*mMax + m-1]), qrdata->vtemp, Q[index]);

    /* Return success */
    return 0;
}
