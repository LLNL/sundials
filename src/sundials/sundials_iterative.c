/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 *                Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
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
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <sundials/sundials.h>

#include "sundials_iterative_impl.h"

#define FACTOR RCONST(1000.0)
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Function : SUNModifiedGS
 * -----------------------------------------------------------------
 * This implementation of SUNModifiedGS is a slight modification of
 * a previous modified Gram-Schmidt routine (called mgs) written by
 * Milo Dorr.
 * -----------------------------------------------------------------
 */

SUNErrCode ModifiedGS(N_Vector* v, realtype** h, int k, int p,
                      realtype* new_vk_norm)
{
  return (SUNModifiedGS(v, h, k, p, new_vk_norm));
}

SUNErrCode SUNModifiedGS(N_Vector* v, realtype** h, int k, int p,
                         realtype* new_vk_norm)
{
  int i, k_minus_1, i0;
  realtype new_norm_2, new_product, vk_norm, temp;
  SUNDeclareContext(v[0]->sunctx);

  vk_norm = SUNCheckCallLastErr((N_VDotProd(v[k], v[k])), sunctx);
  vk_norm = SUNRsqrt(vk_norm);
  k_minus_1 = k - 1;
  i0        = SUNMAX(k - p, 0);

  /* Perform modified Gram-Schmidt */

  for (i = i0; i < k; i++)
  {
    h[i][k_minus_1] = SUNCheckCallLastErr(N_VDotProd(v[i], v[k]), sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, v[k], -h[i][k_minus_1], v[i], v[k]), sunctx);
  }

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = SUNCheckCallLastErr(N_VDotProd(v[k], v[k]), sunctx);
  *new_vk_norm = SUNRsqrt(*new_vk_norm);

  /* If the norm of the new vector at v[k] is less than
     FACTOR (== 1000) times unit roundoff times the norm of the
     input vector v[k], then the vector will be reorthogonalized
     in order to ensure that nonorthogonality is not being masked
     by a very small vector length. */

  temp = FACTOR * vk_norm;
  if ((temp + (*new_vk_norm)) != temp) return (0);

  new_norm_2 = ZERO;

  for (i = i0; i < k; i++)
  {
    new_product = SUNCheckCallLastErr(N_VDotProd(v[i], v[k]), sunctx);
    temp = FACTOR * h[i][k_minus_1];
    if ((temp + new_product) == temp) continue;
    h[i][k_minus_1] += new_product;
    SUNCheckCallLastErr(N_VLinearSum(ONE, v[k], -new_product, v[i], v[k]), sunctx);
    new_norm_2 += SUNSQR(new_product);
  }

  if (new_norm_2 != ZERO)
  {
    new_product  = SUNSQR(*new_vk_norm) - new_norm_2;
    *new_vk_norm = (new_product > ZERO) ? SUNRsqrt(new_product) : ZERO;
  }

  return (0);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNClassicalGS
 * -----------------------------------------------------------------
 * This implementation of SUNClassicalGS was contributed by Homer
 * Walker and Peter Brown.
 * -----------------------------------------------------------------
 */

SUNErrCode ClassicalGS(N_Vector* v, realtype** h, int k, int p,
                       realtype* new_vk_norm, realtype* stemp, N_Vector* vtemp)
{
  return (SUNClassicalGS(v, h, k, p, new_vk_norm, stemp, vtemp));
}

SUNErrCode SUNClassicalGS(N_Vector* v, realtype** h, int k, int p,
                          realtype* new_vk_norm, realtype* stemp, N_Vector* vtemp)
{
  int i, i0, k_minus_1;
  realtype vk_norm;
  SUNDeclareContext(v[0]->sunctx);

  k_minus_1 = k - 1;
  i0        = SUNMAX(k - p, 0);

  /* Perform Classical Gram-Schmidt */

  SUNCheckCall(N_VDotProdMulti(k - i0 + 1, v[k], v + i0, stemp), sunctx);

  vk_norm = SUNRsqrt(stemp[k - i0]);
  for (i = k - i0 - 1; i >= 0; i--)
  {
    h[i][k_minus_1] = stemp[i];
    stemp[i + 1]    = -stemp[i];
    vtemp[i + 1]    = v[i];
  }
  stemp[0] = ONE;
  vtemp[0] = v[k];

  SUNCheckCall(N_VLinearCombination(k - i0 + 1, stemp, vtemp, v[k]),
                     sunctx);

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = SUNCheckCallLastErr(SUNRsqrt(N_VDotProd(v[k], v[k])), sunctx);

  /* Reorthogonalize if necessary */

  if ((FACTOR * (*new_vk_norm)) < vk_norm)
  {
    SUNCheckCall(N_VDotProdMulti(k - i0, v[k], v + i0, stemp + 1), sunctx);

    stemp[0] = ONE;
    vtemp[0] = v[k];
    for (i = i0; i < k; i++)
    {
      h[i][k_minus_1] += stemp[i - i0 + 1];
      stemp[i - i0 + 1] = -stemp[i - i0 + 1];
      vtemp[i - i0 + 1] = v[i - i0];
    }

    SUNCheckCall(N_VLinearCombination(k + 1, stemp, vtemp, v[k]), sunctx);

    *new_vk_norm = SUNCheckCallLastErr(SUNRsqrt(N_VDotProd(v[k], v[k])), sunctx);
  }

  return (0);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRfact
 * -----------------------------------------------------------------
 * This implementation of SUNQRfact is a slight modification of a
 * previous routine (called qrfact) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

SUNLsStatus QRfact(int n, realtype** h, realtype* q, int job)
{
  return (SUNQRfact(n, h, q, job));
}

SUNLsStatus SUNQRfact(int n, realtype** h, realtype* q, int job)
{
  realtype c, s, temp1, temp2, temp3;
  int i, j, k, q_ptr, n_minus_1, code = 0;

  switch (job)
  {
  case 0:

    /* Compute a new factorization of H */

    code = 0;
    for (k = 0; k < n; k++)
    {
      /* Multiply column k by the previous k-1 Givens rotations */

      for (j = 0; j < k - 1; j++)
      {
        i           = 2 * j;
        temp1       = h[j][k];
        temp2       = h[j + 1][k];
        c           = q[i];
        s           = q[i + 1];
        h[j][k]     = c * temp1 - s * temp2;
        h[j + 1][k] = s * temp1 + c * temp2;
      }

      /* Compute the Givens rotation components c and s */

      q_ptr = 2 * k;
      temp1 = h[k][k];
      temp2 = h[k + 1][k];
      if (temp2 == ZERO)
      {
        c = ONE;
        s = ZERO;
      }
      else if (SUNRabs(temp2) >= SUNRabs(temp1))
      {
        temp3 = temp1 / temp2;
        s     = -ONE / SUNRsqrt(ONE + SUNSQR(temp3));
        c     = -s * temp3;
      }
      else
      {
        temp3 = temp2 / temp1;
        c     = ONE / SUNRsqrt(ONE + SUNSQR(temp3));
        s     = -c * temp3;
      }
      q[q_ptr]     = c;
      q[q_ptr + 1] = s;
      if ((h[k][k] = c * temp1 - s * temp2) == ZERO) code = k + 1;
    }
    break;

  default:

    /* Update the factored H to which a new column has been added */

    n_minus_1 = n - 1;
    code      = 0;

    /* Multiply the new column by the previous n-1 Givens rotations */

    for (k = 0; k < n_minus_1; k++)
    {
      i                   = 2 * k;
      temp1               = h[k][n_minus_1];
      temp2               = h[k + 1][n_minus_1];
      c                   = q[i];
      s                   = q[i + 1];
      h[k][n_minus_1]     = c * temp1 - s * temp2;
      h[k + 1][n_minus_1] = s * temp1 + c * temp2;
    }

    /* Compute new Givens rotation and multiply it times the last two
       entries in the new column of H.  Note that the second entry of
       this product will be 0, so it is not necessary to compute it. */

    temp1 = h[n_minus_1][n_minus_1];
    temp2 = h[n][n_minus_1];
    if (temp2 == ZERO)
    {
      c = ONE;
      s = ZERO;
    }
    else if (SUNRabs(temp2) >= SUNRabs(temp1))
    {
      temp3 = temp1 / temp2;
      s     = -ONE / SUNRsqrt(ONE + SUNSQR(temp3));
      c     = -s * temp3;
    }
    else
    {
      temp3 = temp2 / temp1;
      c     = ONE / SUNRsqrt(ONE + SUNSQR(temp3));
      s     = -c * temp3;
    }
    q_ptr        = 2 * n_minus_1;
    q[q_ptr]     = c;
    q[q_ptr + 1] = s;
    if ((h[n_minus_1][n_minus_1] = c * temp1 - s * temp2) == ZERO) code = n;
  }

  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRsol
 * -----------------------------------------------------------------
 * This implementation of SUNQRsol is a slight modification of a
 * previous routine (called qrsol) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

SUNLsStatus QRsol(int n, realtype** h, realtype* q, realtype* b)
{
  return (SUNQRsol(n, h, q, b));
}

SUNLsStatus SUNQRsol(int n, realtype** h, realtype* q, realtype* b)
{
  realtype c, s, temp1, temp2;
  int i, k, q_ptr, code = 0;

  /* Compute Q*b */

  for (k = 0; k < n; k++)
  {
    q_ptr    = 2 * k;
    c        = q[q_ptr];
    s        = q[q_ptr + 1];
    temp1    = b[k];
    temp2    = b[k + 1];
    b[k]     = c * temp1 - s * temp2;
    b[k + 1] = s * temp1 + c * temp2;
  }

  /* Solve  R*x = Q*b */

  for (k = n - 1; k >= 0; k--)
  {
    if (h[k][k] == ZERO)
    {
      code = k + 1;
      break;
    }
    b[k] /= h[k][k];
    for (i = 0; i < k; i++) b[i] -= b[k] * h[i][k];
  }

  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_MGS
 * -----------------------------------------------------------------
 * Implementation of QRAdd to be called in Anderson Acceleration
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_MGS(N_Vector* Q, realtype* R, N_Vector df, int m, int mMax,
                        void* QRdata)
{
  sunindextype j;
  SUNQRData qrdata  = (SUNQRData)QRdata;
  SUNDeclareContext(Q[0]->sunctx);

  SUNCheckCallLastErr(N_VScale(ONE, df, qrdata->vtemp), sunctx);
  for (j = 0; j < m; j++)
  {
    R[m * mMax + j] = SUNCheckCallLastErr(N_VDotProd(Q[j], qrdata->vtemp), sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp, -R[m * mMax + j], Q[j], qrdata->vtemp), sunctx);
  }
  R[m * mMax + m] = SUNCheckCallLastErr(N_VDotProd(qrdata->vtemp, qrdata->vtemp), sunctx);
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  SUNCheckCallLastErr(N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]), sunctx);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_ICWY
 * -----------------------------------------------------------------
 * Low synchronous implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_ICWY(N_Vector* Q, realtype* R, N_Vector df, int m, int mMax,
                         void* QRdata)
{
  sunindextype j, k;
  SUNQRData qrdata  = (SUNQRData)QRdata;
  SUNDeclareContext(Q[0]->sunctx);

  SUNCheckCallLastErr(N_VScale(ONE, df, qrdata->vtemp), sunctx);  /* stores d_fi in temp */

  if (m > 0)
  {
    /* T(1:k-1,k-1)^T = Q(:,1:k-1)^T * Q(:,k-1) */
    SUNCheckCall(N_VDotProdMulti(m, Q[m - 1], Q,
                                       qrdata->temp_array + (m - 1) * mMax),
                       sunctx);

    /* T(k-1,k-1) = 1.0 */
    qrdata->temp_array[(m - 1) * mMax + (m - 1)] = ONE;

    /* R(1:k-1,k) = Q_k-1^T * df */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax),
                       sunctx);

    /* Solve T^T * R(1:k-1,k) = R(1:k-1,k) */
    for (k = 0; k < m; k++)
    {
      /* Skip setting the diagonal element because it doesn't change */
      for (j = k + 1; j < m; j++)
      {
        R[m * mMax + j] -= R[m * mMax + k] * qrdata->temp_array[j * mMax + k];
      }
    }
    /* end */

    /* Q(:,k-1) = df - Q_k-1 R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2),
                       sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp), sunctx);
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = SUNCheckCallLastErr(N_VDotProd(qrdata->vtemp, qrdata->vtemp), sunctx);
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / \| df \| */
  SUNCheckCallLastErr(N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]), sunctx);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_ICWY_SB
 * -----------------------------------------------------------------
 * Low synchronous implementation of QRAdd to be called in
 * Anderson Acceleration which utilizes a single buffer reduction.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_ICWY_SB(N_Vector* Q, realtype* R, N_Vector df, int m,
                            int mMax, void* QRdata)
{
  sunindextype j, k;
  SUNQRData qrdata  = (SUNQRData)QRdata;
  SUNDeclareContext(Q[0]->sunctx);

  SUNCheckCallLastErr(N_VScale(ONE, df, qrdata->vtemp), sunctx);  /* stores d_fi in temp */

  if (m > 0)
  {
    /* T(1:k-1,k-1)^T = Q(:,1:k-1)^T * Q(:,k-1) */
    SUNCheckCall(N_VDotProdMultiLocal(m, Q[m - 1], Q,
                                            qrdata->temp_array + (m - 1) * mMax),
                       sunctx);

    /* R(1:k-1,k) = Q_k-1^T * df */
    /* Put R values at end of temp_array */
    SUNCheckCall(N_VDotProdMultiLocal(m, qrdata->vtemp, Q,
                                            qrdata->temp_array +
                                              (m - 1) * mMax + m),
                       sunctx);

    SUNCheckCall(N_VDotProdMultiAllReduce(m + m, qrdata->vtemp,
                                                qrdata->temp_array +
                                                  (m - 1) * mMax),
                       sunctx);

    /* Move the last values from temp array into R */
    for (k = 0; k < m; k++)
    {
      R[m * mMax + k] = qrdata->temp_array[(m - 1) * mMax + m + k];
    }

    /* T(k-1,k-1) = 1.0 */
    qrdata->temp_array[(m - 1) * mMax + (m - 1)] = ONE;

    /* Solve T^T * R(1:k-1,k) = R(1:k-1,k) */
    for (k = 0; k < m; k++)
    {
      /* Skip setting the diagonal element because it doesn't change */
      for (j = k + 1; j < m; j++)
      {
        R[m * mMax + j] -= R[m * mMax + k] * qrdata->temp_array[j * mMax + k];
      }
    }
    /* end */

    /* Q(:,k-1) = df - Q_k-1 R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2),
                       sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp), sunctx);
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = SUNCheckCallLastErr(N_VDotProd(qrdata->vtemp, qrdata->vtemp), sunctx);
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / \| df \| */
  SUNCheckCallLastErr(N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]), sunctx);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_CGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_CGS2(N_Vector* Q, realtype* R, N_Vector df, int m, int mMax,
                         void* QRdata)
{
  sunindextype j;
  SUNQRData qrdata  = (SUNQRData)QRdata;
  SUNDeclareContext(Q[0]->sunctx);

  SUNCheckCallLastErr(N_VScale(ONE, df, qrdata->vtemp), sunctx);  /* temp = df */

  if (m > 0)
  {
    /* s_k = Q_k-1^T df_aa -- update with sdata as a realtype* array */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax),
                       sunctx);

    /* y = df - Q_k-1 s_k */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2),
                       sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp2), sunctx);

    /* z_k = Q_k-1^T y */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp2, Q, qrdata->temp_array),
                       sunctx);

    /* df = y - Q_k-1 z_k */
    SUNCheckCall(N_VLinearCombination(m, qrdata->temp_array, Q, Q[m]),
                       sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp2, -ONE, Q[m], qrdata->vtemp), sunctx);

    /* R(1:k-1,k) = s_k + z_k */
    for (j = 0; j < m; j++)
    {
      R[m * mMax + j] = R[m * mMax + j] + qrdata->temp_array[j];
    }
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = SUNCheckCallLastErr(N_VDotProd(qrdata->vtemp, qrdata->vtemp), sunctx);
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  SUNCheckCallLastErr(N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]), sunctx);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_DCGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_DCGS2(N_Vector* Q, realtype* R, N_Vector df, int m,
                          int mMax, void* QRdata)
{
  sunindextype j;
  SUNQRData qrdata  = (SUNQRData)QRdata;
  SUNDeclareContext(Q[0]->sunctx);

  SUNCheckCallLastErr(N_VScale(ONE, df, qrdata->vtemp), sunctx);  /* temp = df */

  if (m > 0)
  {
    /* R(1:k-1,k) = Q_k-1^T df_aa */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax),
                       sunctx);
    /* Delayed reorthogonalization */
    if (m > 1)
    {
      /* s = Q_k-2^T Q(:,k-1) */
      SUNCheckCall(N_VDotProdMulti(m - 1, Q[m - 1], Q, qrdata->temp_array),
                         sunctx);

      /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
      SUNCheckCall(N_VLinearCombination(m - 1, qrdata->temp_array, Q,
                                              qrdata->vtemp2),
                         sunctx);
      SUNCheckCallLastErr(N_VLinearSum(ONE, Q[m - 1], -ONE, qrdata->vtemp2, Q[m - 1]),
                                sunctx);

      /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
      for (j = 0; j < m - 1; j++)
      {
        R[(m - 1) * mMax + j] = R[(m - 1) * mMax + j] + qrdata->temp_array[j];
      }
    }

    /* df = df - Q(:,k-1) R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2),
                       sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp), sunctx);
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = SUNCheckCallLastErr(N_VDotProd(qrdata->vtemp, qrdata->vtemp), sunctx);
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  SUNCheckCallLastErr(N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]), sunctx);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_DCGS2_SB
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration which utilizes a single buffer reduction.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_DCGS2_SB(N_Vector* Q, realtype* R, N_Vector df, int m,
                             int mMax, void* QRdata)
{
  sunindextype j;
  SUNQRData qrdata  = (SUNQRData)QRdata;
  SUNDeclareContext(Q[0]->sunctx);

  SUNCheckCallLastErr(N_VScale(ONE, df, qrdata->vtemp), sunctx);  /* temp = df */

  if (m > 0)
  {
    if (m == 1)
    {
      /* R(1:k-1,k) = Q_k-1^T df_aa */
      SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax),
                         sunctx);
    }
    /* Delayed reorthogonalization */
    else if (m > 1)
    {
      /* R(1:k-1,k) = Q_k-1^T df_aa */
      /* Put R values at beginning of temp array */
      SUNCheckCall(N_VDotProdMultiLocal(m, qrdata->vtemp, Q,
                                              qrdata->temp_array),
                         sunctx);

      /* s = Q_k-2^T Q(:,k-1) */
      SUNCheckCall(N_VDotProdMultiLocal(m - 1, Q[m - 1], Q,
                                              qrdata->temp_array + m),
                         sunctx);
      SUNCheckCall(N_VDotProdMultiAllReduce(m + m - 1, qrdata->vtemp,
                                                  qrdata->temp_array),
                         sunctx);

      /* Move R values to R */
      for (j = 0; j < m; j++) { R[m * mMax + j] = qrdata->temp_array[j]; }

      /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
      SUNCheckCall(N_VLinearCombination(m - 1, qrdata->temp_array + m, Q,
                                              qrdata->vtemp2),
                         sunctx);
      SUNCheckCallLastErr(N_VLinearSum(ONE, Q[m - 1], -ONE, qrdata->vtemp2, Q[m - 1]), sunctx);

      /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
      for (j = 0; j < m - 1; j++)
      {
        R[(m - 1) * mMax + j] = R[(m - 1) * mMax + j] + qrdata->temp_array[m + j];
      }
    }

    /* df = df - Q(:,k-1) R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2),
                       sunctx);
    SUNCheckCallLastErr(N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp), sunctx);
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = SUNCheckCallLastErr(N_VDotProd(qrdata->vtemp, qrdata->vtemp), sunctx);
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  SUNCheckCallLastErr(N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]), sunctx);

  return SUN_SUCCESS;
}
