/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This file implements fused stub kernels for CVODE.
 * -----------------------------------------------------------------
 */


#include "cvode_diag_impl.h"
#include "cvode_impl.h"

#define ZERO   RCONST(0.0)
#define PT1    RCONST(0.1)
#define FRACT  RCONST(0.1)
#define ONEPT5 RCONST(1.50)
#define ONE    RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SS.
 * -----------------------------------------------------------------
 */

int cvEwtSetSS_fused(const booleantype atolmin0,
                     const realtype reltol,
                     const realtype Sabstol,
                     const N_Vector ycur,
                     N_Vector tempv,
                     N_Vector weight)
{
  SUNCheckCallLastErrNoRet(N_VAbs(ycur, tempv), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VScale(reltol, tempv, tempv), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VAddConst(tempv, Sabstol, tempv), SUNCTX);
  if (atolmin0) {
    sunrealtype min = SUNCheckCallLastErrNoRet(N_VMin(tempv), SUNCTX);
    if (min <= ZERO) return(-1);
  }
  SUNCheckCallLastErrNoRet(N_VInv(tempv, weight), SUNCTX);
  return 0;
}

/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SV.
 * -----------------------------------------------------------------
 */


int cvEwtSetSV_fused(const booleantype atolmin0,
                     const realtype reltol,
                     const N_Vector Vabstol,
                     const N_Vector ycur,
                     N_Vector tempv,
                     N_Vector weight)
{
  SUNCheckCallLastErrNoRet(N_VAbs(ycur, tempv), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VLinearSum(reltol, tempv, ONE,
               Vabstol, tempv), SUNCTX);
  if (atolmin0) {
    sunrealtype min = SUNCheckCallLastErrNoRet(N_VMin(tempv), SUNCTX);
    if (min <= ZERO) return(-1);
  }
  SUNCheckCallLastErrNoRet(N_VInv(tempv, weight), SUNCTX);
  return 0;
}


/*
 * -----------------------------------------------------------------
 * Determine if the constraints of the problem are satisfied by
 * the proposed step.
 * -----------------------------------------------------------------
 */


int cvCheckConstraints_fused(const N_Vector c,
                             const N_Vector ewt,
                             const N_Vector y,
                             const N_Vector mm,
                             N_Vector tmp)
{
  SUNCheckCallLastErrNoRet(N_VCompare(ONEPT5, c, tmp), SUNCTX);           /* a[i]=1 when |c[i]|=2  */
  SUNCheckCallLastErrNoRet(N_VProd(tmp, c, tmp), SUNCTX);                 /* a * c                 */
  SUNCheckCallLastErrNoRet(N_VDiv(tmp, ewt, tmp), SUNCTX);                /* a * c * wt            */
  SUNCheckCallLastErrNoRet(N_VLinearSum(ONE, y, -PT1, tmp, tmp), SUNCTX); /* y - 0.1 * a * c * wt  */
  SUNCheckCallLastErrNoRet(N_VProd(tmp, mm, tmp), SUNCTX);                /* v = mm*(y-0.1*a*c*wt) */
  return 0;
}


/*
 * -----------------------------------------------------------------
 * Compute the nonlinear residual.
 * -----------------------------------------------------------------
 */


int cvNlsResid_fused(const realtype rl1,
                     const realtype ngamma,
                     const N_Vector zn1,
                     const N_Vector ycor,
                     const N_Vector ftemp,
                     N_Vector res)
{
  SUNCheckCallLastErrNoRet(N_VLinearSum(rl1, zn1, ONE, ycor, res), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VLinearSum(ngamma, ftemp, ONE, res, res), SUNCTX);
  return 0;
}

/*
 * -----------------------------------------------------------------
 * Form y with perturbation = FRACT*(func. iter. correction)
 * -----------------------------------------------------------------
 */

int cvDiagSetup_formY(const realtype h,
                      const realtype r,
                      const N_Vector fpred,
                      const N_Vector zn1,
                      const N_Vector ypred,
                      N_Vector ftemp,
                      N_Vector y)
{
  SUNCheckCallLastErrNoRet(N_VLinearSum(h, fpred, -ONE, zn1, ftemp), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VLinearSum(r, ftemp, ONE, ypred, y), SUNCTX);
  return 0;
}

/*
 * -----------------------------------------------------------------
 * Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i)
 * protecting against deltay_i being at roundoff level.
 * -----------------------------------------------------------------
 */

int cvDiagSetup_buildM(const realtype fract,
                       const realtype uround,
                       const realtype h,
                       const N_Vector ftemp,
                       const N_Vector fpred,
                       const N_Vector ewt,
                       N_Vector bit,
                       N_Vector bitcomp,
                       N_Vector y,
                       N_Vector M)
{
  SUNCheckCallLastErrNoRet(N_VLinearSum(ONE, M, -ONE, fpred, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VLinearSum(FRACT, ftemp, -h, M, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VProd(ftemp, ewt, y), SUNCTX);
  /* Protect against deltay_i being at roundoff level */
  SUNCheckCallLastErrNoRet(N_VCompare(uround, y, bit), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VAddConst(bit, -ONE, bitcomp), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VProd(ftemp, bit, y), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VLinearSum(FRACT, y, -ONE, bitcomp, y), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VDiv(M, y, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VProd(M, bit, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VLinearSum(ONE, M, -ONE, bitcomp, M), SUNCTX);
  return 0;
}


/*
 * -----------------------------------------------------------------
 *  Update M with changed gamma so that M = I - gamma*J.
 * -----------------------------------------------------------------
 */

int cvDiagSolve_updateM(const realtype r, N_Vector M)
{
  SUNCheckCallLastErrNoRet(N_VInv(M, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VAddConst(M, -ONE, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VScale(r, M, M), SUNCTX);
  SUNCheckCallLastErrNoRet(N_VAddConst(M, ONE, M), SUNCTX);
  return 0;
}