/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on sundials_spbcgs.c code, written by Peter Brown and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for the SPBCGS implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sundials/sundials_math.h>

#include "sundials/sundials_context.h"
#include "sundials/sundials_errors.h"
#include "sundials_context_impl.h"
#include "sundials_logger_impl.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * SPBCGS solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SPBCGS_CONTENT(S)  ( (SUNLinearSolverContent_SPBCGS)(S->content) )
#define PRETYPE(S)         ( SPBCGS_CONTENT(S)->pretype )
#define LASTFLAG(S)        ( SPBCGS_CONTENT(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPBCGS linear solver
 */

SUNLinearSolver SUNLinSol_SPBCGS(N_Vector y, int pretype, int maxl, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_SPBCGS content;

  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE)  && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH))
    pretype = SUN_PREC_NONE;
  if (maxl <= 0)
    maxl = SUNSPBCGS_MAXL_DEFAULT;

  /* check that the supplied N_Vector supports all requisite operations */
  SUNAssert((y->ops->nvclone) && (y->ops->nvdestroy) && (y->ops->nvlinearsum) &&
              (y->ops->nvprod) && (y->ops->nvdiv) && (y->ops->nvscale) &&
              (y->ops->nvdotprod),
            SUN_ERR_ARG_OUTOFRANGE, sunctx);

  /* Create linear solver */
  S = NULL;
  S = SUNCheckCallLastErrReturnNull(SUNLinSolNewEmpty(sunctx), sunctx);

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_SPBCGS;
  S->ops->getid             = SUNLinSolGetID_SPBCGS;
  S->ops->setatimes         = SUNLinSolSetATimes_SPBCGS;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_SPBCGS;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_SPBCGS;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_SPBCGS;
  S->ops->initialize        = SUNLinSolInitialize_SPBCGS;
  S->ops->setup             = SUNLinSolSetup_SPBCGS;
  S->ops->solve             = SUNLinSolSolve_SPBCGS;
  S->ops->numiters          = SUNLinSolNumIters_SPBCGS;
  S->ops->resnorm           = SUNLinSolResNorm_SPBCGS;
  S->ops->resid             = SUNLinSolResid_SPBCGS;
  S->ops->lastflag          = SUNLinSolLastFlag_SPBCGS;
  S->ops->space             = SUNLinSolSpace_SPBCGS;
  S->ops->free              = SUNLinSolFree_SPBCGS;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPBCGS) malloc(sizeof *content);
  SUNAssert(content, SUN_ERR_MALLOC_FAIL, sunctx);

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag   = 0;
  content->maxl        = maxl;
  content->pretype     = pretype;
  content->zeroguess   = SUNFALSE;
  content->numiters    = 0;
  content->resnorm     = ZERO;
  content->r_star      = NULL;
  content->r           = NULL;
  content->p           = NULL;
  content->q           = NULL;
  content->u           = NULL;
  content->Ap          = NULL;
  content->vtemp       = NULL;
  content->s1          = NULL;
  content->s2          = NULL;
  content->ATimes      = NULL;
  content->ATData      = NULL;
  content->Psetup      = NULL;
  content->Psolve      = NULL;
  content->PData       = NULL;
  content->print_level = 0;
  content->info_file   = stdout;
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  content->info_file   = (sunctx->logger->info_fp) ? sunctx->logger->info_fp : stdout;
#endif

  /* Allocate content */
  content->r_star = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  content->r = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  content->p = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  content->q = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  content->u = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  content->Ap = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  content->vtemp = SUNCheckCallLastErrReturnNull(N_VClone(y), sunctx);

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the type of preconditioning for SPBCGS to use
 */

SUNErrCode SUNLinSol_SPBCGSSetPrecType(SUNLinearSolver S, int pretype)
{
  /* Check for legal pretype */
  SUNAssert((pretype == SUN_PREC_NONE) || (pretype == SUN_PREC_LEFT) ||
              (pretype == SUN_PREC_RIGHT) || (pretype == SUN_PREC_BOTH),
            SUN_ERR_ARG_CORRUPT, S->sunctx);

  /* Set pretype */
  PRETYPE(S) = pretype;
  return SUN_SUCCESS;
}


/* ----------------------------------------------------------------------------
 * Function to set the maximum number of iterations for SPBCGS to use
 */

SUNErrCode SUNLinSol_SPBCGSSetMaxl(SUNLinearSolver S, int maxl)
{
  /* Check for legal pretype */
  if (maxl <= 0)
    maxl = SUNSPBCGS_MAXL_DEFAULT;

  /* Set pretype */
  SPBCGS_CONTENT(S)->maxl = maxl;
  return SUN_SUCCESS;
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPBCGS(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_ITERATIVE);
}


SUNLinearSolver_ID SUNLinSolGetID_SPBCGS(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_SPBCGS);
}


SUNErrCode SUNLinSolInitialize_SPBCGS(SUNLinearSolver S)
{
  SUNContext sunctx = S->sunctx;

  if (SPBCGS_CONTENT(S)->maxl <= 0)
    SPBCGS_CONTENT(S)->maxl = SUNSPBCGS_MAXL_DEFAULT;

  SUNAssert(SPBCGS_CONTENT(S)->ATimes, SUN_ERR_ARG_CORRUPT, S->sunctx);

  if ( (PRETYPE(S) != SUN_PREC_LEFT) &&
       (PRETYPE(S) != SUN_PREC_RIGHT) &&
       (PRETYPE(S) != SUN_PREC_BOTH) )
    PRETYPE(S) = SUN_PREC_NONE;

  SUNAssert((PRETYPE(S) == SUN_PREC_NONE) || SPBCGS_CONTENT(S)->Psolve, SUN_ERR_ARG_CORRUPT, sunctx);
  
  /* no additional memory to allocate */

  /* return with success */
  return(LASTFLAG(S));
}

SUNErrCode SUNLinSolSetATimes_SPBCGS(SUNLinearSolver S, void* ATData,
                                    SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  SPBCGS_CONTENT(S)->ATimes = ATimes;
  SPBCGS_CONTENT(S)->ATData = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetPreconditioner_SPBCGS(SUNLinearSolver S, void* PData,
                                             SUNPSetupFn Psetup,
                                             SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  SPBCGS_CONTENT(S)->Psetup = Psetup;
  SPBCGS_CONTENT(S)->Psolve = Psolve;
  SPBCGS_CONTENT(S)->PData = PData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_SPBCGS(SUNLinearSolver S, N_Vector s1,
                                             N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  SPBCGS_CONTENT(S)->s1 = s1;
  SPBCGS_CONTENT(S)->s2 = s2;
  return SUN_SUCCESS;
}


SUNErrCode SUNLinSolSetZeroGuess_SPBCGS(SUNLinearSolver S, booleantype onoff)
{
  /* set flag indicating a zero initial guess */
  SPBCGS_CONTENT(S)->zeroguess = onoff;
  return SUN_SUCCESS;
}


SUNLsStatus SUNLinSolSetup_SPBCGS(SUNLinearSolver S, SUNMatrix A)
{
  int ier;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to SPBCGS memory structures */
  Psetup = SPBCGS_CONTENT(S)->Psetup;
  PData = SPBCGS_CONTENT(S)->PData;

  /* no solver-specific setup is required, but if user-supplied
     Psetup routine exists, call that here */
  if (Psetup != NULL) {
    ier = Psetup(PData);
    if (ier != 0) {
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_PSET_FAIL_UNREC : SUNLS_PSET_FAIL_REC;
      return(LASTFLAG(S));
    }
  }

  /* return with success */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

SUNLsStatus SUNLinSolSolve_SPBCGS(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                  N_Vector b, realtype delta)
{
  /* local data and shortcut variables */
  realtype alpha, beta, omega, omega_denom, beta_num, beta_denom, r_norm, rho;
  N_Vector r_star, r, p, q, u, Ap, vtemp;
  booleantype preOnLeft, preOnRight, scale_x, scale_b, converged;
  booleantype *zeroguess;
  int l, l_max, ier;
  void *A_data, *P_data;
  N_Vector sx, sb;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  realtype *res_norm;
  int *nli;

  /* local variables for fused vector operations */
  realtype cv[3];
  N_Vector Xv[3];

  /* Make local shorcuts to solver variables. */
  l_max        = SPBCGS_CONTENT(S)->maxl;
  r_star       = SPBCGS_CONTENT(S)->r_star;
  r            = SPBCGS_CONTENT(S)->r;
  p            = SPBCGS_CONTENT(S)->p;
  q            = SPBCGS_CONTENT(S)->q;
  u            = SPBCGS_CONTENT(S)->u;
  Ap           = SPBCGS_CONTENT(S)->Ap;
  vtemp        = SPBCGS_CONTENT(S)->vtemp;
  sb           = SPBCGS_CONTENT(S)->s1;
  sx           = SPBCGS_CONTENT(S)->s2;
  A_data       = SPBCGS_CONTENT(S)->ATData;
  P_data       = SPBCGS_CONTENT(S)->PData;
  atimes       = SPBCGS_CONTENT(S)->ATimes;
  psolve       = SPBCGS_CONTENT(S)->Psolve;
  zeroguess    = &(SPBCGS_CONTENT(S)->zeroguess);
  nli          = &(SPBCGS_CONTENT(S)->numiters);
  res_norm     = &(SPBCGS_CONTENT(S)->resnorm);

  /* Initialize counters and convergence flag */
  *nli = 0;
  converged = SUNFALSE;

  /* set booleantype flags for internal solver options */
  preOnLeft  = ( (PRETYPE(S) == SUN_PREC_LEFT) ||
                 (PRETYPE(S) == SUN_PREC_BOTH) );
  preOnRight = ( (PRETYPE(S) == SUN_PREC_RIGHT) ||
                 (PRETYPE(S) == SUN_PREC_BOTH) );
  scale_x = (sx != NULL);
  scale_b = (sb != NULL);

  /* Check for unsupported use case */
  if (preOnRight && !(*zeroguess)) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(SUNLS_ILL_INPUT);
  }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  if (SPBCGS_CONTENT(S)->print_level && SPBCGS_CONTENT(S)->info_file
      && (SPBCGS_CONTENT(S)->info_file != S->sunctx->logger->info_fp))
    fprintf(SPBCGS_CONTENT(S)->info_file, "SUNLINSOL_SPBCGS:\n");
#endif

  /* Check if Atimes function has been set */
  SUNAssert(atimes, SUN_ERR_ARG_CORRUPT, S->sunctx);

  /* If preconditioning, check if psolve has been set */
  SUNAssert(!(preOnLeft || preOnRight) || psolve, SUN_ERR_ARG_CORRUPT, S->sunctx);

  /* Set r_star to initial (unscaled) residual r_0 = b - A*x_0 */

  if (*zeroguess) {
    SUNCheckCallLastErr(N_VScale(ONE, b, r_star), S->sunctx);
  } else {
    ier = atimes(A_data, x, r_star);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }
    SUNCheckCallLastErr(N_VLinearSum(ONE, b, -ONE, r_star, r_star), S->sunctx);
  }

  /* Apply left preconditioner and b-scaling to r_star = r_0 */

  if (preOnLeft) {
    ier = psolve(P_data, r_star, r, delta, SUN_PREC_LEFT);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
      return(LASTFLAG(S));
    }
  }
  else {
    SUNCheckCallLastErr(N_VScale(ONE, r_star, r), S->sunctx);
  }

  if (scale_b) {
    SUNCheckCallLastErr(N_VProd(sb, r, r_star), S->sunctx);
  } else {
    SUNCheckCallLastErr(N_VScale(ONE, r, r_star), S->sunctx);
  }

  /* Initialize beta_denom to the dot product of r0 with r0 */

  beta_denom = SUNCheckCallLastErr(N_VDotProd(r_star, r_star), S->sunctx);

  /* Set r_norm to L2 norm of r_star = sb P1_inv r_0, and
     return if small */

  *res_norm = r_norm = rho = SUNRsqrt(beta_denom);

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  /* print the initial residual */
  if (SPBCGS_CONTENT(S)->print_level && SPBCGS_CONTENT(S)->info_file
     && (SPBCGS_CONTENT(S)->info_file != S->sunctx->logger->info_fp))
  {
    fprintf(SPBCGS_CONTENT(S)->info_file,
            SUNLS_MSG_RESIDUAL,
            (long int) 0, *res_norm);
  }
  SUNLogger_QueueMsg(S->sunctx->logger, SUN_LOGLEVEL_INFO,
    "SUNLinSolSolve_SPBCGS", "initial-residual",
    "nli = %li, resnorm = %.16g", (long int) 0, *res_norm);
#endif

  if (r_norm <= delta) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
  }

  /* Copy r_star to r and p */

  SUNCheckCallLastErr(N_VScale(ONE, r_star, r), S->sunctx);
  SUNCheckCallLastErr(N_VScale(ONE, r_star, p), S->sunctx);

  /* Set x = sx x if non-zero guess */
  if (scale_x && !(*zeroguess)) {
    SUNCheckCallLastErr(N_VProd(sx, x, x), S->sunctx);
  }

  /* Begin main iteration loop */

  for(l = 0; l < l_max; l++) {

    (*nli)++;

    /* Generate Ap = A-tilde p, where A-tilde = sb P1_inv A P2_inv sx_inv */

    /*   Apply x-scaling: vtemp = sx_inv p */

    if (scale_x) N_VDiv(p, sx, vtemp);
    else SUNCheckCallLastErr(N_VScale(ONE, p, vtemp), S->sunctx);

    /*   Apply right preconditioner: vtemp = P2_inv sx_inv p */

    if (preOnRight) {
      SUNCheckCallLastErr(N_VScale(ONE, vtemp, Ap), S->sunctx);
      ier = psolve(P_data, Ap, vtemp, delta, SUN_PREC_RIGHT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    }

    /*   Apply A: Ap = A P2_inv sx_inv p */

    ier = atimes(A_data, vtemp, Ap );
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }

    /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */

    if (preOnLeft) {
      ier = psolve(P_data, Ap, vtemp, delta, SUN_PREC_LEFT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    } else {
      SUNCheckCallLastErr(N_VScale(ONE, Ap, vtemp), S->sunctx);
    }

    /*   Apply b-scaling: Ap = sb P1_inv A P2_inv sx_inv p */

    if (scale_b) N_VProd(sb, vtemp, Ap);
    else SUNCheckCallLastErr(N_VScale(ONE, vtemp, Ap), S->sunctx);

    /* Calculate alpha = <r,r_star>/<Ap,r_star> */

    alpha = SUNCheckCallLastErr(N_VDotProd(Ap, r_star), S->sunctx);
    alpha = beta_denom / alpha;

    /* Update q = r - alpha*Ap = r - alpha*(sb P1_inv A P2_inv sx_inv p) */

    SUNCheckCallLastErr(N_VLinearSum(ONE, r, -alpha, Ap, q), S->sunctx);

    /* Generate u = A-tilde q */

    /*   Apply x-scaling: vtemp = sx_inv q */

    if (scale_x) N_VDiv(q, sx, vtemp);
    else SUNCheckCallLastErr(N_VScale(ONE, q, vtemp), S->sunctx);

    /*   Apply right preconditioner: vtemp = P2_inv sx_inv q */

    if (preOnRight) {
      SUNCheckCallLastErr(N_VScale(ONE, vtemp, u), S->sunctx);
      ier = psolve(P_data, u, vtemp, delta, SUN_PREC_RIGHT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    }

    /*   Apply A: u = A P2_inv sx_inv u */

    ier = atimes(A_data, vtemp, u );
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }

    /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */

    if (preOnLeft) {
      ier = psolve(P_data, u, vtemp, delta, SUN_PREC_LEFT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    } else {
      SUNCheckCallLastErr(N_VScale(ONE, u, vtemp), S->sunctx);
    }

    /*   Apply b-scaling: u = sb P1_inv A P2_inv sx_inv u */

    if (scale_b) N_VProd(sb, vtemp, u);
    else SUNCheckCallLastErr(N_VScale(ONE, vtemp, u), S->sunctx);

    /* Calculate omega = <u,q>/<u,u> */

    omega_denom = SUNCheckCallLastErr(N_VDotProd(u, u), S->sunctx);
    if (omega_denom == ZERO) omega_denom = ONE;
    omega = SUNCheckCallLastErr(N_VDotProd(u, q), S->sunctx);
    omega /= omega_denom;

    /* Update x = x + alpha*p + omega*q */
    if (l == 0 && *zeroguess) {
      SUNCheckCallLastErr(N_VLinearSum(alpha, p, omega, q, x), S->sunctx);
    } else {
      cv[0] = ONE;
      Xv[0] = x;

      cv[1] = alpha;
      Xv[1] = p;

      cv[2] = omega;
      Xv[2] = q;

      ier = N_VLinearCombination(3, cv, Xv, x);
      SUNCheck(ier == SUN_SUCCESS, ier, S->sunctx);

    }

    /* Update the residual r = q - omega*u */

    SUNCheckCallLastErr(N_VLinearSum(ONE, q, -omega, u, r), S->sunctx);

    /* Set rho = norm(r) and check convergence */

    *res_norm = rho = SUNCheckCallLastErr(SUNRsqrt(N_VDotProd(r, r)), S->sunctx);

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    /* print current iteration number and the residual */
    if (SPBCGS_CONTENT(S)->print_level && SPBCGS_CONTENT(S)->info_file
        && (SPBCGS_CONTENT(S)->info_file != S->sunctx->logger->info_fp))
    {
      fprintf(SPBCGS_CONTENT(S)->info_file,
              SUNLS_MSG_RESIDUAL,
              (long int) *nli, *res_norm);
    }
    SUNLogger_QueueMsg(S->sunctx->logger, SUN_LOGLEVEL_INFO,
      "SUNLinSolSolve_SPBCGS", "iterate-residual",
      "nli = %li, resnorm = %.16g", (long int) 0, *res_norm);
#endif

    if (rho <= delta) {
      converged = SUNTRUE;
      break;
    }

    /* Not yet converged, continue iteration */
    /* Update beta = <rnew,r_star> / <rold,r_start> * alpha / omega */

    beta_num = SUNCheckCallLastErr(N_VDotProd(r, r_star), S->sunctx);
    beta = ((beta_num / beta_denom) * (alpha / omega));

    /* Update p = r + beta*(p - omega*Ap) = beta*p - beta*omega*Ap + r */
    cv[0] = beta;
    Xv[0] = p;

    cv[1] = -alpha*(beta_num / beta_denom);
    Xv[1] = Ap;

    cv[2] = ONE;
    Xv[2] = r;

    ier = N_VLinearCombination(3, cv, Xv, p);
    SUNCheck(ier == SUN_SUCCESS, ier, S->sunctx);

    /* udpate beta_denom for next iteration */
    beta_denom = beta_num;
  }

  /* Main loop finished */

  if ((converged == SUNTRUE) || (rho < r_norm)) {

    /* Apply the x-scaling and right preconditioner: x = P2_inv sx_inv x */

    if (scale_x) {
      SUNCheckCallLastErr(N_VDiv(x, sx, x), S->sunctx);
    }
    if (preOnRight) {
      ier = psolve(P_data, x, vtemp, delta, SUN_PREC_RIGHT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
      SUNCheckCallLastErr(N_VScale(ONE, vtemp, x), S->sunctx);
    }

    *zeroguess = SUNFALSE;
    if (converged == SUNTRUE)
      LASTFLAG(S) = SUNLS_SUCCESS;
    else
      LASTFLAG(S) = SUNLS_RES_REDUCED;
    return(LASTFLAG(S));

  }
  else {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_CONV_FAIL;
    return(LASTFLAG(S));
  }
}


int SUNLinSolNumIters_SPBCGS(SUNLinearSolver S) SUNDIALS_NOEXCEPT
{
  /* return the stored 'numiters' value */
  return (SPBCGS_CONTENT(S)->numiters);
}


realtype SUNLinSolResNorm_SPBCGS(SUNLinearSolver S) SUNDIALS_NOEXCEPT
{
  /* return the stored 'resnorm' value */
  return (SPBCGS_CONTENT(S)->resnorm);
}


N_Vector SUNLinSolResid_SPBCGS(SUNLinearSolver S) SUNDIALS_NOEXCEPT
{
  /* return the stored 'r' vector */
  return (SPBCGS_CONTENT(S)->r);
}


sunindextype SUNLinSolLastFlag_SPBCGS(SUNLinearSolver S) SUNDIALS_NOEXCEPT
{
  /* return the stored 'last_flag' value */
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_SPBCGS(SUNLinearSolver S, long int* lenrwLS,
                                 long int* leniwLS)
{
  sunindextype liw1, lrw1;
  if (SPBCGS_CONTENT(S)->vtemp->ops->nvspace) {
    SUNCheckCallLastErr(N_VSpace(SPBCGS_CONTENT(S)->vtemp, &lrw1, &liw1), S->sunctx);
  } else {
    lrw1 = liw1 = 0;
  }
  *lenrwLS = lrw1*9;
  *leniwLS = liw1*9;
  return SUN_SUCCESS;
}


SUNErrCode SUNLinSolFree_SPBCGS(SUNLinearSolver S)
{
  if (S->content) {
    /* delete items from within the content structure */
    if (SPBCGS_CONTENT(S)->r_star) {
      N_VDestroy(SPBCGS_CONTENT(S)->r_star);
      SPBCGS_CONTENT(S)->r_star = NULL;
    }
    if (SPBCGS_CONTENT(S)->r) {
      N_VDestroy(SPBCGS_CONTENT(S)->r);
      SPBCGS_CONTENT(S)->r = NULL;
    }
    if (SPBCGS_CONTENT(S)->p) {
      N_VDestroy(SPBCGS_CONTENT(S)->p);
      SPBCGS_CONTENT(S)->p = NULL;
    }
    if (SPBCGS_CONTENT(S)->q) {
      N_VDestroy(SPBCGS_CONTENT(S)->q);
      SPBCGS_CONTENT(S)->q = NULL;
    }
    if (SPBCGS_CONTENT(S)->u) {
      N_VDestroy(SPBCGS_CONTENT(S)->u);
      SPBCGS_CONTENT(S)->u = NULL;
    }
    if (SPBCGS_CONTENT(S)->Ap) {
      N_VDestroy(SPBCGS_CONTENT(S)->Ap);
      SPBCGS_CONTENT(S)->Ap = NULL;
    }
    if (SPBCGS_CONTENT(S)->vtemp) {
      N_VDestroy(SPBCGS_CONTENT(S)->vtemp);
      SPBCGS_CONTENT(S)->vtemp = NULL;
    }
    free(S->content); S->content = NULL;
  }
  if (S->ops) { free(S->ops); S->ops = NULL; }
  free(S); S = NULL;
  return SUN_SUCCESS;
}


SUNErrCode SUNLinSolSetInfoFile_SPBCGS(SUNLinearSolver S,
                                FILE* info_file)
{
  SPBCGS_CONTENT(S)->info_file = info_file;
  return(SUNLS_SUCCESS);
}


SUNErrCode SUNLinSolSetPrintLevel_SPBCGS(SUNLinearSolver S,
                                  int print_level)
{
  /* check for valid print level */
  SUNAssert(print_level >= 0 && print_level <= 1, SUN_ERR_ARG_OUTOFRANGE,
            S->sunctx);
  SPBCGS_CONTENT(S)->print_level = print_level;
  return(SUNLS_SUCCESS);
}
