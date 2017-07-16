/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * Based on sundials_spbcgs.c code, written by Peter Brown and 
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * This is the implementation file for the SPBCGS implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPBCGS linear solver
 */

SUNLinearSolver SUNSPBCGS(N_Vector y, int pretype, int maxl)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_SPBCGS content;
  
  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH))
    pretype = PREC_NONE;
  if (maxl <= 0)
    maxl = SUNSPBCGS_MAXL_DEFAULT;

  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_SPBCGS;
  ops->setatimes         = SUNLinSolSetATimes_SPBCGS;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_SPBCGS;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_SPBCGS;
  ops->initialize        = SUNLinSolInitialize_SPBCGS;
  ops->setup             = SUNLinSolSetup_SPBCGS;
  ops->solve             = SUNLinSolSolve_SPBCGS;
  ops->numiters          = SUNLinSolNumIters_SPBCGS;
  ops->resnorm           = SUNLinSolResNorm_SPBCGS;
  ops->numpsolves        = SUNLinSolNumPSolves_SPBCGS;
  ops->lastflag          = SUNLinSolLastFlag_SPBCGS;  
  ops->free              = SUNLinSolFree_SPBCGS;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPBCGS) malloc(sizeof(struct _SUNLinearSolverContent_SPBCGS));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  content->maxl = maxl;
  content->pretype = pretype;
  content->numiters = 0;
  content->numpsolves = 0;
  content->resnorm = ZERO;
  content->r_star = N_VClone(y);
  if (content->r_star == NULL)  return NULL;
  content->r = N_VClone(y);
  if (content->r == NULL)  return NULL;
  content->p = N_VClone(y);
  if (content->p == NULL)  return NULL;
  content->q = N_VClone(y);
  if (content->q == NULL)  return NULL;
  content->u = N_VClone(y);
  if (content->u == NULL)  return NULL;
  content->Ap = N_VClone(y);
  if (content->Ap == NULL)  return NULL;
  content->vtemp = N_VClone(y);
  if (content->vtemp == NULL)  return NULL;
  content->s1 = NULL;
  content->s2 = NULL;
  content->ATSetup = NULL;
  content->ATimes = NULL;
  content->ATData = NULL;
  content->Psetup = NULL;
  content->Psolve = NULL;
  content->PData = NULL;

  /* Attach content and ops */
  S->content = content;
  S->ops     = ops;

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the type of preconditioning for SPBCGS to use 
 */

SUNDIALS_EXPORT int SUNSPBCGSSetPrecType(SUNLinearSolver S, int pretype) 
{
  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    return(SPBCG_ILL_INPUT);
  }

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPBCG_MEM_NULL);

  /* Set pretype */
  SLS_CONTENT_SPBCGS(S)->pretype = pretype;
  return(SPBCG_SUCCESS);
}


/* ----------------------------------------------------------------------------
 * Function to return the SPBCGS workspace sizes
 */

SUNDIALS_EXPORT int SUNSPBCGSGetWorkspace(SUNLinearSolver S, 
                                         long int *lenrwLS, 
                                         long int *leniwLS)
{
  long int liw1, lrw1;

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPBCG_MEM_NULL);

  N_VSpace(SLS_CONTENT_SPBCGS(S)->vtemp, &lrw1, &liw1);
  *lenrwLS = lrw1*9;
  *leniwLS = liw1*9;
  return(SPBCG_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPBCGS(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_ITERATIVE;
}


int SUNLinSolInitialize_SPBCGS(SUNLinearSolver S)
{
  int k;

  /* set shortcut to SPBCGS memory structure */
  if (S == NULL) return(SPBCG_MEM_NULL);  
  SUNLinearSolverContent_SPBCGS content = SLS_CONTENT_SPBCGS(S);

  /* ensure valid options */
  if ( (content->pretype != PREC_LEFT) && 
       (content->pretype != PREC_RIGHT) && 
       (content->pretype != PREC_BOTH) )
    content->pretype = PREC_NONE;
  if (content->maxl <= 0) 
    content->maxl = SUNSPBCGS_MAXL_DEFAULT;

  /* no additional memory to allocate */
  
  /* return with success */
  content->last_flag = 0;
  return 0;
}


int SUNLinSolSetATimes_SPBCGS(SUNLinearSolver S, void* ATData, 
                            ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATSetup 
     and ATimes routines and data, and return with success */
  if (S == NULL) return(SPBCG_MEM_NULL);
  SLS_CONTENT_SPBCGS(S)->ATSetup = ATSetup;
  SLS_CONTENT_SPBCGS(S)->ATimes = ATimes;
  SLS_CONTENT_SPBCGS(S)->ATData = ATData;
  SLS_CONTENT_SPBCGS(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetPreconditioner_SPBCGS(SUNLinearSolver S, void* PData,
                                    PSetupFn Psetup, PSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  if (S == NULL) return(SPBCG_MEM_NULL);
  SLS_CONTENT_SPBCGS(S)->Psetup = Psetup;
  SLS_CONTENT_SPBCGS(S)->Psolve = Psolve;
  SLS_CONTENT_SPBCGS(S)->PData = PData;
  SLS_CONTENT_SPBCGS(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetScalingVectors_SPBCGS(SUNLinearSolver S, N_Vector s1,
                                     N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors, 
     and return with success */
  if (S == NULL) return(SPBCG_MEM_NULL);
  SLS_CONTENT_SPBCGS(S)->s1 = s1;
  SLS_CONTENT_SPBCGS(S)->s2 = s2;
  SLS_CONTENT_SPBCGS(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetup_SPBCGS(SUNLinearSolver S, SUNMatrix A)
{
  /* Set shortcuts to SPBCGS memory structures */
  if (S == NULL) return(SPBCG_MEM_NULL);
  ATSetupFn ATSetup = SLS_CONTENT_SPBCGS(S)->ATSetup;
  PSetupFn Psetup = SLS_CONTENT_SPBCGS(S)->Psetup;
  void* ATData = SLS_CONTENT_SPBCGS(S)->ATData;
  void* PData = SLS_CONTENT_SPBCGS(S)->PData;
  
  /* no solver-specific setup is required, but if user-supplied 
     ATSetup or Psetup routines exist, call those here */
  if (ATSetup != NULL) {
    SLS_CONTENT_SPBCGS(S)->last_flag = ATSetup(ATData);
    if (SLS_CONTENT_SPBCGS(S)->last_flag != 0)  return 1;
  }
  if (Psetup != NULL) {
    SLS_CONTENT_SPBCGS(S)->last_flag = Psetup(PData);
    if (SLS_CONTENT_SPBCGS(S)->last_flag != 0)  return 1;
  }
  
  /* return with success */ 
  return 0;
}


int SUNLinSolSolve_SPBCGS(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                         N_Vector b, realtype delta)
{
  /* local data and shortcut variables */
  realtype alpha, beta, omega, omega_denom, beta_num, beta_denom, r_norm, rho;
  N_Vector r_star, r, p, q, u, Ap, vtemp;
  booleantype preOnLeft, preOnRight, scale_x, scale_b, converged;
  int l, l_max, ier, pretype;
  void *A_data, *P_data;
  N_Vector sx, sb;
  ATimesFn atimes;
  PSolveFn psolve;
  realtype *res_norm;
  int *nli, *nps;
  
  /* Make local shorcuts to solver variables. */
  if (S == NULL) return(SPBCG_MEM_NULL);
  l_max        = SLS_CONTENT_SPBCGS(S)->maxl;
  r_star       = SLS_CONTENT_SPBCGS(S)->r_star;
  r            = SLS_CONTENT_SPBCGS(S)->r;
  p            = SLS_CONTENT_SPBCGS(S)->p;
  q            = SLS_CONTENT_SPBCGS(S)->q;
  u            = SLS_CONTENT_SPBCGS(S)->u;
  Ap           = SLS_CONTENT_SPBCGS(S)->Ap;
  vtemp        = SLS_CONTENT_SPBCGS(S)->vtemp;
  sb           = SLS_CONTENT_SPBCGS(S)->s1;
  sx           = SLS_CONTENT_SPBCGS(S)->s2;
  A_data       = SLS_CONTENT_SPBCGS(S)->ATData;
  P_data       = SLS_CONTENT_SPBCGS(S)->PData;
  atimes       = SLS_CONTENT_SPBCGS(S)->ATimes;
  psolve       = SLS_CONTENT_SPBCGS(S)->Psolve;
  nli          = &(SLS_CONTENT_SPBCGS(S)->numiters);
  nps          = &(SLS_CONTENT_SPBCGS(S)->numpsolves);
  res_norm     = &(SLS_CONTENT_SPBCGS(S)->resnorm);

  /* Initialize counters and convergence flag */
  *nli = *nps = 0;
  converged = FALSE;

  /* set booleantype flags for internal solver options */
  preOnLeft  = ( (SLS_CONTENT_SPBCGS(S)->pretype == PREC_LEFT) || 
                 (SLS_CONTENT_SPBCGS(S)->pretype == PREC_BOTH) );
  preOnRight = ( (SLS_CONTENT_SPBCGS(S)->pretype == PREC_RIGHT) || 
                 (SLS_CONTENT_SPBCGS(S)->pretype == PREC_BOTH) );
  scale_x = (sx != NULL);
  scale_b = (sb != NULL);

  /* Set r_star to initial (unscaled) residual r_0 = b - A*x_0 */

  if (N_VDotProd(x, x) == ZERO) N_VScale(ONE, b, r_star);
  else {
    ier = atimes(A_data, x, r_star);
    if (ier != 0) {
      SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
        SPBCG_ATIMES_FAIL_UNREC : SPBCG_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPBCGS(S)->last_flag);
    }
    N_VLinearSum(ONE, b, -ONE, r_star, r_star);
  }

  /* Apply left preconditioner and b-scaling to r_star = r_0 */

  if (preOnLeft) {
    ier = psolve(P_data, r_star, r, delta, PREC_LEFT);
    (*nps)++;
    if (ier != 0) {
      SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
        SPBCG_PSOLVE_FAIL_UNREC : SPBCG_PSOLVE_FAIL_REC;
      return(SLS_CONTENT_SPBCGS(S)->last_flag);
    }
  }
  else N_VScale(ONE, r_star, r);

  if (scale_b) N_VProd(sb, r, r_star);
  else N_VScale(ONE, r, r_star);

  /* Initialize beta_denom to the dot product of r0 with r0 */

  beta_denom = N_VDotProd(r_star, r_star);

  /* Set r_norm to L2 norm of r_star = sb P1_inv r_0, and
     return if small */

  *res_norm = r_norm = rho = SUNRsqrt(beta_denom);
  if (r_norm <= delta) {
    SLS_CONTENT_SPBCGS(S)->last_flag = SPBCG_SUCCESS;
    return(SLS_CONTENT_SPBCGS(S)->last_flag);
  }

  /* Copy r_star to r and p */

  N_VScale(ONE, r_star, r);
  N_VScale(ONE, r_star, p);

  /* Begin main iteration loop */

  for(l = 0; l < l_max; l++) {

    (*nli)++;

    /* Generate Ap = A-tilde p, where A-tilde = sb P1_inv A P2_inv sx_inv */

    /*   Apply x-scaling: vtemp = sx_inv p */

    if (scale_x) N_VDiv(p, sx, vtemp);
    else N_VScale(ONE, p, vtemp);

    /*   Apply right preconditioner: vtemp = P2_inv sx_inv p */

    if (preOnRight) {
      N_VScale(ONE, vtemp, Ap);
      ier = psolve(P_data, Ap, vtemp, delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
          SPBCG_PSOLVE_FAIL_UNREC : SPBCG_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPBCGS(S)->last_flag);
      }
    }

    /*   Apply A: Ap = A P2_inv sx_inv p */

    ier = atimes(A_data, vtemp, Ap );
    if (ier != 0) {
      SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
        SPBCG_ATIMES_FAIL_UNREC : SPBCG_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPBCGS(S)->last_flag);
    }

    /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */

    if (preOnLeft) {
      ier = psolve(P_data, Ap, vtemp, delta, PREC_LEFT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
          SPBCG_PSOLVE_FAIL_UNREC : SPBCG_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPBCGS(S)->last_flag);
      }
    }
    else N_VScale(ONE, Ap, vtemp);

    /*   Apply b-scaling: Ap = sb P1_inv A P2_inv sx_inv p */

    if (scale_b) N_VProd(sb, vtemp, Ap);
    else N_VScale(ONE, vtemp, Ap);


    /* Calculate alpha = <r,r_star>/<Ap,r_star> */

    alpha = ((beta_denom / N_VDotProd(Ap, r_star)));

    /* Update q = r - alpha*Ap = r - alpha*(sb P1_inv A P2_inv sx_inv p) */

    N_VLinearSum(ONE, r, -alpha, Ap, q);

    /* Generate u = A-tilde q */

    /*   Apply x-scaling: vtemp = sx_inv q */

    if (scale_x) N_VDiv(q, sx, vtemp);
    else N_VScale(ONE, q, vtemp);

    /*   Apply right preconditioner: vtemp = P2_inv sx_inv q */

    if (preOnRight) {
      N_VScale(ONE, vtemp, u);
      ier = psolve(P_data, u, vtemp, delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
          SPBCG_PSOLVE_FAIL_UNREC : SPBCG_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPBCGS(S)->last_flag);
      }
    }

    /*   Apply A: u = A P2_inv sx_inv u */

    ier = atimes(A_data, vtemp, u );
    if (ier != 0) {
      SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
        SPBCG_ATIMES_FAIL_UNREC : SPBCG_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPBCGS(S)->last_flag);
    }

    /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */

    if (preOnLeft) {
      ier = psolve(P_data, u, vtemp, delta, PREC_LEFT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
          SPBCG_PSOLVE_FAIL_UNREC : SPBCG_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPBCGS(S)->last_flag);
      }
    }
    else N_VScale(ONE, u, vtemp);

    /*   Apply b-scaling: u = sb P1_inv A P2_inv sx_inv u */

    if (scale_b) N_VProd(sb, vtemp, u);
    else N_VScale(ONE, vtemp, u);


    /* Calculate omega = <u,q>/<u,u> */

    omega_denom = N_VDotProd(u, u);
    if (omega_denom == ZERO) omega_denom = ONE;
    omega = (N_VDotProd(u, q) / omega_denom);

    /* Update x = x + alpha*p + omega*q */

    N_VLinearSum(alpha, p, omega, q, vtemp);
    N_VLinearSum(ONE, x, ONE, vtemp, x);

    /* Update the residual r = q - omega*u */

    N_VLinearSum(ONE, q, -omega, u, r);

    /* Set rho = norm(r) and check convergence */

    *res_norm = rho = SUNRsqrt(N_VDotProd(r, r));
    if (rho <= delta) {
      converged = TRUE;
      break;
    }

    /* Not yet converged, continue iteration */
    /* Update beta = <rnew,r_star> / <rold,r_start> * alpha / omega */

    beta_num = N_VDotProd(r, r_star);
    beta = ((beta_num / beta_denom) * (alpha / omega));
    beta_denom = beta_num;

    /* Update p = r + beta*(p - omega*Ap) */

    N_VLinearSum(ONE, p, -omega, Ap, vtemp);
    N_VLinearSum(ONE, r, beta, vtemp, p);

  }

  /* Main loop finished */

  if ((converged == TRUE) || (rho < r_norm)) {

    /* Apply the x-scaling and right preconditioner: x = P2_inv sx_inv x */

    if (scale_x) N_VDiv(x, sx, x);
    if (preOnRight) {
      ier = psolve(P_data, x, vtemp, delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPBCGS(S)->last_flag = (ier < 0) ?
          SPBCG_PSOLVE_FAIL_UNREC : SPBCG_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPBCGS(S)->last_flag);
      }
      N_VScale(ONE, vtemp, x);
    }

    if (converged == TRUE) 
      SLS_CONTENT_SPBCGS(S)->last_flag = SPBCG_SUCCESS;
    else 
      SLS_CONTENT_SPBCGS(S)->last_flag = SPBCG_RES_REDUCED;
    return(SLS_CONTENT_SPBCGS(S)->last_flag);
    
  }
  else {
    SLS_CONTENT_SPBCGS(S)->last_flag = SPBCG_CONV_FAIL;
    return(SLS_CONTENT_SPBCGS(S)->last_flag);
  }
}


int SUNLinSolNumIters_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'numiters' value */
  if (S == NULL) return(SPBCG_MEM_NULL);
  return (SLS_CONTENT_SPBCGS(S)->numiters);
}


realtype SUNLinSolResNorm_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'resnorm' value */
  if (S == NULL) return(SPBCG_MEM_NULL);
  return (SLS_CONTENT_SPBCGS(S)->resnorm);
}


int SUNLinSolNumPSolves_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'numpsolves' value */
  if (S == NULL) return(SPBCG_MEM_NULL);
  return (SLS_CONTENT_SPBCGS(S)->numpsolves);
}


long int SUNLinSolLastFlag_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(SPBCG_MEM_NULL);
  return (SLS_CONTENT_SPBCGS(S)->last_flag);
}


int SUNLinSolFree_SPBCGS(SUNLinearSolver S)
{
  int k;

  if (S == NULL) return(SPBCG_MEM_NULL);

  /* delete items from within the content structure */
  if (SLS_CONTENT_SPBCGS(S)->r_star)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->r_star);
  if (SLS_CONTENT_SPBCGS(S)->r)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->r);
  if (SLS_CONTENT_SPBCGS(S)->p)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->p);
  if (SLS_CONTENT_SPBCGS(S)->q)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->q);
  if (SLS_CONTENT_SPBCGS(S)->u)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->u);
  if (SLS_CONTENT_SPBCGS(S)->Ap)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->Ap);
  if (SLS_CONTENT_SPBCGS(S)->vtemp)
    N_VDestroy(SLS_CONTENT_SPBCGS(S)->vtemp);

  /* delete generic structures */
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}
