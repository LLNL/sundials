/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * Based on sundials_sptfqmr.c code, written by Aaron Collier @ LLNL
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
 * This is the implementation file for the SPTFQMR implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_sptfqmr.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPTFQMR linear solver
 */

SUNLinearSolver SUNSPTFQMR(N_Vector y, int pretype, int maxl)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_SPTFQMR content;
  
  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH))
    pretype = PREC_NONE;
  if (maxl <= 0)
    maxl = SUNSPTFQMR_MAXL_DEFAULT;

  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_SPTFQMR;
  ops->setatimes         = SUNLinSolSetATimes_SPTFQMR;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_SPTFQMR;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_SPTFQMR;
  ops->initialize        = SUNLinSolInitialize_SPTFQMR;
  ops->setup             = SUNLinSolSetup_SPTFQMR;
  ops->solve             = SUNLinSolSolve_SPTFQMR;
  ops->numiters          = SUNLinSolNumIters_SPTFQMR;
  ops->resnorm           = SUNLinSolResNorm_SPTFQMR;
  ops->numpsolves        = SUNLinSolNumPSolves_SPTFQMR;
  ops->lastflag          = SUNLinSolLastFlag_SPTFQMR;  
  ops->free              = SUNLinSolFree_SPTFQMR;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPTFQMR) malloc(sizeof(struct _SUNLinearSolverContent_SPTFQMR));
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
  content->q = N_VClone(y);
  if (content->q == NULL)  return NULL;
  content->d = N_VClone(y);
  if (content->d == NULL)  return NULL;
  content->v = N_VClone(y);
  if (content->v == NULL)  return NULL;
  content->p = N_VClone(y);
  if (content->p == NULL)  return NULL;
  content->r = N_VCloneVectorArray(2, y);
  if (content->r == NULL)  return NULL;
  content->u = N_VClone(y);
  if (content->u == NULL)  return NULL;
  content->vtemp1 = N_VClone(y);
  if (content->vtemp1 == NULL)  return NULL;
  content->vtemp2 = N_VClone(y);
  if (content->vtemp2 == NULL)  return NULL;
  content->vtemp3 = N_VClone(y);
  if (content->vtemp3 == NULL)  return NULL;
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
 * Function to set the type of preconditioning for SPTFQMR to use 
 */

SUNDIALS_EXPORT int SUNSPTFQMRSetPrecType(SUNLinearSolver S, int pretype) 
{
  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    return(SPTFQMR_ILL_INPUT);
  }

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPTFQMR_MEM_NULL);

  /* Set pretype */
  SLS_CONTENT_SPTFQMR(S)->pretype = pretype;
  return(SPTFQMR_SUCCESS);
}


/* ----------------------------------------------------------------------------
 * Function to return the SPTFQMR workspace sizes
 */

SUNDIALS_EXPORT int SUNSPTFQMRGetWorkspace(SUNLinearSolver S, 
                                           long int *lenrwLS, 
                                           long int *leniwLS)
{
  long int liw1, lrw1;

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPTFQMR_MEM_NULL);

  N_VSpace(SLS_CONTENT_SPTFQMR(S)->vtemp1, &lrw1, &liw1);
  *lenrwLS = lrw1*11;
  *leniwLS = liw1*11;
  return(SPTFQMR_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPTFQMR(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_ITERATIVE;
}


int SUNLinSolInitialize_SPTFQMR(SUNLinearSolver S)
{
  /* set shortcut to SPTFQMR memory structure */
  if (S == NULL) return(SPTFQMR_MEM_NULL);  
  SUNLinearSolverContent_SPTFQMR content = SLS_CONTENT_SPTFQMR(S);

  /* ensure valid options */
  if ( (content->pretype != PREC_LEFT) && 
       (content->pretype != PREC_RIGHT) && 
       (content->pretype != PREC_BOTH) )
    content->pretype = PREC_NONE;
  if (content->maxl <= 0) 
    content->maxl = SUNSPTFQMR_MAXL_DEFAULT;

  /* no additional memory to allocate */
  
  /* return with success */
  content->last_flag = 0;
  return 0;
}


int SUNLinSolSetATimes_SPTFQMR(SUNLinearSolver S, void* ATData, 
                               ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATSetup 
     and ATimes routines and data, and return with success */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  SLS_CONTENT_SPTFQMR(S)->ATSetup = ATSetup;
  SLS_CONTENT_SPTFQMR(S)->ATimes = ATimes;
  SLS_CONTENT_SPTFQMR(S)->ATData = ATData;
  SLS_CONTENT_SPTFQMR(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetPreconditioner_SPTFQMR(SUNLinearSolver S, void* PData,
                                       PSetupFn Psetup, PSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  SLS_CONTENT_SPTFQMR(S)->Psetup = Psetup;
  SLS_CONTENT_SPTFQMR(S)->Psolve = Psolve;
  SLS_CONTENT_SPTFQMR(S)->PData = PData;
  SLS_CONTENT_SPTFQMR(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetScalingVectors_SPTFQMR(SUNLinearSolver S,
                                       N_Vector s1,
                                       N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors, 
     and return with success */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  SLS_CONTENT_SPTFQMR(S)->s1 = s1;
  SLS_CONTENT_SPTFQMR(S)->s2 = s2;
  SLS_CONTENT_SPTFQMR(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetup_SPTFQMR(SUNLinearSolver S, SUNMatrix A)
{
  /* Set shortcuts to SPTFQMR memory structures */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  ATSetupFn ATSetup = SLS_CONTENT_SPTFQMR(S)->ATSetup;
  PSetupFn Psetup = SLS_CONTENT_SPTFQMR(S)->Psetup;
  void* ATData = SLS_CONTENT_SPTFQMR(S)->ATData;
  void* PData = SLS_CONTENT_SPTFQMR(S)->PData;
  
  /* no solver-specific setup is required, but if user-supplied 
     ATSetup or Psetup routines exist, call those here */
  if (ATSetup != NULL) {
    SLS_CONTENT_SPTFQMR(S)->last_flag = ATSetup(ATData);
    if (SLS_CONTENT_SPTFQMR(S)->last_flag != 0)  return 1;
  }
  if (Psetup != NULL) {
    SLS_CONTENT_SPTFQMR(S)->last_flag = Psetup(PData);
    if (SLS_CONTENT_SPTFQMR(S)->last_flag != 0)  return 1;
  }
  
  /* return with success */ 
  return 0;
}


int SUNLinSolSolve_SPTFQMR(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                         N_Vector b, realtype delta)
{
  /* local data and shortcut variables */
  realtype alpha, tau, eta, beta, c, sigma, v_bar, omega;
  realtype rho[2];
  realtype r_init_norm, r_curr_norm;
  realtype temp_val;
  booleantype preOnLeft, preOnRight, scale_x, scale_b, converged;
  booleantype b_ok;
  int n, m, ier, l_max;
  void *A_data, *P_data;
  ATimesFn atimes;
  PSolveFn psolve;
  realtype *res_norm;
  int *nli, *nps;
  N_Vector sx, sb, r_star, q, d, v, p, *r, u, vtemp1, vtemp2, vtemp3;
  
  /* Make local shorcuts to solver variables. */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  l_max        = SLS_CONTENT_SPTFQMR(S)->maxl;
  r_star       = SLS_CONTENT_SPTFQMR(S)->r_star;
  q            = SLS_CONTENT_SPTFQMR(S)->q;
  d            = SLS_CONTENT_SPTFQMR(S)->d;
  v            = SLS_CONTENT_SPTFQMR(S)->v;
  p            = SLS_CONTENT_SPTFQMR(S)->p;
  r            = SLS_CONTENT_SPTFQMR(S)->r;
  u            = SLS_CONTENT_SPTFQMR(S)->u;
  vtemp1       = SLS_CONTENT_SPTFQMR(S)->vtemp1;
  vtemp2       = SLS_CONTENT_SPTFQMR(S)->vtemp2;
  vtemp3       = SLS_CONTENT_SPTFQMR(S)->vtemp3;
  sb           = SLS_CONTENT_SPTFQMR(S)->s1;
  sx           = SLS_CONTENT_SPTFQMR(S)->s2;
  A_data       = SLS_CONTENT_SPTFQMR(S)->ATData;
  P_data       = SLS_CONTENT_SPTFQMR(S)->PData;
  atimes       = SLS_CONTENT_SPTFQMR(S)->ATimes;
  psolve       = SLS_CONTENT_SPTFQMR(S)->Psolve;
  nli          = &(SLS_CONTENT_SPTFQMR(S)->numiters);
  nps          = &(SLS_CONTENT_SPTFQMR(S)->numpsolves);
  res_norm     = &(SLS_CONTENT_SPTFQMR(S)->resnorm);

  /* Initialize counters and convergence flag */
  temp_val = r_curr_norm = -ONE;
  *nli = *nps = 0;
  converged = FALSE;
  b_ok = FALSE;

  /* set booleantype flags for internal solver options */
  preOnLeft  = ( (SLS_CONTENT_SPTFQMR(S)->pretype == PREC_LEFT) || 
                 (SLS_CONTENT_SPTFQMR(S)->pretype == PREC_BOTH) );
  preOnRight = ( (SLS_CONTENT_SPTFQMR(S)->pretype == PREC_RIGHT) || 
                 (SLS_CONTENT_SPTFQMR(S)->pretype == PREC_BOTH) );
  scale_x = (sx != NULL);
  scale_b = (sb != NULL);

  /* Set r_star to initial (unscaled) residual r_star = r_0 = b - A*x_0 */
  /* NOTE: if x == 0 then just set residual to b and continue */
  if (N_VDotProd(x, x) == ZERO) N_VScale(ONE, b, r_star);
  else {
    ier = atimes(A_data, x, r_star);
    if (ier != 0) {
      SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
        SPTFQMR_ATIMES_FAIL_UNREC : SPTFQMR_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPTFQMR(S)->last_flag);
    }
    N_VLinearSum(ONE, b, -ONE, r_star, r_star);
  }

  /* Apply left preconditioner and b-scaling to r_star (or really just r_0) */
  if (preOnLeft) {
    ier = psolve(P_data, r_star, vtemp1, delta, PREC_LEFT);
    (*nps)++;
    if (ier != 0) {
      SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
        SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
      return(SLS_CONTENT_SPTFQMR(S)->last_flag);
    }
  }
  else N_VScale(ONE, r_star, vtemp1);
  if (scale_b) N_VProd(sb, vtemp1, r_star);
  else N_VScale(ONE, vtemp1, r_star);

  /* Initialize rho[0] */
  /* NOTE: initialized here to reduce number of computations - avoid need
           to compute r_star^T*r_star twice, and avoid needlessly squaring
           values */
  rho[0] = N_VDotProd(r_star, r_star);

  /* Compute norm of initial residual (r_0) to see if we really need
     to do anything */
  *res_norm = r_init_norm = SUNRsqrt(rho[0]);
  if (r_init_norm <= delta) {
    SLS_CONTENT_SPTFQMR(S)->last_flag = SPTFQMR_SUCCESS;
    return(SLS_CONTENT_SPTFQMR(S)->last_flag);
  }

  /* Set v = A*r_0 (preconditioned and scaled) */
  if (scale_x) N_VDiv(r_star, sx, vtemp1);
  else N_VScale(ONE, r_star, vtemp1);
  if (preOnRight) {
    N_VScale(ONE, vtemp1, v);
    ier = psolve(P_data, v, vtemp1, delta, PREC_RIGHT);
    (*nps)++;
    if (ier != 0) {
      SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
        SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
      return(SLS_CONTENT_SPTFQMR(S)->last_flag);
    }
  }
  ier = atimes(A_data, vtemp1, v);
  if (ier != 0) {
    SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
      SPTFQMR_ATIMES_FAIL_UNREC : SPTFQMR_ATIMES_FAIL_REC;
    return(SLS_CONTENT_SPTFQMR(S)->last_flag);
  }
  if (preOnLeft) {
    ier = psolve(P_data, v, vtemp1, delta, PREC_LEFT);
    (*nps)++;
    if (ier != 0) {
      SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
        SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
      return(SLS_CONTENT_SPTFQMR(S)->last_flag);
    }
  }
  else N_VScale(ONE, v, vtemp1);
  if (scale_b) N_VProd(sb, vtemp1, v);
  else N_VScale(ONE, vtemp1, v);

  /* Initialize remaining variables */
  N_VScale(ONE, r_star, r[0]);
  N_VScale(ONE, r_star, u);
  N_VScale(ONE, r_star, p);
  N_VConst(ZERO, d);

  tau = r_init_norm;
  v_bar = eta = ZERO;

  /* START outer loop */
  for (n = 0; n < l_max; ++n) {

    /* Increment linear iteration counter */
    (*nli)++;

    /* sigma = r_star^T*v */
    sigma = N_VDotProd(r_star, v);

    /* alpha = rho[0]/sigma */
    alpha = rho[0]/sigma;

    /* q = u-alpha*v */
    N_VLinearSum(ONE, u, -alpha, v, q);

    /* r[1] = r[0]-alpha*A*(u+q) */
    N_VLinearSum(ONE, u, ONE, q, r[1]);
    if (scale_x) N_VDiv(r[1], sx, r[1]);
    if (preOnRight) {
      N_VScale(ONE, r[1], vtemp1);
      ier = psolve(P_data, vtemp1, r[1], delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
          SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPTFQMR(S)->last_flag);
      }
    }
    ier = atimes(A_data, r[1], vtemp1);
    if (ier != 0) {
      SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
        SPTFQMR_ATIMES_FAIL_UNREC : SPTFQMR_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPTFQMR(S)->last_flag);
    }
    if (preOnLeft) {
      ier = psolve(P_data, vtemp1, r[1], delta, PREC_LEFT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
          SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPTFQMR(S)->last_flag);
      }
    }
    else N_VScale(ONE, vtemp1, r[1]);
    if (scale_b) N_VProd(sb, r[1], vtemp1);
    else N_VScale(ONE, r[1], vtemp1);
    N_VLinearSum(ONE, r[0], -alpha, vtemp1, r[1]);

    /* START inner loop */
    for (m = 0; m < 2; ++m) {

      /* d = [*]+(v_bar^2*eta/alpha)*d */
      /* NOTES:
       *   (1) [*] = u if m == 0, and q if m == 1
       *   (2) using temp_val reduces the number of required computations
       *       if the inner loop is executed twice
       */
      if (m == 0) {
	temp_val = SUNRsqrt(N_VDotProd(r[1], r[1]));
	omega = SUNRsqrt(SUNRsqrt(N_VDotProd(r[0], r[0]))*temp_val);
	N_VLinearSum(ONE, u, SUNSQR(v_bar)*eta/alpha, d, d);
      }
      else {
	omega = temp_val;
	N_VLinearSum(ONE, q, SUNSQR(v_bar)*eta/alpha, d, d);
      }

      /* v_bar = omega/tau */
      v_bar = omega/tau;

      /* c = (1+v_bar^2)^(-1/2) */
      c = ONE / SUNRsqrt(ONE+SUNSQR(v_bar));

      /* tau = tau*v_bar*c */
      tau = tau*v_bar*c;

      /* eta = c^2*alpha */
      eta = SUNSQR(c)*alpha;

      /* x = x+eta*d */
      N_VLinearSum(ONE, x, eta, d, x);

      /* Check for convergence... */
      /* NOTE: just use approximation to norm of residual, if possible */
      *res_norm = r_curr_norm = tau*SUNRsqrt(m+1);

      /* Exit inner loop if iteration has converged based upon approximation
	 to norm of current residual */
      if (r_curr_norm <= delta) {
	converged = TRUE;
	break;
      }

      /* Decide if actual norm of residual vector should be computed */
      /* NOTES:
       *   (1) if r_curr_norm > delta, then check if actual residual norm
       *       is OK (recall we first compute an approximation)
       *   (2) if r_curr_norm >= r_init_norm and m == 1 and n == l_max, then
       *       compute actual residual norm to see if the iteration can be
       *       saved
       *   (3) the scaled and preconditioned right-hand side of the given
       *       linear system (denoted by b) is only computed once, and the
       *       result is stored in vtemp3 so it can be reused - reduces the
       *       number of psovles if using left preconditioning
       */
      if ((r_curr_norm > delta) ||
	  (r_curr_norm >= r_init_norm && m == 1 && n == l_max)) {

	/* Compute norm of residual ||b-A*x||_2 (preconditioned and scaled) */
	if (scale_x) N_VDiv(x, sx, vtemp1);
	else N_VScale(ONE, x, vtemp1);
	if (preOnRight) {
	  ier = psolve(P_data, vtemp1, vtemp2, delta, PREC_RIGHT);
	  (*nps)++;
	  if (ier != 0) {
            SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
              SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_UNREC;
            return(SLS_CONTENT_SPTFQMR(S)->last_flag);
          }
	  N_VScale(ONE, vtemp2, vtemp1);
	}
	ier = atimes(A_data, vtemp1, vtemp2);
        if (ier != 0) {
          SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
            SPTFQMR_ATIMES_FAIL_UNREC : SPTFQMR_ATIMES_FAIL_REC;
          return(SLS_CONTENT_SPTFQMR(S)->last_flag);
        }
	if (preOnLeft) {
	  ier = psolve(P_data, vtemp2, vtemp1, delta, PREC_LEFT);
	  (*nps)++;
	  if (ier != 0) {
            SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
              SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
            return(SLS_CONTENT_SPTFQMR(S)->last_flag);
          }
	}
	else N_VScale(ONE, vtemp2, vtemp1);
	if (scale_b) N_VProd(sb, vtemp1, vtemp2);
	else N_VScale(ONE, vtemp1, vtemp2);
	/* Only precondition and scale b once (result saved for reuse) */
	if (!b_ok) {
	  b_ok = TRUE;
	  if (preOnLeft) {
	    ier = psolve(P_data, b, vtemp3, delta, PREC_LEFT);
	    (*nps)++;
	    if (ier != 0) {
              SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
                SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
              return(SLS_CONTENT_SPTFQMR(S)->last_flag);
            }
	  }
	  else N_VScale(ONE, b, vtemp3);
	  if (scale_b) N_VProd(sb, vtemp3, vtemp3);
	}
	N_VLinearSum(ONE, vtemp3, -ONE, vtemp2, vtemp1);
	*res_norm = r_curr_norm = SUNRsqrt(N_VDotProd(vtemp1, vtemp1));

	/* Exit inner loop if inequality condition is satisfied 
	   (meaning exit if we have converged) */
	if (r_curr_norm <= delta) {
	  converged = TRUE;
	  break;
	}

      }

    }  /* END inner loop */

    /* If converged, then exit outer loop as well */
    if (converged == TRUE) break;

    /* rho[1] = r_star^T*r_[1] */
    rho[1] = N_VDotProd(r_star, r[1]);

    /* beta = rho[1]/rho[0] */
    beta = rho[1]/rho[0];

    /* u = r[1]+beta*q */
    N_VLinearSum(ONE, r[1], beta, q, u);

    /* p = u+beta*(q+beta*p) */
    N_VLinearSum(beta, q, SUNSQR(beta), p, p);
    N_VLinearSum(ONE, u, ONE, p, p);

    /* v = A*p */
    if (scale_x) N_VDiv(p, sx, vtemp1);
    else N_VScale(ONE, p, vtemp1);
    if (preOnRight) {
      N_VScale(ONE, vtemp1, v);
      ier = psolve(P_data, v, vtemp1, delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
          SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPTFQMR(S)->last_flag);
      }
    }
    ier = atimes(A_data, vtemp1, v);
    if (ier != 0) {
      SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
        SPTFQMR_ATIMES_FAIL_UNREC : SPTFQMR_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPTFQMR(S)->last_flag);
    }
    if (preOnLeft) {
      ier = psolve(P_data, v, vtemp1, delta, PREC_LEFT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
          SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPTFQMR(S)->last_flag);
      }
    }
    else N_VScale(ONE, v, vtemp1);
    if (scale_b) N_VProd(sb, vtemp1, v);
    else N_VScale(ONE, vtemp1, v);

    /* Shift variable values */
    /* NOTE: reduces storage requirements */
    N_VScale(ONE, r[1], r[0]);
    rho[0] = rho[1];

  }  /* END outer loop */

  /* Determine return value */
  /* If iteration converged or residual was reduced, then return current iterate (x) */
  if ((converged == TRUE) || (r_curr_norm < r_init_norm)) {
    if (scale_x) N_VDiv(x, sx, x);
    if (preOnRight) {
      ier = psolve(P_data, x, vtemp1, delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPTFQMR(S)->last_flag = (ier < 0) ?
          SPTFQMR_PSOLVE_FAIL_UNREC : SPTFQMR_PSOLVE_FAIL_UNREC;
        return(SLS_CONTENT_SPTFQMR(S)->last_flag);
      }
      N_VScale(ONE, vtemp1, x);
    }
    if (converged == TRUE) 
      SLS_CONTENT_SPTFQMR(S)->last_flag = SPTFQMR_SUCCESS;
    else 
      SLS_CONTENT_SPTFQMR(S)->last_flag = SPTFQMR_RES_REDUCED;
    return(SLS_CONTENT_SPTFQMR(S)->last_flag);
  }
  /* Otherwise, return error code */
  else {
    SLS_CONTENT_SPTFQMR(S)->last_flag = SPTFQMR_CONV_FAIL;
    return(SLS_CONTENT_SPTFQMR(S)->last_flag);
  }
}


int SUNLinSolNumIters_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'numiters' value */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  return (SLS_CONTENT_SPTFQMR(S)->numiters);
}


realtype SUNLinSolResNorm_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'resnorm' value */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  return (SLS_CONTENT_SPTFQMR(S)->resnorm);
}


int SUNLinSolNumPSolves_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'numpsolves' value */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  return (SLS_CONTENT_SPTFQMR(S)->numpsolves);
}


long int SUNLinSolLastFlag_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(SPTFQMR_MEM_NULL);
  return (SLS_CONTENT_SPTFQMR(S)->last_flag);
}


int SUNLinSolFree_SPTFQMR(SUNLinearSolver S)
{
  if (S == NULL) return(SPTFQMR_MEM_NULL);

  /* delete items from within the content structure */
  if (SLS_CONTENT_SPTFQMR(S)->r_star)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->r_star);
  if (SLS_CONTENT_SPTFQMR(S)->q)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->q);
  if (SLS_CONTENT_SPTFQMR(S)->d)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->d);
  if (SLS_CONTENT_SPTFQMR(S)->v)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->v);
  if (SLS_CONTENT_SPTFQMR(S)->p)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->p);
  if (SLS_CONTENT_SPTFQMR(S)->r)
    N_VDestroyVectorArray(SLS_CONTENT_SPTFQMR(S)->r, 2);
  if (SLS_CONTENT_SPTFQMR(S)->u)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->u);
  if (SLS_CONTENT_SPTFQMR(S)->vtemp1)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->vtemp1);
  if (SLS_CONTENT_SPTFQMR(S)->vtemp2)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->vtemp2);
  if (SLS_CONTENT_SPTFQMR(S)->vtemp3)
    N_VDestroy(SLS_CONTENT_SPTFQMR(S)->vtemp3);

  /* delete generic structures */
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}
