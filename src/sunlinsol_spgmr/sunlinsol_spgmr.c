/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * Based on sundials_spgmr.c code, written by Scott D. Cohen, 
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the implementation file for the SPGMR implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPGMR linear solver
 */

SUNLinearSolver SUNSPGMR(N_Vector y, int pretype, int maxl)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_SPGMR content;
  
  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH))
    pretype = PREC_NONE;
  if (maxl <= 0)
    maxl = SUNSPGMR_MAXL_DEFAULT;

  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_SPGMR;
  ops->setatimes         = SUNLinSolSetATimes_SPGMR;
  ops->setpreconditioner = SUNLinSolSetPreconditioner_SPGMR;
  ops->setscalingvectors = SUNLinSolSetScalingVectors_SPGMR;
  ops->initialize        = SUNLinSolInitialize_SPGMR;
  ops->setup             = SUNLinSolSetup_SPGMR;
  ops->solve             = SUNLinSolSolve_SPGMR;
  ops->numiters          = SUNLinSolNumIters_SPGMR;
  ops->resnorm           = SUNLinSolResNorm_SPGMR;
  ops->numpsolves        = SUNLinSolNumPSolves_SPGMR;
  ops->lastflag          = SUNLinSolLastFlag_SPGMR;  
  ops->free              = SUNLinSolFree_SPGMR;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPGMR) malloc(sizeof(struct _SUNLinearSolverContent_SPGMR));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->last_flag = 0;
  content->maxl = maxl;
  content->pretype = pretype;
  content->gstype = SUNSPGMR_GSTYPE_DEFAULT;
  content->max_restarts = 0;
  content->numiters = 0;
  content->numpsolves = 0;
  content->resnorm = ZERO;
  content->xcor = N_VClone(y);
  content->vtemp = N_VClone(y);
  content->s1 = NULL;
  content->s2 = NULL;
  content->ATSetup = NULL;
  content->ATimes = NULL;
  content->ATData = NULL;
  content->Psetup = NULL;
  content->Psolve = NULL;
  content->PData = NULL;
  content->V = NULL;
  content->Hes = NULL;
  content->givens = NULL;
  content->yg = NULL;

  /* Attach content and ops */
  S->content = content;
  S->ops     = ops;

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the type of preconditioning for SPGMR to use 
 */

SUNDIALS_EXPORT int SUNSPGMRSetPrecType(SUNLinearSolver S, int pretype) 
{
  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    return(SPGMR_ILL_INPUT);
  }

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPGMR_MEM_NULL);

  /* Set pretype */
  SLS_CONTENT_SPGMR(S)->pretype = pretype;
  return(SPGMR_SUCCESS);
}


/* ----------------------------------------------------------------------------
 * Function to set the type of Gram-Schmidt orthogonalization for SPGMR to use
 */

SUNDIALS_EXPORT int SUNSPGMRSetGSType(SUNLinearSolver S, int gstype)
{
  /* Check for legal gstype */ 
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    return(SPGMR_ILL_INPUT);
  }

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPGMR_MEM_NULL);

  /* Set pretype */
  SLS_CONTENT_SPGMR(S)->gstype = gstype;
  return(SPGMR_SUCCESS);
}


/* ----------------------------------------------------------------------------
 * Function to return the SPGMR workspace sizes
 */

SUNDIALS_EXPORT int SUNSPGMRGetWorkspace(SUNLinearSolver S, 
                                         long int *lenrwLS, 
                                         long int *leniwLS)
{
  int maxl;
  long int liw1, lrw1;

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SPGMR_MEM_NULL);

  maxl = SLS_CONTENT_SPGMR(S)->maxl;
  N_VSpace(SLS_CONTENT_SPGMR(S)->vtemp, &lrw1, &liw1);
  *lenrwLS = lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;
  *leniwLS = liw1*(maxl + 5);
  return(SPGMR_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPGMR(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_ITERATIVE;
}


int SUNLinSolInitialize_SPGMR(SUNLinearSolver S)
{
  int k;

  /* set shortcut to SPGMR memory structure */
  if (S == NULL) return(SPGMR_MEM_NULL);  
  SUNLinearSolverContent_SPGMR content = SLS_CONTENT_SPGMR(S);

  /* ensure valid options */
  if (content->max_restarts < 0) 
    content->max_restarts = 0;
  if ( (content->pretype != PREC_LEFT) && 
       (content->pretype != PREC_RIGHT) && 
       (content->pretype != PREC_BOTH) )
    content->pretype = PREC_NONE;


  /* allocate solver-specific memory (where the size depends on the
     choice of maxl) here */

  /*   Krylov subspace vectors */
  content->V = N_VCloneVectorArray(content->maxl+1, content->vtemp);
  if (content->V == NULL) {
    SUNLinSolFree(S);
    content->last_flag = 1;
    return(1);
  }

  /*   Hessenberg matrix Hes */
  content->Hes = (realtype **) malloc((content->maxl+1)*sizeof(realtype *)); 
  if (content->Hes == NULL) {
    SUNLinSolFree(S);
    content->last_flag = 1;
    return(1);
  }

  for (k=0; k<=content->maxl; k++) {
    content->Hes[k] = NULL;
    content->Hes[k] = (realtype *) malloc(content->maxl*sizeof(realtype));
    if (content->Hes[k] == NULL) {
      SUNLinSolFree(S);
      content->last_flag = 1;
      return(1);
    }
  }
  
  /*   Givens rotation components */
  content->givens = NULL;
  content->givens = (realtype *) malloc(2*content->maxl*sizeof(realtype));
  if (content->givens == NULL) {
    SUNLinSolFree(S);
    content->last_flag = 1;
    return(1);
  }

  /*    y and g vectors */
  content->yg = (realtype *) malloc((content->maxl+1)*sizeof(realtype));
  if (content->yg == NULL) {
    SUNLinSolFree(S);
    content->last_flag = 1;
    return(1);
  }

  /* return with success */
  content->last_flag = 0;
  return 0;
}


int SUNLinSolSetATimes_SPGMR(SUNLinearSolver S, void* ATData, 
                            ATSetupFn ATSetup, ATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATSetup 
     and ATimes routines and data, and return with success */
  if (S == NULL) return(SPGMR_MEM_NULL);
  SLS_CONTENT_SPGMR(S)->ATSetup = ATSetup;
  SLS_CONTENT_SPGMR(S)->ATimes = ATimes;
  SLS_CONTENT_SPGMR(S)->ATData = ATData;
  SLS_CONTENT_SPGMR(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetPreconditioner_SPGMR(SUNLinearSolver S, void* PData,
                                    PSetupFn Psetup, PSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  if (S == NULL) return(SPGMR_MEM_NULL);
  SLS_CONTENT_SPGMR(S)->Psetup = Psetup;
  SLS_CONTENT_SPGMR(S)->Psolve = Psolve;
  SLS_CONTENT_SPGMR(S)->PData = PData;
  SLS_CONTENT_SPGMR(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetScalingVectors_SPGMR(SUNLinearSolver S, N_Vector s1,
                                     N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors, 
     and return with success */
  if (S == NULL) return(SPGMR_MEM_NULL);
  SLS_CONTENT_SPGMR(S)->s1 = s1;
  SLS_CONTENT_SPGMR(S)->s2 = s2;
  SLS_CONTENT_SPGMR(S)->last_flag = 0;
  return 0;
}


int SUNLinSolSetup_SPGMR(SUNLinearSolver S, SUNMatrix A)
{
  /* Set shortcuts to SPGMR memory structures */
  if (S == NULL) return(SPGMR_MEM_NULL);
  ATSetupFn ATSetup = SLS_CONTENT_SPGMR(S)->ATSetup;
  PSetupFn Psetup = SLS_CONTENT_SPGMR(S)->Psetup;
  void* ATData = SLS_CONTENT_SPGMR(S)->ATData;
  void* PData = SLS_CONTENT_SPGMR(S)->PData;
  
  /* no solver-specific setup is required, but if user-supplied 
     ATSetup or Psetup routines exist, call those here */
  if (ATSetup != NULL) {
    SLS_CONTENT_SPGMR(S)->last_flag = ATSetup(ATData);
    if (SLS_CONTENT_SPGMR(S)->last_flag != 0)  return 1;
  }
  if (Psetup != NULL) {
    SLS_CONTENT_SPGMR(S)->last_flag = Psetup(PData);
    if (SLS_CONTENT_SPGMR(S)->last_flag != 0)  return 1;
  }
  
  /* return with success */ 
  return 0;
}


int SUNLinSolSolve_SPGMR(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                         N_Vector b, realtype delta)
{
  /* local data and shortcut variables */
  N_Vector *V, xcor, vtemp, s1, s2;
  realtype **Hes, *givens, *yg, *res_norm;
  realtype beta, rotation_product, r_norm, s_product, rho;
  booleantype preOnLeft, preOnRight, scale2, scale1, converged;
  int i, j, k, l, l_plus_1, l_max, krydim, ier, ntries, max_restarts, gstype;
  int *nli, *nps;
  void *A_data, *P_data;
  ATimesFn atimes;
  PSolveFn psolve;

  /* Initialize some variables */
  l_plus_1 = 0;
  krydim = 0;

  /* Make local shorcuts to solver variables. */
  if (S == NULL) return(SPGMR_MEM_NULL);
  l_max        = SLS_CONTENT_SPGMR(S)->maxl;
  max_restarts = SLS_CONTENT_SPGMR(S)->max_restarts;
  gstype       = SLS_CONTENT_SPGMR(S)->gstype;
  V            = SLS_CONTENT_SPGMR(S)->V;
  Hes          = SLS_CONTENT_SPGMR(S)->Hes;
  givens       = SLS_CONTENT_SPGMR(S)->givens;
  xcor         = SLS_CONTENT_SPGMR(S)->xcor;
  yg           = SLS_CONTENT_SPGMR(S)->yg;
  vtemp        = SLS_CONTENT_SPGMR(S)->vtemp;
  s1           = SLS_CONTENT_SPGMR(S)->s1;
  s2           = SLS_CONTENT_SPGMR(S)->s2;
  A_data       = SLS_CONTENT_SPGMR(S)->ATData;
  P_data       = SLS_CONTENT_SPGMR(S)->PData;
  atimes       = SLS_CONTENT_SPGMR(S)->ATimes;
  psolve       = SLS_CONTENT_SPGMR(S)->Psolve;
  nli          = &(SLS_CONTENT_SPGMR(S)->numiters);
  nps          = &(SLS_CONTENT_SPGMR(S)->numpsolves);
  res_norm     = &(SLS_CONTENT_SPGMR(S)->resnorm);

  /* Initialize counters and convergence flag */
  *nli = *nps = 0;
  converged = FALSE;

  /* set booleantype flags for internal solver options */
  preOnLeft  = ( (SLS_CONTENT_SPGMR(S)->pretype == PREC_LEFT) || 
                 (SLS_CONTENT_SPGMR(S)->pretype == PREC_BOTH) );
  preOnRight = ( (SLS_CONTENT_SPGMR(S)->pretype == PREC_RIGHT) || 
                 (SLS_CONTENT_SPGMR(S)->pretype == PREC_BOTH) );
  scale1 = (s1 != NULL);
  scale2 = (s2 != NULL);

  /* Set vtemp and V[0] to initial (unscaled) residual r_0 = b - A*x_0 */
  if (N_VDotProd(x, x) == ZERO) {
    N_VScale(ONE, b, vtemp);
  } else {
    ier = atimes(A_data, x, vtemp);
    if (ier != 0) {
      SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
          SPGMR_ATIMES_FAIL_UNREC : SPGMR_ATIMES_FAIL_REC;
      return(SLS_CONTENT_SPGMR(S)->last_flag);
    }
    N_VLinearSum(ONE, b, -ONE, vtemp, vtemp);
  }
  N_VScale(ONE, vtemp, V[0]);

  /* Apply left preconditioner and left scaling to V[0] = r_0 */
  if (preOnLeft) {
    ier = psolve(P_data, V[0], vtemp, delta, PREC_LEFT);
    (*nps)++;
    if (ier != 0) {
      SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
          SPGMR_PSOLVE_FAIL_UNREC : SPGMR_PSOLVE_FAIL_REC;
      return(SLS_CONTENT_SPGMR(S)->last_flag);
    }
  } else {
    N_VScale(ONE, V[0], vtemp);
  }
  
  if (scale1) {
    N_VProd(s1, vtemp, V[0]);   
  } else {
    N_VScale(ONE, vtemp, V[0]);
  }

  /* Set r_norm = beta to L2 norm of V[0] = s1 P1_inv r_0, and
     return if small  */
  *res_norm = r_norm = beta = SUNRsqrt(N_VDotProd(V[0], V[0])); 
  if (r_norm <= delta) {
    SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_SUCCESS;
    return(SLS_CONTENT_SPGMR(S)->last_flag);
  }

  /* Initialize rho to avoid compiler warning message */
  rho = beta;

  /* Set xcor = 0 */
  N_VConst(ZERO, xcor);

  /* Begin outer iterations: up to (max_restarts + 1) attempts */  
  for (ntries=0; ntries<=max_restarts; ntries++) {
    
    /* Initialize the Hessenberg matrix Hes and Givens rotation
       product.  Normalize the initial vector V[0] */
    for (i=0; i<=l_max; i++)
      for (j=0; j<l_max; j++)
        Hes[i][j] = ZERO;

    rotation_product = ONE;    
    N_VScale(ONE/r_norm, V[0], V[0]);
    
    /* Inner loop: generate Krylov sequence and Arnoldi basis */    
    for (l=0; l<l_max; l++) {
      (*nli)++;      
      krydim = l_plus_1 = l + 1;
      
      /* Generate A-tilde V[l], where A-tilde = s1 P1_inv A P2_inv s2_inv */
      
      /*   Apply right scaling: vtemp = s2_inv V[l] */
      if (scale2) N_VDiv(V[l], s2, vtemp);
      else N_VScale(ONE, V[l], vtemp);
      
      /*   Apply right preconditioner: vtemp = P2_inv s2_inv V[l] */ 
      if (preOnRight) {
        N_VScale(ONE, vtemp, V[l_plus_1]);
        ier = psolve(P_data, V[l_plus_1], vtemp, delta, PREC_RIGHT);
        (*nps)++;
        if (ier != 0) {
          SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
              SPGMR_PSOLVE_FAIL_UNREC : SPGMR_PSOLVE_FAIL_REC;
          return(SLS_CONTENT_SPGMR(S)->last_flag);
        }
      }
      
      /* Apply A: V[l+1] = A P2_inv s2_inv V[l] */
      ier = atimes( A_data, vtemp, V[l_plus_1] );
      if (ier != 0) {
        SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
            SPGMR_ATIMES_FAIL_UNREC : SPGMR_ATIMES_FAIL_REC;
        return(SLS_CONTENT_SPGMR(S)->last_flag);
      }
      
      /* Apply left preconditioning: vtemp = P1_inv A P2_inv s2_inv V[l] */
      if (preOnLeft) {
        ier = psolve(P_data, V[l_plus_1], vtemp, delta, PREC_LEFT);
        (*nps)++;
        if (ier != 0) {
          SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
              SPGMR_PSOLVE_FAIL_UNREC : SPGMR_PSOLVE_FAIL_REC;
          return(SLS_CONTENT_SPGMR(S)->last_flag);
        }
      } else {
        N_VScale(ONE, V[l_plus_1], vtemp);
      }
      
      /* Apply left scaling: V[l+1] = s1 P1_inv A P2_inv s2_inv V[l] */
      if (scale1) {
        N_VProd(s1, vtemp, V[l_plus_1]);
      } else {
        N_VScale(ONE, vtemp, V[l_plus_1]);
      }
      
      /*  Orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde */
      if (gstype == CLASSICAL_GS) {
        if (ClassicalGS(V, Hes, l_plus_1, l_max, &(Hes[l_plus_1][l]),
                        vtemp, yg) != 0) {
          SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_GS_FAIL;
          return(SLS_CONTENT_SPGMR(S)->last_flag);
        }
      } else {
        if (ModifiedGS(V, Hes, l_plus_1, l_max, &(Hes[l_plus_1][l])) != 0) {
          SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_GS_FAIL;
          return(SLS_CONTENT_SPGMR(S)->last_flag);
        }
      }
      
      /*  Update the QR factorization of Hes */
      if(QRfact(krydim, Hes, givens, l) != 0 ) {
        SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_QRFACT_FAIL;
        return(SLS_CONTENT_SPGMR(S)->last_flag);
      }
      
      /*  Update residual norm estimate; break if convergence test passes */
      rotation_product *= givens[2*l+1];
      *res_norm = rho = SUNRabs(rotation_product*r_norm);
      
      if (rho <= delta) { converged = TRUE; break; }
      
      /* Normalize V[l+1] with norm value from the Gram-Schmidt routine */
      N_VScale(ONE/Hes[l_plus_1][l], V[l_plus_1], V[l_plus_1]);
    }
    
    /* Inner loop is done.  Compute the new correction vector xcor */
    
    /*   Construct g, then solve for y */
    yg[0] = r_norm;
    for (i=1; i<=krydim; i++) yg[i]=ZERO;
    if (QRsol(krydim, Hes, givens, yg) != 0) {
      SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_QRSOL_FAIL;
      return(SLS_CONTENT_SPGMR(S)->last_flag);
    }
    
    /*   Add correction vector V_l y to xcor */
    for (k=0; k<krydim; k++)
      N_VLinearSum(yg[k], V[k], ONE, xcor, xcor);
    
    /* If converged, construct the final solution vector x and return */
    if (converged) {
      
      /* Apply right scaling and right precond.: vtemp = P2_inv s2_inv xcor */
      if (scale2) N_VDiv(xcor, s2, xcor);
      if (preOnRight) {
        ier = psolve(P_data, xcor, vtemp, delta, PREC_RIGHT);
        (*nps)++;
        if (ier != 0) {
          SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
              SPGMR_PSOLVE_FAIL_UNREC : SPGMR_PSOLVE_FAIL_REC;
          return(SLS_CONTENT_SPGMR(S)->last_flag);
        }
      } else {
        N_VScale(ONE, xcor, vtemp);
      }
      
      /* Add vtemp to initial x to get final solution x, and return */
      N_VLinearSum(ONE, x, ONE, vtemp, x);
      
      SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_SUCCESS;
      return(SLS_CONTENT_SPGMR(S)->last_flag);
    }
    
    /* Not yet converged; if allowed, prepare for restart */    
    if (ntries == max_restarts) break;
    
    /* Construct last column of Q in yg */
    s_product = ONE;
    for (i=krydim; i>0; i--) {
      yg[i] = s_product*givens[2*i-2];
      s_product *= givens[2*i-1];
    }
    yg[0] = s_product;
    
    /* Scale r_norm and yg */
    r_norm *= s_product;
    for (i=0; i<=krydim; i++)
      yg[i] *= r_norm;
    r_norm = SUNRabs(r_norm);
    
    /* Multiply yg by V_(krydim+1) to get last residual vector; restart */
    N_VScale(yg[0], V[0], V[0]);
    for (k=1; k<=krydim; k++)
      N_VLinearSum(yg[k], V[k], ONE, V[0], V[0]);
    
  }
  
  /* Failed to converge, even after allowed restarts.
     If the residual norm was reduced below its initial value, compute
     and return x anyway.  Otherwise return failure flag. */
  if (rho < beta) {
    
    /* Apply right scaling and right precond.: vtemp = P2_inv s2_inv xcor */    
    if (scale2) N_VDiv(xcor, s2, xcor);
    if (preOnRight) {
      ier = psolve(P_data, xcor, vtemp, delta, PREC_RIGHT);
      (*nps)++;
      if (ier != 0) {
        SLS_CONTENT_SPGMR(S)->last_flag = (ier < 0) ? 
            SPGMR_PSOLVE_FAIL_UNREC : SPGMR_PSOLVE_FAIL_REC;
        return(SLS_CONTENT_SPGMR(S)->last_flag);
      }
      } else {
      N_VScale(ONE, xcor, vtemp);
    }

    /* Add vtemp to initial x to get final solution x, and return */
    N_VLinearSum(ONE, x, ONE, vtemp, x);
    
    SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_RES_REDUCED;
    return(SLS_CONTENT_SPGMR(S)->last_flag);
  }

  SLS_CONTENT_SPGMR(S)->last_flag = SPGMR_CONV_FAIL;
  return(SLS_CONTENT_SPGMR(S)->last_flag);
}


int SUNLinSolNumIters_SPGMR(SUNLinearSolver S)
{
  /* return the stored 'numiters' value */
  if (S == NULL) return(SPGMR_MEM_NULL);
  return (SLS_CONTENT_SPGMR(S)->numiters);
}


realtype SUNLinSolResNorm_SPGMR(SUNLinearSolver S)
{
  /* return the stored 'resnorm' value */
  if (S == NULL) return(SPGMR_MEM_NULL);
  return (SLS_CONTENT_SPGMR(S)->resnorm);
}


int SUNLinSolNumPSolves_SPGMR(SUNLinearSolver S)
{
  /* return the stored 'numpsolves' value */
  if (S == NULL) return(SPGMR_MEM_NULL);
  return (SLS_CONTENT_SPGMR(S)->numpsolves);
}


long int SUNLinSolLastFlag_SPGMR(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(SPGMR_MEM_NULL);
  return (SLS_CONTENT_SPGMR(S)->last_flag);
}


int SUNLinSolFree_SPGMR(SUNLinearSolver S)
{
  int k;

  if (S == NULL) return(SPGMR_MEM_NULL);

  /* delete items from within the content structure */
  if (SLS_CONTENT_SPGMR(S)->xcor)
    N_VDestroy(SLS_CONTENT_SPGMR(S)->xcor);
  if (SLS_CONTENT_SPGMR(S)->vtemp)
    N_VDestroy(SLS_CONTENT_SPGMR(S)->vtemp);
  if (SLS_CONTENT_SPGMR(S)->V)
    N_VDestroyVectorArray(SLS_CONTENT_SPGMR(S)->V, 
                          SLS_CONTENT_SPGMR(S)->maxl+1);
  if (SLS_CONTENT_SPGMR(S)->Hes) {
    for (k=0; k<=SLS_CONTENT_SPGMR(S)->maxl; k++)
      if (SLS_CONTENT_SPGMR(S)->Hes[k])
        free(SLS_CONTENT_SPGMR(S)->Hes[k]);
    free(SLS_CONTENT_SPGMR(S)->Hes);
  }
  if (SLS_CONTENT_SPGMR(S)->givens)
    free(SLS_CONTENT_SPGMR(S)->givens);
  if (SLS_CONTENT_SPGMR(S)->yg)
    free(SLS_CONTENT_SPGMR(S)->yg);

  /* delete generic structures */
  free(S->content);  S->content = NULL;
  free(S->ops);  S->ops = NULL;
  free(S); S = NULL;
  return 0;
}
