/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-31 18:30:46 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the KINSOL interface to the
 * scaled, preconditioned TFQMR (SPTFQMR) iterative linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include "kinsol_sptfqmr_impl.h"
#include "sundials_math.h"

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/*
 * -----------------------------------------------------------------
 * keys for KINSptfqmrPrintInfo
 * -----------------------------------------------------------------
 */

#define PRNT_NLI 1
#define PRNT_EPS 2

/*
 * -----------------------------------------------------------------
 * function prototypes
 * -----------------------------------------------------------------
 */

/* KINSptfqmr linit, lsetup, lsolve, and lfree routines */

static int KINSptfqmrInit(KINMem kin_mem);
static int KINSptfqmrSetup(KINMem kin_mem);
static int KINSptfqmrSolve(KINMem kin_mem, N_Vector xx,
			   N_Vector bb, realtype *res_norm);
static int KINSptfqmrFree(KINMem kin_mem);

/* KINSptfqmr Atimes and PSolve routines called by generic SPTFQMR solver */

static int KINSptfqmrAtimes(void *kinsol_mem, N_Vector v, N_Vector z);
static int KINSptfqmrPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lr);

/* difference quotient approximation for Jacobian-vector product */

static int KINSptfqmrDQJtimes(N_Vector v, N_Vector Jv, N_Vector u,
			      booleantype *new_u, void *jac_data);

static void KINSptfqmrPrintInfo(KINMem kin_mem, char *funcname, int key,...);

/*
 * -----------------------------------------------------------------
 * readability replacements
 * -----------------------------------------------------------------
 */

#define lrw1         (kin_mem->kin_lrw1)
#define liw1         (kin_mem->kin_liw1)
#define nni          (kin_mem->kin_nni)
#define nnilset      (kin_mem->kin_nnilset)
#define func         (kin_mem->kin_func)
#define f_data       (kin_mem->kin_f_data)
#define printfl      (kin_mem->kin_printfl)
#define linit        (kin_mem->kin_linit)
#define lsetup       (kin_mem->kin_lsetup)
#define lsolve       (kin_mem->kin_lsolve)
#define lfree        (kin_mem->kin_lfree)
#define lmem         (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu           (kin_mem->kin_uu)
#define fval         (kin_mem->kin_fval)
#define uscale       (kin_mem->kin_uscale)
#define fscale       (kin_mem->kin_fscale)
#define sqrt_relfunc (kin_mem->kin_sqrt_relfunc)
#define eps          (kin_mem->kin_eps)
#define sJpnorm      (kin_mem->kin_sJpnorm)
#define sfdotJp      (kin_mem->kin_sfdotJp)
#define errfp        (kin_mem->kin_errfp)
#define infofp       (kin_mem->kin_infofp)
#define setupNonNull (kin_mem->kin_setupNonNull)
#define vtemp1       (kin_mem->kin_vtemp1)
#define vec_tmpl     (kin_mem->kin_vtemp1)
#define vtemp2       (kin_mem->kin_vtemp2)

#define pretype     (kinsptfqmr_mem->q_pretype)
#define nli         (kinsptfqmr_mem->q_nli)
#define npe         (kinsptfqmr_mem->q_npe)
#define nps         (kinsptfqmr_mem->q_nps)
#define ncfl        (kinsptfqmr_mem->q_ncfl)
#define njtimes     (kinsptfqmr_mem->q_njtimes)
#define nfeSG       (kinsptfqmr_mem->q_nfeSG)
#define new_uu      (kinsptfqmr_mem->q_new_uu)
#define sptfqmr_mem (kinsptfqmr_mem->q_sptfqmr_mem)
#define last_flag   (kinsptfqmr_mem->q_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmr
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPTFQMR linear solver module.
 * KINSptfqmr sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in *kinmem to be KINSptfqmrInit, KINSptfqmrSetup,
 * KINSptfqmrSolve, and KINSptfqmrFree, respectively. It allocates
 * memory for a structure of type KINSptfqmrMemRec and sets the
 * kin_lmem field in *kinmem to the address of this structure. It
 * also calls SptfqmrMalloc to allocate memory for the module
 * SPTFQMR. In summary, KINSptfqmr sets the following fields in the
 * KINSptfqmrMemRec structure:
 *
 *  q_pretype   = PREC_NONE
 *  q_maxl      = KINSPTFQMR_MAXL  if maxl <= 0
 *              = maxl             if maxl >  0
 *  q_pset      = NULL
 *  q_psolve    = NULL
 *  q_P_data    = NULL
 *  q_jtimes    = NULL
 *  q_J_data    = NULL
 *  q_last_flag = KINSPTFQMR_SUCCESS
 * -----------------------------------------------------------------
 */

int KINSptfqmr(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;
  int maxl1;

  if (kinmem == NULL){
    fprintf(stderr, MSGQ_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale, N_VDiv, 
     or N_VWL2Norm because they are required by KINSOL */

  if ((vec_tmpl->ops->nvconst == NULL) ||
      (vec_tmpl->ops->nvdotprod == NULL) ||
      (vec_tmpl->ops->nvl1norm == NULL)) {
    if (errfp != NULL) fprintf(errfp, MSGQ_BAD_NVECTOR);
    return(KINSPTFQMR_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* set four main function fields in kin_mem */

  linit  = KINSptfqmrInit; 
  lsetup = KINSptfqmrSetup;
  lsolve = KINSptfqmrSolve;
  lfree  = KINSptfqmrFree;

  /* get memory for KINSptfqmrMemRec */
  kinsptfqmr_mem = NULL;
  kinsptfqmr_mem = (KINSptfqmrMem) malloc(sizeof(KINSptfqmrMemRec));
  if (kinsptfqmr_mem == NULL){
    fprintf(errfp, MSGQ_MEM_FAIL);
    return(KINSPTFQMR_MEM_FAIL);  
  }

  /* set SPTFQMR parameters that were passed in call sequence */

  maxl1 = (maxl <= 0) ? KINSPTFQMR_MAXL : maxl;
  kinsptfqmr_mem->q_maxl = maxl1;  

  /* set default values for the rest of the SPTFQMR parameters */

  kinsptfqmr_mem->q_pretype   = PREC_NONE;
  kinsptfqmr_mem->q_last_flag = KINSPTFQMR_SUCCESS;
  kinsptfqmr_mem->q_pset      = NULL;
  kinsptfqmr_mem->q_psolve    = NULL;
  kinsptfqmr_mem->q_P_data    = NULL;
  kinsptfqmr_mem->q_jtimes    = NULL;
  kinsptfqmr_mem->q_J_data    = NULL;

  /* call SptfqmrMalloc to allocate workspace for SPTFQMR */

  /* vec_tmpl passed as template vector */

  sptfqmr_mem = NULL;
  sptfqmr_mem = SptfqmrMalloc(maxl1, vec_tmpl);
  if (sptfqmr_mem == NULL) {
    fprintf(errfp, MSGQ_MEM_FAIL);
    free(kinsptfqmr_mem); kinsptfqmr_mem = NULL;
    return(KINSPTFQMR_MEM_FAIL);
  }

  /* this is an iterative linear solver */

  inexact_ls = TRUE;

  /* attach linear solver memory to KINSOL memory */

  lmem = kinsptfqmr_mem;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrSetPrecSetupFn
 * -----------------------------------------------------------------
 */

int KINSptfqmrSetPreconditioner(void *kinmem,
				KINSpilsPrecSetupFn pset,
				KINSpilsPrecSolveFn psolve,
                                void *P_data)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  kinsptfqmr_mem->q_pset   = pset;
  kinsptfqmr_mem->q_psolve = psolve;
  kinsptfqmr_mem->q_P_data = P_data;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int KINSptfqmrSetJacTimesVecFn(void *kinmem,
			       KINSpilsJacTimesVecFn jtimes,
			       void *J_data)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  kinsptfqmr_mem->q_jtimes = jtimes;
  kinsptfqmr_mem->q_J_data = J_data;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;
  int maxl;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  maxl = kinsptfqmr_mem->q_maxl;

#ifdef DEBUG
  *lenrwSG = lrw1 * 12;
  *leniwSG = liw1 * 12;
#else
  *lenrwSG = lrw1 * 11;
  *leniwSG = liw1 * 11;
#endif

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetNumPrecEvals(void *kinmem, long int *npevals)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;
  *npevals = npe;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetNumPrecSolves(void *kinmem, long int *npsolves)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;
  *npsolves = nps;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetNumLinIters
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetNumLinIters(void *kinmem, long int *nliters)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;
  *nliters = nli;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetNumConvFails
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetNumConvFails(void *kinmem, long int *nlcfails)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;
  *nlcfails = ncfl;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetNumJtimesEvals(void *kinmem, long int *njvevals)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;
  *njvevals = njtimes;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetNumFuncEvals(void *kinmem, long int *nfevalsSG)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;
  *nfevalsSG = nfeSG;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrGetLastFlag
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetLastFlag(void *kinmem, int *flag)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_KINMEM_NULL);
    return(KINSPTFQMR_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(stderr, MSGQ_SETGET_LMEM_NULL);
    return(KINSPTFQMR_LMEM_NULL);
  }
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  *flag = last_flag;

  return(KINSPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

#define maxl   (kinsptfqmr_mem->q_maxl)
#define pset   (kinsptfqmr_mem->q_pset)
#define psolve (kinsptfqmr_mem->q_psolve)
#define P_data (kinsptfqmr_mem->q_P_data)
#define jtimes (kinsptfqmr_mem->q_jtimes)
#define J_data (kinsptfqmr_mem->q_J_data)

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrInit
 * -----------------------------------------------------------------
 * This routine initializes variables associated with the SPTFQMR
 * iterative linear solver. Memory allocation was done previously
 * in KINSptfqmr.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrInit(KINMem kin_mem)
{
  KINSptfqmrMem kinsptfqmr_mem;

  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  /* initialize counters */

  npe = nli = nps = ncfl = 0;
  njtimes = nfeSG = 0;

  /* set preconditioner type */

  if (psolve != NULL) {
    pretype = PREC_RIGHT;
  } else {
    pretype = PREC_NONE;
  }
  
  /* set setupNonNull to TRUE iff there is preconditioning with setup */

  setupNonNull = ((psolve != NULL) && (pset != NULL));

  /* if jtimes is NULL at this time, set it to private DQ routine */

  if (jtimes == NULL) {
    jtimes = KINSptfqmrDQJtimes;
    J_data = kin_mem;
  }

  last_flag = KINSPTFQMR_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the SPTFQMR linear
 * solver, that is, it is an interface to the user-supplied
 * routine pset.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrSetup(KINMem kin_mem)
{
  KINSptfqmrMem kinsptfqmr_mem;
  int ret;

  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  /* call pset routine */

  ret = pset(uu, uscale, fval, fscale, P_data, vtemp1, vtemp2);

  last_flag = ret;

  if (ret != 0) return(1);

  npe++;
  nnilset = nni;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic SPTFQMR solver routine
 * called SptfqmrSolve for the solution of the linear system Ax = b.
 *
 * Appropriate variables are passed to SptfqmrSolve and the counters
 * nli, nps, and ncfl are incremented, and the return value is set
 * according to the success of SptfqmrSolve. The success flag is
 * returned if SptfqmrSolve converged, or if the residual was reduced.
 * Of the other error conditions, only preconditioner solver
 * failure is specifically returned. Otherwise a generic flag is
 * returned to denote failure of this routine.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrSolve(KINMem kin_mem, N_Vector xx, N_Vector bb, 
			   realtype *res_norm)
{
  KINSptfqmrMem kinsptfqmr_mem;
  int ret, nli_inc, nps_inc;
  
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling KINSptfqmrSolve, to the RHS vector for the system
     to be solved. */ 
 
  N_VConst(ZERO, xx);

  new_uu = TRUE;  /* set flag required for user Jacobian routine */

  /* call SptfqmrSolve */

  ret = SptfqmrSolve(sptfqmr_mem, kin_mem, xx, bb, pretype, eps,
		     kin_mem, fscale, fscale, KINSptfqmrAtimes,
		     KINSptfqmrPSolve, res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */

  nli = nli + (long int) nli_inc;
  nps = nps + (long int) nps_inc;

  if (printfl > 2) KINSptfqmrPrintInfo(kin_mem, "KINSptfqmrSolve", PRNT_NLI, &nli_inc);

  if (ret != 0) ncfl++;

  /* Compute the terms sJpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.

     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */

  KINSptfqmrAtimes(kin_mem, xx, bb);
  sJpnorm = N_VWL2Norm(bb,fscale);
  N_VProd(bb, fscale, bb);
  N_VProd(bb, fscale, bb);
  sfdotJp = N_VDotProd(fval, bb);

  if (printfl > 2) KINSptfqmrPrintInfo(kin_mem, "KINSptfqmrSolve", PRNT_EPS, res_norm, &eps);

  /* set return value to appropriate value */

  last_flag = ret;

  if ((ret == SPTFQMR_SUCCESS) || (ret == SPTFQMR_RES_REDUCED)) return(0);
  else if (ret == SPTFQMR_PSOLVE_FAIL_REC) return(1);
  else return(-1);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrFree
 * -----------------------------------------------------------------
 * Frees memory specific to the SPTFQMR linear solver module.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrFree(KINMem kin_mem)
{
  KINSptfqmrMem kinsptfqmr_mem;

  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  SptfqmrFree(sptfqmr_mem);
  free(lmem); lmem = NULL;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrAtimes
 * -----------------------------------------------------------------
 * This routine coordinates the generation of the matrix-vector
 * product z = J*v by calling either KINSptfqmrDQJtimes, which uses a
 * difference quotient approximation for J*v, or by calling the
 * user-supplied routine KINSptfqmrJacTimesVecFn if it is non-NULL.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrAtimes(void *kinsol_mem, N_Vector v, N_Vector z)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;
  int ret;

  kin_mem = (KINMem) kinsol_mem;
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  ret = jtimes(v, z, uu, &new_uu, J_data);
  njtimes++;

  return(ret);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SptfqmrSolve routine
 * and the user's psolve routine. It passes to psolve all required
 * state information from kinsol_mem. Its return value is the same
 * as that returned by psolve. Note that the generic SPTFQMR solver
 * guarantees that KINSptfqmrPSolve will not be called in the case in
 * which preconditioning is not done. This is the only case in which
 * the user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lrdummy)
{
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;
  int ret;

  kin_mem = (KINMem) kinsol_mem;
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  /* copy the rhs into z before the psolve call */   
  /* Note: z returns with the solution */

  N_VScale(ONE, r, z);

  /* this call is counted in nps within the KINSptfqmrSolve routine */

  ret = psolve(uu, uscale, fval, fscale, z, P_data, vtemp1);

  return(ret);     
}

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmrDQJtimes
 * -----------------------------------------------------------------
 * This routine computes the matrix-vector product z = J*v using a
 * difference quotient approximation. The approximation is 
 * J*v = [func(uu + sigma*v) - func(uu)]/sigma. Here sigma is based
 * on the dot products (uscale*uu, uscale*v)
 * and (uscale*v, uscale*v), ||uscale*v||_L1, and on sqrt_relfunc
 * (the square root of the relative error in the function). Note
 * that v in the argument list has already been both preconditioned
 * and unscaled.
 * -----------------------------------------------------------------
 */

static int KINSptfqmrDQJtimes(N_Vector v, N_Vector Jv,
			      N_Vector u, booleantype *new_u, 
			      void *jac_data)
{
  realtype sigma, sigma_inv, sutsv, sq1norm, sign, vtv;
  KINMem kin_mem;
  KINSptfqmrMem kinsptfqmr_mem;

  /* jac_data is kin_mem */

  kin_mem = (KINMem) jac_data;
  kinsptfqmr_mem = (KINSptfqmrMem) lmem;

  /* scale the vector v and put Du*v into vtemp1 */

  N_VProd(v, uscale, vtemp1);

  /* scale u and put into Jv (used as a temporary storage) */

  N_VProd(u, uscale, Jv);

  /* compute dot product (Du*u).(Du*v) */

  sutsv = N_VDotProd(Jv, vtemp1);

  /* compute dot product (Du*v).(Du*v) */

  vtv = N_VDotProd(vtemp1, vtemp1);

  sq1norm = N_VL1Norm(vtemp1);

  sign = (sutsv >= ZERO) ? ONE : -ONE ;
 
  /*  this expression for sigma is from p. 469, Brown and Saad paper */

  sigma = sign*sqrt_relfunc*MAX(ABS(sutsv),sq1norm)/vtv; 

  sigma_inv = ONE/sigma;

  /* compute the u-prime at which to evaluate the function func */

  N_VLinearSum(ONE, u, sigma, v, vtemp1);
 
  /* call the system function to calculate func(u+sigma*v) */

  func(vtemp1, vtemp2, f_data);    
  nfeSG++;

  /* finish the computation of the difference quotient */

  N_VLinearSum(sigma_inv, vtemp2, -sigma_inv, fval, Jv);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINSptfqmrPrintInfo
 * -----------------------------------------------------------------
 */

static void KINSptfqmrPrintInfo(KINMem kin_mem, char *funcname, int key,...)
{
  va_list ap;
  realtype rnum1, rnum2;
  int inum1;

  fprintf(infofp, "---%s\n   ", funcname);

  /* initialize argument processing */

  va_start(ap, key); 

  switch(key) {

  case PRNT_NLI:
    inum1 = *(va_arg(ap, int *));
    fprintf(infofp, "nli_inc = %d\n", inum1);
    break;
    
  case PRNT_EPS:
    rnum1 = *(va_arg(ap, realtype *));
    rnum2 = *(va_arg(ap, realtype *));
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "residual norm = %12.3Lg  eps = %12.3Lg\n", rnum1, rnum2);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "residual norm = %12.3lg  eps = %12.3lg\n", rnum1, rnum2);
#else
    fprintf(infofp, "residual norm = %12.3g  eps = %12.3g\n", rnum1, rnum2);
#endif
      break;

  }

  /* finalize argument processing */

  va_end(ap);

  return;
}
