/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-31 18:30:46 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the KINSOL interface to the
 * scaled, preconditioned Bi-CGSTAB (SPBCG) iterative linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include "kinsol_spbcgs_impl.h"
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
 * keys for KINSpbcgPrintInfo
 * -----------------------------------------------------------------
 */

#define PRNT_NLI 1
#define PRNT_EPS 2

/*
 * -----------------------------------------------------------------
 * function prototypes
 * -----------------------------------------------------------------
 */

/* KINSpbcg linit, lsetup, lsolve, and lfree routines */

static int KINSpbcgInit(KINMem kin_mem);
static int KINSpbcgSetup(KINMem kin_mem);
static int KINSpbcgSolve(KINMem kin_mem, N_Vector xx,
			 N_Vector bb, realtype *res_norm);
static int KINSpbcgFree(KINMem kin_mem);

/* KINSpbcg Atimes and PSolve routines called by generic SPBCG solver */

static int KINSpbcgAtimes(void *kinsol_mem, N_Vector v, N_Vector z);
static int KINSpbcgPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lr);

/* difference quotient approximation for Jacobian-vector product */

static int KINSpbcgDQJtimes(N_Vector v, N_Vector Jv, N_Vector u,
			    booleantype *new_u, void *jac_data);

static void KINSpbcgPrintInfo(KINMem kin_mem, char *funcname, int key,...);

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
#define inexact_ls   (kin_mem->kin_inexact_ls)
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

#define pretype   (kinspbcg_mem->b_pretype)
#define nli       (kinspbcg_mem->b_nli)
#define npe       (kinspbcg_mem->b_npe)
#define nps       (kinspbcg_mem->b_nps)
#define ncfl      (kinspbcg_mem->b_ncfl)
#define njtimes   (kinspbcg_mem->b_njtimes)
#define nfeSG     (kinspbcg_mem->b_nfeSG)
#define new_uu    (kinspbcg_mem->b_new_uu)
#define spbcg_mem (kinspbcg_mem->b_spbcg_mem)
#define last_flag (kinspbcg_mem->b_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcg
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPBCG linear solver module.
 * KINSpbcg sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in *kinmem to be KINSpbcgInit, KINSpbcgSetup,
 * KINSpbcgSolve, and KINSpbcgFree, respectively. It allocates
 * memory for a structure of type KINSpbcgMemRec and sets the
 * kin_lmem field in *kinmem to the address of this structure. It
 * also calls SpbcgMalloc to allocate memory for the module
 * SPBCG. In summary, KINSpbcg sets the following fields in the
 * KINSpbcgMemRec structure:
 *
 *  pretype     = PREC_NONE
 *  b_maxl      = KINSPBCG_MAXL  if maxl <= 0
 *              = maxl           if maxl >  0
 *  b_pset      = NULL
 *  b_psolve    = NULL
 *  b_P_data    = NULL
 *  b_jtimes    = NULL
 *  b_J_data    = NULL
 *  b_last_flag = KINSPBCG_SUCCESS
 * -----------------------------------------------------------------
 */

int KINSpbcg(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;
  int maxl1;

  if (kinmem == NULL){
    fprintf(stderr, MSGB_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale, N_VDiv, 
     or N_VWL2Norm because they are required by KINSOL */

  if ((vec_tmpl->ops->nvconst == NULL) ||
      (vec_tmpl->ops->nvdotprod == NULL) ||
      (vec_tmpl->ops->nvl1norm == NULL)) {
    if (errfp != NULL) fprintf(errfp, MSGB_BAD_NVECTOR);
    return(KINSPBCG_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* set four main function fields in kin_mem */

  linit  = KINSpbcgInit; 
  lsetup = KINSpbcgSetup;
  lsolve = KINSpbcgSolve;
  lfree  = KINSpbcgFree;

  /* get memory for KINSpbcgMemRec */
  kinspbcg_mem = NULL;
  kinspbcg_mem = (KINSpbcgMem) malloc(sizeof(KINSpbcgMemRec));
  if (kinspbcg_mem == NULL){
    fprintf(errfp, MSGB_MEM_FAIL);
    return(KINSPBCG_MEM_FAIL);  
  }

  /* set SPBCG parameters that were passed in call sequence */

  maxl1 = (maxl <= 0) ? KINSPBCG_MAXL : maxl;
  kinspbcg_mem->b_maxl = maxl1;  

  /* set default values for the rest of the SPBCG parameters */

  kinspbcg_mem->b_pretype   = PREC_NONE;
  kinspbcg_mem->b_last_flag = KINSPBCG_SUCCESS;
  kinspbcg_mem->b_pset      = NULL;
  kinspbcg_mem->b_psolve    = NULL;
  kinspbcg_mem->b_P_data    = NULL;
  kinspbcg_mem->b_jtimes    = NULL;
  kinspbcg_mem->b_J_data    = NULL;

  /* call SpbcgMalloc to allocate workspace for SPBCG */

  /* vec_tmpl passed as template vector */
  spbcg_mem = NULL;
  spbcg_mem = SpbcgMalloc(maxl1, vec_tmpl);
  if (spbcg_mem == NULL) {
    fprintf(errfp, MSGB_MEM_FAIL);
    free(kinspbcg_mem); kinspbcg_mem = NULL;
    return(KINSPBCG_MEM_FAIL);
  }

  /* This is an iterative linear solver */

  inexact_ls = TRUE;

  /* attach linear solver memory to KINSOL memory */

  lmem = kinspbcg_mem;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgSetPreconditioner
 * -----------------------------------------------------------------
 */

int KINSpbcgSetPreconditioner(void *kinmem,
			      KINSpilsPrecSetupFn pset,
			      KINSpilsPrecSolveFn psolve,
			      void *P_data)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;

  kinspbcg_mem->b_pset   = pset;
  kinspbcg_mem->b_psolve = psolve;
  kinspbcg_mem->b_P_data = P_data;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int KINSpbcgSetJacTimesVecFn(void *kinmem,
			     KINSpilsJacTimesVecFn jtimes,
			     void *J_data)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;

  kinspbcg_mem->b_jtimes = jtimes;
  kinspbcg_mem->b_J_data = J_data;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINSpbcgGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG)
{
  KINMem kin_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }

  *lenrwSG = lrw1 * 7;
  *leniwSG = liw1 * 7;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int KINSpbcgGetNumPrecEvals(void *kinmem, long int *npevals)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;
  *npevals = npe;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int KINSpbcgGetNumPrecSolves(void *kinmem, long int *npsolves)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;
  *npsolves = nps;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetNumLinIters
 * -----------------------------------------------------------------
 */

int KINSpbcgGetNumLinIters(void *kinmem, long int *nliters)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;
  *nliters = nli;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetNumConvFails
 * -----------------------------------------------------------------
 */

int KINSpbcgGetNumConvFails(void *kinmem, long int *nlcfails)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;
  *nlcfails = ncfl;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int KINSpbcgGetNumJtimesEvals(void *kinmem, long int *njvevals)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;
  *njvevals = njtimes;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINSpbcgGetNumFuncEvals(void *kinmem, long int *nfevalsSG)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;
  *nfevalsSG = nfeSG;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgGetLastFlag
 * -----------------------------------------------------------------
 */

int KINSpbcgGetLastFlag(void *kinmem, int *flag)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* return immediately if kinmem is NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSGB_SETGET_KINMEM_NULL);
    return(KINSPBCG_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(stderr, MSGB_SETGET_LMEM_NULL);
    return(KINSPBCG_LMEM_NULL);
  }
  kinspbcg_mem = (KINSpbcgMem) lmem;

  *flag = last_flag;

  return(KINSPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

#define maxl   (kinspbcg_mem->b_maxl)
#define pset   (kinspbcg_mem->b_pset)
#define psolve (kinspbcg_mem->b_psolve)
#define P_data (kinspbcg_mem->b_P_data)
#define jtimes (kinspbcg_mem->b_jtimes)
#define J_data (kinspbcg_mem->b_J_data)

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgInit
 * -----------------------------------------------------------------
 * This routine initializes variables associated with the SPBCG
 * iterative linear solver. Mmemory allocation was done previously
 * in KINSpbcg.
 * -----------------------------------------------------------------
 */

static int KINSpbcgInit(KINMem kin_mem)
{
  KINSpbcgMem kinspbcg_mem;

  kinspbcg_mem = (KINSpbcgMem) lmem;

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
    jtimes = KINSpbcgDQJtimes;
    J_data = kin_mem;
  }

  last_flag = KINSPBCG_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the SPBCG linear
 * solver, that is, it is an interface to the user-supplied
 * routine pset.
 * -----------------------------------------------------------------
 */

static int KINSpbcgSetup(KINMem kin_mem)
{
  KINSpbcgMem kinspbcg_mem;
  int ret;

  kinspbcg_mem = (KINSpbcgMem) lmem;

  /* call pset routine */

  ret = pset(uu, uscale, fval, fscale, P_data, vtemp1, vtemp2);

  last_flag = ret;

  if (ret != 0) return(1);

  npe++;
  nnilset = nni; 

  /* return the same value ret that pset returned */

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic SPBCG solver routine
 * called SpbcgSolve for the solution of the linear system Ax = b.
 *
 * Appropriate variables are passed to SpbcgSolve and the counters
 * nli, nps, and ncfl are incremented, and the return value is set
 * according to the success of SpbcgSolve. The success flag is
 * returned if SpbcgSolve converged, or if the residual was reduced.
 * Of the other error conditions, only preconditioner solver
 * failure is specifically returned. Otherwise a generic flag is
 * returned to denote failure of this routine.
 * -----------------------------------------------------------------
 */

static int KINSpbcgSolve(KINMem kin_mem, N_Vector xx, N_Vector bb, 
                         realtype *res_norm)
{
  KINSpbcgMem kinspbcg_mem;
  int ret, nli_inc, nps_inc;
  
  kinspbcg_mem = (KINSpbcgMem) lmem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling KINSpbcgSolve, to the RHS vector for the system
     to be solved. */ 
 
  N_VConst(ZERO, xx);

  new_uu = TRUE;  /* set flag required for user Jacobian routine */

  /* call SpbcgSolve */

  ret = SpbcgSolve(spbcg_mem, kin_mem, xx, bb, pretype, eps,
                   kin_mem, fscale, fscale, KINSpbcgAtimes,
                   KINSpbcgPSolve, res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */

  nli = nli + (long int) nli_inc;
  nps = nps + (long int) nps_inc;

  if (printfl > 2) KINSpbcgPrintInfo(kin_mem, "KINSpbcgSolve", PRNT_NLI, &nli_inc);

  if (ret != 0) ncfl++;

  /* Compute the terms sJpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.

     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */

  KINSpbcgAtimes(kin_mem, xx, bb);
  sJpnorm = N_VWL2Norm(bb,fscale);
  N_VProd(bb, fscale, bb);
  N_VProd(bb, fscale, bb);
  sfdotJp = N_VDotProd(fval, bb);

  if (printfl > 2) KINSpbcgPrintInfo(kin_mem, "KINSpbcgSolve", PRNT_EPS, res_norm, &eps);

  /* set return value to appropriate value */

  last_flag = ret;

  if ((ret == SPBCG_SUCCESS) || (ret == SPBCG_RES_REDUCED)) return(0);
  else if (ret == SPBCG_PSOLVE_FAIL_REC) return(1);
  else return(-1);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgFree
 * -----------------------------------------------------------------
 * Frees memory specific to the SPBCG linear solver module.
 * -----------------------------------------------------------------
 */

static int KINSpbcgFree(KINMem kin_mem)
{
  KINSpbcgMem kinspbcg_mem;

  kinspbcg_mem = (KINSpbcgMem) lmem;

  SpbcgFree(spbcg_mem);
  free(lmem); lmem = NULL;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgAtimes
 * -----------------------------------------------------------------
 * This routine coordinates the generation of the matrix-vector
 * product z = J*v by calling either KINSpbcgDQJtimes, which uses a
 * difference quotient approximation for J*v, or by calling the
 * user-supplied routine KINSpbcgJacTimesVecFn if it is non-NULL.
 * -----------------------------------------------------------------
 */

static int KINSpbcgAtimes(void *kinsol_mem, N_Vector v, N_Vector z)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;
  int ret;

  kin_mem = (KINMem) kinsol_mem;
  kinspbcg_mem = (KINSpbcgMem) lmem;

  ret = jtimes(v, z, uu, &new_uu, J_data);
  njtimes++;

  return(ret);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SpbcgSolve routine
 * and the user's psolve routine. It passes to psolve all required
 * state information from kinsol_mem. Its return value is the same
 * as that returned by psolve. Note that the generic SPBCG solver
 * guarantees that KINSpbcgPSolve will not be called in the case in
 * which preconditioning is not done. This is the only case in which
 * the user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

static int KINSpbcgPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lrdummy)
{
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;
  int ret;

  kin_mem = (KINMem) kinsol_mem;
  kinspbcg_mem = (KINSpbcgMem) lmem;

  /* copy the rhs into z before the psolve call */   
  /* Note: z returns with the solution */

  N_VScale(ONE, r, z);

  /* this call is counted in nps within the KINSpbcgSolve routine */

  ret = psolve(uu, uscale, fval, fscale, z, P_data, vtemp1);

  return(ret);     
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcgDQJtimes
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

static int KINSpbcgDQJtimes(N_Vector v, N_Vector Jv,
                            N_Vector u, booleantype *new_u, 
                            void *jac_data)
{
  realtype sigma, sigma_inv, sutsv, sq1norm, sign, vtv;
  KINMem kin_mem;
  KINSpbcgMem kinspbcg_mem;

  /* jac_data is kin_mem */

  kin_mem = (KINMem) jac_data;
  kinspbcg_mem = (KINSpbcgMem) lmem;

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
 * KINSpbcgPrintInfo
 * -----------------------------------------------------------------
 */

static void KINSpbcgPrintInfo(KINMem kin_mem, char *funcname, int key,...)
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
