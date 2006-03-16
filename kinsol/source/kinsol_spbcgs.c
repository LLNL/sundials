/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2006-03-16 20:34:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
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
#include "kinsol_spbcgs.h"
#include "kinsol_spils_impl.h"

#include "sundials_spbcgs.h"
#include "sundials_math.h"

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)

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

/*
 * -----------------------------------------------------------------
 * readability replacements
 * -----------------------------------------------------------------
 */

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

#define pretype   (kinspils_mem->s_pretype)
#define nli       (kinspils_mem->s_nli)
#define npe       (kinspils_mem->s_npe)
#define nps       (kinspils_mem->s_nps)
#define ncfl      (kinspils_mem->s_ncfl)
#define njtimes   (kinspils_mem->s_njtimes)
#define nfes      (kinspils_mem->s_nfes)
#define new_uu    (kinspils_mem->s_new_uu)
#define spils_mem (kinspils_mem->s_spils_mem)
#define last_flag (kinspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcg
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPBCG linear solver module.
 * KINSpbcg sets the kin_linit, kin_lsetup, kin_lsolve, and
 * kin_lfree fields in *kinmem to be KINSpbcgInit, KINSpbcgSetup,
 * KINSpbcgSolve, and KINSpbcgFree, respectively. It allocates
 * memory for a structure of type KINSpilsMemRec and sets the
 * kin_lmem field in *kinmem to the address of this structure. It
 * also calls SpbcgMalloc to allocate memory for the module
 * SPBCG. In summary, KINSpbcg sets the following fields in the
 * KINSpilsMemRec structure:
 *
 *  pretype     = PREC_NONE
 *  s_maxl      = KINSPILS_MAXL  if maxl <= 0
 *              = maxl           if maxl >  0
 *  s_pset      = NULL
 *  s_psolve    = NULL
 *  s_P_data    = NULL
 *  s_jtimes    = NULL
 *  s_J_data    = NULL
 *  s_last_flag = KINSPILS_SUCCESS
 * -----------------------------------------------------------------
 */

int KINSpbcg(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;
  int maxl1;

  if (kinmem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", "KINSpbcg", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale, N_VDiv, 
     or N_VWL2Norm because they are required by KINSOL */

  if ((vec_tmpl->ops->nvconst == NULL) ||
      (vec_tmpl->ops->nvdotprod == NULL) ||
      (vec_tmpl->ops->nvl1norm == NULL)) {
    KINProcessError(NULL, KINSPILS_ILL_INPUT, "KINSPILS", "KINSpbcg", MSGS_BAD_NVECTOR);
    return(KINSPILS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* set four main function fields in kin_mem */

  linit  = KINSpbcgInit; 
  lsetup = KINSpbcgSetup;
  lsolve = KINSpbcgSolve;
  lfree  = KINSpbcgFree;

  /* get memory for KINSpilsMemRec */
  kinspils_mem = NULL;
  kinspils_mem = (KINSpilsMem) malloc(sizeof(KINSpilsMemRec));
  if (kinspils_mem == NULL){
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpbcg", MSGS_MEM_FAIL);
    return(KINSPILS_MEM_FAIL);  
  }

  /* Set ILS type */
  kinspils_mem->s_type = SPILS_SPBCG;

  /* set SPBCG parameters that were passed in call sequence */

  maxl1 = (maxl <= 0) ? KINSPILS_MAXL : maxl;
  kinspils_mem->s_maxl = maxl1;  

  /* set default values for the rest of the SPBCG parameters */

  kinspils_mem->s_pretype   = PREC_NONE;
  kinspils_mem->s_last_flag = KINSPILS_SUCCESS;
  kinspils_mem->s_pset      = NULL;
  kinspils_mem->s_psolve    = NULL;
  kinspils_mem->s_P_data    = NULL;
  kinspils_mem->s_jtimes    = NULL;
  kinspils_mem->s_J_data    = NULL;

  /* call SpbcgMalloc to allocate workspace for SPBCG */

  /* vec_tmpl passed as template vector */
  spbcg_mem = NULL;
  spbcg_mem = SpbcgMalloc(maxl1, vec_tmpl);
  if (spbcg_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_FAIL, "KINSPILS", "KINSpbcg", MSGS_MEM_FAIL);
    free(kinspils_mem); kinspils_mem = NULL;
    return(KINSPILS_MEM_FAIL);
  }

  /* This is an iterative linear solver */

  inexact_ls = TRUE;

  /* Attach SPBCG memory to spils memory structure */
  spils_mem = (void *) spbcg_mem;

  /* attach linear solver memory to KINSOL memory */
  lmem = kinspils_mem;

  return(KINSPILS_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

#define maxl   (kinspils_mem->s_maxl)
#define pset   (kinspils_mem->s_pset)
#define psolve (kinspils_mem->s_psolve)
#define P_data (kinspils_mem->s_P_data)
#define jtimes (kinspils_mem->s_jtimes)
#define J_data (kinspils_mem->s_J_data)

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
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;

  kinspils_mem = (KINSpilsMem) lmem;
  spbcg_mem = (SpbcgMem) spils_mem;

  /* initialize counters */

  npe = nli = nps = ncfl = 0;
  njtimes = nfes = 0;

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
    jtimes = KINSpilsDQJtimes;
    J_data = kin_mem;
  }

  /*  Set maxl in the SPBCG memory in case it was changed by the user */
  spbcg_mem->l_max  = maxl;

  last_flag = KINSPILS_SUCCESS;
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
  KINSpilsMem kinspils_mem;
  int ret;

  kinspils_mem = (KINSpilsMem) lmem;

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
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;
  int ret, nli_inc, nps_inc;
  
  kinspils_mem = (KINSpilsMem) lmem;
  spbcg_mem = (SpbcgMem) spils_mem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling KINSpbcgSolve, to the RHS vector for the system
     to be solved. */ 
 
  N_VConst(ZERO, xx);

  new_uu = TRUE;  /* set flag required for user Jacobian routine */

  /* call SpbcgSolve */

  ret = SpbcgSolve(spbcg_mem, kin_mem, xx, bb, pretype, eps,
                   kin_mem, fscale, fscale, KINSpilsAtimes,
                   KINSpilsPSolve, res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */

  nli = nli + (long int) nli_inc;
  nps = nps + (long int) nps_inc;

  if (printfl > 2) 
    KINPrintInfo(kin_mem, PRNT_NLI, "KINSPBCG", "KINSpbcgSolve", INFO_NLI, nli_inc);

  if (ret != 0) ncfl++;

  /* Compute the terms sJpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.

     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */

  KINSpilsAtimes(kin_mem, xx, bb);
  sJpnorm = N_VWL2Norm(bb,fscale);
  N_VProd(bb, fscale, bb);
  N_VProd(bb, fscale, bb);
  sfdotJp = N_VDotProd(fval, bb);

  if (printfl > 2) 
    KINPrintInfo(kin_mem, PRNT_EPS, "KINSPBCG", "KINSpbcgSolve", INFO_EPS, *res_norm, eps);

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
  KINSpilsMem kinspils_mem;
  SpbcgMem spbcg_mem;

  kinspils_mem = (KINSpilsMem) lmem;
  spbcg_mem = (SpbcgMem) spils_mem;

  SpbcgFree(spbcg_mem);
  free(lmem); lmem = NULL;

  return(0);
}
