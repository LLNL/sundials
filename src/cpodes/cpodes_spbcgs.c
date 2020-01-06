/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the CPSPBCG linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_spbcgs.h>
#include "cpodes_spils_impl.h"
#include "cpodes_private.h"

#include <sundials/sundials_spbcgs.h>
#include <sundials/sundials_math.h>

/* CPSPBCG linit, lsetup, lsolve, and lfree routines */

static int cpSpbcgInit(CPodeMem cp_mem);
static int cpSpbcgSetup(CPodeMem cp_mem, int convfail, 
                        N_Vector yP, N_Vector ypP, N_Vector fctP, 
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpSpbcgSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                        N_Vector yC, N_Vector ypC, N_Vector fctC);
static int cpSpbcgFree(CPodeMem cp_mem);


/* Readability Replacements */

#define ode_type      (cp_mem->cp_ode_type)

#define tq            (cp_mem->cp_tq)
#define nst           (cp_mem->cp_nst)
#define tn            (cp_mem->cp_tn)
#define gamma         (cp_mem->cp_gamma)
#define gammap        (cp_mem->cp_gammap)
#define f             (cp_mem->cp_f)
#define ewt           (cp_mem->cp_ewt)
#define errfp         (cp_mem->cp_errfp)
#define mnewt         (cp_mem->cp_mnewt)
#define linit         (cp_mem->cp_linit)
#define lsetup        (cp_mem->cp_lsetup)
#define lsolve        (cp_mem->cp_lsolve)
#define lfree         (cp_mem->cp_lfree)
#define lmem          (cp_mem->cp_lmem)
#define vec_tmpl      (cp_mem->cp_tempv)
#define lsetup_exists (cp_mem->cp_lsetup_exists)

#define sqrtN     (cpspils_mem->s_sqrtN)   
#define ytemp     (cpspils_mem->s_ytemp)
#define yptemp    (cpspils_mem->s_yptemp)
#define x         (cpspils_mem->s_x)
#define ycur      (cpspils_mem->s_ycur)
#define ypcur     (cpspils_mem->s_ypcur)
#define fcur      (cpspils_mem->s_fcur)
#define delta     (cpspils_mem->s_delta)
#define deltar    (cpspils_mem->s_deltar)

#define npe       (cpspils_mem->s_npe)
#define nli       (cpspils_mem->s_nli)
#define nps       (cpspils_mem->s_nps)
#define ncfl      (cpspils_mem->s_ncfl)
#define nstlpre   (cpspils_mem->s_nstlpre)
#define njtimes   (cpspils_mem->s_njtimes)
#define nfes      (cpspils_mem->s_nfes)

#define spils_mem (cpspils_mem->s_spils_mem)

#define last_flag (cpspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : CPSpbcg
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the Spbcg linear solver module. CPSpbcg first
 * calls the existing lfree routine if this is not NULL. It then sets
 * the cp_linit, cp_lsetup, cp_lsolve, cp_lfree fields in (*cpode_mem)
 * to be CPSpbcgInit, CPSpbcgSetup, CPSpbcgSolve, and CPSpbcgFree,
 * respectively. It allocates memory for a structure of type
 * CPSpilsMemRec and sets the cp_lmem field in (*cpode_mem) to the
 * address of this structure. It sets lsetup_exists in (*cpode_mem),
 * and sets various fields in the CPSpilsMemRec structure.
 * Finally, CPSpbcg allocates memory for ytemp and x, and calls
 * SpbcgMalloc to allocate memory for the Spbcg solver.
 * -----------------------------------------------------------------
 */

int CPSpbcg(void *cpode_mem, int pretype, int maxl)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  SpbcgMem spbcg_mem;
  int mxl;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPSPBCG", "CPSpbcg", MSGS_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPBCG", "CPSpbcg", MSGS_BAD_NVECTOR);
    return(CPSPILS_ILL_INPUT);
  }

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPSPBCG", "CPSpbcg", MSGS_BAD_PRETYPE);
    return(CPSPILS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */
  linit  = cpSpbcgInit;
  lsetup = cpSpbcgSetup;
  lsolve = cpSpbcgSolve;
  lfree  = cpSpbcgFree;

  /* Get memory for CPSpilsMemRec */
  cpspils_mem = NULL;
  cpspils_mem = (CPSpilsMem) malloc(sizeof(struct CPSpilsMemRec));
  if (cpspils_mem == NULL) {
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPSPBCG", "CPSpbcg", MSGS_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Set ILS type */
  cpspils_mem->s_type = SPILS_SPBCG;

  /* Set Spbcg parameters that have been passed in call sequence */
  cpspils_mem->s_pretype = pretype;
  mxl = cpspils_mem->s_maxl = (maxl <= 0) ? CPSPILS_MAXL : maxl;


  /* Set defaults for Jacobian-related fileds */
  cpspils_mem->s_jtimesDQ = SUNTRUE;
  cpspils_mem->s_jtvE     = NULL;
  cpspils_mem->s_jtvI     = NULL;
  cpspils_mem->s_j_data   = NULL;

  /* Set defaults for preconditioner-related fields */
  cpspils_mem->s_psetE  = NULL;
  cpspils_mem->s_psetI  = NULL;
  cpspils_mem->s_pslvE  = NULL;
  cpspils_mem->s_pslvI  = NULL;
  cpspils_mem->s_P_data = NULL;

  /* Set default values for the rest of the Spbcg parameters */
  cpspils_mem->s_delt      = CPSPILS_DELT;

  cpspils_mem->s_last_flag = CPSPILS_SUCCESS;

  lsetup_exists = SUNFALSE;

  /* Alocate memory */
  spbcg_mem = NULL;
  ytemp = NULL;
  yptemp = NULL;
  x = NULL;

  /* Call SpbcgMalloc to allocate workspace for Spbcg */
  spbcg_mem = SpbcgMalloc(mxl, vec_tmpl);
  if (spbcg_mem == NULL) {
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPSPBCG", "CPSpbcg", MSGS_MEM_FAIL);
    free(cpspils_mem);
    return(CPSPILS_MEM_FAIL);
  }

  /* Allocate memory for x, ytemp and (if needed) yptemp */ 
  x = N_VClone(vec_tmpl);
  if (x == NULL) {
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPSPBCG", "CPSpbcg", MSGS_MEM_FAIL);
    SpbcgFree(spbcg_mem);
    free(cpspils_mem);
    return(CPSPILS_MEM_FAIL);
  }
  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPSPBCG", "CPSpbcg", MSGS_MEM_FAIL);
    SpbcgFree(spbcg_mem);
    N_VDestroy(x);
    free(cpspils_mem);
    return(CPSPILS_MEM_FAIL);
  }
  if (ode_type == CP_IMPL) {
    yptemp = N_VClone(vec_tmpl);
    if (yptemp == NULL) {
      cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPSPBCG", "CPSpbcg", MSGS_MEM_FAIL);
      SpbcgFree(spbcg_mem);
      N_VDestroy(x);
      N_VDestroy(ytemp);
      free(cpspils_mem);
      return(CPSPILS_MEM_FAIL);
    }
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = SUNRsqrt(N_VDotProd(ytemp, ytemp));
  
  /* Attach SPBCG memory to spils memory structure */
  spils_mem = (void *) spbcg_mem;

  /* Attach linear solver memory to integrator memory */
  lmem = cpspils_mem;

  return(CPSPILS_SUCCESS);
}

/* Additional readability replacements */

#define pretype  (cpspils_mem->s_pretype)
#define delt     (cpspils_mem->s_delt)

#define maxl     (cpspils_mem->s_maxl)
#define psetE    (cpspils_mem->s_psetE)
#define psetI    (cpspils_mem->s_psetI)
#define pslvE    (cpspils_mem->s_pslvE)
#define pslvI    (cpspils_mem->s_pslvI)
#define P_data   (cpspils_mem->s_P_data)

#define jtimesDQ (cpspils_mem->s_jtimesDQ)
#define jtvE     (cpspils_mem->s_jtvE)
#define jtvI     (cpspils_mem->s_jtvI)
#define j_data   (cpspils_mem->s_j_data)

/*
 * -----------------------------------------------------------------
 * Function : cpSpbcgInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Spbcg
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpSpbcgInit(CPodeMem cp_mem)
{
  CPSpilsMem cpspils_mem;
  SpbcgMem spbcg_mem;

  cpspils_mem = (CPSpilsMem) lmem;
  spbcg_mem = (SpbcgMem) spils_mem;

  /* Initialize counters */
  npe = nli = nps = ncfl = nstlpre = 0;
  njtimes = nfes = 0;

  /* 
   * Check for legal combination pretype - psolve
   *
   * Set lsetup_exists = SUNTRUE iff there is preconditioning (pretype != PREC_NONE)
   * and there is a preconditioning setup phase (pset != NULL)             
   *
   * If jtimes is NULL at this time, set it to DQ 
   */

  if (ode_type == CP_EXPL) {

    if ((pretype != PREC_NONE) && (pslvE == NULL)) {
      cpProcessError(cp_mem, -1, "CPSPBCG", "CPSpbcgInit", MSGS_PSOLVE_REQ);
      last_flag = CPSPILS_ILL_INPUT;
      return(-1);
    }
    lsetup_exists = (pretype != PREC_NONE) && (psetE != NULL);
    if (jtimesDQ) {
      jtvE = cpSpilsDQjtvExpl;
      j_data = cp_mem;
    } else {
      j_data = cp_mem->cp_user_data;
    }

  } else {

    if ((pretype != PREC_NONE) && (pslvI == NULL)) {
      cpProcessError(cp_mem, -1, "CPSPBCG", "CPSpbcgInit", MSGS_PSOLVE_REQ);
      last_flag = CPSPILS_ILL_INPUT;
      return(-1);
    }
    lsetup_exists = (pretype != PREC_NONE) && (psetI != NULL);
    if (jtimesDQ) {
      jtvI = cpSpilsDQjtvImpl;
      j_data = cp_mem;
    } else {
      j_data = cp_mem->cp_user_data;
    }

  }

  /*  Set maxl in the SPBCG memory in case it was changed by the user */
  spbcg_mem->l_max  = maxl;

  last_flag = CPSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : cpSpbcgSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the Spbcg linear solver.
 * It makes a decision as to whether or not to signal for reevaluation
 * of Jacobian data in the pset routine, based on various state
 * variables, then it calls pset. If we signal for reevaluation,
 * then we reset jcur = *jcurPtr to SUNTRUE, regardless of the pset output.
 * In any case, if jcur == SUNTRUE, we increment npe and save nst in nstlpre.
 * -----------------------------------------------------------------
 */
static int cpSpbcgSetup(CPodeMem cp_mem, int convfail, 
                        N_Vector yP, N_Vector ypP, N_Vector fctP, 
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  retval=0;
  CPSpilsMem cpspils_mem;

  cpspils_mem = (CPSpilsMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = SUNRabs((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlpre + CPSPILS_MSBPRE) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPSPILS_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    *jcurPtr = jbad;
    jok = !jbad;
    
    /* Call pset routine and possibly reset jcur */
    retval = psetE(tn, yP, fctP, jok, jcurPtr, gamma, P_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      if (jbad) *jcurPtr = SUNTRUE;
      /* If jcur = SUNTRUE, increment npe and save nst value */
      if (*jcurPtr) {
        npe++;
        nstlpre = nst;
      }
      last_flag = SPBCG_SUCCESS;
    } else if (retval < 0) {
      cpProcessError(cp_mem, SPBCG_PSET_FAIL_UNREC, "CPSPBCG", "CPSpbcgSetup", MSGS_PSET_FAILED);
      last_flag = SPBCG_PSET_FAIL_UNREC;
    } else if (retval > 0) {
      last_flag = SPBCG_PSET_FAIL_REC;
    }
    
    break;

  case CP_IMPL:

    /* Call psetI routine */
    retval = psetI(tn, yP, ypP, fctP, gamma, P_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      last_flag = SPBCG_SUCCESS;
    } else if (retval < 0) {
      cpProcessError(cp_mem, SPBCG_PSET_FAIL_UNREC, "CPSPBCG", "CPSpbcgSetup", MSGS_PSET_FAILED);
      last_flag = SPBCG_PSET_FAIL_UNREC;
    }else if (retval > 0) {
      last_flag = SPBCG_PSET_FAIL_REC;
    }

    npe++;
    nstlpre = nst;

    break;

  }

  /* Return the same value that pset returned */
  return(retval);
}

/*
 * -----------------------------------------------------------------
 * Function : cpSpbcgSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic solver SpbcgSolve
 * for the solution of the linear system Ax = b with the SPBCG method.
 * The solution x is returned in the vector b.
 *
 * If the WRMS norm of b is small, we return x = b (if this is the first
 * Newton iteration) or x = 0 (if a later Newton iteration).
 *
 * Otherwise, we set the tolerance parameter and initial guess (x = 0),
 * call SpbcgSolve, and copy the solution x into b. The x-scaling and
 * b-scaling arrays are both equal to weight.
 *
 * The counters nli, nps, and ncfl are incremented, and the return value
 * is set according to the success of SpbcgSolve. The success flag is
 * returned if SpbcgSolve converged, or if this is the first Newton
 * iteration and the residual norm was reduced below its initial value.
 * -----------------------------------------------------------------
 */
static int cpSpbcgSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                        N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  realtype bnorm, res_norm;
  CPSpilsMem cpspils_mem;
  SpbcgMem spbcg_mem;
  int nli_inc, nps_inc, retval;
  
  cpspils_mem = (CPSpilsMem) lmem;

  spbcg_mem = (SpbcgMem) spils_mem;

  /* Test norm(b); if small, return x = 0 or x = b */
  deltar = delt * tq[4]; 

  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= deltar) {
    if (mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  ycur  = yC;
  ypcur = ypC;
  fcur  = fctC;

  /* Set inputs delta and initial guess x = 0 to SpbcgSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SpbcgSolve and copy x to b */
  retval = SpbcgSolve(spbcg_mem, cp_mem, x, b, pretype, delta,
                      cp_mem, weight, weight, cpSpilsAtimes, cpSpilsPSolve,
                      &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (retval != SPBCG_SUCCESS) ncfl++;

  /* Interpret return value from SpbcgSolve */

  last_flag = retval;

  switch(retval) {

  case SPBCG_SUCCESS:
    return(0);
    break;
  case SPBCG_RES_REDUCED:
    if (mnewt == 0) return(0);
    else            return(1);
    break;
  case SPBCG_CONV_FAIL:
    return(1);
    break;
  case SPBCG_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPBCG_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPBCG_MEM_NULL:
    return(-1);
    break;
  case SPBCG_ATIMES_FAIL_UNREC:
    cpProcessError(cp_mem, SPBCG_ATIMES_FAIL_UNREC, "CPSPBCG", "CPSpbcgSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPBCG_PSOLVE_FAIL_UNREC:
    cpProcessError(cp_mem, SPBCG_PSOLVE_FAIL_UNREC, "CPSPBCG", "CPSpbcgSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : cpSpbcgFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the Spbcg linear solver.
 * -----------------------------------------------------------------
 */

static int cpSpbcgFree(CPodeMem cp_mem)
{
  CPSpilsMem cpspils_mem;
  SpbcgMem spbcg_mem;

  cpspils_mem = (CPSpilsMem) lmem;

  spbcg_mem = (SpbcgMem) spils_mem;

  SpbcgFree(spbcg_mem);
  N_VDestroy(x);
  N_VDestroy(ytemp);
  if (ode_type == CP_EXPL) N_VDestroy(yptemp);
  if (cpspils_mem->s_pfree != NULL) (cpspils_mem->s_pfree)(cp_mem);
  free(cpspils_mem); cpspils_mem = NULL;

  return(0);
}

