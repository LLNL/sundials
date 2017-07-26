/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSPBCG linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes_spbcgs.h>
#include "cvodes_spils_impl.h"
#include "cvodes_impl.h"

#include <sundials/sundials_spbcgs.h>
#include <sundials/sundials_math.h>

/* Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* CVSPBCG linit, lsetup, lsolve, and lfree routines */

static int CVSpbcgInit(CVodeMem cv_mem);

static int CVSpbcgSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3);

static int CVSpbcgSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ynow, N_Vector fnow);

static int CVSpbcgFree(CVodeMem cv_mem);

/* CVSPBCG lfreeB function */

static int CVSpbcgFreeB(CVodeBMem cvB_mem);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcg
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the Spbcg linear solver module. CVSpbcg first
 * calls the existing lfree routine if this is not NULL. It then sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be CVSpbcgInit, CVSpbcgSetup, CVSpbcgSolve, and CVSpbcgFree,
 * respectively. It allocates memory for a structure of type
 * CVSpilsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure. It sets setupNonNull in (*cvode_mem),
 * and sets various fields in the CVSpilsMemRec structure.
 * Finally, CVSpbcg allocates memory for ytemp and x, and calls
 * SpbcgMalloc to allocate memory for the Spbcg solver.
 * -----------------------------------------------------------------
 */

int CVSpbcg(void *cvode_mem, int pretype, int maxl)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  SpbcgMem spbcg_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPBCG", "CVSpbcg", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if N_VDotProd is present */
  if (cv_mem->cv_tempv->ops->nvdotprod == NULL) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPBCG", "CVSpbcg", MSGS_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }

  if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  cv_mem->cv_linit  = CVSpbcgInit;
  cv_mem->cv_lsetup = CVSpbcgSetup;
  cv_mem->cv_lsolve = CVSpbcgSolve;
  cv_mem->cv_lfree  = CVSpbcgFree;

  /* Get memory for CVSpilsMemRec */
  cvspils_mem = NULL;
  cvspils_mem = (CVSpilsMem) malloc(sizeof(struct CVSpilsMemRec));
  if (cvspils_mem == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPBCG", "CVSpbcg", MSGS_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Set ILS type */
  cvspils_mem->s_type = SPILS_SPBCG;

  /* Set Spbcg parameters that have been passed in call sequence */
  cvspils_mem->s_pretype = pretype;
  mxl = cvspils_mem->s_maxl = (maxl <= 0) ? CVSPILS_MAXL : maxl;

  /* Set defaults for Jacobian-related fileds */
  cvspils_mem->s_jtimesDQ = TRUE;
  cvspils_mem->s_jtimes   = NULL;
  cvspils_mem->s_j_data   = NULL;

  /* Set defaults for preconditioner-related fields */
  cvspils_mem->s_pset   = NULL;
  cvspils_mem->s_psolve = NULL;
  cvspils_mem->s_pfree  = NULL;
  cvspils_mem->s_P_data = cv_mem->cv_user_data;

  /* Set default values for the rest of the Spbcg parameters */
  cvspils_mem->s_eplifac = CVSPILS_EPLIN;

  cvspils_mem->s_last_flag = CVSPILS_SUCCESS;

  cvSpilsInitializeCounters(cvspils_mem);

  cv_mem->cv_setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPBCG", "CVSpbcg", MSGS_BAD_PRETYPE);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  cvspils_mem->s_ytemp = N_VClone(cv_mem->cv_tempv);
  if (cvspils_mem->s_ytemp == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPBCG", "CVSpbcg", MSGS_MEM_FAIL);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }
  cvspils_mem->s_x = N_VClone(cv_mem->cv_tempv);
  if (cvspils_mem->s_x == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPBCG", "CVSpbcg", MSGS_MEM_FAIL);
    N_VDestroy(cvspils_mem->s_ytemp);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, cvspils_mem->s_ytemp);
  cvspils_mem->s_sqrtN = SUNRsqrt(N_VDotProd(cvspils_mem->s_ytemp, cvspils_mem->s_ytemp));

  /* Call SpbcgMalloc to allocate workspace for Spbcg */
  spbcg_mem = NULL;
  spbcg_mem = SpbcgMalloc(mxl, cv_mem->cv_tempv);
  if (spbcg_mem == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPBCG", "CVSpbcg", MSGS_MEM_FAIL);
    N_VDestroy(cvspils_mem->s_ytemp);
    N_VDestroy(cvspils_mem->s_x);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }
  
  /* Attach SPBCG memory to spils memory structure */
  cvspils_mem->s_spils_mem = (void *) spbcg_mem;

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvspils_mem;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Spbcg
 * linear solver.
 * -----------------------------------------------------------------
 */

static int CVSpbcgInit(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  SpbcgMem spbcg_mem;

  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  spbcg_mem = (SpbcgMem) cvspils_mem->s_spils_mem;


  /* Initialize counters */
  cvSpilsInitializeCounters(cvspils_mem);

  /* Check for legal combination pretype - psolve */
  if ((cvspils_mem->s_pretype != PREC_NONE) && (cvspils_mem->s_psolve == NULL)) {
    cvProcessError(cv_mem, -1, "CVSPBCG", "CVSpbcgInit", MSGS_PSOLVE_REQ);
    cvspils_mem->s_last_flag = CVSPILS_ILL_INPUT;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE)  and there is a preconditioning
     setup phase (pset != NULL) */
  cv_mem->cv_setupNonNull = (cvspils_mem->s_pretype != PREC_NONE) &&
    (cvspils_mem->s_pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (cvspils_mem->s_jtimesDQ) {
    cvspils_mem->s_jtimes = CVSpilsDQJtimes;
    cvspils_mem->s_j_data = cv_mem;
  } else {
    cvspils_mem->s_j_data = cv_mem->cv_user_data;
  }

  /*  Set maxl in the SPBCG memory in case it was changed by the user */
  spbcg_mem->l_max = cvspils_mem->s_maxl;

  cvspils_mem->s_last_flag = CVSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the Spbcg linear solver.
 * It makes a decision as to whether or not to signal for reevaluation
 * of Jacobian data in the pset routine, based on various state
 * variables, then it calls pset. If we signal for reevaluation,
 * then we reset jcur = *jcurPtr to TRUE, regardless of the pset output.
 * In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.
 * -----------------------------------------------------------------
 */

static int CVSpbcgSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  retval;
  CVSpilsMem cvspils_mem;

  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  jbad = (cv_mem->cv_nst == 0) ||
    (cv_mem->cv_nst > cvspils_mem->s_nstlpre + CVSPILS_MSBPRE) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVSPILS_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  retval = cvspils_mem->s_pset(cv_mem->cv_tn, ypred, fpred, jok,
                               jcurPtr, cv_mem->cv_gamma, cvspils_mem->s_P_data, 
                               vtemp1, vtemp2, vtemp3);
  if (retval < 0) {
    cvProcessError(cv_mem, SPBCG_PSET_FAIL_UNREC, "CVSPBCG", "CVSpbcgSetup", MSGS_PSET_FAILED);
    cvspils_mem->s_last_flag = SPBCG_PSET_FAIL_UNREC;
  }
  if (retval > 0) {
    cvspils_mem->s_last_flag = SPBCG_PSET_FAIL_REC;
  }

  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    cvspils_mem->s_npe++;
    cvspils_mem->s_nstlpre = cv_mem->cv_nst;
  }

  cvspils_mem->s_last_flag = SPBCG_SUCCESS;

  /* Return the same value that pset returned */
  return(retval);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgSolve
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

static int CVSpbcgSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ynow, N_Vector fnow)
{
  realtype bnorm, res_norm;
  CVSpilsMem cvspils_mem;
  SpbcgMem spbcg_mem;
  int nli_inc, nps_inc, retval;
  
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  spbcg_mem = (SpbcgMem) cvspils_mem->s_spils_mem;

  /* Test norm(b); if small, return x = 0 or x = b */
  cvspils_mem->s_deltar = cvspils_mem->s_eplifac * cv_mem->cv_tq[4]; 

  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= cvspils_mem->s_deltar) {
    if (cv_mem->cv_mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  cvspils_mem->s_ycur = ynow;
  cvspils_mem->s_fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SpbcgSolve */  
  cvspils_mem->s_delta = cvspils_mem->s_deltar * cvspils_mem->s_sqrtN;
  N_VConst(ZERO, cvspils_mem->s_x);
  
  /* Call SpbcgSolve and copy x to b */
  retval = SpbcgSolve(spbcg_mem, cv_mem, cvspils_mem->s_x, b, cvspils_mem->s_pretype, cvspils_mem->s_delta,
                      cv_mem, weight, weight, CVSpilsAtimes, CVSpilsPSolve,
                      &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, cvspils_mem->s_x, b);
  
  /* Increment counters nli, nps, and ncfl */
  cvspils_mem->s_nli += nli_inc;
  cvspils_mem->s_nps += nps_inc;
  if (retval != SPBCG_SUCCESS) cvspils_mem->s_ncfl++;

  /* Interpret return value from SpbcgSolve */

  cvspils_mem->s_last_flag = retval;

  switch(retval) {

  case SPBCG_SUCCESS:
    return(0);
    break;
  case SPBCG_RES_REDUCED:
    if (cv_mem->cv_mnewt == 0) return(0);
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
    cvProcessError(cv_mem, SPBCG_ATIMES_FAIL_UNREC, "CVSPBCG", "CVSpbcgSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPBCG_PSOLVE_FAIL_UNREC:
    cvProcessError(cv_mem, SPBCG_PSOLVE_FAIL_UNREC, "CVSPBCG", "CVSpbcgSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);  

}

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the Spbcg linear solver.
 * -----------------------------------------------------------------
 */

static int CVSpbcgFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  SpbcgMem spbcg_mem;

  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  N_VDestroy(cvspils_mem->s_ytemp);
  N_VDestroy(cvspils_mem->s_x);

  spbcg_mem = (SpbcgMem) cvspils_mem->s_spils_mem;
  SpbcgFree(spbcg_mem);

  if (cvspils_mem->s_pfree != NULL) (cvspils_mem->s_pfree)(cv_mem);

  free(cvspils_mem);
  cv_mem->cv_lmem = NULL;
  
  return(0);
}


/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * CVSpbcgB
 *
 * Wrapper for the backward phase
 */

int CVSpbcgB(void *cvode_mem, int which, int pretypeB, int maxlB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  CVSpilsMemB cvspilsB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPBCG", "CVSpbcgB", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVSPILS_NO_ADJ, "CVSPBCG", "CVSpbcgB", MSGS_NO_ADJ);
    return(CVSPILS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPBCG", "CVSpbcgB", MSGS_BAD_WHICH);
    return(CVSPILS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Get memory for CVSpilsMemRecB */
  cvspilsB_mem = NULL;
  cvspilsB_mem = (CVSpilsMemB) malloc(sizeof(struct CVSpilsMemRecB));
  if (cvspilsB_mem == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPBCG", "CVSpbcgB", MSGS_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  cvspilsB_mem->s_psetB = NULL;
  cvspilsB_mem->s_psolveB = NULL;
  cvspilsB_mem->s_P_dataB = NULL;

  /* initialize Jacobian function */
  cvspilsB_mem->s_jtimesB = NULL;

  /* attach lmemB and lfree */
  cvB_mem->cv_lmem = cvspilsB_mem;
  cvB_mem->cv_lfree = CVSpbcgFreeB;

  flag = CVSpbcg(cvodeB_mem, pretypeB, maxlB);

  if (flag != CVSPILS_SUCCESS) {
    free(cvspilsB_mem);
    cvspilsB_mem = NULL;
  }

  return(flag);
}

/*
 * CVSpbcgFreeB 
 */


static int CVSpbcgFreeB(CVodeBMem cvB_mem)
{
  CVSpilsMemB cvspilsB_mem;

  cvspilsB_mem = (CVSpilsMemB) (cvB_mem->cv_lmem);

  free(cvspilsB_mem);

  return(0);
}
