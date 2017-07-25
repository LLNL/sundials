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
 * This is the implementation file for the CVSPTFQMR linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvode/cvode_sptfqmr.h>
#include "cvode_spils_impl.h"
#include "cvode_impl.h"

#include <sundials/sundials_sptfqmr.h>
#include <sundials/sundials_math.h>

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* CVSPTFQMR linit, lsetup, lsolve, and lfree routines */

static int CVSptfqmrInit(CVodeMem cv_mem);

static int CVSptfqmrSetup(CVodeMem cv_mem, N_Vector vtemp1,
                          N_Vector vtemp2, N_Vector vtemp3);

static int CVSptfqmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                          N_Vector ynow, N_Vector fnow);

static int CVSptfqmrFree(CVodeMem cv_mem);


/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmr
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the Sptfqmr linear solver module. CVSptfqmr first
 * calls the existing lfree routine if this is not NULL. It then sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be CVSptfqmrInit, CVSptfqmrSetup, CVSptfqmrSolve, and CVSptfqmrFree,
 * respectively. It allocates memory for a structure of type
 * CVSpilsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure. It sets setupNonNull in (*cvode_mem),
 * and sets various fields in the CVSpilsMemRec structure.
 * Finally, CVSptfqmr allocates memory for ytemp and x, and calls
 * SptfqmrMalloc to allocate memory for the Sptfqmr solver.
 * -----------------------------------------------------------------
 */

int CVSptfqmr(void *cvode_mem, int pretype, int maxl)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  SptfqmrMem sptfqmr_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPTFQMR", "CVSptfqmr", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if N_VDotProd is present */
  if (cv_mem->cv_tempv->ops->nvdotprod == NULL) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPTFQMR", "CVSptfqmr", MSGS_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }

  if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  cv_mem->cv_linit  = CVSptfqmrInit;
  cv_mem->cv_lsetup = CVSptfqmrSetup;
  cv_mem->cv_lsolve = CVSptfqmrSolve;
  cv_mem->cv_lfree  = CVSptfqmrFree;

  /* Get memory for CVSpilsMemRec */
  cvspils_mem = NULL;
  cvspils_mem = (CVSpilsMem) malloc(sizeof(struct CVSpilsMemRec));
  if (cvspils_mem == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPTFQMR", "CVSptfqmr", MSGS_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Set ILS type */
  cvspils_mem->s_type = SPILS_SPTFQMR;

  /* Set Sptfqmr parameters that have been passed in call sequence */
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

  /* Set default values for the rest of the Sptfqmr parameters */
  cvspils_mem->s_eplifac = CVSPILS_EPLIN;

  cvspils_mem->s_last_flag = CVSPILS_SUCCESS;

  cvSpilsInitializeCounters(cvspils_mem);

  cv_mem->cv_setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPTFQMR", "CVSptfqmr", MSGS_BAD_PRETYPE);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */

  cvspils_mem->s_ytemp = N_VClone(cv_mem->cv_tempv);
  if (cvspils_mem->s_ytemp == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPTFQMR", "CVSptfqmr", MSGS_MEM_FAIL);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }

  cvspils_mem->s_x = N_VClone(cv_mem->cv_tempv);
  if (cvspils_mem->s_x == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPTFQMR", "CVSptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(cvspils_mem->s_ytemp);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, cvspils_mem->s_ytemp);
  cvspils_mem->s_sqrtN = SUNRsqrt(N_VDotProd(cvspils_mem->s_ytemp, cvspils_mem->s_ytemp));

  /* Call SptfqmrMalloc to allocate workspace for Sptfqmr */
  sptfqmr_mem = NULL;
  sptfqmr_mem = SptfqmrMalloc(mxl, cv_mem->cv_tempv);
  if (sptfqmr_mem == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPTFQMR", "CVSptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(cvspils_mem->s_ytemp);
    N_VDestroy(cvspils_mem->s_x);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }
  
  /* Attach SPTFQMR memory to spils memory structure */
  cvspils_mem->s_spils_mem = (void *) sptfqmr_mem;

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvspils_mem;

  return(CVSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Sptfqmr
 * linear solver.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrInit(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  SptfqmrMem sptfqmr_mem;

  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  sptfqmr_mem = (SptfqmrMem) cvspils_mem->s_spils_mem;

  /* Initialize counters */
  cvSpilsInitializeCounters(cvspils_mem);

  /* Check for legal combination pretype - psolve */
  if ((cvspils_mem->s_pretype != PREC_NONE) && (cvspils_mem->s_psolve == NULL)) {
    cvProcessError(cv_mem, -1, "CVSPTFQMR", "CVSptfqmrInit", MSGS_PSOLVE_REQ);
    cvspils_mem->s_last_flag = CVSPILS_ILL_INPUT;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE)  and there is a preconditioning
     setup phase (pset != NULL) */
  cv_mem->cv_setupNonNull = (cvspils_mem->s_pretype != PREC_NONE) && (cvspils_mem->s_pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (cvspils_mem->s_jtimesDQ) {
    cvspils_mem->s_jtimes = CVSpilsDQJtimes;
    cvspils_mem->s_j_data = cv_mem;
  } else {
    cvspils_mem->s_j_data = cv_mem->cv_user_data;
  }

  /*  Set maxl in the SPTFQMR memory in case it was changed by the user */
  sptfqmr_mem->l_max  = cvspils_mem->s_maxl;

  cvspils_mem->s_last_flag = CVSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the Sptfqmr linear solver.
 * It makes a decision as to whether or not to signal for reevaluation
 * of Jacobian data in the pset routine, based on various state
 * variables, then it calls pset. If we signal for reevaluation,
 * then we reset jcur = *jcurPtr to TRUE, regardless of the pset output.
 * In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrSetup(CVodeMem cv_mem, N_Vector vtemp1,
			  N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  retval;
  CVSpilsMem cvspils_mem;
  N_Vector ypred, fpred;
  
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  ypred = cv_mem->cv_zn[0];
  fpred = cv_mem->cv_ftemp;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  jbad = (cv_mem->cv_nst == 0) || (cv_mem->cv_nst > cvspils_mem->s_nstlpre + CVSPILS_MSBPRE) ||
      ((cv_mem->cv_convfail == CV_FAIL_BAD_J) && (dgamma < CVSPILS_DGMAX)) ||
      (cv_mem->cv_convfail == CV_FAIL_OTHER);
  cv_mem->cv_jcur = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  retval = cvspils_mem->s_pset(cv_mem->cv_tn, ypred, fpred, jok, &cv_mem->cv_jcur,
                               cv_mem->cv_gamma, cvspils_mem->s_P_data,
                               vtemp1, vtemp2, vtemp3);
  if (retval < 0) {
    cvProcessError(cv_mem, SPTFQMR_PSET_FAIL_UNREC, "CVSPTFQMR", "CVSptfqmrSetup", MSGS_PSET_FAILED);
    cvspils_mem->s_last_flag = SPTFQMR_PSET_FAIL_UNREC;
  }
  if (retval > 0) {
    cvspils_mem->s_last_flag = SPTFQMR_PSET_FAIL_REC;
  }

  if (jbad) cv_mem->cv_jcur = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (cv_mem->cv_jcur) {
    cvspils_mem->s_npe++;
    cvspils_mem->s_nstlpre = cv_mem->cv_nst;
  }

  cvspils_mem->s_last_flag = SPTFQMR_SUCCESS;

  /* Return the same value that pset returned */
  return(retval);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic solver SptfqmrSolve
 * for the solution of the linear system Ax = b with the SPTFQMR method.
 * The solution x is returned in the vector b.
 *
 * If the WRMS norm of b is small, we return x = b (if this is the first
 * Newton iteration) or x = 0 (if a later Newton iteration).
 *
 * Otherwise, we set the tolerance parameter and initial guess (x = 0),
 * call SptfqmrSolve, and copy the solution x into b. The x-scaling and
 * b-scaling arrays are both equal to weight.
 *
 * The counters nli, nps, and ncfl are incremented, and the return value
 * is set according to the success of SptfqmrSolve. The success flag is
 * returned if SptfqmrSolve converged, or if this is the first Newton
 * iteration and the residual norm was reduced below its initial value.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
			  N_Vector ynow, N_Vector fnow)
{
  realtype bnorm, res_norm;
  CVSpilsMem cvspils_mem;
  SptfqmrMem sptfqmr_mem;
  int nli_inc, nps_inc, retval;
  
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  sptfqmr_mem = (SptfqmrMem) cvspils_mem->s_spils_mem;

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

  /* Set inputs delta and initial guess x = 0 to SptfqmrSolve */  
  cvspils_mem->s_delta = cvspils_mem->s_deltar * cvspils_mem->s_sqrtN;
  N_VConst(ZERO, cvspils_mem->s_x);
  
  /* Call SptfqmrSolve and copy x to b */
  retval = SptfqmrSolve(sptfqmr_mem, cv_mem, cvspils_mem->s_x, b, cvspils_mem->s_pretype, cvspils_mem->s_delta,
                        cv_mem, weight, weight, CVSpilsAtimes, CVSpilsPSolve,
                        &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, cvspils_mem->s_x, b);
  
  /* Increment counters nli, nps, and ncfl */
  cvspils_mem->s_nli += nli_inc;
  cvspils_mem->s_nps += nps_inc;
  if (retval != SPTFQMR_SUCCESS) cvspils_mem->s_ncfl++;

  /* Interpret return value from SpgmrSolve */

  cvspils_mem->s_last_flag = retval;

  switch(retval) {

  case SPTFQMR_SUCCESS:
    return(0);
    break;
  case SPTFQMR_RES_REDUCED:
    if (cv_mem->cv_mnewt == 0) return(0);
    else            return(1);
    break;
  case SPTFQMR_CONV_FAIL:
    return(1);
    break;
  case SPTFQMR_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPTFQMR_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPTFQMR_MEM_NULL:
    return(-1);
    break;
  case SPTFQMR_ATIMES_FAIL_UNREC:
    cvProcessError(cv_mem, SPTFQMR_ATIMES_FAIL_UNREC, "CVSPTFQMR", "CVSptfqmrSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPTFQMR_PSOLVE_FAIL_UNREC:
    cvProcessError(cv_mem, SPTFQMR_PSOLVE_FAIL_UNREC, "CVSPTFQMR", "CVSptfqmrSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmrFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the Sptfqmr linear solver.
 * -----------------------------------------------------------------
 */

static int CVSptfqmrFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  SptfqmrMem sptfqmr_mem;
    
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  N_VDestroy(cvspils_mem->s_ytemp);
  N_VDestroy(cvspils_mem->s_x);

  sptfqmr_mem = (SptfqmrMem) cvspils_mem->s_spils_mem;
  SptfqmrFree(sptfqmr_mem);

  if (cvspils_mem->s_pfree != NULL) (cvspils_mem->s_pfree)(cv_mem);

  free(cvspils_mem);
  cv_mem->cv_lmem = NULL;

  return(0);
}

