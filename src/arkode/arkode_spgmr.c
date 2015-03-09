/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 *---------------------------------------------------------------
 * This is the implementation file for the ARKSPGMR linear solver.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_spgmr.h>
#include "arkode_spils_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_spgmr.h>
#include <sundials/sundials_math.h>

/* Constants */
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* ARKSPGMR linit, lsetup, lsolve, and lfree routines */
static int ARKSpgmrInit(ARKodeMem ark_mem);
static int ARKSpgmrSetup(ARKodeMem ark_mem, int convfail, 
			 N_Vector ypred, N_Vector fpred, 
			 booleantype *jcurPtr, N_Vector vtemp1,
			 N_Vector vtemp2, N_Vector vtemp3);
static int ARKSpgmrSolve(ARKodeMem ark_mem, N_Vector b, 
			 N_Vector weight, N_Vector ynow, 
			 N_Vector fnow);
static void ARKSpgmrFree(ARKodeMem ark_mem);

/* ARKSPGMR minit, msetup, msolve, and mfree routines */
static int ARKMassSpgmrInit(ARKodeMem ark_mem);
static int ARKMassSpgmrSetup(ARKodeMem ark_mem, N_Vector vtemp1, 
			     N_Vector vtemp2, N_Vector vtemp3);
static int ARKMassSpgmrSolve(ARKodeMem ark_mem, N_Vector b, 
			     N_Vector weight);
static void ARKMassSpgmrFree(ARKodeMem ark_mem);


/*---------------------------------------------------------------
 ARKSpgmr:

 This routine initializes the memory record and sets various 
 function fields specific to the Spgmr linear solver module. 
 ARKSpgmr first calls the existing lfree routine if this is not 
 NULL.  It then sets the ark_linit, ark_lsetup, ark_lsolve, 
 ark_lfree fields in (*arkode_mem) to be ARKSpgmrInit, 
 ARKSpgmrSetup, ARKSpgmrSolve, and ARKSpgmrFree, respectively.  
 It allocates memory for a structure of type ARKSpilsMemRec and
 sets the ark_lmem field in (*arkode_mem) to the address of this
 structure.  It sets setupNonNull in (*arkode_mem), and sets 
 various fields in the ARKSpilsMemRec structure. Finally, 
 ARKSpgmr allocates memory for ytemp and x, and calls
 SpgmrMalloc to allocate memory for the Spgmr solver.
---------------------------------------------------------------*/
int ARKSpgmr(void *arkode_mem, int pretype, int maxl)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  SpgmrMem spgmr_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if N_VDotProd and N_VProd are present */
  if ((ark_mem->ark_tempv->ops->nvdotprod == NULL) ||
      (ark_mem->ark_tempv->ops->nvprod == NULL)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  ark_mem->ark_linit  = ARKSpgmrInit;
  ark_mem->ark_lsetup = ARKSpgmrSetup;
  ark_mem->ark_lsolve = ARKSpgmrSolve;
  ark_mem->ark_lfree  = ARKSpgmrFree;
  ark_mem->ark_lsolve_type = 0;

  /* Get memory for ARKSpilsMemRec */
  arkspils_mem = NULL;
  arkspils_mem = (ARKSpilsMem) malloc(sizeof(struct ARKSpilsMemRec));
  if (arkspils_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set ILS type */
  arkspils_mem->s_type = SPILS_SPGMR;

  /* Set Spgmr parameters that have been passed in call sequence */
  arkspils_mem->s_pretype    = pretype;
  mxl = arkspils_mem->s_maxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;

  /* Set defaults for Jacobian-related fields */
  arkspils_mem->s_jtimesDQ = TRUE;
  arkspils_mem->s_jtimes   = NULL;
  arkspils_mem->s_j_data   = NULL;

  /* Set defaults for preconditioner-related fields */
  arkspils_mem->s_pset   = NULL;
  arkspils_mem->s_psolve = NULL;
  arkspils_mem->s_pfree  = NULL;
  arkspils_mem->s_P_data = ark_mem->ark_user_data;

  /* Initialize counters */
  arkspils_mem->s_npe = arkspils_mem->s_nli = 0;
  arkspils_mem->s_nps = arkspils_mem->s_ncfl = 0;
  arkspils_mem->s_nstlpre = arkspils_mem->s_njtimes = 0;
  arkspils_mem->s_nfes = 0;

  /* Set default values for the rest of the Spgmr parameters */
  arkspils_mem->s_gstype = MODIFIED_GS;
  arkspils_mem->s_eplifac = ARKSPILS_EPLIN;
  arkspils_mem->s_last_flag  = ARKSPILS_SUCCESS;
  ark_mem->ark_setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_BAD_PRETYPE);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  arkspils_mem->s_ytemp = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->s_ytemp == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_MEM_FAIL);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  arkspils_mem->s_x = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->s_x == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->s_ytemp);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, arkspils_mem->s_ytemp);
  arkspils_mem->s_sqrtN = SUNRsqrt( N_VDotProd(arkspils_mem->s_ytemp, 
					    arkspils_mem->s_ytemp) );

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = NULL;
  spgmr_mem = SpgmrMalloc(mxl, ark_mem->ark_tempv);
  if (spgmr_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKSpgmr", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->s_ytemp);
    N_VDestroy(arkspils_mem->s_x);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }
  
  /* Attach SPGMR memory to spils memory structure */
  arkspils_mem->s_spils_mem = (void *) spgmr_mem;

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkspils_mem;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpgmrInit:

 This routine does remaining initializations specific to the Spgmr 
 linear solver.
---------------------------------------------------------------*/
static int ARKSpgmrInit(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Initialize counters */
  arkspils_mem->s_npe = arkspils_mem->s_nli = 0;
  arkspils_mem->s_nps = arkspils_mem->s_ncfl = 0;
  arkspils_mem->s_nstlpre = arkspils_mem->s_njtimes = 0;
  arkspils_mem->s_nfes = 0;

  /* Check for legal combination pretype - psolve */
  if ((arkspils_mem->s_pretype != PREC_NONE) 
      && (arkspils_mem->s_psolve == NULL)) {
    arkProcessError(ark_mem, -1, "ARKSPGMR", "ARKSpgmrInit", 
		    MSGS_PSOLVE_REQ);
    arkspils_mem->s_last_flag = ARKSPILS_ILL_INPUT;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning (pretype != PREC_NONE)
     and there is a preconditioning setup phase (pset != NULL)             */
  ark_mem->ark_setupNonNull = (arkspils_mem->s_pretype != PREC_NONE) 
    && (arkspils_mem->s_pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (arkspils_mem->s_jtimesDQ) {
    arkspils_mem->s_jtimes = ARKSpilsDQJtimes;
    arkspils_mem->s_j_data = ark_mem;
  } else {
    arkspils_mem->s_j_data = ark_mem->ark_user_data;
  }

  arkspils_mem->s_last_flag = ARKSPILS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 ARKSpgmrSetup:

 This routine does the setup operations for the Spgmr linear 
 solver. It makes a decision as to whether or not to signal for 
 re-evaluation of Jacobian data in the pset routine, based on 
 various state variables, then it calls pset.  If we signal for 
 re-evaluation, then we reset jcur = *jcurPtr to TRUE, regardless 
 of the pset output. In any case, if jcur == TRUE, we increment 
 npe and save nst in nstlpre.
---------------------------------------------------------------*/
static int ARKSpgmrSetup(ARKodeMem ark_mem, int convfail, 
			 N_Vector ypred, N_Vector fpred, 
			 booleantype *jcurPtr, N_Vector vtemp1,
			 N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  retval;
  ARKSpilsMem arkspils_mem;

  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkspils_mem->s_nstlpre + ARKSPILS_MSBPRE) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKSPILS_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  retval = arkspils_mem->s_pset(ark_mem->ark_tn, ypred, fpred, jok, 
				jcurPtr, ark_mem->ark_gamma, 
				arkspils_mem->s_P_data, vtemp1, 
				vtemp2, vtemp3);
  if (retval < 0) {
    arkProcessError(ark_mem, SPGMR_PSET_FAIL_UNREC, "ARKSPGMR", 
		    "ARKSpgmrSetup", MSGS_PSET_FAILED);
    arkspils_mem->s_last_flag = SPGMR_PSET_FAIL_UNREC;
  }
  if (retval > 0) {
    arkspils_mem->s_last_flag = SPGMR_PSET_FAIL_REC;
  }

  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    arkspils_mem->s_npe++;
    arkspils_mem->s_nstlpre = ark_mem->ark_nst;
  }

  arkspils_mem->s_last_flag = SPGMR_SUCCESS;

  /* Return the same value that pset returned */
  return(retval);
}


/*---------------------------------------------------------------
 ARKSpgmrSolve:

 This routine handles the call to the generic solver SpgmrSolve
 for the solution of the linear system Ax = b with the SPGMR 
 method, without restarts.  The solution x is returned in the 
 vector b.

 If the WRMS norm of b is small, we return x = b (if this is the 
 first Newton iteration) or x = 0 (if a later Newton iteration).

 Otherwise, we set the tolerance parameter and initial guess 
 (x = 0), call SpgmrSolve, and copy the solution x into b.  The 
 x-scaling and b-scaling arrays are both equal to weight, and no 
 restarts are allowed.

 The counters nli, nps, and ncfl are incremented, and the return
 value is set according to the success of SpgmrSolve.  The 
 success flag is returned if SpgmrSolve converged, or if this is 
 the first Newton iteration and the residual norm was reduced 
 below its initial value.
---------------------------------------------------------------*/
static int ARKSpgmrSolve(ARKodeMem ark_mem, N_Vector b, 
			 N_Vector weight, N_Vector ynow, 
			 N_Vector fnow)
{
  realtype bnorm, res_norm;
  ARKSpilsMem arkspils_mem;
  SpgmrMem spgmr_mem;
  int nli_inc, nps_inc, retval;
  
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  spgmr_mem = (SpgmrMem) arkspils_mem->s_spils_mem;

  /* Test norm(b); if small, return x = 0 or x = b */
  arkspils_mem->s_deltar = arkspils_mem->s_eplifac * ark_mem->ark_eLTE; 

  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= arkspils_mem->s_deltar) {
    if (ark_mem->ark_mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  arkspils_mem->s_ycur = ynow;
  arkspils_mem->s_fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SpgmrSolve */  
  arkspils_mem->s_delta = arkspils_mem->s_deltar * arkspils_mem->s_sqrtN;
  N_VConst(ZERO, arkspils_mem->s_x);
  /* N_VConst(ark_mem->ark_uround, arkspils_mem->s_x); */
  
  /* Call SpgmrSolve and copy x to b */
  retval = SpgmrSolve(spgmr_mem, ark_mem, arkspils_mem->s_x, b, 
		      arkspils_mem->s_pretype, arkspils_mem->s_gstype, 
		      arkspils_mem->s_delta, 0, ark_mem, weight, weight, 
		      ARKSpilsAtimes, ARKSpilsPSolve, &res_norm, 
		      &nli_inc, &nps_inc);
  N_VScale(ONE, arkspils_mem->s_x, b);
  
  /* Increment counters nli, nps, and ncfl */
  arkspils_mem->s_nli += nli_inc;
  arkspils_mem->s_nps += nps_inc;
  if (retval != SPGMR_SUCCESS) arkspils_mem->s_ncfl++;

  /* Log solver statistics to diagnostics file (if requested) */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "      kry  %19.16g  %19.16g  %i  %i\n", 
	    bnorm, res_norm, nli_inc, nps_inc);
  
  /* Interpret return value from SpgmrSolve */
  arkspils_mem->s_last_flag = retval;

  switch(retval) {

  case SPGMR_SUCCESS:
    return(0);
    break;
  case SPGMR_RES_REDUCED:
    /* allow reduction but not solution on first Newton iteration, 
       otherwise return with a recoverable failure */
    if (ark_mem->ark_mnewt == 0) return(0);
    else                         return(1);
    break;
  case SPGMR_CONV_FAIL:
    return(1);
    break;
  case SPGMR_QRFACT_FAIL:
    return(1);
    break;
  case SPGMR_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPGMR_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPGMR_MEM_NULL:
    return(-1);
    break;
  case SPGMR_ATIMES_FAIL_UNREC:
    arkProcessError(ark_mem, SPGMR_ATIMES_FAIL_UNREC, "ARKSPGMR", 
		    "ARKSpgmrSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPGMR_PSOLVE_FAIL_UNREC:
    arkProcessError(ark_mem, SPGMR_PSOLVE_FAIL_UNREC, "ARKSPGMR", 
		    "ARKSpgmrSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  case SPGMR_GS_FAIL:
    return(-1);
    break;
  case SPGMR_QRSOL_FAIL:
    return(-1);
    break;
  }

  return(0);
}


/*---------------------------------------------------------------
 ARKSpgmrFree:

 This routine frees memory specific to the Spgmr linear solver.
---------------------------------------------------------------*/
static void ARKSpgmrFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  SpgmrMem spgmr_mem;

  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  
  N_VDestroy(arkspils_mem->s_ytemp);
  N_VDestroy(arkspils_mem->s_x);

  spgmr_mem = (SpgmrMem) arkspils_mem->s_spils_mem;
  SpgmrFree(spgmr_mem);

  if (arkspils_mem->s_pfree != NULL) (arkspils_mem->s_pfree)(ark_mem);

  free(arkspils_mem);
  ark_mem->ark_lmem = NULL;
}




/*---------------------------------------------------------------
 ARKMassSpgmr:

 This routine initializes the memory record and sets various 
 function fields specific to the Spgmr mass matrix solver module. 
 ARKMassSpgmr first calls the existing mfree routine if this is 
 not NULL.  It then sets the ark_minit, ark_msetup, ark_msolve, 
 ark_mfree fields in (*arkode_mem) to be ARKMassSpgmrInit, 
 ARKMassSpgmrSetup, ARKMassSpgmrSolve, and ARKMassSpgmrFree, 
 respectively.  It allocates memory for a structure of type 
 ARKSpilsMassMemRec and sets the ark_mass_mem field in 
 (*arkode_mem) to the address of this structure.  It sets 
 MassSetupNonNull in (*arkode_mem), and sets various fields in
 the ARKSpilsMassMemRec structure, allocates memory for ytemp
 and x, and calls SpgmrMalloc to allocate memory for the Spgmr 
 solver.
---------------------------------------------------------------*/
int ARKMassSpgmr(void *arkode_mem, int pretype, int maxl, 
		 ARKSpilsMassTimesVecFn mtimes, void *mtimes_data)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;
  SpgmrMem spgmr_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if N_VDotProd is present */
  if(ark_mem->ark_tempv->ops->nvdotprod == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  if (ark_mem->ark_mfree != NULL) ark_mem->ark_mfree(ark_mem);

  /* Set four main function fields in ark_mem, enable mass matrix */
  ark_mem->ark_mass_matrix = TRUE;
  ark_mem->ark_minit  = ARKMassSpgmrInit;
  ark_mem->ark_msetup = ARKMassSpgmrSetup;
  ark_mem->ark_msolve = ARKMassSpgmrSolve;
  ark_mem->ark_mfree  = ARKMassSpgmrFree;
  ark_mem->ark_msolve_type = 0;

  /* Get memory for ARKSpilsMassMemRec */
  arkspils_mem = NULL;
  arkspils_mem = (ARKSpilsMassMem) malloc(sizeof(struct ARKSpilsMassMemRec));
  if (arkspils_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set mass-matrix-vector product routine */
  ark_mem->ark_mtimes      = mtimes;
  ark_mem->ark_mtimes_data = mtimes_data;

  /* Set ILS type */
  arkspils_mem->s_type = SPILS_SPGMR;

  /* Set Spgmr parameters that have been passed in call sequence */
  arkspils_mem->s_pretype    = pretype;
  mxl = arkspils_mem->s_maxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;

  /* Set defaults for preconditioner-related fields */
  arkspils_mem->s_pset   = NULL;
  arkspils_mem->s_psolve = NULL;
  arkspils_mem->s_pfree  = NULL;
  arkspils_mem->s_P_data = ark_mem->ark_user_data;

  /* Initialize counters */
  arkspils_mem->s_npe = arkspils_mem->s_nli  = 0;
  arkspils_mem->s_nps = arkspils_mem->s_ncfl = 0;

  /* Set default values for the rest of the Spgmr parameters */
  arkspils_mem->s_gstype        = MODIFIED_GS;
  arkspils_mem->s_eplifac       = ARKSPILS_EPLIN;
  arkspils_mem->s_last_flag     = ARKSPILS_SUCCESS;
  ark_mem->ark_MassSetupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_BAD_PRETYPE);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  arkspils_mem->s_ytemp = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->s_ytemp == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_MEM_FAIL);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  arkspils_mem->s_x = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->s_x == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->s_ytemp);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, arkspils_mem->s_ytemp);
  arkspils_mem->s_sqrtN = SUNRsqrt( N_VDotProd(arkspils_mem->s_ytemp, 
					    arkspils_mem->s_ytemp) );

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = NULL;
  spgmr_mem = SpgmrMalloc(mxl, ark_mem->ark_tempv);
  if (spgmr_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPGMR", 
		    "ARKMassSpgmr", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->s_ytemp);
    N_VDestroy(arkspils_mem->s_x);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }
  
  /* Attach SPGMR memory to spils memory structure */
  arkspils_mem->s_spils_mem = (void *) spgmr_mem;

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_mass_mem = arkspils_mem;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKMassSpgmrInit:

 This routine does remaining initializations specific to the Spgmr 
 linear solver.
---------------------------------------------------------------*/
static int ARKMassSpgmrInit(ARKodeMem ark_mem)
{
  ARKSpilsMassMem arkspils_mem;
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Initialize counters */
  arkspils_mem->s_npe = arkspils_mem->s_nli  = 0;
  arkspils_mem->s_nps = arkspils_mem->s_ncfl = 0;

  /* Check for legal combination pretype - psolve */
  if ((arkspils_mem->s_pretype != PREC_NONE) 
      && (arkspils_mem->s_psolve == NULL)) {
    arkProcessError(ark_mem, -1, "ARKSPGMR", 
		    "ARKMassSpgmrInit", MSGS_PSOLVE_REQ);
    arkspils_mem->s_last_flag = ARKSPILS_ILL_INPUT;
    return(-1);
  }

  /* Set MassSetupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE) and there is a preconditioning setup 
     phase (pset != NULL)             */
  ark_mem->ark_MassSetupNonNull = (arkspils_mem->s_pretype != PREC_NONE) 
    && (arkspils_mem->s_pset != NULL);

  arkspils_mem->s_last_flag = ARKSPILS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 ARKMassSpgmrSetup:

 This routine does the setup operations for the Spgmr mass matrix
 solver. It calls pset and increments npe.
---------------------------------------------------------------*/
static int ARKMassSpgmrSetup(ARKodeMem ark_mem, N_Vector vtemp1, 
			     N_Vector vtemp2, N_Vector vtemp3)
{
  int  retval;
  ARKSpilsMassMem arkspils_mem;
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Call pset routine */
  retval = arkspils_mem->s_pset(ark_mem->ark_tn, 
				arkspils_mem->s_P_data, 
				vtemp1, vtemp2, vtemp3);
  arkspils_mem->s_npe++;
  if (retval < 0) {
    arkProcessError(ark_mem, SPGMR_PSET_FAIL_UNREC, "ARKSPGMR", 
		    "ARKMassSpgmrSetup", MSGS_PSET_FAILED);
    arkspils_mem->s_last_flag = SPGMR_PSET_FAIL_UNREC;
  }
  if (retval > 0) {
    arkspils_mem->s_last_flag = SPGMR_PSET_FAIL_REC;
  }
  if (retval == 0) {
    arkspils_mem->s_last_flag = SPGMR_SUCCESS;
  }

  /* Return the same value that pset returned */
  return(retval);
}


/*---------------------------------------------------------------
 ARKMassSpgmrSolve:

 This routine handles the call to the generic solver SpgmrSolve
 for the solution of the mass matrix system Mx = b with the SPGMR 
 method, without restarts.  The solution x is returned in the 
 vector b.

 We set the tolerance parameter and initial guess (x = 0), call 
 SpgmrSolve, and copy the solution x into b.  The x-scaling and 
 b-scaling arrays are both equal to weight, and no restarts are
 allowed.

 The counters nli, nps, and ncfl are incremented, and the return
 value is set according to the success of SpgmrSolve.
---------------------------------------------------------------*/
static int ARKMassSpgmrSolve(ARKodeMem ark_mem, N_Vector b, 
			     N_Vector weight)
{
  realtype res_norm;
  ARKSpilsMassMem arkspils_mem;
  SpgmrMem spgmr_mem;
  int nli_inc, nps_inc, retval;
  
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;
  spgmr_mem = (SpgmrMem) arkspils_mem->s_spils_mem;

  /* Set inputs delta and initial guess x = 0 to SpgmrSolve */  
  arkspils_mem->s_deltar = arkspils_mem->s_eplifac * ark_mem->ark_eLTE; 
  arkspils_mem->s_delta  = arkspils_mem->s_deltar * arkspils_mem->s_sqrtN;
  N_VConst(ZERO, arkspils_mem->s_x);
  
  /* Call SpgmrSolve and copy x to b */
  retval = SpgmrSolve(spgmr_mem, ark_mem, arkspils_mem->s_x, b, 
		      arkspils_mem->s_pretype, arkspils_mem->s_gstype, 
		      arkspils_mem->s_delta, 0, ark_mem, weight, weight, 
		      ARKSpilsMtimes, ARKSpilsMPSolve, &res_norm, 
		      &nli_inc, &nps_inc);
  N_VScale(ONE, arkspils_mem->s_x, b);
  
  /* Increment counters nli, nps, and ncfl */
  arkspils_mem->s_nli += nli_inc;
  arkspils_mem->s_nps += nps_inc;
  if (retval != SPGMR_SUCCESS) arkspils_mem->s_ncfl++;

  /* Log solver statistics to diagnostics file (if requested) */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "      mass  %19.16g  %i  %i\n", 
	    res_norm, nli_inc, nps_inc);
  
  /* Interpret return value from SpgmrSolve */
  arkspils_mem->s_last_flag = retval;

  switch(retval) {

  case SPGMR_SUCCESS:
    return(0);
    break;
  case SPGMR_RES_REDUCED:
    return(1);
    break;
  case SPGMR_CONV_FAIL:
    return(1);
    break;
  case SPGMR_QRFACT_FAIL:
    return(1);
    break;
  case SPGMR_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPGMR_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPGMR_MEM_NULL:
    return(-1);
    break;
  case SPGMR_ATIMES_FAIL_UNREC:
    arkProcessError(ark_mem, SPGMR_ATIMES_FAIL_UNREC, "ARKSPGMR", 
		    "ARKMassSpgmrSolve", MSGS_MTIMES_FAILED);    
    return(-1);
    break;
  case SPGMR_PSOLVE_FAIL_UNREC:
    arkProcessError(ark_mem, SPGMR_PSOLVE_FAIL_UNREC, "ARKSPGMR", 
		    "ARKMassSpgmrSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  case SPGMR_GS_FAIL:
    return(-1);
    break;
  case SPGMR_QRSOL_FAIL:
    return(-1);
    break;
  }

  return(0);
}


/*---------------------------------------------------------------
 ARKMassSpgmrFree:

 This routine frees memory specific to the Spgmr linear solver.
---------------------------------------------------------------*/
static void ARKMassSpgmrFree(ARKodeMem ark_mem)
{
  ARKSpilsMassMem arkspils_mem;
  SpgmrMem spgmr_mem;
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;
  
  N_VDestroy(arkspils_mem->s_ytemp);
  N_VDestroy(arkspils_mem->s_x);

  spgmr_mem = (SpgmrMem) arkspils_mem->s_spils_mem;
  SpgmrFree(spgmr_mem);

  if (arkspils_mem->s_pfree != NULL) (arkspils_mem->s_pfree)(ark_mem);

  free(arkspils_mem);
  ark_mem->ark_mass_mem = NULL;
}


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
