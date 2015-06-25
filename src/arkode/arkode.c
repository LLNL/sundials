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
 * This is the implementation file for the main ARKODE integrator.
 * It is independent of the ARKODE linear solver in use.
 *--------------------------------------------------------------*/

/*===============================================================
             Import Header Files                                 
===============================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#define NO_DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
#endif

#ifdef __GNUC__
#define SUNDIALS_UNUSED __attribute__ ((unused))
#else
#define SUNDIALS_UNUSED
#endif

/*===============================================================
             Private Functions Prototypes
===============================================================*/
static void arkPrintMem(ARKodeMem ark_mem) SUNDIALS_UNUSED;
static booleantype arkCheckNvector(N_Vector tmpl);
static booleantype arkAllocVectors(ARKodeMem ark_mem, 
				   N_Vector tmpl);
static int arkAllocRKVectors(ARKodeMem ark_mem);
static void arkFreeVectors(ARKodeMem ark_mem);

static int arkAllocFPData(ARKodeMem ark_mem);
static int arkResizeFPData(ARKodeMem ark_mem, 
			   ARKVecResizeFn resize,
			   void *resize_data,
			   long int lrw_diff,
			   long int liw_diff);
static void arkFreeFPData(ARKodeMem ark_mem);

static int arkInitialSetup(ARKodeMem ark_mem);
static int arkHin(ARKodeMem ark_mem, realtype tout);
static realtype arkUpperBoundH0(ARKodeMem ark_mem, 
				realtype tdist);
static int arkYddNorm(ARKodeMem ark_mem, realtype hg, 
		      realtype *yddnrm);
static int arkSetButcherTables(ARKodeMem ark_mem);
static int arkCheckButcherTables(ARKodeMem ark_mem);

static int arkStep(ARKodeMem ark_mem);
static void arkPredict(ARKodeMem ark_mem, int istage);
static int arkSet(ARKodeMem ark_mem);
static int arkComputeSolutions(ARKodeMem ark_mem, realtype *dsm);
static int arkDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
			  realtype saved_t, int *nefPtr, 
			  realtype dsm);
static int arkCompleteStep(ARKodeMem ark_mem, realtype dsm);
static int arkPrepareNextStep(ARKodeMem ark_mem);

static int arkNls(ARKodeMem ark_mem, int nflag);
static int arkNlsResid(ARKodeMem ark_mem, N_Vector y, 
			N_Vector fy, N_Vector r);
static int arkNlsNewton(ARKodeMem ark_mem, int nflag);
static int arkNlsAccelFP(ARKodeMem ark_mem, int nflag);
static int arkAndersonAcc(ARKodeMem ark_mem, N_Vector gval, 
			  N_Vector fv, N_Vector x, N_Vector xold, 
			  int iter, realtype *R, realtype *gamma);
static int arkLs(ARKodeMem ark_mem, int nflag);
static int arkHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, 
			  realtype saved_t, int *ncfPtr);

static int arkHandleFailure(ARKodeMem ark_mem,int flag);

static int arkDenseEval(ARKodeMem ark_mem, realtype tau,
			int d, int order, N_Vector yout);

static int arkFullRHS(ARKodeMem ark_mem, realtype t, 
		      N_Vector y, N_Vector tmp, N_Vector f);

static int arkEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, 
		       N_Vector weight);
static int arkEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, 
		       N_Vector weight);
static int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, 
		       N_Vector weight);
static int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, 
		       N_Vector weight);

static int arkAdapt(ARKodeMem ark_mem);
static int arkAdaptPID(ARKodeMem ark_mem, realtype *hnew);
static int arkAdaptPI(ARKodeMem ark_mem, realtype *hnew);
static int arkAdaptI(ARKodeMem ark_mem, realtype *hnew);
static int arkAdaptExpGus(ARKodeMem ark_mem, realtype *hnew);
static int arkAdaptImpGus(ARKodeMem ark_mem, realtype *hnew);
static int arkAdaptImExGus(ARKodeMem ark_mem, realtype *hnew);

static int arkRootCheck1(ARKodeMem ark_mem);
static int arkRootCheck2(ARKodeMem ark_mem);
static int arkRootCheck3(ARKodeMem ark_mem);
static int arkRootfind(ARKodeMem ark_mem);




/*===============================================================
  EXPORTED FUNCTIONS
===============================================================*/

/*---------------------------------------------------------------
 ARKodeCreate:

 ARKodeCreate creates an internal memory block for a problem to 
 be solved by ARKODE.  If successful, ARKodeCreate returns a 
 pointer to the problem memory. This pointer should be passed to
 ARKodeInit. If an initialization error occurs, ARKodeCreate 
 prints an error message to standard err and returns NULL. 
---------------------------------------------------------------*/
void *ARKodeCreate()
{
  int i, j, iret;
  ARKodeMem ark_mem;

  ark_mem = NULL;
  ark_mem = (ARKodeMem) malloc(sizeof(struct ARKodeMemRec));
  if (ark_mem == NULL) {
    arkProcessError(NULL, 0, "ARKODE", "ARKodeCreate", 
		    MSGARK_ARKMEM_FAIL);
    return(NULL);
  }

  /* Zero out ark_mem */
  memset(ark_mem, 0, sizeof(struct ARKodeMemRec));

  /* Set uround */
  ark_mem->ark_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  iret = ARKodeSetDefaults((void *)ark_mem);
  if (iret != ARK_SUCCESS) {
    arkProcessError(NULL, 0, "ARKODE", "ARKodeCreate", 
		    "Error setting default solver options");
    return(NULL);
  }

  /* Initialize internal RK parameters and coefficients */
  ark_mem->ark_stages = 0;
  ark_mem->ark_istage = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ae,i,j) = ZERO;
      ARK_A(ark_mem->ark_Ai,i,j) = ZERO;
    }
    ark_mem->ark_c[i]   = ZERO;
    ark_mem->ark_b[i]   = ZERO;
    ark_mem->ark_b2[i]  = ZERO;
    ark_mem->ark_Fi[i]  = NULL;
    ark_mem->ark_Fe[i]  = NULL;
  }

  /* Initialize root finding variables */
  ark_mem->ark_glo     = NULL;
  ark_mem->ark_ghi     = NULL;
  ark_mem->ark_grout   = NULL;
  ark_mem->ark_iroots  = NULL;
  ark_mem->ark_rootdir = NULL;
  ark_mem->ark_gfun    = NULL;
  ark_mem->ark_nrtfn   = 0;
  ark_mem->ark_gactive = NULL;
  ark_mem->ark_mxgnull = 1;

  /* Set default nonlinear solver choice to Newton,
     initialize fixed-point solver variables */
  ark_mem->ark_use_fp   = FALSE;
  ark_mem->ark_fp_R     = NULL;
  ark_mem->ark_fp_gamma = NULL;
  ark_mem->ark_fp_df    = NULL;
  ark_mem->ark_fp_dg    = NULL;
  ark_mem->ark_fp_q     = NULL;
  ark_mem->ark_fp_fval  = NULL;
  ark_mem->ark_fp_fold  = NULL;
  ark_mem->ark_fp_gold  = NULL;

  /* Initialize diagnostics reporting variables */
  ark_mem->ark_report  = FALSE;
  ark_mem->ark_diagfp  = NULL;

  /* Initialize lrw and liw */
  ark_mem->ark_lrw = 58;   /* to be updated */
  ark_mem->ark_liw = 40;   /* to be updated */

  /* No mallocs have been done yet */
  ark_mem->ark_VabstolMallocDone  = FALSE;
  ark_mem->ark_VRabstolMallocDone = FALSE;
  ark_mem->ark_MallocDone         = FALSE;

  /* Return pointer to ARKODE memory block */
  return((void *)ark_mem);
}


/*---------------------------------------------------------------
 ARKodeInit:
 
 ARKodeInit allocates and initializes memory for a problem. All 
 inputs are checked for errors. If any error occurs during 
 initialization, it is reported to the file whose file pointer 
 is errfp and an error flag is returned. Otherwise, it returns 
 ARK_SUCCESS.  This routine must be called prior to calling 
 ARKode to evolve the problem.
---------------------------------------------------------------*/
int ARKodeInit(void *arkode_mem, ARKRhsFn fe, ARKRhsFn fi, 
	       realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  booleantype nvectorOK, allocOK;
  long int lrw1, liw1;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeInit", 
		    MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeInit", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Set implicit/explicit problem based on function pointers */
  if (fe == NULL) ark_mem->ark_implicit = TRUE;
  if (fi == NULL) ark_mem->ark_explicit = TRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (ark_mem->ark_implicit && ark_mem->ark_explicit) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
  		    "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = arkCheckNvector(y0);
  if (!nvectorOK) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeInit", MSGARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Set space requirements for one N_Vector */
  if (y0->ops->nvspace != NULL) {
    N_VSpace(y0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  ark_mem->ark_lrw1 = lrw1;
  ark_mem->ark_liw1 = liw1;


  /* Allocate the solver vectors (using y0 as a template) */
  allocOK = arkAllocVectors(ark_mem, y0);
  if (!allocOK) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* All error checking is complete at this point */

  /* Copy the input parameters into ARKODE state */
  ark_mem->ark_fe = fe;
  ark_mem->ark_fi = fi;
  ark_mem->ark_tn = t0;
  ark_mem->ark_tnew = t0;

  /* Set step parameters */
  ark_mem->ark_hold  = ZERO;
  ark_mem->ark_tolsf = ONE;

  /* Set the linear solver addresses to NULL.
     (We check != NULL later, in ARKode.) */
  ark_mem->ark_linit  = NULL;
  ark_mem->ark_lsetup = NULL;
  ark_mem->ark_lsolve = NULL;
  ark_mem->ark_lfree  = NULL;
  ark_mem->ark_lmem   = NULL;
  ark_mem->ark_lsolve_type = -1;

  /* Set the mass matrix to identity, and mass matrix solver 
     addresses to NULL. (We check != NULL later, in ARKode.) */
  ark_mem->ark_mass_matrix = FALSE;
  ark_mem->ark_minit       = NULL;
  ark_mem->ark_msetup      = NULL;
  ark_mem->ark_msolve      = NULL;
  ark_mem->ark_mfree       = NULL;
  ark_mem->ark_mtimes      = NULL;
  ark_mem->ark_mtimes_data = NULL;
  ark_mem->ark_mass_mem    = NULL;
  ark_mem->ark_msolve_type = -1;

  /* Initialize ycur */
  N_VScale(ONE, y0, ark_mem->ark_ycur);

  /* Initialize error history (1.0 => met target accuracy) */
  ark_mem->ark_hadapt_ehist[0] = ONE;
  ark_mem->ark_hadapt_ehist[1] = ONE;
  ark_mem->ark_hadapt_ehist[2] = ONE;
  ark_mem->ark_eLTE = 1.0;

  /* Initialize step history */
  ark_mem->ark_hadapt_hhist[0] = ZERO;
  ark_mem->ark_hadapt_hhist[1] = ZERO;
  ark_mem->ark_hadapt_hhist[2] = ZERO;

  /* Initialize all the counters */
  ark_mem->ark_nst          = 0;
  ark_mem->ark_nst_acc      = 0;
  ark_mem->ark_nst_exp      = 0;
  ark_mem->ark_nst_attempts = 0;
  ark_mem->ark_nfe          = 0;
  ark_mem->ark_nfi          = 0;
  ark_mem->ark_ncfn         = 0;
  ark_mem->ark_nmassfails   = 0;
  ark_mem->ark_netf         = 0;
  ark_mem->ark_nni          = 0;
  ark_mem->ark_nsetups      = 0;
  ark_mem->ark_nhnil        = 0;
  ark_mem->ark_nstlp        = 0;
  ark_mem->ark_nge          = 0;
  ark_mem->ark_irfnd        = 0;
  ark_mem->ark_mass_solves  = 0;
  ark_mem->ark_mass_mult    = 0;

  /* Initialize other integrator optional outputs */
  ark_mem->ark_h0u    = ZERO;
  ark_mem->ark_next_h = ZERO;

  /* Initially, rwt should point to ewt */
  ark_mem->ark_rwt_is_ewt = TRUE;

  /* Indicate that problem size is new */
  ark_mem->ark_resized = TRUE;
  ark_mem->ark_firststage = TRUE;

  /* Problem has been successfully initialized */
  ark_mem->ark_MallocDone = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeReInit:

 ARKodeReInit re-initializes ARKODE's memory for a problem, 
 assuming it has already been allocated in a prior ARKodeInit 
 call.  All problem specification inputs are checked for errors.
 If any error occurs during initialization, it is reported to 
 the file whose file pointer is errfp.  This routine should only
 be called after ARKodeInit, and only when the problem dynamics 
 or desired solvers have changed dramatically, so that the
 problem integration should resume as if started from scratch.

 The return value is ARK_SUCCESS = 0 if no errors occurred, or
 a negative value otherwise.
---------------------------------------------------------------*/
int ARKodeReInit(void *arkode_mem, ARKRhsFn fe, ARKRhsFn fi, 
		 realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
 
  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeReInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */
  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeReInit", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeReInit", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }
  
  /* Set implicit/explicit problem based on function pointers */
  ark_mem->ark_implicit = ark_mem->ark_explicit = FALSE;
  if (fe == NULL) ark_mem->ark_implicit = TRUE;
  if (fi == NULL) ark_mem->ark_explicit = TRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (ark_mem->ark_implicit && ark_mem->ark_explicit) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
  		    "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  ark_mem->ark_fe = fe;
  ark_mem->ark_fi = fi;
  ark_mem->ark_tn = t0;
  ark_mem->ark_tnew = t0;
  
  /* Set step parameters */
  ark_mem->ark_hold  = ZERO;
  ark_mem->ark_tolsf = ONE;

  /* Do not reset the linear solver addresses to NULL.  This means 
     that if the user does not re-set these manually, we'll re-use 
     the linear solver routines that were set during ARKodeInit. */

  /* Initialize ycur */
  N_VScale(ONE, y0, ark_mem->ark_ycur);
 
  /* Initialize error history (1.0 => met target accuracy) */
  ark_mem->ark_hadapt_ehist[0] = ONE;
  ark_mem->ark_hadapt_ehist[1] = ONE;
  ark_mem->ark_hadapt_ehist[2] = ONE;
  ark_mem->ark_eLTE = 1.0;

  /* Initialize step history */
  ark_mem->ark_hadapt_hhist[0] = ZERO;
  ark_mem->ark_hadapt_hhist[1] = ZERO;
  ark_mem->ark_hadapt_hhist[2] = ZERO;

  /* Initialize all the counters */
  ark_mem->ark_nst          = 0;
  ark_mem->ark_nst_acc      = 0;
  ark_mem->ark_nst_exp      = 0;
  ark_mem->ark_nst_attempts = 0;
  ark_mem->ark_nfe          = 0;
  ark_mem->ark_nfi          = 0;
  ark_mem->ark_ncfn         = 0;
  ark_mem->ark_nmassfails   = 0;
  ark_mem->ark_netf         = 0;
  ark_mem->ark_nni          = 0;
  ark_mem->ark_nsetups      = 0;
  ark_mem->ark_nhnil        = 0;
  ark_mem->ark_nstlp        = 0;
  ark_mem->ark_nge          = 0;
  ark_mem->ark_irfnd        = 0;
  ark_mem->ark_mass_solves  = 0;
  ark_mem->ark_mass_mult    = 0;

  /* Indicate that problem size is new */
  ark_mem->ark_resized = TRUE;
  ark_mem->ark_firststage = TRUE;

  /* Initialize other integrator optional outputs */
  ark_mem->ark_h0u    = ZERO;
  ark_mem->ark_next_h = ZERO;

  /* Problem has been successfully re-initialized */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeResize:

 ARKodeResize re-initializes ARKODE's memory for a problem with a
 changing vector size.  It is assumed that the problem dynamics 
 before and after the vector resize will be comparable, so that 
 all time-stepping heuristics prior to calling ARKodeResize 
 remain valid after the call.  If instead the dynamics should be 
 re-calibrated, the ARKode memory structure should be deleted 
 with a call to ARKodeFree, and re-created with calls to 
 ARKodeCreate and ARKodeInit.

 To aid in the vector-resize operation, the user can supply a 
 vector resize function, that will take as input an N_Vector with
 the previous size, and return as output a corresponding vector 
 of the new size.  If this function (of type ARKVecResizeFn) is 
 not supplied (i.e. is set to NULL), then all existing N_Vectors 
 will be destroyed and re-cloned from the input vector.

 In the case that the dynamical time scale should be modified 
 slightly from the previous time scale, an input "hscale" is 
 allowed, that will re-scale the upcoming time step by the 
 specified factor.  If a value <= 0 is specified, the default of 
 1.0 will be used.

 Other arguments:
   arkode_mem       Existing ARKode memory data structure.
   ynew             The newly-sized solution vector, holding 
                    the current dependent variable values.
   t0               The current value of the independent 
                    variable.
   resize_data      User-supplied data structure that will be 
                    passed to the supplied resize function.

 The return value is ARK_SUCCESS = 0 if no errors occurred, or
 a negative value otherwise.
---------------------------------------------------------------*/
int ARKodeResize(void *arkode_mem, N_Vector y0, 
		 realtype hscale, realtype t0, 
		 ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  long int lrw1=0, liw1=0;
  long int lrw_diff, liw_diff;
  int ier, i;
 
  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeResize", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */
  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeResize", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeResize", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }
  
  /* Copy the input parameters into ARKODE state */
  ark_mem->ark_tn = t0;
  ark_mem->ark_tnew = t0;

  /* Update time-stepping parameters */
  /*   adjust upcoming step size depending on hscale */
  if (hscale < 0.0)  hscale = 1.0;
  if (hscale != 1.0) {

    /* Encode hscale into ark_mem structure */
    ark_mem->ark_eta = hscale;
    ark_mem->ark_hprime *= hscale;

    /* If next step would overtake tstop, adjust stepsize */
    if ( ark_mem->ark_tstopset ) 
      if ( (ark_mem->ark_tn+ark_mem->ark_hprime-ark_mem->ark_tstop)*ark_mem->ark_hprime > ZERO ) {
	ark_mem->ark_hprime = (ark_mem->ark_tstop-ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
	ark_mem->ark_eta = ark_mem->ark_hprime/ark_mem->ark_h;
      }

  }

  /* Determing change in vector sizes */
  if (y0->ops->nvspace != NULL) 
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->ark_lrw1;
  liw_diff = liw1 - ark_mem->ark_liw1;
  ark_mem->ark_lrw1 = lrw1;
  ark_mem->ark_liw1 = liw1;

  /* Resize the ARKode vectors */
  /*     Vabstol */
  if (ark_mem->ark_Vabstol != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_Vabstol);
      ark_mem->ark_Vabstol = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_Vabstol, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     VRabstol */
  if (ark_mem->ark_VRabstol != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_VRabstol);
      ark_mem->ark_VRabstol = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_VRabstol, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     ark_Fe */
  for (i=0; i<ARK_S_MAX; i++) {
    if (ark_mem->ark_Fe[i] != NULL) {
      if (resize == NULL) {
	N_VDestroy(ark_mem->ark_Fe[i]);
	ark_mem->ark_Fe[i] = N_VClone(y0);
      } else {
	if (resize(ark_mem->ark_Fe[i], y0, resize_data)) {
	  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			  "ARKodeResize", MSGARK_RESIZE_FAIL);
	  return(ARK_ILL_INPUT);
	}
      }
      ark_mem->ark_lrw += lrw_diff;
      ark_mem->ark_liw += liw_diff;
    }
  }
  /*     ark_Fi */
  for (i=0; i<ARK_S_MAX; i++) {
    if (ark_mem->ark_Fi[i] != NULL) {
      if (resize == NULL) {
	N_VDestroy(ark_mem->ark_Fi[i]);
	ark_mem->ark_Fi[i] = N_VClone(y0);
      } else {
	if (resize(ark_mem->ark_Fi[i], y0, resize_data)) {
	  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			  "ARKodeResize", MSGARK_RESIZE_FAIL);
	  return(ARK_ILL_INPUT);
	}
      }
      ark_mem->ark_lrw += lrw_diff;
      ark_mem->ark_liw += liw_diff;
    }
  }
  /*     ewt */
  if (ark_mem->ark_ewt != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_ewt);
      ark_mem->ark_ewt = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_ewt, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     rwt  */
  if (ark_mem->ark_rwt_is_ewt) {      /* update pointer to ewt */
    ark_mem->ark_rwt = ark_mem->ark_ewt;
  } else {                            /* resize if distinct from ewt */
    if (ark_mem->ark_rwt != NULL) {
      if (resize == NULL) {
	N_VDestroy(ark_mem->ark_rwt);
	ark_mem->ark_rwt = N_VClone(y0);
      } else {
	if (resize(ark_mem->ark_rwt, y0, resize_data)) {
	  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			  "ARKodeResize", MSGARK_RESIZE_FAIL);
	  return(ARK_ILL_INPUT);
	}
      }
      ark_mem->ark_lrw += lrw_diff;
      ark_mem->ark_liw += liw_diff;
    }
  }
  /*     ycur */
  if (ark_mem->ark_ycur != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_ycur);
      ark_mem->ark_ycur = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_ycur, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     acor */
  if (ark_mem->ark_acor != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_acor);
      ark_mem->ark_acor = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_acor, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     sdata */
  if (ark_mem->ark_sdata != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_sdata);
      ark_mem->ark_sdata = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_sdata, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     tempv */
  if (ark_mem->ark_tempv != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_tempv);
      ark_mem->ark_tempv = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_tempv, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     ftemp */
  if (ark_mem->ark_ftemp != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_ftemp);
      ark_mem->ark_ftemp = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_ftemp, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     fold */
  if (ark_mem->ark_fold != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_fold);
      ark_mem->ark_fold = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_fold, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     fnew */
  if (ark_mem->ark_fnew != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_fnew);
      ark_mem->ark_fnew = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_fnew, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     yold */
  if (ark_mem->ark_yold != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_yold);
      ark_mem->ark_yold = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_yold, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }
  /*     ynew */
  if (ark_mem->ark_ynew != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_ynew);
      ark_mem->ark_ynew = N_VClone(y0);
    } else {
      if (resize(ark_mem->ark_ynew, y0, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }


  /* Resize fixed-point solver memory */
  if (ark_mem->ark_use_fp) {
    ier = arkResizeFPData(ark_mem, resize, resize_data, 
			  lrw_diff, liw_diff);
    if (ier != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		      "ARKodeResize", MSGARK_MEM_FAIL);
      return(ARK_MEM_FAIL);
    }
  }


  /* Copy y0 into ark_ycur to set the current solution */
  N_VScale(ONE, y0, ark_mem->ark_ycur);

  /* Load updated error weights */
  ier = ark_mem->ark_efun(ark_mem->ark_ycur,
			  ark_mem->ark_ewt, 
			  ark_mem->ark_e_data);
  if (ier != 0) {
    if (ark_mem->ark_itol == ARK_WF) 
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeResize", MSGARK_EWT_FAIL);
    else 
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeResize", MSGARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }

  /* Load updated residual weights */
  if (!ark_mem->ark_rwt_is_ewt) {
    ier = ark_mem->ark_rfun(ark_mem->ark_ycur,
			    ark_mem->ark_rwt, 
			    ark_mem->ark_r_data);
    if (ier != 0) {
      if (ark_mem->ark_itol == ARK_WF) 
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_RWT_FAIL);
      else 
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeResize", MSGARK_BAD_RWT);
      return(ARK_ILL_INPUT);
    }
  }

  /* Indicate that problem size is new */
  ark_mem->ark_resized = TRUE;
  ark_mem->ark_firststage = TRUE;
  
  /* Problem has been successfully re-sized */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSStolerances:
 ARKodeSVtolerances:
 ARKodeWFtolerances:

 These functions specify the integration tolerances. One of them
 SHOULD be called before the first call to ARKode; otherwise 
 default values of reltol=1e-4 and abstol=1e-9 will be used, 
 which may be entirely incorrect for a specific problem.

 ARKodeSStolerances specifies scalar relative and absolute 
   tolerances.

 ARKodeSVtolerances specifies scalar relative tolerance and a 
   vector absolute tolerance (a potentially different absolute 
   tolerance for each vector component).

 ARKodeWFtolerances specifies a user-provides function (of type
   ARKEwtFn) which will be called to set the error weight vector.
---------------------------------------------------------------*/
int ARKodeSStolerances(void *arkode_mem, realtype reltol, 
		       realtype abstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */
  if (reltol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }
  if (abstol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  ark_mem->ark_reltol    = reltol;
  ark_mem->ark_Sabstol   = abstol;
  ark_mem->ark_itol      = ARK_SS;
  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun      = arkEwtSet;
  ark_mem->ark_e_data    = NULL; /* set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}

int ARKodeSVtolerances(void *arkode_mem, realtype reltol, 
		       N_Vector abstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */
  if (reltol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }
  if (N_VMin(abstol) < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  if ( !(ark_mem->ark_VabstolMallocDone) ) {
    ark_mem->ark_Vabstol = N_VClone(ark_mem->ark_ewt);
    ark_mem->ark_lrw += ark_mem->ark_lrw1;
    ark_mem->ark_liw += ark_mem->ark_liw1;
    ark_mem->ark_VabstolMallocDone = TRUE;
  }
  N_VScale(ONE, abstol, ark_mem->ark_Vabstol);
  ark_mem->ark_reltol    = reltol;
  ark_mem->ark_itol      = ARK_SV;
  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun      = arkEwtSet;
  ark_mem->ark_e_data    = NULL; /* set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}

int ARKodeWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeWFtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeWFtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Copy tolerance data into memory */
  ark_mem->ark_itol      = ARK_WF;
  ark_mem->ark_user_efun = TRUE;
  ark_mem->ark_efun      = efun;
  ark_mem->ark_e_data    = NULL; /* set to user_data in InitialSetup */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeResStolerance:
 ARKodeResVtolerance:
 ARKodeResFtolerance:

 These functions specify the absolute residual tolerance. 
 Specification of the absolute residual tolerance is only 
 necessary for problems with non-identity mass matrices in which
 the units of the solution vector y dramatically differ from the 
 units of the ODE right-hand side f(t,y).  If this occurs, one 
 of these routines SHOULD be called before the first call to 
 ARKode; otherwise the default value of rabstol=1e-9 will be 
 used, which may be entirely incorrect for a specific problem.

 ARKodeResStolerances specifies a scalar residual tolerance.

 ARKodeResVtolerances specifies a vector residual tolerance 
  (a potentially different absolute residual tolerance for 
   each vector component).

 ARKodeResFtolerances specifies a user-provides function (of 
   type ARKRwtFn) which will be called to set the residual 
   weight vector.
---------------------------------------------------------------*/
int ARKodeResStolerance(void *arkode_mem, realtype rabstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeResStolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeResStolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */
  if (rabstol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeResStolerances", MSGARK_BAD_RABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->ark_rwt_is_ewt) {
    ark_mem->ark_rwt_is_ewt = FALSE;
    ark_mem->ark_rwt = N_VClone(ark_mem->ark_ewt);
    ark_mem->ark_lrw += ark_mem->ark_lrw1;
    ark_mem->ark_liw += ark_mem->ark_liw1;
  }

  /* Copy tolerances into memory */
  ark_mem->ark_SRabstol  = rabstol;
  ark_mem->ark_ritol     = ARK_SS;
  ark_mem->ark_user_rfun = FALSE;
  ark_mem->ark_rfun      = arkRwtSet;
  ark_mem->ark_r_data    = NULL; /* set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}

int ARKodeResVtolerance(void *arkode_mem, N_Vector rabstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeResVtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeResVtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */
  if (N_VMin(rabstol) < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeResVtolerances", MSGARK_BAD_RABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->ark_rwt_is_ewt) {
    ark_mem->ark_rwt_is_ewt = FALSE;
    ark_mem->ark_rwt = N_VClone(ark_mem->ark_ewt);
    ark_mem->ark_lrw += ark_mem->ark_lrw1;
    ark_mem->ark_liw += ark_mem->ark_liw1;
  }

  /* Copy tolerances into memory */
  if ( !(ark_mem->ark_VRabstolMallocDone) ) {
    ark_mem->ark_VRabstol = N_VClone(ark_mem->ark_rwt);
    ark_mem->ark_lrw += ark_mem->ark_lrw1;
    ark_mem->ark_liw += ark_mem->ark_liw1;
    ark_mem->ark_VRabstolMallocDone = TRUE;
  }
  N_VScale(ONE, rabstol, ark_mem->ark_VRabstol);
  ark_mem->ark_ritol     = ARK_SV;
  ark_mem->ark_user_rfun = FALSE;
  ark_mem->ark_rfun      = arkRwtSet;
  ark_mem->ark_r_data    = NULL; /* set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}

int ARKodeResFtolerance(void *arkode_mem, ARKRwtFn rfun)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeResFtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeResFtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->ark_rwt_is_ewt) {
    ark_mem->ark_rwt_is_ewt = FALSE;
    ark_mem->ark_rwt = N_VClone(ark_mem->ark_ewt);
    ark_mem->ark_lrw += ark_mem->ark_lrw1;
    ark_mem->ark_liw += ark_mem->ark_liw1;
  }

  /* Copy tolerance data into memory */
  ark_mem->ark_ritol     = ARK_WF;
  ark_mem->ark_user_rfun = TRUE;
  ark_mem->ark_rfun      = rfun;
  ark_mem->ark_r_data    = NULL; /* set to user_data in InitialSetup */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKode:

 This routine is the main driver of the ARKODE package. 

 It integrates over a time interval defined by the user, by 
 calling arkStep to do internal time steps.

 The first time that ARKode is called for a successfully 
 initialized problem, it computes a tentative initial step size.

 ARKode supports two modes as specified by itask: ARK_NORMAL and
 ARK_ONE_STEP.  In the ARK_NORMAL mode, the solver steps until 
 it reaches or passes tout and then interpolates to obtain 
 y(tout).  In the ARK_ONE_STEP mode, it takes one internal step 
 and returns.  The behavior of both modes can be over-rided 
 through user-specification of ark_tstop (through the
 ARKodeSetStopTime function), in which case if a solver step 
 would pass tstop, the step is shortened so that it stops at
 exactly the specified stop time, and hence interpolation of 
 y(tout) is not required.
---------------------------------------------------------------*/
int ARKode(void *arkode_mem, realtype tout, N_Vector yout, 
	   realtype *tret, int itask)
{
  ARKodeMem ark_mem;
  long int nstloc;
  int retval, hflag, kflag, istate, ir, ier, irfndp;
  int ewtsetOK;
  realtype troundoff, tout_hin, rh, nrm;
  booleantype inactive_roots;


  /*-------------------------------------
    1. Check and process inputs
  -------------------------------------*/

  /* Check if arkode_mem exists */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKode", 
		    MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */
  if (ark_mem->ark_MallocDone == FALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKode", 
		    MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  
  /* Check for yout != NULL */
  if ((ark_mem->ark_y = yout) == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		    MSGARK_YOUT_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for tret != NULL */
  if (tret == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		    MSGARK_TRET_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for valid itask */
  if ( (itask != ARK_NORMAL) && (itask != ARK_ONE_STEP) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		    MSGARK_BAD_ITASK);
    return(ARK_ILL_INPUT);
  }

  if (itask == ARK_NORMAL) ark_mem->ark_toutc = tout;
  ark_mem->ark_taskc = itask;


  /*----------------------------------------
    2. Initializations performed only at
       the first step (nst=0):
       - initial setup
       - compute initial step size
       - check for approach to tstop
       - check for approach to a root
  ----------------------------------------*/
  if (ark_mem->ark_nst == 0) {

    /* initialize tret values to initialization time */
    ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;

    /* Temporarily set ark_h and perform initial integrator setup */
    ark_mem->ark_h = SUNRabs(tout - ark_mem->ark_tn);
    if (ark_mem->ark_h == ZERO)  ark_mem->ark_h = ONE;
    ier = arkInitialSetup(ark_mem);
    if (ier!= ARK_SUCCESS) return(ier);
    
    /* Copy f(t0,y0) into ark_fold */
    /*   if (!ark_mem->ark_explicit)*/
    N_VScale(ONE, ark_mem->ark_fnew, ark_mem->ark_fold);
    
    /* Test input tstop for legality. */
    if ( ark_mem->ark_tstopset ) {
      if ( (ark_mem->ark_tstop - ark_mem->ark_tn)*(tout - ark_mem->ark_tn) <= ZERO ) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKVODE", "ARKode", 
			MSGARK_BAD_TSTOP, ark_mem->ark_tstop, ark_mem->ark_tn);
        return(ARK_ILL_INPUT);
      }
    }

    /* Set initial h (from H0 or arkHin). */
    ark_mem->ark_h = ark_mem->ark_hin;
    if ( (ark_mem->ark_h != ZERO) && 
	 ((tout-ark_mem->ark_tn)*ark_mem->ark_h < ZERO) ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		      MSGARK_BAD_H0);
      return(ARK_ILL_INPUT);
    }
    if (ark_mem->ark_h == ZERO) {
      /* Again, temporarily set ark_h for estimating an optimal value */
      ark_mem->ark_h = SUNRabs(tout - ark_mem->ark_tn);
      if (ark_mem->ark_h == ZERO)  ark_mem->ark_h = ONE;
      /* Estimate the first step size */
      tout_hin = tout;
      if ( ark_mem->ark_tstopset && 
	   (tout-ark_mem->ark_tn)*(tout-ark_mem->ark_tstop) > ZERO ) 
	tout_hin = ark_mem->ark_tstop; 
      hflag = arkHin(ark_mem, tout_hin);
      if (hflag != ARK_SUCCESS) {
        istate = arkHandleFailure(ark_mem, hflag);
        return(istate);
      }
    }
    rh = SUNRabs(ark_mem->ark_h)*ark_mem->ark_hmax_inv;
    if (rh > ONE) ark_mem->ark_h /= rh;
    if (SUNRabs(ark_mem->ark_h) < ark_mem->ark_hmin)
      ark_mem->ark_h *= ark_mem->ark_hmin/SUNRabs(ark_mem->ark_h);

    /* Check for approach to tstop */
    if (ark_mem->ark_tstopset) {
      if ( (ark_mem->ark_tn + ark_mem->ark_h - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO ) {
        ark_mem->ark_h = (ark_mem->ark_tstop - ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
      }
    }

    /* Set initial time step factors */
    ark_mem->ark_h0u    = ark_mem->ark_h;
    ark_mem->ark_hprime = ark_mem->ark_h;

    /* Check for zeros of root function g at and near t0. */
    if (ark_mem->ark_nrtfn > 0) {
      retval = arkRootCheck1(ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "arkRootCheck1", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tn);
        return(ARK_RTFUNC_FAIL);
      }
    }

  } /* end of first call block */


  /*----------------------------------------
    2b. Initializations performed only in 
        the first step after the problem has
        been resized:
       - re-initialize the linear solver
       - fills ynew, fnew and fold arrays
       - checks for approach to tstop
       - checks for root near t0
  ----------------------------------------*/
  if (ark_mem->ark_nst > 0 && ark_mem->ark_resized) {

    /* Re-initialize the linear solver (assumes that solver 
       memory has already been resized appropriately) */
    if (!ark_mem->ark_explicit) {
      if (ark_mem->ark_linit != NULL) {
	ier = ark_mem->ark_linit(ark_mem);
	if (ier != 0) {
	  arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE", 
			  "ARKode", MSGARK_LINIT_FAIL);
	  return(ARK_LINIT_FAIL);
	}
      }
    }

    /* Re-initialize the mass matrix solver (assumes that solver 
       memory has already been resized appropriately) */
    if (ark_mem->ark_minit != NULL) {
      ier = ark_mem->ark_minit(ark_mem);
      if (ier != 0) {
	arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKODE", 
			"ARKode", MSGARK_MASSINIT_FAIL);
	return(ARK_MASSINIT_FAIL);
      }
    }

    /* Fill initial ynew and fnew arrays */
    N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_ynew);
    ier = arkFullRHS(ark_mem, ark_mem->ark_tn, ark_mem->ark_ycur,
		     ark_mem->ark_ftemp, ark_mem->ark_fnew);

    /* if the problem involves a non-identity mass matrix, update fnew here */
    if (ark_mem->ark_mass_matrix) {
      N_VScale(ark_mem->ark_h, ark_mem->ark_fnew, ark_mem->ark_fnew);   /* scale RHS */
      retval = ark_mem->ark_msolve(ark_mem, ark_mem->ark_fnew, ark_mem->ark_rwt); 
      /* retval = ark_mem->ark_msolve(ark_mem, ark_mem->ark_fnew, ark_mem->ark_ewt);  */
      N_VScale(ONE/ark_mem->ark_h, ark_mem->ark_fnew, ark_mem->ark_fnew);   /* scale result */
      ark_mem->ark_mass_solves++;
      if (retval != ARK_SUCCESS) {
	ark_mem->ark_nmassfails++;
	arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE", 
			"ARKode", "Mass matrix solver failure");
	return(ARK_MASSSOLVE_FAIL);
      }
    }
  
    /* Copy f(t0,y0) into ark_fold */
    /* if (!ark_mem->ark_explicit) */
    N_VScale(ONE, ark_mem->ark_fnew, ark_mem->ark_fold);

    /* Check for legal tstop (correct direction of integration) */
    if (ark_mem->ark_tstopset) {
      if ( (ark_mem->ark_tstop - ark_mem->ark_tn)*ark_mem->ark_h < ZERO ) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			MSGARK_BAD_TSTOP, ark_mem->ark_tstop, ark_mem->ark_tn);
        return(ARK_ILL_INPUT);
      }
    }

    /* Check for zeros of root function g at and near t0. */
    if (ark_mem->ark_nrtfn > 0) {
      retval = arkRootCheck1(ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "arkRootCheck1", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tn);
        return(ARK_RTFUNC_FAIL);
      }
    }

  } /* end of first-step-after-resize call block */


  /*------------------------------------------------------
    3. At following steps, perform stop tests:
       - check for root in last step
       - check if we passed tstop
       - check if we passed tout (NORMAL mode)
       - check if current tn was returned (ONE_STEP mode)
       - check if we are close to tstop
         (adjust step size if needed)
  -------------------------------------------------------*/
  if (ark_mem->ark_nst > 0 && !ark_mem->ark_resized) {

    /* Estimate an infinitesimal time interval to be used as
       a roundoff for time quantities (based on current time 
       and step size) */
    troundoff = FUZZ_FACTOR*ark_mem->ark_uround *
      (SUNRabs(ark_mem->ark_tn) + SUNRabs(ark_mem->ark_h));

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */
    if (ark_mem->ark_nrtfn > 0) {

      irfndp = ark_mem->ark_irfnd;
      
      retval = arkRootCheck2(ark_mem);

      if (retval == CLOSERT) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkRootCheck2", 
			MSGARK_CLOSE_ROOTS, ark_mem->ark_tlo);
        return(ARK_ILL_INPUT);
      } else if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "arkRootCheck2", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tlo);
        return(ARK_RTFUNC_FAIL);
      } else if (retval == RTFOUND) {
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tlo;
        return(ARK_ROOT_RETURN);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if ( SUNRabs(ark_mem->ark_tn - ark_mem->ark_tretlast) > troundoff ) {

        retval = arkRootCheck3(ark_mem);

        if (retval == ARK_SUCCESS) {     /* no root found */
          ark_mem->ark_irfnd = 0;
          if ((irfndp == 1) && (itask == ARK_ONE_STEP)) {
            ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
            N_VScale(ONE, ark_mem->ark_ycur, yout);
            return(ARK_SUCCESS);
          }
        } else if (retval == RTFOUND) {  /* a new root was found */
          ark_mem->ark_irfnd = 1;
          ark_mem->ark_tretlast = *tret = ark_mem->ark_tlo;
          return(ARK_ROOT_RETURN);
        } else if (retval == ARK_RTFUNC_FAIL) {  /* g failed */
          arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "arkRootCheck3", 
			  MSGARK_RTFUNC_FAILED, ark_mem->ark_tlo);
          return(ARK_RTFUNC_FAIL);
        }

      }

    } /* end of root stop check */

    /* In ARK_NORMAL mode, test if tout was reached */
    if ( (itask == ARK_NORMAL) && 
	 ((ark_mem->ark_tn-tout)*ark_mem->ark_h >= ZERO) ) {
      ark_mem->ark_tretlast = *tret = tout;
      ier =  ARKodeGetDky(ark_mem, tout, 0, yout);
      if (ier != ARK_SUCCESS) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKode", MSGARK_BAD_TOUT, tout);
        return(ARK_ILL_INPUT);
      }
      return(ARK_SUCCESS);
    }

    /* In ARK_ONE_STEP mode, test if tn was returned */
    if ( itask == ARK_ONE_STEP && 
	 SUNRabs(ark_mem->ark_tn - ark_mem->ark_tretlast) > troundoff ) {
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_ycur, yout);
      return(ARK_SUCCESS);
    }

    /* Test for tn at tstop or near tstop */
    if ( ark_mem->ark_tstopset ) {

      if ( SUNRabs(ark_mem->ark_tn - ark_mem->ark_tstop) <= troundoff) {
        ier =  ARKodeGetDky(ark_mem, ark_mem->ark_tstop, 0, yout);
        if (ier != ARK_SUCCESS) {
          arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			  MSGARK_BAD_TSTOP, ark_mem->ark_tstop, ark_mem->ark_tn);
          return(ARK_ILL_INPUT);
        }
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tstop;
        ark_mem->ark_tstopset = FALSE;
        return(ARK_TSTOP_RETURN);
      }
      
      /* If next step would overtake tstop, adjust stepsize */
      if ( (ark_mem->ark_tn + ark_mem->ark_hprime - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO ) {
        ark_mem->ark_hprime = (ark_mem->ark_tstop - ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
        ark_mem->ark_eta = ark_mem->ark_hprime/ark_mem->ark_h;
      }
    }
  } /* end stopping tests block */  


  /*--------------------------------------------------
    4. Looping point for internal steps
 
       - update the ewt vector for the next step
       - check for errors (too many steps, too much
         accuracy requested, step size too small)
       - take a new step (via arkStep); stop on error 
       - perform stop tests:
         - check for root in last step taken
         - check if tout was passed
         - check if close to tstop
         - check if in ONE_STEP mode (must return)
  --------------------------------------------------*/
  nstloc = 0;
  for(;;) {
   
    ark_mem->ark_next_h = ark_mem->ark_h;
    
    /* Reset and check ewt */
    if (ark_mem->ark_nst > 0 && !ark_mem->ark_resized) {
      ewtsetOK = ark_mem->ark_efun(ark_mem->ark_ycur,
				   ark_mem->ark_ewt, 
				   ark_mem->ark_e_data);
      if (ewtsetOK != 0) {
        if (ark_mem->ark_itol == ARK_WF) 
          arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			  MSGARK_EWT_NOW_FAIL, ark_mem->ark_tn);
        else 
          arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			  MSGARK_EWT_NOW_BAD, ark_mem->ark_tn);
	
        istate = ARK_ILL_INPUT;
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
        N_VScale(ONE, ark_mem->ark_ycur, yout);
        break;
      }
    }
    
    /* Reset and check rwt */
    if (!ark_mem->ark_rwt_is_ewt) {
      if (ark_mem->ark_nst > 0 && !ark_mem->ark_resized) {
	ewtsetOK = ark_mem->ark_rfun(ark_mem->ark_ycur,
				     ark_mem->ark_rwt,
				     ark_mem->ark_r_data);
	if (ewtsetOK != 0) {
	  if (ark_mem->ark_itol == ARK_WF) 
	    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			    MSGARK_RWT_NOW_FAIL, ark_mem->ark_tn);
	  else 
	    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			    MSGARK_RWT_NOW_BAD, ark_mem->ark_tn);
	
	  istate = ARK_ILL_INPUT;
	  ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
	  N_VScale(ONE, ark_mem->ark_ycur, yout);
	  break;
	}
      }
    }
    
    /* Check for too many steps */
    if ( (ark_mem->ark_mxstep>0) && (nstloc >= ark_mem->ark_mxstep) ) {
      arkProcessError(ark_mem, ARK_TOO_MUCH_WORK, "ARKODE", "ARKode", 
		      MSGARK_MAX_STEPS, ark_mem->ark_tn);
      istate = ARK_TOO_MUCH_WORK;
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_ycur, yout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm = N_VWrmsNorm(ark_mem->ark_ycur, ark_mem->ark_ewt);
    ark_mem->ark_tolsf = ark_mem->ark_uround * nrm;
    if (ark_mem->ark_tolsf > ONE) {
      arkProcessError(ark_mem, ARK_TOO_MUCH_ACC, "ARKODE", "ARKode", 
		      MSGARK_TOO_MUCH_ACC, ark_mem->ark_tn);
      istate = ARK_TOO_MUCH_ACC;
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_ycur, yout);
      ark_mem->ark_tolsf *= TWO;
      break;
    } else {
      ark_mem->ark_tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */
    if (ark_mem->ark_tn + ark_mem->ark_h == ark_mem->ark_tn) {
      ark_mem->ark_nhnil++;
      if (ark_mem->ark_nhnil <= ark_mem->ark_mxhnil) 
        arkProcessError(ark_mem, ARK_WARNING, "ARKODE", "ARKode", 
			MSGARK_HNIL, ark_mem->ark_tn, ark_mem->ark_h);
      if (ark_mem->ark_nhnil == ark_mem->ark_mxhnil) 
        arkProcessError(ark_mem, ARK_WARNING, "ARKODE", "ARKode", 
			MSGARK_HNIL_DONE);
    }

    /* Call arkStep to take a step */
    kflag = arkStep(ark_mem);

    /* Process failed step cases, and exit loop */
    if (kflag != ARK_SUCCESS) {
      istate = arkHandleFailure(ark_mem, kflag);
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_ycur, yout);
      break;
    }
    
    nstloc++;

    /* Check for root in last step taken. */
    if (ark_mem->ark_nrtfn > 0) {

      retval = arkRootCheck3(ark_mem);
      if (retval == RTFOUND) {  /* A new root was found */
        ark_mem->ark_irfnd = 1;
        istate = ARK_ROOT_RETURN;
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tlo;
        break;
      } else if (retval == ARK_RTFUNC_FAIL) { /* g failed */
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "arkRootCheck3", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tlo);
        istate = ARK_RTFUNC_FAIL;
        break;
      }

      /* If we are at the end of the first step and we still have
       * some event functions that are inactive, issue a warning
       * as this may indicate a user error in the implementation
       * of the root function. */
      if (ark_mem->ark_nst==1) {
        inactive_roots = FALSE;
        for (ir=0; ir<ark_mem->ark_nrtfn; ir++) { 
          if (!ark_mem->ark_gactive[ir]) {
            inactive_roots = TRUE;
            break;
          }
        }
        if ((ark_mem->ark_mxgnull > 0) && inactive_roots) {
          arkProcessError(ark_mem, ARK_WARNING, "ARKODES", "ARKode", 
			  MSGARK_INACTIVE_ROOTS);
        }
      }
    }

    /* In NORMAL mode, check if tout reached */
    if ( (itask == ARK_NORMAL) &&
	 (ark_mem->ark_tn-tout)*ark_mem->ark_h >= ZERO ) {
      istate = ARK_SUCCESS;
      ark_mem->ark_tretlast = *tret = tout;
      (void) ARKodeGetDky(ark_mem, tout, 0, yout);
      ark_mem->ark_next_h = ark_mem->ark_hprime;
      break;
    }

    /* Check if tn is at tstop or near tstop */
    if ( ark_mem->ark_tstopset ) {
      troundoff = FUZZ_FACTOR*ark_mem->ark_uround * 
	(SUNRabs(ark_mem->ark_tn) + SUNRabs(ark_mem->ark_h));
      if ( SUNRabs(ark_mem->ark_tn - ark_mem->ark_tstop) <= troundoff) {
        (void) ARKodeGetDky(ark_mem, ark_mem->ark_tstop, 0, yout);
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tstop;
        ark_mem->ark_tstopset = FALSE;
        istate = ARK_TSTOP_RETURN;
        break;
      }
      if ( (ark_mem->ark_tn + ark_mem->ark_hprime - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO ) {
        ark_mem->ark_hprime = (ark_mem->ark_tstop - ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
        ark_mem->ark_eta = ark_mem->ark_hprime/ark_mem->ark_h;
      }
    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (itask == ARK_ONE_STEP) {
      istate = ARK_SUCCESS;
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_ycur, yout);
      ark_mem->ark_next_h = ark_mem->ark_hprime;
      break;
    }

  } /* end looping for internal steps */

  return(istate);
}


/*---------------------------------------------------------------
 ARKodeGetDky:

 This routine computes the k-th derivative of the interpolating
 polynomial at the time t and stores the result in the vector 
 dky. This routine internally calls arkDenseEval to perform the
 interpolation.  We have the restriction that 0 <= k <= 3.  This 
 routine uses an interpolating polynomial of degree 
 max(ark_dense_q, k), i.e. it will form a polynomial of the
 degree requested by the user through ark_dense_q, unless 
 higher-order derivatives are requested.

 This function is called by ARKode with k=0 and t=tout to perform
 interpolation of outputs, but may also be called directly by the
 user.  Note: in all cases it will be called after ark_tn has 
 been updated to correspond with the end time of the last 
 successful step.
---------------------------------------------------------------*/
int ARKodeGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, tfuzz, tp, tn1;
  int retval, degree;
  ARKodeMem ark_mem;
  
  /* Check all inputs for legality */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeGetDky", 
		    MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (dky == NULL) {
    arkProcessError(ark_mem, ARK_BAD_DKY, "ARKODE", "ARKodeGetDky", 
		    MSGARK_NULL_DKY);
    return(ARK_BAD_DKY);
  }
  if ((k < 0) || (k > 3)) {
    arkProcessError(ark_mem, ARK_BAD_K, "ARKODE", "ARKodeGetDky", 
		    MSGARK_BAD_K);
    return(ARK_BAD_K);
  }
  
  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * ark_mem->ark_uround * 
    (SUNRabs(ark_mem->ark_tn) + SUNRabs(ark_mem->ark_hold));
  if (ark_mem->ark_hold < ZERO) tfuzz = -tfuzz;
  tp = ark_mem->ark_tn - ark_mem->ark_hold - tfuzz;
  tn1 = ark_mem->ark_tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    arkProcessError(ark_mem, ARK_BAD_T, "ARKODE", "ARKodeGetDky", 
		    MSGARK_BAD_T, t, ark_mem->ark_tn-ark_mem->ark_hold, 
		    ark_mem->ark_tn);
    return(ARK_BAD_T);
  }

  /* call arkDenseEval to evaluate result */
  s = (t - ark_mem->ark_tn) / ark_mem->ark_h;
  degree = (k > ark_mem->ark_dense_q) ? k : ark_mem->ark_dense_q;
  retval = arkDenseEval(ark_mem, s, k, degree, dky);
  if (retval != ARK_SUCCESS) 
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "arkDenseEval", 
		    MSGARK_RHSFUNC_FAILED, t);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeFree:

 This routine frees the problem memory allocated by ARKodeInit.
 Such memory includes all the vectors allocated by 
 arkAllocVectors, and arkAllocRKVectors, and the memory lmem for
 the linear solver (deallocated by a call to lfree).
---------------------------------------------------------------*/
void ARKodeFree(void **arkode_mem)
{
  ARKodeMem ark_mem;

  if (*arkode_mem == NULL) return;

  ark_mem = (ARKodeMem) (*arkode_mem);
  
  arkFreeVectors(ark_mem);

  arkFreeFPData(ark_mem);

  if (ark_mem->ark_lfree != NULL) 
    ark_mem->ark_lfree(ark_mem);

  if (ark_mem->ark_mfree != NULL) 
    ark_mem->ark_mfree(ark_mem);

  if (ark_mem->ark_nrtfn > 0) {
    free(ark_mem->ark_glo);     ark_mem->ark_glo     = NULL;
    free(ark_mem->ark_ghi);     ark_mem->ark_ghi     = NULL;
    free(ark_mem->ark_grout);   ark_mem->ark_grout   = NULL;
    free(ark_mem->ark_iroots);  ark_mem->ark_iroots  = NULL;
    free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
    free(ark_mem->ark_gactive); ark_mem->ark_gactive = NULL;
  }

  free(*arkode_mem);
  *arkode_mem = NULL;
}


/*---------------------------------------------------------------
 ARKodeRootInit:

 ARKodeRootInit initializes a rootfinding problem to be solved
 during the integration of the ODE system.  It loads the root
 function pointer and the number of root functions, and allocates
 workspace memory.  The return value is ARK_SUCCESS = 0 if no 
 errors occurred, or a negative value otherwise.
---------------------------------------------------------------*/
int ARKodeRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  ARKodeMem ark_mem;
  int i, nrt;

  /* Check arkode_mem pointer */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  nrt = (nrtfn < 0) ? 0 : nrtfn;

  /* If rerunning ARKodeRootInit() with a different number of root
     functions (changing number of gfun components), then free
     currently held memory resources */
  if ((nrt != ark_mem->ark_nrtfn) && (ark_mem->ark_nrtfn > 0)) {
    free(ark_mem->ark_glo);     ark_mem->ark_glo     = NULL;
    free(ark_mem->ark_ghi);     ark_mem->ark_ghi     = NULL;
    free(ark_mem->ark_grout);   ark_mem->ark_grout   = NULL;
    free(ark_mem->ark_iroots);  ark_mem->ark_iroots  = NULL;
    free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
    free(ark_mem->ark_gactive); ark_mem->ark_gactive = NULL;

    ark_mem->ark_lrw -= 3 * (ark_mem->ark_nrtfn);
    ark_mem->ark_liw -= 3 * (ark_mem->ark_nrtfn);
  }

  /* If ARKodeRootInit() was called with nrtfn == 0, then set 
     ark_nrtfn to zero and ark_gfun to NULL before returning */
  if (nrt == 0) {
    ark_mem->ark_nrtfn = nrt;
    ark_mem->ark_gfun = NULL;
    return(ARK_SUCCESS);
  }

  /* If rerunning ARKodeRootInit() with the same number of root 
     functions (not changing number of gfun components), then 
     check if the root function argument has changed */
  /* If g != NULL then return as currently reserved memory 
     resources will suffice */
  if (nrt == ark_mem->ark_nrtfn) {
    if (g != ark_mem->ark_gfun) {
      if (g == NULL) {
        free(ark_mem->ark_glo);     ark_mem->ark_glo     = NULL;
        free(ark_mem->ark_ghi);     ark_mem->ark_ghi     = NULL;
        free(ark_mem->ark_grout);   ark_mem->ark_grout   = NULL;
        free(ark_mem->ark_iroots);  ark_mem->ark_iroots  = NULL;
        free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
        free(ark_mem->ark_gactive); ark_mem->ark_gactive = NULL;

        ark_mem->ark_lrw -= 3*nrt;
        ark_mem->ark_liw -= 3*nrt;

        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeRootInit", MSGARK_NULL_G);
        return(ARK_ILL_INPUT);
      }
      else {
        ark_mem->ark_gfun = g;
        return(ARK_SUCCESS);
      }
    }
    else return(ARK_SUCCESS);
  }

  /* Set variable values in ARKode memory block */
  ark_mem->ark_nrtfn = nrt;
  if (g == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeRootInit", MSGARK_NULL_G);
    return(ARK_ILL_INPUT);
  }
  else ark_mem->ark_gfun = g;

  /* Allocate necessary memory and return */
  ark_mem->ark_glo = NULL;
  ark_mem->ark_glo = (realtype *) malloc(nrt*sizeof(realtype));
  if (ark_mem->ark_glo == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_ghi = NULL;
  ark_mem->ark_ghi = (realtype *) malloc(nrt*sizeof(realtype));
  if (ark_mem->ark_ghi == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_grout = NULL;
  ark_mem->ark_grout = (realtype *) malloc(nrt*sizeof(realtype));
  if (ark_mem->ark_grout == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL;
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_iroots = NULL;
  ark_mem->ark_iroots = (int *) malloc(nrt*sizeof(int));
  if (ark_mem->ark_iroots == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL; 
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    free(ark_mem->ark_grout); ark_mem->ark_grout = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_rootdir = NULL;
  ark_mem->ark_rootdir = (int *) malloc(nrt*sizeof(int));
  if (ark_mem->ark_rootdir == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL; 
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    free(ark_mem->ark_grout); ark_mem->ark_grout = NULL;
    free(ark_mem->ark_iroots); ark_mem->ark_iroots = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_gactive = NULL;
  ark_mem->ark_gactive = (booleantype *) malloc(nrt*sizeof(booleantype));
  if (ark_mem->ark_gactive == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL; 
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    free(ark_mem->ark_grout); ark_mem->ark_grout = NULL;
    free(ark_mem->ark_iroots); ark_mem->ark_iroots = NULL;
    free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODES", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* Set default values for rootdir (both directions) */
  for(i=0; i<nrt; i++) ark_mem->ark_rootdir[i] = 0;

  /* Set default values for gactive (all active) */
  for(i=0; i<nrt; i++) ark_mem->ark_gactive[i] = TRUE;

  ark_mem->ark_lrw += 3*nrt;
  ark_mem->ark_liw += 3*nrt;

  return(ARK_SUCCESS);
}



/*===============================================================
  Internal functions that may be replaced by the user
===============================================================*/

/*---------------------------------------------------------------
 arkEwtSet

 This routine is responsible for setting the error weight vector ewt,
 according to tol_type, as follows:

 (1) ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol), i=0,...,neq-1
     if tol_type = ARK_SS
 (2) ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol[i]), i=0,...,neq-1
     if tol_type = ARK_SV

 arkEwtSet returns 0 if ewt is successfully set as above to a
 positive vector and -1 otherwise. In the latter case, ewt is
 considered undefined.

 All the real work is done in the routines arkEwtSetSS, arkEwtSetSV.
---------------------------------------------------------------*/
int arkEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  ARKodeMem ark_mem;
  int flag = 0;

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  switch(ark_mem->ark_itol) {
  case ARK_SS: 
    flag = arkEwtSetSS(ark_mem, ycur, weight);
    break;
  case ARK_SV: 
    flag = arkEwtSetSV(ark_mem, ycur, weight);
    break;
  }
  
  return(flag);
}


/*---------------------------------------------------------------
 arkRwtSet

 This routine is responsible for setting the residual weight 
 vector rwt, according to tol_type, as follows:

 (1) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol), i=0,...,neq-1
     if tol_type = ARK_SS
 (2) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol[i]), i=0,...,neq-1
     if tol_type = ARK_SV
 (3) unset if tol_type is any other value (occurs rwt=ewt)

 arkRwtSet returns 0 if rwt is successfully set as above to a
 positive vector and -1 otherwise. In the latter case, rwt is
 considered undefined.

 All the real work is done in the routines arkRwtSetSS, arkRwtSetSV.
---------------------------------------------------------------*/
int arkRwtSet(N_Vector y, N_Vector weight, void *data)
{
  ARKodeMem ark_mem;
  N_Vector My;
  int flag = 0;

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  /* return if rwt is just ewt */
  if (ark_mem->ark_rwt_is_ewt)  return(0);

  /* put M*y into ark_ftemp */
  My = ark_mem->ark_ftemp;
  if (ark_mem->ark_mass_matrix) {
    flag = ark_mem->ark_mtimes(y, My, ark_mem->ark_tn, 
			       ark_mem->ark_mtimes_data);
    ark_mem->ark_mass_mult++;
    if (flag != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
  } else {  /* this condition should not apply, but just in case */
    N_VScale(ONE, y, My);
  }

  /* call appropriate routine to fill rwt */
  switch(ark_mem->ark_ritol) {
  case ARK_SS: 
    flag = arkRwtSetSS(ark_mem, My, weight);
    break;
  case ARK_SV: 
    flag = arkRwtSetSV(ark_mem, My, weight);
    break;
  }
  
  return(flag);
}


/*---------------------------------------------------------------
 arkErrHandler is the default error handling function.
   It sends the error message to the stream pointed to by ark_errfp 
---------------------------------------------------------------*/
void arkErrHandler(int error_code, const char *module,
		   const char *function, char *msg, void *data)
{
  ARKodeMem ark_mem;
  char err_type[10];

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  if (error_code == ARK_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

#ifndef NO_FPRINTF_OUTPUT
  if (ark_mem->ark_errfp!=NULL) {
    fprintf(ark_mem->ark_errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(ark_mem->ark_errfp,"  %s\n\n",msg);
  }
#endif

  return;
}



/*===============================================================
   Private Helper Functions
===============================================================*/

/*---------------------------------------------------------------
 arkPrintMem:

 This routine outputs the ark_mem structure to stdout.
---------------------------------------------------------------*/
static void arkPrintMem(ARKodeMem ark_mem)
{
  int i, j;

  /* output integer quantities */
  printf("ark_itol = %i\n", ark_mem->ark_itol);
  printf("ark_ritol = %i\n", ark_mem->ark_ritol);
  printf("ark_q = %i\n", ark_mem->ark_q);
  printf("ark_p = %i\n", ark_mem->ark_p);
  printf("ark_istage = %i\n", ark_mem->ark_istage);
  printf("ark_stages = %i\n", ark_mem->ark_stages);
  printf("ark_dense_q = %i\n", ark_mem->ark_dense_q);
  printf("ark_mnewt = %i\n", ark_mem->ark_mnewt);
  printf("ark_hadapt_imethod = %i\n", ark_mem->ark_hadapt_imethod);
  printf("ark_maxcor = %i\n", ark_mem->ark_maxcor);
  printf("ark_mxhnil = %i\n", ark_mem->ark_mxhnil);
  printf("ark_maxnef = %i\n", ark_mem->ark_maxnef);
  printf("ark_maxncf = %i\n", ark_mem->ark_maxncf);
  printf("ark_small_nef = %i\n", ark_mem->ark_small_nef);
  printf("ark_msbp = %i\n", ark_mem->ark_msbp);
  printf("ark_predictor = %i\n", ark_mem->ark_predictor);
  printf("ark_nhnil = %i\n", ark_mem->ark_nhnil);
  printf("ark_lsolve_type = %i\n", ark_mem->ark_lsolve_type);
  printf("ark_msolve_type = %i\n", ark_mem->ark_msolve_type);
  printf("ark_nrtfn = %i\n", ark_mem->ark_nrtfn);
  if (ark_mem->ark_iroots != NULL) 
    for (i=0; i<ark_mem->ark_nrtfn; i++)
      printf("ark_iroots[%i] = %i\n", i, ark_mem->ark_iroots[i]);
  if (ark_mem->ark_rootdir != NULL) 
    for (i=0; i<ark_mem->ark_nrtfn; i++)
      printf("ark_rootdir[%i] = %i\n", i, ark_mem->ark_rootdir[i]);
  printf("ark_taskc = %i\n", ark_mem->ark_taskc);
  printf("ark_irfnd = %i\n", ark_mem->ark_irfnd);
  printf("ark_mxgnull = %i\n", ark_mem->ark_mxgnull);

  /* output long integer quantities */
  printf("ark_mxstep = %li\n", ark_mem->ark_mxstep);
  printf("ark_nst = %li\n", ark_mem->ark_nst);
  printf("ark_nst_acc = %li\n", ark_mem->ark_nst_acc);
  printf("ark_nst_exp = %li\n", ark_mem->ark_nst_exp);
  printf("ark_nst_attempts = %li\n", ark_mem->ark_nst_attempts);
  printf("ark_nfe = %li\n", ark_mem->ark_nfe);
  printf("ark_nfi = %li\n", ark_mem->ark_nfi);
  printf("ark_ncfn = %li\n", ark_mem->ark_ncfn);
  printf("ark_nmassfails = %li\n", ark_mem->ark_nmassfails);
  printf("ark_netf = %li\n", ark_mem->ark_netf);
  printf("ark_nni = %li\n", ark_mem->ark_nni);
  printf("ark_nsetups = %li\n", ark_mem->ark_nsetups);
  printf("ark_lrw1 = %li\n", ark_mem->ark_lrw1);
  printf("ark_liw1 = %li\n", ark_mem->ark_liw1);
  printf("ark_lrw = %li\n", ark_mem->ark_lrw);
  printf("ark_liw = %li\n", ark_mem->ark_liw);
  printf("ark_fp_m = %li\n", ark_mem->ark_fp_m);
  if (ark_mem->ark_fp_imap != NULL)
    for (i=0; i<ark_mem->ark_fp_m; i++)
      printf("ark_fp_imap[%i] = %li\n", i, ark_mem->ark_fp_imap[i]);
  printf("ark_mass_solves = %li\n", ark_mem->ark_mass_solves);
  printf("ark_mass_mult = %li\n", ark_mem->ark_mass_mult);
  printf("ark_nstlp = %li\n", ark_mem->ark_nstlp);
  printf("ark_nge = %li\n", ark_mem->ark_nge);

  /* output boolean quantities */
  printf("ark_user_efun = %i\n", ark_mem->ark_user_efun);
  printf("ark_user_linear = %i\n", ark_mem->ark_linear);
  printf("ark_user_linear_timedep = %i\n", ark_mem->ark_linear_timedep);
  printf("ark_user_explicit = %i\n", ark_mem->ark_explicit);
  printf("ark_user_implicit = %i\n", ark_mem->ark_implicit);
  printf("ark_tstopset = %i\n", ark_mem->ark_tstopset);
  printf("ark_hadapt_pq = %i\n", ark_mem->ark_hadapt_pq);
  printf("ark_report = %i\n", ark_mem->ark_report);
  printf("ark_use_fp = %i\n", ark_mem->ark_use_fp);
  printf("ark_mass_matrix = %i\n", ark_mem->ark_mass_matrix);
  printf("ark_jcur = %i\n", ark_mem->ark_jcur);
  printf("ark_setupNonNull = %i\n", ark_mem->ark_setupNonNull);
  printf("ark_MassSetupNonNull = %i\n", ark_mem->ark_MassSetupNonNull);
  printf("ark_VabstolMallocDone = %i\n", ark_mem->ark_VabstolMallocDone);
  printf("ark_MallocDone = %i\n", ark_mem->ark_MallocDone);
  printf("ark_resized = %i\n", ark_mem->ark_resized);
  printf("ark_firststage = %i\n", ark_mem->ark_firststage);
  if (ark_mem->ark_gactive != NULL)
    for (i=0; i<ark_mem->ark_nrtfn; i++)
      printf("ark_gactive[%i] = %i\n", i, ark_mem->ark_gactive[i]);

  /* output realtype quantities */
  printf("ark_uround = %.16g\n", ark_mem->ark_uround);
  printf("ark_reltol = %.16g\n", ark_mem->ark_reltol);
  printf("ark_Sabstol = %.16g\n", ark_mem->ark_Sabstol);
  printf("ark_tstop = %.16g\n", ark_mem->ark_tstop);
  printf("ark_Ae = \n");
  for (i=0; i<ARK_S_MAX; i++) {
    printf("    ");
    for (j=0; j<ARK_S_MAX; j++)
      printf("%.16g  ", ARK_A(ark_mem->ark_Ae,i,j));
    printf("\n");
  }
  printf("ark_Ai = \n");
  for (i=0; i<ARK_S_MAX; i++) {
    printf("    ");
    for (j=0; j<ARK_S_MAX; j++)
      printf("%.16g  ", ARK_A(ark_mem->ark_Ai,i,j));
    printf("\n");
  }
  printf("ark_c = ");
  for (i=0; i<ARK_S_MAX; i++) 
    printf("%.16g  ", ark_mem->ark_c[i]);
  printf("\n");
  printf("ark_b = ");
  for (i=0; i<ARK_S_MAX; i++) 
    printf("%.16g  ", ark_mem->ark_b[i]);
  printf("\n");
  printf("ark_b2 = ");
  for (i=0; i<ARK_S_MAX; i++) 
    printf("%.16g  ", ark_mem->ark_b2[i]);
  printf("\n");
  printf("ark_hin = %.16g\n", ark_mem->ark_hin);
  printf("ark_h = %.16g\n", ark_mem->ark_h);
  printf("ark_hprime = %.16g\n", ark_mem->ark_hprime);
  printf("ark_next_h = %.16g\n", ark_mem->ark_next_h);
  printf("ark_eta = %.16g\n", ark_mem->ark_eta);
  printf("ark_tn = %.16g\n", ark_mem->ark_tn);
  printf("ark_tretlast = %.16g\n", ark_mem->ark_tretlast);
  printf("ark_gamma = %.16g\n", ark_mem->ark_gamma);
  printf("ark_gammap = %.16g\n", ark_mem->ark_gammap);
  printf("ark_gamrat = %.16g\n", ark_mem->ark_gamrat);
  printf("ark_crate = %.16g\n", ark_mem->ark_crate);
  printf("ark_eLTE = %.16g\n", ark_mem->ark_eLTE);
  printf("ark_nlscoef = %.16g\n", ark_mem->ark_nlscoef);
  printf("ark_fixedstep = %i\n", ark_mem->ark_fixedstep);
  printf("ark_hadapt_ehist =  %.16g  %.16g  %.16g\n",
	 ark_mem->ark_hadapt_ehist[0], ark_mem->ark_hadapt_ehist[1], ark_mem->ark_hadapt_ehist[2]);
  printf("ark_hadapt_hhist =  %.16g  %.16g  %.16g\n",
	 ark_mem->ark_hadapt_hhist[0], ark_mem->ark_hadapt_hhist[1], ark_mem->ark_hadapt_hhist[2]);
  printf("ark_hadapt_cfl = %.16g\n", ark_mem->ark_hadapt_cfl);
  printf("ark_hadapt_safety = %.16g\n", ark_mem->ark_hadapt_safety);
  printf("ark_hadapt_bias = %.16g\n", ark_mem->ark_hadapt_bias);
  printf("ark_hadapt_growth = %.16g\n", ark_mem->ark_hadapt_growth);
  printf("ark_hadapt_lbound = %.16g\n", ark_mem->ark_hadapt_lbound);
  printf("ark_hadapt_ubound = %.16g\n", ark_mem->ark_hadapt_ubound);
  printf("ark_hadapt_k1 = %.16g\n", ark_mem->ark_hadapt_k1);
  printf("ark_hadapt_k2 = %.16g\n", ark_mem->ark_hadapt_k2);
  printf("ark_hadapt_k3 = %.16g\n", ark_mem->ark_hadapt_k3);
  printf("ark_hmin = %.16g\n", ark_mem->ark_hmin);
  printf("ark_hmax_inv = %.16g\n", ark_mem->ark_hmax_inv);
  printf("ark_etamax = %.16g\n", ark_mem->ark_etamax);
  printf("ark_etamx1 = %.16g\n", ark_mem->ark_etamx1);
  printf("ark_etamxf = %.16g\n", ark_mem->ark_etamxf);
  printf("ark_etacf = %.16g\n", ark_mem->ark_etacf);
  printf("ark_crdown = %.16g\n", ark_mem->ark_crdown);
  printf("ark_rdiv = %.16g\n", ark_mem->ark_rdiv);
  printf("ark_dgmax = %.16g\n", ark_mem->ark_dgmax);
  if (ark_mem->ark_fp_R != NULL) {
    printf("ark_fp_R =  ");
    for (i=0; i<ark_mem->ark_fp_m*ark_mem->ark_fp_m; i++)
      printf("%.16g  ", ark_mem->ark_fp_R[i]);
    printf("\n");
  }
  if (ark_mem->ark_fp_gamma != NULL) {
    printf("ark_fp_gamma =  ");
    for (i=0; i<ark_mem->ark_fp_m; i++)
      printf("%.16g  ", ark_mem->ark_fp_gamma[i]);
    printf("\n");
  }
  printf("ark_h0u = %.16g\n", ark_mem->ark_h0u);
  printf("ark_tnew = %.16g\n", ark_mem->ark_tnew);
  printf("ark_hold = %.16g\n", ark_mem->ark_hold);
  printf("ark_tolsf = %.16g\n", ark_mem->ark_tolsf);
  printf("ark_tlo = %.16g\n", ark_mem->ark_tlo);
  printf("ark_thi = %.16g\n", ark_mem->ark_thi);
  printf("ark_trout = %.16g\n", ark_mem->ark_trout);
  if (ark_mem->ark_glo != NULL) 
    for (i=0; i<ark_mem->ark_nrtfn; i++)
      printf("ark_glo[%i] = %.16g\n", i, ark_mem->ark_glo[i]);
  if (ark_mem->ark_ghi != NULL) 
    for (i=0; i<ark_mem->ark_nrtfn; i++)
      printf("ark_ghi[%i] = %.16g\n", i, ark_mem->ark_ghi[i]);
  if (ark_mem->ark_grout != NULL) 
    for (i=0; i<ark_mem->ark_nrtfn; i++)
      printf("ark_grout[%i] = %.16g\n", i, ark_mem->ark_grout[i]);
  printf("ark_toutc = %.16g\n", ark_mem->ark_toutc);
  printf("ark_ttol = %.16g\n", ark_mem->ark_ttol);

#ifdef DEBUG_OUTPUT
  /* output vector quantities */  
  if (ark_mem->ark_Vabstol != NULL) {
    printf("ark_Vapbsol:\n");
    N_VPrint_Serial(ark_mem->ark_Vabstol);
  }
  for (i=0; i<ARK_S_MAX; i++) 
    if (ark_mem->ark_Fe[i] != NULL) {
      printf("ark_Fe[%i]:\n", i);
      N_VPrint_Serial(ark_mem->ark_Fe[i]);
    }
  for (i=0; i<ARK_S_MAX; i++) 
    if (ark_mem->ark_Fi[i] != NULL) {
      printf("ark_Fi[%i]:\n", i);
      N_VPrint_Serial(ark_mem->ark_Fi[i]);
    }
  if (ark_mem->ark_ewt != NULL) {
    printf("ark_ewt:\n");
    N_VPrint_Serial(ark_mem->ark_ewt);
  }
  if (!ark_mem->ark_rwt_is_ewt && ark_mem->ark_rwt != NULL) {
    printf("ark_rwt:\n");
    N_VPrint_Serial(ark_mem->ark_rwt);
  }
  if (ark_mem->ark_y != NULL) {
    printf("ark_y:\n");
    N_VPrint_Serial(ark_mem->ark_y);
  }
  if (ark_mem->ark_ycur != NULL) {
    printf("ark_ycur:\n");
    N_VPrint_Serial(ark_mem->ark_ycur);
  }
  if (ark_mem->ark_sdata != NULL) {
    printf("ark_sdata:\n");
    N_VPrint_Serial(ark_mem->ark_sdata);
  }
  if (ark_mem->ark_tempv != NULL) {
    printf("ark_tempv:\n");
    N_VPrint_Serial(ark_mem->ark_tempv);
  }
  if (ark_mem->ark_acor != NULL) {
    printf("ark_acor:\n");
    N_VPrint_Serial(ark_mem->ark_acor);
  }
  if (ark_mem->ark_ftemp != NULL) {
    printf("ark_ftemp:\n");
    N_VPrint_Serial(ark_mem->ark_ftemp);
  }
  if (ark_mem->ark_fold != NULL) {
    printf("ark_fold:\n");
    N_VPrint_Serial(ark_mem->ark_fold);
  }
  if (ark_mem->ark_fnew != NULL) {
    printf("ark_fnew:\n");
    N_VPrint_Serial(ark_mem->ark_fnew);
  }
  if (ark_mem->ark_yold != NULL) {
    printf("ark_yold:\n");
    N_VPrint_Serial(ark_mem->ark_yold);
  }
  if (ark_mem->ark_ynew != NULL) {
    printf("ark_ynew:\n");
    N_VPrint_Serial(ark_mem->ark_ynew);
  }
  if (ark_mem->ark_fp_df != NULL) 
    for (i=0; i<ark_mem->ark_fp_m; i++)
      if (ark_mem->ark_fp_df[i] != NULL) {
	printf("ark_fp_df[%i]:\n", i);
	N_VPrint_Serial(ark_mem->ark_fp_df[i]);
      }
  if (ark_mem->ark_fp_dg != NULL) 
    for (i=0; i<ark_mem->ark_fp_m; i++)
      if (ark_mem->ark_fp_dg[i] != NULL) {
	printf("ark_fp_dg[%i]:\n", i);
	N_VPrint_Serial(ark_mem->ark_fp_dg[i]);
      }
  if (ark_mem->ark_fp_q != NULL) 
    for (i=0; i<ark_mem->ark_fp_m; i++)
      if (ark_mem->ark_fp_q[i] != NULL) {
	printf("ark_fp_q[%i]:\n", i);
	N_VPrint_Serial(ark_mem->ark_fp_q[i]);
      }
  if (ark_mem->ark_fp_fval != NULL) {
    printf("ark_fp_fval:\n");
    N_VPrint_Serial(ark_mem->ark_fp_fval);
  }
  if (ark_mem->ark_fp_fold != NULL) {
    printf("ark_fp_fold:\n");
    N_VPrint_Serial(ark_mem->ark_fp_fold);
  }
  if (ark_mem->ark_fp_gold != NULL) {
    printf("ark_fp_gold:\n");
    N_VPrint_Serial(ark_mem->ark_fp_gold);
  }
#endif

}


/*---------------------------------------------------------------
 arkCheckNvector:

 This routine checks if all required vector operations are 
 present.  If any of them is missing it returns FALSE.
---------------------------------------------------------------*/
static booleantype arkCheckNvector(N_Vector tmpl)  /* to be updated?? */
{
  if ((tmpl->ops->nvclone     == NULL) ||
      (tmpl->ops->nvdestroy   == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) ||
      (tmpl->ops->nvconst     == NULL) ||
      (tmpl->ops->nvdiv       == NULL) ||
      (tmpl->ops->nvscale     == NULL) ||
      (tmpl->ops->nvabs       == NULL) ||
      (tmpl->ops->nvinv       == NULL) ||
      (tmpl->ops->nvaddconst  == NULL) ||
      (tmpl->ops->nvmaxnorm   == NULL) ||
      (tmpl->ops->nvwrmsnorm  == NULL) ||
      (tmpl->ops->nvmin       == NULL))
    return(FALSE);
  else
    return(TRUE);
}


/*---------------------------------------------------------------
 arkAllocVectors:

 This routine allocates the ARKODE vectors ewt, acor, ycur, 
 sdata, tempv, ftemp, fold, fnew, yold, and ynew.  If any of 
 these vectors already exist, they are left alone.  Otherwise, it
 will allocate each vector by cloning the input vector. If all 
 memory allocations are successful, arkAllocVectors returns TRUE. 
 Otherwise all vector memory is freed and arkAllocVectors 
 returns FALSE. This routine also updates the optional outputs 
 lrw and liw, which are (respectively) the lengths of the real 
 and integer work spaces.
---------------------------------------------------------------*/
static booleantype arkAllocVectors(ARKodeMem ark_mem, N_Vector tmpl)
{
  /* Allocate ewt if needed */
  if (ark_mem->ark_ewt == NULL) {
    ark_mem->ark_ewt = N_VClone(tmpl);
    if (ark_mem->ark_ewt == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Set rwt to point at ewt */
  if (ark_mem->ark_rwt_is_ewt) 
    ark_mem->ark_rwt = ark_mem->ark_ewt;

  /* Allocate acor if needed */
  if (ark_mem->ark_acor == NULL) {
    ark_mem->ark_acor = N_VClone(tmpl);
    if (ark_mem->ark_acor == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ycur if needed */
  if (ark_mem->ark_ycur == NULL) {
    ark_mem->ark_ycur = N_VClone(tmpl);
    if (ark_mem->ark_ycur == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate sdata if needed */
  if (ark_mem->ark_sdata == NULL) {
    ark_mem->ark_sdata = N_VClone(tmpl);
    if (ark_mem->ark_sdata == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate tempv if needed */
  if (ark_mem->ark_tempv == NULL) {
    ark_mem->ark_tempv = N_VClone(tmpl);
    if (ark_mem->ark_tempv == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ftemp if needed */
  if (ark_mem->ark_ftemp == NULL) {
    ark_mem->ark_ftemp = N_VClone(tmpl);
    if (ark_mem->ark_ftemp == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate fold if needed */
  if (ark_mem->ark_fold == NULL) {
    ark_mem->ark_fold = N_VClone(tmpl);
    if (ark_mem->ark_fold == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate fnew if needed */
  if (ark_mem->ark_fnew == NULL) {
    ark_mem->ark_fnew = N_VClone(tmpl);
    if (ark_mem->ark_fnew == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate yold if needed */
  if (ark_mem->ark_yold == NULL) {
    ark_mem->ark_yold = N_VClone(tmpl);
    if (ark_mem->ark_yold == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ynew if needed */
  if (ark_mem->ark_ynew == NULL) {
    ark_mem->ark_ynew = N_VClone(tmpl);
    if (ark_mem->ark_ynew == NULL) {
      arkFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  return(TRUE);
}


/*---------------------------------------------------------------
 arkAllocRKVectors

 This routine allocates the ARKODE vectors 
 Fe[0], ..., Fe[stages-1], Fi[0], ..., Fi[stages-1], and fa and 
 fb (if necessary).

 If any of these vectors already exist, they are left alone.  
 Otherwise, it will allocate each vector by cloning the input 
 vector. If all memory allocations are successful, 
 arkAllocRKVectors returns TRUE.  Otherwise all vector memory is 
 freed and arkAllocRKVectors returns FALSE. This routine also 
 updates the optional outputs lrw and liw, which are 
 (respectively) the lengths of the real and integer work spaces.

 NOTE: this routine is separate from arkAllocVectors, since that
 routine is called in ARKodeInit, which may be called by a user
 prior to setting their desired method order, 
 implicit/explicit/imex problem, and number of RK stages.  Hence
 this routine is called during arkInitialStep, which is called 
 on the first call to the ARKode() integrator.
---------------------------------------------------------------*/
static int arkAllocRKVectors(ARKodeMem ark_mem)
{
  int j;

  /* Allocate Fe[0] ... Fe[stages-1] if needed */
  if (!ark_mem->ark_implicit) {
    for (j=0; j<ark_mem->ark_stages; j++) {
      if (ark_mem->ark_Fe[j] == NULL) {
	ark_mem->ark_Fe[j] = N_VClone(ark_mem->ark_ewt);
	if (ark_mem->ark_Fe[j] == NULL) {
	  arkFreeVectors(ark_mem);
	  return(ARK_MEM_FAIL);
	} else {
	  ark_mem->ark_lrw += ark_mem->ark_lrw1;
	  ark_mem->ark_liw += ark_mem->ark_liw1;
	}
      }
    }
  }

  /* Allocate Fi[0] ... Fi[stages-1] if needed */
  if (!ark_mem->ark_explicit) {
    for (j=0; j<ark_mem->ark_stages; j++) {
      if (ark_mem->ark_Fi[j] == NULL) {
	ark_mem->ark_Fi[j] = N_VClone(ark_mem->ark_ewt);
	if (ark_mem->ark_Fi[j] == NULL) {
	  arkFreeVectors(ark_mem);
	  return(ARK_MEM_FAIL);
	} else {
	  ark_mem->ark_lrw += ark_mem->ark_lrw1;
	  ark_mem->ark_liw += ark_mem->ark_liw1;
	}
      }
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkFreeVectors

 This routine frees the ARKODE vectors allocated in both
 arkAllocVectors and arkAllocRKVectors.
---------------------------------------------------------------*/
static void arkFreeVectors(ARKodeMem ark_mem)
{
  int j;
  if (ark_mem->ark_ewt != NULL) {
    N_VDestroy(ark_mem->ark_ewt);
    ark_mem->ark_ewt = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if ((!ark_mem->ark_rwt_is_ewt) && (ark_mem->ark_rwt != NULL)) {
    N_VDestroy(ark_mem->ark_rwt);
    ark_mem->ark_rwt = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_acor != NULL) {
    N_VDestroy(ark_mem->ark_acor);
    ark_mem->ark_acor = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_ycur != NULL) {
    N_VDestroy(ark_mem->ark_ycur);
    ark_mem->ark_ycur = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_sdata != NULL) {
    N_VDestroy(ark_mem->ark_sdata);
    ark_mem->ark_sdata = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_tempv != NULL) {
    N_VDestroy(ark_mem->ark_tempv);
    ark_mem->ark_tempv = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_ftemp != NULL) {
    N_VDestroy(ark_mem->ark_ftemp);
    ark_mem->ark_ftemp = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_fold != NULL) {
    N_VDestroy(ark_mem->ark_fold);
    ark_mem->ark_fold = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_fnew != NULL) {
    N_VDestroy(ark_mem->ark_fnew);
    ark_mem->ark_fnew = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_yold != NULL) {
    N_VDestroy(ark_mem->ark_yold);
    ark_mem->ark_yold = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_ynew != NULL) {
    N_VDestroy(ark_mem->ark_ynew);
    ark_mem->ark_ynew = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  /* if (ark_mem->ark_fa != NULL) { */
  /*   N_VDestroy(ark_mem->ark_fa); */
  /*   ark_mem->ark_fa = NULL; */
  /*   ark_mem->ark_lrw -= ark_mem->ark_lrw1; */
  /*   ark_mem->ark_liw -= ark_mem->ark_liw1; */
  /* } */
  /* if (ark_mem->ark_fb != NULL) { */
  /*   N_VDestroy(ark_mem->ark_fb); */
  /*   ark_mem->ark_fb = NULL; */
  /*   ark_mem->ark_lrw -= ark_mem->ark_lrw1; */
  /*   ark_mem->ark_liw -= ark_mem->ark_liw1; */
  /* } */
  for(j=0; j<ARK_S_MAX; j++) {
    if (ark_mem->ark_Fe[j] != NULL) {
      N_VDestroy(ark_mem->ark_Fe[j]);
      ark_mem->ark_Fe[j] = NULL;
      ark_mem->ark_lrw -= ark_mem->ark_lrw1;
      ark_mem->ark_liw -= ark_mem->ark_liw1;
    }
    if (ark_mem->ark_Fi[j] != NULL) {
      N_VDestroy(ark_mem->ark_Fi[j]);
      ark_mem->ark_Fi[j] = NULL;
      ark_mem->ark_lrw -= ark_mem->ark_lrw1;
      ark_mem->ark_liw -= ark_mem->ark_liw1;
    }
  }
  if (ark_mem->ark_Vabstol != NULL) {
    N_VDestroy(ark_mem->ark_Vabstol);
    ark_mem->ark_Vabstol = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
}


/*---------------------------------------------------------------
 arkAllocFPData

 This routine allocates all required memory for performing the 
 accelerated fixed-point solver.
---------------------------------------------------------------*/
static int arkAllocFPData(ARKodeMem ark_mem)
{
  long int maa = ark_mem->ark_fp_m;

  /* ensure that DotProd function is defined */
  if (ark_mem->ark_ewt->ops->nvdotprod == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "arkAllocFPData", MSGARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Allocate ark_fp_fval if needed */
  if (ark_mem->ark_fp_fval == NULL) {
    ark_mem->ark_fp_fval = N_VClone(ark_mem->ark_ewt);
    if (ark_mem->ark_fp_fval == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ark_fp_fold if needed */
  if (ark_mem->ark_fp_fold == NULL) {
    ark_mem->ark_fp_fold = N_VClone(ark_mem->ark_ewt);
    if (ark_mem->ark_fp_fold == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ark_fp_gold if needed */
  if (ark_mem->ark_fp_gold == NULL) {
    ark_mem->ark_fp_gold = N_VClone(ark_mem->ark_ewt);
    if (ark_mem->ark_fp_gold == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ark_fp_df if needed */
  if ((ark_mem->ark_fp_df == NULL) && (maa > 0)) {
    ark_mem->ark_fp_df = N_VCloneVectorArray(maa, ark_mem->ark_ewt);
    if (ark_mem->ark_fp_df == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += maa*ark_mem->ark_lrw1;
      ark_mem->ark_liw += maa*ark_mem->ark_liw1;
    }
  }

  /* Allocate ark_fp_dg if needed */
  if ((ark_mem->ark_fp_dg == NULL) && (maa > 0)) {
    ark_mem->ark_fp_dg = N_VCloneVectorArray(maa, ark_mem->ark_ewt);
    if (ark_mem->ark_fp_dg == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += maa*ark_mem->ark_lrw1;
      ark_mem->ark_liw += maa*ark_mem->ark_liw1;
    }
  }

  /* Allocate ark_fp_q if needed */
  if ((ark_mem->ark_fp_q == NULL) && (maa > 0)) {
    ark_mem->ark_fp_q = N_VCloneVectorArray(maa, ark_mem->ark_ewt);
    if (ark_mem->ark_fp_q == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += maa*ark_mem->ark_lrw1;
      ark_mem->ark_liw += maa*ark_mem->ark_liw1;
    }
  }

  /* Allocate ark_fp_R if needed */
  if ((ark_mem->ark_fp_R == NULL) && (maa > 0)) {
    ark_mem->ark_fp_R = (realtype *) malloc((maa*maa) * sizeof(realtype));
    if (ark_mem->ark_fp_R == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += maa*maa;
    }
  }

  /* Allocate ark_fp_gamma if needed */
  if ((ark_mem->ark_fp_gamma == NULL) && (maa > 0)) {
    ark_mem->ark_fp_gamma = (realtype *) malloc(maa * sizeof(realtype));
    if (ark_mem->ark_fp_gamma == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_lrw += maa;
    }
  }

  /* Allocate ark_fp_imap if needed */
  if ((ark_mem->ark_fp_imap == NULL) && (maa > 0)) {
    ark_mem->ark_fp_imap = (long int *) malloc(maa * sizeof(long int));
    if (ark_mem->ark_fp_imap == NULL) {
      arkFreeFPData(ark_mem);
      return(ARK_MEM_FAIL);
    } else {
      ark_mem->ark_liw += maa;
    }
  }

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
 arkResizeFPData

 This routine resizes all required memory for the accelerated
 fixed-point solver (called from ARKodeResize()).
---------------------------------------------------------------*/
static int arkResizeFPData(ARKodeMem ark_mem, ARKVecResizeFn resize,
			   void *resize_data, long int lrw_diff, 
			   long int liw_diff)
{
  long int i;
  long int maa = ark_mem->ark_fp_m;

  /* Resize ark_fp_fval if needed */
  if (ark_mem->ark_fp_fval != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_fp_fval);
      ark_mem->ark_fp_fval = N_VClone(ark_mem->ark_ewt);
    } else {
      if (resize(ark_mem->ark_fp_fval, ark_mem->ark_ewt, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"arkResizeFPData", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }

  /* Resize ark_fp_fold if needed */
  if (ark_mem->ark_fp_fold != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_fp_fold);
      ark_mem->ark_fp_fold = N_VClone(ark_mem->ark_ewt);
    } else {
      if (resize(ark_mem->ark_fp_fold, ark_mem->ark_ewt, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"arkResizeFPData", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }

  /* Resize ark_fp_gold if needed */
  if (ark_mem->ark_fp_gold != NULL) {
    if (resize == NULL) {
      N_VDestroy(ark_mem->ark_fp_gold);
      ark_mem->ark_fp_gold = N_VClone(ark_mem->ark_ewt);
    } else {
      if (resize(ark_mem->ark_fp_gold, ark_mem->ark_ewt, resize_data)) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"arkResizeFPData", MSGARK_RESIZE_FAIL);
	return(ARK_ILL_INPUT);
      }
    }
    ark_mem->ark_lrw += lrw_diff;
    ark_mem->ark_liw += liw_diff;
  }

  /* Resize ark_fp_df if needed */
  if ((ark_mem->ark_fp_df != NULL) && (maa > 0)) {
    for (i=0; i<maa; i++) {
      if (resize == NULL) {
	N_VDestroy(ark_mem->ark_fp_df[i]);
	ark_mem->ark_fp_df[i] = N_VClone(ark_mem->ark_ewt);
      } else {
	if (resize(ark_mem->ark_fp_df[i], ark_mem->ark_ewt, resize_data)) {
	  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			  "arkResizeFPData", MSGARK_RESIZE_FAIL);
	  return(ARK_ILL_INPUT);
	}
      }
    }
    ark_mem->ark_lrw += maa*lrw_diff;
    ark_mem->ark_liw += maa*liw_diff;
  }

  /* Resize ark_fp_dg if needed */
  if ((ark_mem->ark_fp_dg != NULL) && (maa > 0)) {
    for (i=0; i<maa; i++) {
      if (resize == NULL) {
	N_VDestroy(ark_mem->ark_fp_dg[i]);
	ark_mem->ark_fp_dg[i] = N_VClone(ark_mem->ark_ewt);
      } else {
	if (resize(ark_mem->ark_fp_dg[i], ark_mem->ark_ewt, resize_data)) {
	  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			  "arkResizeFPData", MSGARK_RESIZE_FAIL);
	  return(ARK_ILL_INPUT);
	}
      }
    }
    ark_mem->ark_lrw += maa*lrw_diff;
    ark_mem->ark_liw += maa*liw_diff;
  }

  /* Resize ark_fp_q if needed */
  if ((ark_mem->ark_fp_q != NULL) && (maa > 0)) {
    for (i=0; i<maa; i++) {
      if (resize == NULL) {
	N_VDestroy(ark_mem->ark_fp_q[i]);
	ark_mem->ark_fp_q[i] = N_VClone(ark_mem->ark_ewt);
      } else {
	if (resize(ark_mem->ark_fp_q[i], ark_mem->ark_ewt, resize_data)) {
	  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			  "arkResizeFPData", MSGARK_RESIZE_FAIL);
	  return(ARK_ILL_INPUT);
	}
      }
    }
    ark_mem->ark_lrw += maa*lrw_diff;
    ark_mem->ark_liw += maa*liw_diff;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkFreeFPData

 This routine frees all memory allocated by arkAllocFPData.
---------------------------------------------------------------*/
static void arkFreeFPData(ARKodeMem ark_mem)
{
  long int maa = ark_mem->ark_fp_m;

  /* free ark_fp_fval if needed */
  if (ark_mem->ark_fp_fval != NULL) {
    N_VDestroy(ark_mem->ark_fp_fval);
    ark_mem->ark_fp_fval = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }

  /* free ark_fp_fold if needed */
  if (ark_mem->ark_fp_fold != NULL) {
    N_VDestroy(ark_mem->ark_fp_fold);
    ark_mem->ark_fp_fold = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }

  /* free ark_fp_gold if needed */
  if (ark_mem->ark_fp_gold != NULL) {
    N_VDestroy(ark_mem->ark_fp_gold);
    ark_mem->ark_fp_gold = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }

  /* free ark_fp_df if needed */
  if ((ark_mem->ark_fp_df != NULL) && (maa>0)) {
    N_VDestroyVectorArray(ark_mem->ark_fp_df, maa);
    ark_mem->ark_fp_df = NULL;
    ark_mem->ark_lrw -= maa*ark_mem->ark_lrw1;
    ark_mem->ark_liw -= maa*ark_mem->ark_liw1;
  }

  /* free ark_fp_dg if needed */
  if ((ark_mem->ark_fp_dg != NULL) && (maa>0)) {
    N_VDestroyVectorArray(ark_mem->ark_fp_dg, maa);
    ark_mem->ark_fp_dg = NULL;
    ark_mem->ark_lrw -= maa*ark_mem->ark_lrw1;
    ark_mem->ark_liw -= maa*ark_mem->ark_liw1;
  }

  /* free ark_fp_q if needed */
  if ((ark_mem->ark_fp_q != NULL) && (maa>0)) {
    N_VDestroyVectorArray(ark_mem->ark_fp_q, maa);
    ark_mem->ark_fp_q = NULL;
    ark_mem->ark_lrw -= maa*ark_mem->ark_lrw1;
    ark_mem->ark_liw -= maa*ark_mem->ark_liw1;
  }

  /* free ark_fp_R if needed */
  if (ark_mem->ark_fp_R != NULL) {
    free(ark_mem->ark_fp_R);
    ark_mem->ark_fp_R = NULL;
    ark_mem->ark_lrw -= maa*maa;
  }

  /* free ark_fp_gamma if needed */
  if (ark_mem->ark_fp_gamma != NULL) {
    free(ark_mem->ark_fp_gamma);
    ark_mem->ark_fp_gamma = NULL;
    ark_mem->ark_lrw -= maa;
  }

  /* free ark_fp_imap if needed */
  if (ark_mem->ark_fp_imap != NULL) {
    free(ark_mem->ark_fp_imap);
    ark_mem->ark_fp_imap = NULL;
    ark_mem->ark_liw -= maa;
  }
}


/*---------------------------------------------------------------
 arkInitialSetup

 This routine performs input consistency checks at the first step.
 If needed, it also checks the linear solver module and calls the
 linear solver initialization routine.
---------------------------------------------------------------*/
static int arkInitialSetup(ARKodeMem ark_mem)
{
  int ier;

  /* Set first step growth factor */
  ark_mem->ark_etamax = ark_mem->ark_etamx1;

  /* Create Butcher tables (if not already set) */
  ier = arkSetButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "arkInitialStep", "Could not create Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher tables are OK */
  ier = arkCheckButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "arkInitialStep", "Error in Butcher tables");
    return(ARK_ILL_INPUT);
  }

  /* Allocate ARK RHS vector memory */
  ier = arkAllocRKVectors(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "arkInitialSetup", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* Check if lsolve function exists and call linit (if it exists) */
  if (!ark_mem->ark_explicit && !ark_mem->ark_use_fp) {
    if (ark_mem->ark_lsolve == NULL) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "arkInitialSetup", MSGARK_LSOLVE_NULL);
      return(ARK_ILL_INPUT);
    }
    if (ark_mem->ark_linit != NULL) {
      ier = ark_mem->ark_linit(ark_mem);
      if (ier != 0) {
	arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE", 
			"arkInitialSetup", MSGARK_LINIT_FAIL);
	return(ARK_LINIT_FAIL);
      }
    }
  }

  /* Call minit (if it exists) */
  if (ark_mem->ark_minit != NULL) {
    ier = ark_mem->ark_minit(ark_mem);
    if (ier != 0) {
      arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKODE", 
		      "arkInitialSetup", MSGARK_MASSINIT_FAIL);
      return(ARK_MASSINIT_FAIL);
    }
  }
  
  /* Call msetup (it necessary) */
  if (ark_mem->ark_mass_matrix && ark_mem->ark_MassSetupNonNull) {
    ier = ark_mem->ark_msetup(ark_mem, ark_mem->ark_ewt, 
			      ark_mem->ark_acor, ark_mem->ark_sdata);
    if (ier != 0) {
      arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKODE", 
		      "arkInitialSetup", MSGARK_MASSSETUP_FAIL);
      return(ARK_MASSSETUP_FAIL);
    }
  }

  /* Check for consistency between linear system modules 
     (if lsolve is direct, msolve needs to match) */
  if (ark_mem->ark_mass_matrix) {  /* M != I */
    if (ark_mem->ark_lsolve_type != 0) { /* non-iterative */
      if (ark_mem->ark_lsolve_type != ark_mem->ark_msolve_type) {
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"arkInitialSetup", 
			"Incompatible linear and mass matrix solvers");
	return(ARK_MEM_FAIL);
      }
    }
  }

  /* Set data for efun (if left unspecified) */
  if (ark_mem->ark_user_efun) 
    ark_mem->ark_e_data = ark_mem->ark_user_data;
  else                        
    ark_mem->ark_e_data = ark_mem;

  /* Load initial error weights */
  ier = ark_mem->ark_efun(ark_mem->ark_ycur,
			  ark_mem->ark_ewt, 
			  ark_mem->ark_e_data);
  if (ier != 0) {
    if (ark_mem->ark_itol == ARK_WF) 
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "arkInitialSetup", MSGARK_EWT_FAIL);
    else 
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "arkInitialSetup", MSGARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }

  /* Set data for rfun (if left unspecified) */
  if (ark_mem->ark_user_rfun) 
    ark_mem->ark_r_data = ark_mem->ark_user_data;
  else                        
    ark_mem->ark_r_data = ark_mem;

  /* Load initial residual weights */
  if (ark_mem->ark_rwt_is_ewt) {      /* update pointer to ewt */
    ark_mem->ark_rwt = ark_mem->ark_ewt;
  } else { 
    ier = ark_mem->ark_rfun(ark_mem->ark_ycur,
			    ark_mem->ark_rwt, 
			    ark_mem->ark_r_data);
    if (ier != 0) {
      if (ark_mem->ark_itol == ARK_WF) 
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"arkInitialSetup", MSGARK_RWT_FAIL);
      else 
	arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"arkInitialSetup", MSGARK_BAD_RWT);
      return(ARK_ILL_INPUT);
    }
  }

  /* Allocate fixed-point solver memory */
  if (ark_mem->ark_use_fp) {
    ier = arkAllocFPData(ark_mem);
    if (ier != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		      "arkInitialSetup", MSGARK_MEM_FAIL);
      return(ARK_MEM_FAIL);
    }
  }

  /* Fill initial ynew and fnew arrays */
  N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_ynew);
  ier = arkFullRHS(ark_mem, ark_mem->ark_tn, ark_mem->ark_ycur,
		   ark_mem->ark_ftemp, ark_mem->ark_fnew);
  if (ier != 0)  return(ARK_RHSFUNC_FAIL);

  /* if the problem involves a non-identity mass matrix, update fnew here */
  if (ark_mem->ark_mass_matrix) {
    N_VScale(ark_mem->ark_h, ark_mem->ark_fnew, ark_mem->ark_fnew);   /* scale RHS */
    ier = ark_mem->ark_msolve(ark_mem, ark_mem->ark_fnew, ark_mem->ark_rwt); 
    N_VScale(ONE/ark_mem->ark_h, ark_mem->ark_fnew, ark_mem->ark_fnew);   /* scale result */
    ark_mem->ark_mass_solves++;
    if (ier != ARK_SUCCESS) {
      ark_mem->ark_nmassfails++;
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE", 
		      "arkInitialSetup", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkHin

 This routine computes a tentative initial step size h0. 
 If tout is too close to tn (= t0), then arkHin returns 
 ARK_TOO_CLOSE and h remains uninitialized. Note that here tout 
 is either the value passed to ARKode at the first call or the 
 value of tstop (if tstop is enabled and it is closer to t0=tn 
 than tout). If the RHS function fails unrecoverably, arkHin 
 returns ARK_RHSFUNC_FAIL. If the RHS function fails recoverably 
 too many times and recovery is not possible, arkHin returns 
 ARK_REPTD_RHSFUNC_ERR. Otherwise, arkHin sets h to the chosen 
 value h0 and returns ARK_SUCCESS.

 The algorithm used seeks to find h0 as a solution of
       (WRMS norm of (h0^2 ydd / 2)) = 1, 
 where ydd = estimated second derivative of y.  Although this 
 choice is based on an error expansion of the Backward Euler
 method, and hence results in an overly-conservative time step 
 for our higher-order ARK methods, it does find an order-of-
 magnitude estimate of the initial time scale of the solution.
 Since this method is only used on the first time step, the
 additional caution will not overly hinder solver efficiency.

 We start with an initial estimate equal to the geometric mean 
 of the lower and upper bounds on the step size.

 Loop up to H0_ITERS times to find h0.
 Stop if new and previous values differ by a factor < 2.
 Stop if hnew/hg > 2 after one iteration, as this probably 
 means that the ydd value is bad because of cancellation error.        
  
 For each new proposed hg, we allow H0_ITERS attempts to
 resolve a possible recoverable failure from f() by reducing
 the proposed stepsize by a factor of 0.2. If a legal stepsize
 still cannot be found, fall back on a previous value if 
 possible, or else return ARK_REPTD_RHSFUNC_ERR.

 Finally, we apply a bias (0.5) and verify that h0 is within 
 bounds.
---------------------------------------------------------------*/
static int arkHin(ARKodeMem ark_mem, realtype tout)
{
  int retval, sign, count1, count2;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  booleantype hgOK, hnewOK;

  /* If tout is too close to tn, give up */
  if ((tdiff = tout-ark_mem->ark_tn) == ZERO) return(ARK_TOO_CLOSE);
  
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = SUNRabs(tdiff);
  tround = ark_mem->ark_uround * SUNMAX(SUNRabs(ark_mem->ark_tn), SUNRabs(tout));

  if (tdist < TWO*tround) return(ARK_TOO_CLOSE);
  
  /* Set lower and upper bounds on h0, and take geometric mean 
     as first trial value.
     Exit with this value if the bounds cross each other. */
  hlb = H0_LBFACTOR * tround;
  hub = arkUpperBoundH0(ark_mem, tdist);

  hg  = SUNRsqrt(hlb*hub);

  if (hub < hlb) {
    if (sign == -1) ark_mem->ark_h = -hg;
    else            ark_mem->ark_h =  hg;
    return(ARK_SUCCESS);
  }
  
  /* Outer loop */
  hnewOK = FALSE;
  hs = hg;     /* safeguard against 'uninitialized variable' warning */
  for(count1 = 1; count1 <= H0_ITERS; count1++) {

    /* Attempts to estimate ydd */
    hgOK = FALSE;

    for (count2 = 1; count2 <= H0_ITERS; count2++) {
      hgs = hg*sign;
      retval = arkYddNorm(ark_mem, hgs, &yddnrm);
      /* If f() failed unrecoverably, give up */
      if (retval < 0) return(ARK_RHSFUNC_FAIL);
      /* If successful, we can use ydd */
      if (retval == ARK_SUCCESS) {hgOK = TRUE; break;}
      /* f() failed recoverably; cut step size and test it again */
      hg *= POINT2;
    }

    /* If f() failed recoverably H0_ITERS times */
    if (!hgOK) {
      /* Exit if this is the first or second pass. No recovery possible */
      if (count1 <= 2) return(ARK_REPTD_RHSFUNC_ERR);
      /* We have a fall-back option. The value hs is a previous hnew which
         passed through f(). Use it and break */
      hnew = hs;
      break;
    }

    /* The proposed step size is feasible. Save it. */
    hs = hg;

    /* If stopping criteria was met, or if this is the last pass, stop */
    if ( (hnewOK) || (count1 == H0_ITERS))  {hnew = hg; break;}

    /* Propose new step size */
    hnew = (yddnrm*hub*hub > TWO) ? SUNRsqrt(TWO/yddnrm) : SUNRsqrt(hg*hub);
    hrat = hnew/hg;
    
    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO))  hnewOK = TRUE;

    /* After one pass, if ydd seems to be bad, use fall-back value. */
    if ((count1 > 1) && (hrat > TWO)) {
      hnew = hg;
      hnewOK = TRUE;
    }

    /* Send this value back through f() */
    hg = hnew;
  }

  /* Apply bounds, bias factor, and attach sign */
  h0 = H0_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  ark_mem->ark_h = h0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkUpperBoundH0

 This routine sets an upper bound on abs(h0) based on
 tdist = tn - t0 and the values of y[i]/y'[i].
---------------------------------------------------------------*/
static realtype arkUpperBoundH0(ARKodeMem ark_mem, realtype tdist)
{
  realtype hub_inv, hub;
  N_Vector temp1, temp2;

  /* Bound based on |y0|/|y0'| -- allow at most an increase of
   * H0_UBFACTOR in y0 (based on a forward Euler step). The weight 
   * factor is used as a safeguard against zero components in y0. */
  temp1 = ark_mem->ark_tempv;
  temp2 = ark_mem->ark_acor;

  N_VAbs(ark_mem->ark_ynew, temp2);
  ark_mem->ark_efun(ark_mem->ark_ynew, temp1, ark_mem->ark_e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(H0_UBFACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(ark_mem->ark_fnew, temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /* bound based on tdist -- allow at most a step of magnitude
   * H0_UBFACTOR * tdist */
  hub = H0_UBFACTOR*tdist;

  /* Use the smaller of the two */
  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}


/*---------------------------------------------------------------
 arkYddNorm

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.
---------------------------------------------------------------*/
static int arkYddNorm(ARKodeMem ark_mem, realtype hg, realtype *yddnrm)
{
  int retval;

  /* increment y with a multiple of f */
  N_VLinearSum(hg, ark_mem->ark_fnew, ONE, 
	       ark_mem->ark_ynew, ark_mem->ark_y);

  /* compute y', via the ODE RHS routine */
  retval = arkFullRHS(ark_mem, ark_mem->ark_tn+hg, ark_mem->ark_y, 
		      ark_mem->ark_ftemp, ark_mem->ark_tempv);
  if (retval != 0) return(ARK_RHSFUNC_FAIL);

  /* if using a non-identity mass matrix, update fnew here to get y' */
  if (ark_mem->ark_mass_matrix) {
    N_VScale(ark_mem->ark_h, ark_mem->ark_tempv, ark_mem->ark_tempv);
    retval = ark_mem->ark_msolve(ark_mem, ark_mem->ark_tempv, ark_mem->ark_rwt); 
    /* retval = ark_mem->ark_msolve(ark_mem, ark_mem->ark_tempv, ark_mem->ark_ewt);  */
    N_VScale(ONE/ark_mem->ark_h, ark_mem->ark_tempv, ark_mem->ark_tempv);
    ark_mem->ark_mass_solves++;
    if (retval != ARK_SUCCESS) {
      ark_mem->ark_nmassfails++;
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE", 
		      "arkYddNorm", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }

  /* difference new f and original f to estimate y'' */
  N_VLinearSum(ONE, ark_mem->ark_tempv, -ONE, 
	       ark_mem->ark_fnew, ark_mem->ark_tempv);
  N_VScale(ONE/hg, ark_mem->ark_tempv, ark_mem->ark_tempv);

  /* compute norm of y'' */
  *yddnrm = N_VWrmsNorm(ark_mem->ark_tempv, ark_mem->ark_ewt);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkSetButcherTables

 This routine determines the ERK/DIRK/ARK method to use, based 
 on the desired accuracy and information on whether the problem
 is explicit, implicit or imex.
---------------------------------------------------------------*/
static int arkSetButcherTables(ARKodeMem ark_mem)
{

  /* if tables have already been specified, just return */
  int i,j,q;
  booleantype A_set = FALSE;
  for (i=0; i<ARK_S_MAX; i++)
    for (j=0; j<ARK_S_MAX; j++) {
      if (SUNRabs(ARK_A(ark_mem->ark_Ae,i,j)) > TINY)  A_set = TRUE;
      if (SUNRabs(ARK_A(ark_mem->ark_Ai,i,j)) > TINY)  A_set = TRUE;
    }
  if (A_set)  return (ARK_SUCCESS);

  /**** explicit methods ****/
  if (ark_mem->ark_explicit) {

    switch (ark_mem->ark_q) {
    case(2):    /* Heun-Euler: q=2, p=1, s=2 */
      ARKodeLoadButcherTable(0, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(3):    /* Bogacki-Shampine: q=3, p=2, s=4 */
      ARKodeLoadButcherTable(1, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(4):    /* Zonneveld: q=4, p=3, s=5 */
      ARKodeLoadButcherTable(3, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(5):    /* Cash-Karp: q=5, p=4, s=6 */
      ARKodeLoadButcherTable(6, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(6):    /* Verner: q=6, p=5, s=8 */
      ARKodeLoadButcherTable(10, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "arkSetButcherTables", 
		      "No explicit method at requested order, using q=6.");
      ARKodeLoadButcherTable(12, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;
    }

    /**** implicit methods ****/
  } else if (ark_mem->ark_implicit) {

    switch (ark_mem->ark_q) {

    case(2):    /* SDIRK-2-1: q=2, p=1, s=2 */
      ARKodeLoadButcherTable(11, &ark_mem->ark_stages,
    			     &q,
    			     &ark_mem->ark_p,
    			     ark_mem->ark_Ai,
    			     ark_mem->ark_b,
    			     ark_mem->ark_c,
    			     ark_mem->ark_b2);
      break;

    /* case(3):    /\* Billington: q=3, p=2, s=3 *\/ */
    /*   ARKodeLoadButcherTable(12, &ark_mem->ark_stages, */
    /* 			     &q, */
    /* 			     &ark_mem->ark_p, */
    /* 			     ark_mem->ark_Ai, */
    /* 			     ark_mem->ark_b, */
    /* 			     ark_mem->ark_c, */
    /* 			     ark_mem->ark_b2); */
    /* case(3):    /\* Kvaerno(4,2,3): q=3, p=2, s=4 *\/ */
    /*   ARKodeLoadButcherTable(14, &ark_mem->ark_stages, */
    /* 			     &q, */
    /* 			     &ark_mem->ark_p, */
    /* 			     ark_mem->ark_Ai, */
    /* 			     ark_mem->ark_b, */
    /* 			     ark_mem->ark_c, */
    /* 			     ark_mem->ark_b2); */
    case(3):    /* ARK3(2)4L[2]SA: q=3, p=2, s=4 */
      ARKodeLoadButcherTable(15, &ark_mem->ark_stages,
    			     &q,
    			     &ark_mem->ark_p,
    			     ark_mem->ark_Ai,
    			     ark_mem->ark_b,
    			     ark_mem->ark_c,
    			     ark_mem->ark_b2);
      break;

    case(4):    /* SDIRK-5-4: q=4, p=3, s=5, A/L stable */
      ARKodeLoadButcherTable(18, &ark_mem->ark_stages,
    			     &q,
    			     &ark_mem->ark_p,
    			     ark_mem->ark_Ai,
    			     ark_mem->ark_b,
    			     ark_mem->ark_c,
    			     ark_mem->ark_b2);
      break;

    /* case(5):    /\* Kvaerno(7,4,5): q=5, p=4, s=7, A/L stable *\/ */
    /*   ARKodeLoadButcherTable(21, &ark_mem->ark_stages, */
    /* 			     &q, */
    /* 			     &ark_mem->ark_p, */
    /* 			     ark_mem->ark_Ai, */
    /* 			     ark_mem->ark_b, */
    /* 			     ark_mem->ark_c, */
    /* 			     ark_mem->ark_b2); */
    case(5):    /* ARK5(4)8L[2]SA: q=5, p=4, s=8, A stable */
      ARKodeLoadButcherTable(22, &ark_mem->ark_stages,
    			     &q,
    			     &ark_mem->ark_p,
    			     ark_mem->ark_Ai,
    			     ark_mem->ark_b,
    			     ark_mem->ark_c,
    			     ark_mem->ark_b2);
      break;

    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "arkSetButcherTables", 
		      "No implicit method at requested order, using q=5.");
      ARKodeLoadButcherTable(26, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;
    }
 
    /**** ImEx methods ****/
  } else {

    switch (ark_mem->ark_q) {

    case(2):
    case(3):    /* ARK3(2)4L[2]SA: q=3, p=2, s=4, DIRK is A/L stable */
      ARKodeLoadButcherTable(2, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      ARKodeLoadButcherTable(15, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(4):    /* ARK4(3)6L[2]SA: q=4, p=3, s=6, DIRK is A/L stable */
      ARKodeLoadButcherTable(4, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      ARKodeLoadButcherTable(20, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(5):    /* ARK5(4)8L[2]SA: q=5, p=4, s=8, DIRK is A stable */
      ARKodeLoadButcherTable(9, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      ARKodeLoadButcherTable(22, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "arkSetButcherTables", 
		      "No ImEx method at requested order, using q=5.");
      ARKodeLoadButcherTable(11, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      ARKodeLoadButcherTable(27, &ark_mem->ark_stages, 
			     &q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;
    }

  }

  /* if input method order does not match actual, update here */
  ark_mem->ark_q = q;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkCheckButcherTables

 This routine runs through the explicit and/or implicit Butcher 
 tables to ensure that they meet all necessary requirements, 
 including:
     strictly lower-triangular (ERK)
     lower-triangular with some nonzeros on diagonal (IRK)
     method order q > 0 (all)
     embedding order q > 0 (all)
     stages > 0 (all)
     matching stages (ARK)
     matching p (ARK)
     matching q (ARK)
     matching c (ARK)
     matching b (ARK)
     matching b2 (ARK)

Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
---------------------------------------------------------------*/
static int arkCheckButcherTables(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  realtype tol = RCONST(1.0e-12);

  /* check that ERK table is strictly lower triangular */
  if (!ark_mem->ark_implicit) {
    okay = TRUE;
    for (i=0; i<ark_mem->ark_stages; i++)
      for (j=i; j<ark_mem->ark_stages; j++)
	if (SUNRabs(ARK_A(ark_mem->ark_Ae,i,j)) > tol)
	  okay = FALSE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "arkCheckButcherTables",
		      "Ae Butcher table is implicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is implicit */
  if (!ark_mem->ark_explicit) {
    okay = FALSE;
    for (i=0; i<ark_mem->ark_stages; i++)
      if (SUNRabs(ARK_A(ark_mem->ark_Ai,i,i)) > tol)
	okay = TRUE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "arkCheckButcherTables",
		      "Ai Butcher table is explicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is lower triangular */
  if (!ark_mem->ark_explicit) {
    okay = TRUE;
    for (i=0; i<ark_mem->ark_stages; i++)
      for (j=i+1; j<ark_mem->ark_stages; j++)
	if (SUNRabs(ARK_A(ark_mem->ark_Ai,i,j)) > tol)
	  okay = FALSE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "arkCheckButcherTables",
		      "Ai Butcher table has entries above diagonal!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that method order q > 0 */
  if (ark_mem->ark_q < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "arkCheckButcherTables",
		    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if (ark_mem->ark_p < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "arkCheckButcherTables",
		    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that stages > 0 */
  if (ark_mem->ark_stages < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "arkCheckButcherTables",
		    "stages < 1!");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep

 This routine performs one internal arkode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
 - preliminary adjustments if a new step size was chosen;
 - loop over internal time step stages:
   - predict stage solution (if implicit);
   - set up RHS data using old solutions;
   - solve the implicit system, or evaluate explicit update;
 - testing the local error;
 - reset stepsize for the next step.
 On a failure in the linear/nonlinear solvers or error test, the
 step may be reattempted, depending on the nature of the failure.
---------------------------------------------------------------*/
static int arkStep(ARKodeMem ark_mem)
{
  realtype saved_t, dsm;
  int retval, ncf, nef, is, nflag, kflag, eflag;
  
  saved_t = ark_mem->ark_tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;
  eflag = ARK_SUCCESS;
  kflag = SOLVE_SUCCESS;

  if ((ark_mem->ark_nst > 0) && (ark_mem->ark_hprime != ark_mem->ark_h)) {
    ark_mem->ark_h = ark_mem->ark_h * ark_mem->ark_eta;
    ark_mem->ark_next_h = ark_mem->ark_h;
  }

  /* Looping point for attempts to take a step */
  for(;;) {  

    /* increment attempt counter */
    ark_mem->ark_nst_attempts++;

    /* Loop over internal stages to the step */
    for (is=0; is<ark_mem->ark_stages; is++) {

      /* store current stage index */
      ark_mem->ark_istage = is;

      /* Update current time with new time step.  If tstop is enabled, it is
	 possible for tn + h to be past tstop by roundoff, and in that case, 
	 we reset tn (after incrementing by h) to tstop. 
	 NOTE: we have already guarded against larger overshoots of tstop. */
      ark_mem->ark_tn = ark_mem->ark_tnew + ark_mem->ark_c[is]*ark_mem->ark_h;
      if (ark_mem->ark_tstopset)
	if ((ark_mem->ark_tn - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO) 
	  ark_mem->ark_tn = ark_mem->ark_tstop;

#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %19.16g,  t_n = %19.16g\n", 
	 ark_mem->ark_nst, is, ark_mem->ark_h, ark_mem->ark_tn);
#endif
      
      /* Call predictor for current stage solution (result placed in ycur) */
      arkPredict(ark_mem, is);  

#ifdef DEBUG_OUTPUT
 printf("predictor:\n");
 N_VPrint_Serial(ark_mem->ark_ycur);
#endif
      
      /* Set up data for evaluation of ARK stage residual (data stored in ark_sdata) */
      if (arkSet(ark_mem) != ARK_SUCCESS)  
	return (ARK_MASSMULT_FAIL);

#ifdef DEBUG_OUTPUT
 printf("rhs data:\n");
 N_VPrint_Serial(ark_mem->ark_sdata);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->ark_report) 	
	fprintf(ark_mem->ark_diagfp, "step  %li  %19.16g  %i  %19.16g\n",
		ark_mem->ark_nst, ark_mem->ark_h, is, ark_mem->ark_tn);

      /* solve implicit problem (if required) */
      if (SUNRabs(ARK_A(ark_mem->ark_Ai,is,is)) > TINY) {

	/* perform implicit solve */
	nflag = arkNls(ark_mem, nflag);

#ifdef DEBUG_OUTPUT
 printf("nonlinear solution:\n");
 N_VPrint_Serial(ark_mem->ark_y);
#endif

	/* check for convergence (on failure, h will have been modified) */
	kflag = arkHandleNFlag(ark_mem, &nflag, saved_t, &ncf);

	/* If fixed time-stepping is used, then anything other than a 
	   successful solve must result in an error */
	if (ark_mem->ark_fixedstep && (kflag != SOLVE_SUCCESS)) 
	  return(kflag);

	/* If h reduced and step needs to be retried, break loop */
	if (kflag == PREDICT_AGAIN) break;

	/* Return if nonlinear solve failed and recovery not possible. */
	if (kflag != SOLVE_SUCCESS) return(kflag);

      /* otherwise no nonlinear solve is needed */
      } else {

	/* if M!=I, solve with M to compute update (place back in sdata) */
	if (ark_mem->ark_mass_matrix) {

	  /* perform mass matrix solve */
	  /* nflag = ark_mem->ark_msolve(ark_mem, ark_mem->ark_sdata,  */
	  /* 			      ark_mem->ark_ewt);  */
	  nflag = ark_mem->ark_msolve(ark_mem, ark_mem->ark_sdata, 
				      ark_mem->ark_rwt); 
	  ark_mem->ark_mass_solves++;
	  
	  /* check for convergence (on failure, h will have been modified) */
	  kflag = arkHandleNFlag(ark_mem, &nflag, saved_t, &ncf);

	  /* increment failure counter if necessary */
	  if (kflag != SOLVE_SUCCESS)
	    ark_mem->ark_nmassfails++;

	  /* If fixed time-stepping is used, then anything other than a 
	     successful solve must result in an error */
	  if (ark_mem->ark_fixedstep && (kflag != SOLVE_SUCCESS)) 
	    return(kflag);

	  /* If h reduced and step needs to be retried, break loop */
	  if (kflag == PREDICT_AGAIN) break;

	  /* Return if solve failed and recovery not possible. */
	  if (kflag != SOLVE_SUCCESS) return(kflag);

      	/* if M==I, set y to be ycur + RHS data computed in arkSet */
	} else {
	  N_VLinearSum(ONE, ark_mem->ark_sdata, ONE, 
		       ark_mem->ark_ycur, ark_mem->ark_y);
	}

#ifdef DEBUG_OUTPUT
 printf("explicit solution:\n");
 N_VPrint_Serial(ark_mem->ark_y);
#endif

      }

      /* successful stage solve, fill in implicit and explicit RHS for this stage */
      if (!ark_mem->ark_explicit) {
	retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_y,
				 ark_mem->ark_Fi[is], ark_mem->ark_user_data);
	ark_mem->ark_nfi++; 
	if (retval < 0)  return(ARK_RHSFUNC_FAIL);
	if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }
      if (!ark_mem->ark_implicit) {
	retval = ark_mem->ark_fe(ark_mem->ark_tn, ark_mem->ark_y, 
				 ark_mem->ark_Fe[is], ark_mem->ark_user_data);
	ark_mem->ark_nfe++;
	if (retval < 0)  return(ARK_RHSFUNC_FAIL);
	if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

    } /* loop over stages */

    /* if h has changed due to convergence failure and a new 
       prediction is needed, continue to next attempt at step
       (cannot occur if fixed time stepping is enabled) */
    if (kflag == PREDICT_AGAIN)  continue;
      
    /* compute time-evolved solution (in ark_y), error estimate (in dsm) */
    retval = arkComputeSolutions(ark_mem, &dsm);
    if (retval < 0)  return(retval);    /* msetup failure */

#ifdef DEBUG_OUTPUT
 printf("error estimate = %19.16g\n", dsm);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->ark_report) 
      fprintf(ark_mem->ark_diagfp, "  etest  %li  %19.16g  %19.16g\n", 
	      ark_mem->ark_nst, ark_mem->ark_h, dsm);

    /* Perform time accuracy error test (if failure, updates h for next try) */
    if (!ark_mem->ark_fixedstep) 
      eflag = arkDoErrorTest(ark_mem, &nflag, saved_t, &nef, dsm);
	
#ifdef DEBUG_OUTPUT
 printf("error test flag = %i\n", eflag);
#endif

    /* Restart step attempt (recompute all stages) if error test fails recoverably */
    if (eflag == TRY_AGAIN)  continue;
	
    /* Return if error test failed and recovery not possible. */
    if (eflag != ARK_SUCCESS)  return(eflag);
	
    /* Error test passed (eflag=ARK_SUCCESS), break from loop */
    break;

  } /* loop over step attempts */

  /* Nonlinear system solves and error test were all successful.
     Update data, and consider change of step. */
  eflag = arkCompleteStep(ark_mem, dsm); 
  if (eflag != ARK_SUCCESS)  return(eflag);
  retval = arkPrepareNextStep(ark_mem); 
  if (retval != ARK_SUCCESS)  return(retval);

  /* Reset growth factor for subsequent time step */
  ark_mem->ark_etamax = ark_mem->ark_hadapt_growth;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkPredict

 This routine computes the prediction for a specific internal 
 stage solution, storing the result in ark_mem->ark_ycur.  The 
 prediction is done using the dense output routine in 
 extrapolation mode, hence stages "far" from the previous time 
 interval are predicted using lower order polynomials than the
 "nearby" stages.
---------------------------------------------------------------*/
static void arkPredict(ARKodeMem ark_mem, int istage)
{
  int i, retval, ord, jstage;
  realtype tau, tau_tol = 0.5;
  realtype h, a0, a1, a2;
  N_Vector yguess = ark_mem->ark_ycur;

  /* if the first step (or if resized), use initial condition as guess */
  if (ark_mem->ark_nst == 0 || ark_mem->ark_resized) {
    N_VScale(ONE, ark_mem->ark_ynew, yguess);
    return;
  }

  /* set evaluation time tau relative shift from previous successful time */
  tau = ark_mem->ark_c[istage]*ark_mem->ark_h/ark_mem->ark_hold;

  /* use requested predictor formula */
  switch (ark_mem->ark_predictor) {

  case 1:

    /***** Dense Output Predictor 1 -- all to max order *****/
    retval = arkDenseEval(ark_mem, tau, 0, ark_mem->ark_dense_q, yguess);
    if (retval == ARK_SUCCESS) return;
    break;

  case 2:

    /***** Dense Output Predictor 2 -- decrease order w/ increasing stage *****/
    ord = SUNMAX(ark_mem->ark_dense_q - istage, 1);
    retval = arkDenseEval(ark_mem, tau, 0, ord, yguess);
    if (retval == ARK_SUCCESS)  return;
    break;

  case 3:

    /***** Max order dense output for stages "close" to previous step, 
	   first-order dense output predictor for subsequent stages *****/
    if (tau < tau_tol) {
      retval = arkDenseEval(ark_mem, tau, 0, ark_mem->ark_dense_q, yguess);
    } else {
      retval = arkDenseEval(ark_mem, tau, 0, 1, yguess);
    }
    if (retval == ARK_SUCCESS)  return;
    break;

  case 4:

    /***** Bootstrap predictor: if any previous stage in step has nonzero c_i, 
	   construct a quadratic Hermite interpolant for prediction; otherwise 
	   use the trivial predictor. *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (ark_mem->ark_mass_matrix)  break;

    /* determine if any previous stages in step meet criteria */
    jstage = -1;
    for (i=0; i<istage; i++)
      jstage = (ark_mem->ark_c[i] != ZERO) ? i : jstage;

    /* if using the trivial predictor, break */
    if (jstage == -1)  break;
    
    /* find the "optimal" previous stage to use */
    for (i=0; i<istage; i++) 
      if ((ark_mem->ark_c[i] > ark_mem->ark_c[jstage]) && (ark_mem->ark_c[i] != ZERO))
	jstage = i;

    /* evaluate the quadratic Hermite interpolant for the prediction */
    h = ark_mem->ark_h * ark_mem->ark_c[jstage];
    tau = ark_mem->ark_h * ark_mem->ark_c[istage];
    a0 = ONE;
    a2 = tau*tau/TWO/h;
    a1 = tau - a2;

    N_VLinearSum(a0, ark_mem->ark_ynew, a1, ark_mem->ark_fnew, yguess);
    if (!ark_mem->ark_explicit) 
      N_VLinearSum(a2, ark_mem->ark_Fi[jstage], ONE, yguess, yguess);
    if (!ark_mem->ark_implicit) 
      N_VLinearSum(a2, ark_mem->ark_Fe[jstage], ONE, yguess, yguess);
    break;

  }

  /* if we made it here, use the trivial predictor (previous step solution) */
  N_VScale(ONE, ark_mem->ark_ynew, yguess);

}


/*---------------------------------------------------------------
 arkSet

 This routine sets up the stage data for computing the RK 
 residual, along with the step- and method-related factors 
 gamma, gammap and gamrat.
---------------------------------------------------------------*/
static int arkSet(ARKodeMem ark_mem)
{
  /* local data */
  realtype hA;
  int retval, j, i = ark_mem->ark_istage;
  N_Vector tmp = ark_mem->ark_tempv;

  /* Initialize sdata to yn - ycur (ycur holds guess) */
  N_VLinearSum(ONE, ark_mem->ark_ynew, -ONE, ark_mem->ark_ycur, 
	       ark_mem->ark_sdata);


  /* If M!=I, replace sdata with M*sdata, so that sdata = M*(yn-ycur) */
  if (ark_mem->ark_mass_matrix) {
    N_VScale(ONE, ark_mem->ark_sdata, tmp);
    retval = ark_mem->ark_mtimes(tmp, ark_mem->ark_sdata, ark_mem->ark_tn, 
				 ark_mem->ark_mtimes_data);
    ark_mem->ark_mass_mult++;
    if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
  }

  /* Iterate over each prior stage updating rhs */
  /*    Explicit pieces */
  if (!ark_mem->ark_implicit)
    for (j=0; j<i; j++) {
      hA = ark_mem->ark_h * ARK_A(ark_mem->ark_Ae,i,j);
      N_VLinearSum(hA, ark_mem->ark_Fe[j], ONE, 
		   ark_mem->ark_sdata, ark_mem->ark_sdata);
    }
  /*    Implicit pieces */
  if (!ark_mem->ark_explicit)
    for (j=0; j<i; j++) {
      hA = ark_mem->ark_h * ARK_A(ark_mem->ark_Ai,i,j);
      N_VLinearSum(hA, ark_mem->ark_Fi[j], ONE, 
		   ark_mem->ark_sdata, ark_mem->ark_sdata);
    }

  /* Update gamma */
  ark_mem->ark_gamma = ark_mem->ark_h * ARK_A(ark_mem->ark_Ai,i,i);
  if (ark_mem->ark_firststage)  
    ark_mem->ark_gammap = ark_mem->ark_gamma;
  ark_mem->ark_gamrat = (ark_mem->ark_firststage) ? 
    ONE : ark_mem->ark_gamma / ark_mem->ark_gammap;  /* protect x/x != 1.0 */

  /* return with success */
  return (ARK_SUCCESS);

}


/*---------------------------------------------------------------
 arkComputeSolutions

 This routine calculates the final RK solution using the existing 
 data.  This solution is placed directly in ark_y.  This routine 
 also computes the error estimate ||y-ytilde||_WRMS, where ytilde 
 is the embedded solution, and the norm weights come from 
 ark_ewt.  This norm value is returned.  The vector form of this 
 estimated error (y-ytilde) is stored in ark_tempv, in case the 
 calling routine wishes to examine the error locations.

 Note: at this point in the step, the vectors ark_tempv, 
 ark_sdata and ark_ycur may all be used as temporary vectors.

 For now, we assume that the ODE is of the form 
    y' = fe(t,y) + fi(t,y)
 so that both y and ytilde can be formed from existing data.  
 However, once we extend the solver to allow for non-identity
 mass matrices, this routine will perform two additional solves
 to compute the solution and embedding.
---------------------------------------------------------------*/
static int arkComputeSolutions(ARKodeMem ark_mem, realtype *dsm)
{
  /* local data */
  realtype hb;
  int ier, j;
  N_Vector y    = ark_mem->ark_y;
  N_Vector yerr = ark_mem->ark_tempv;
  N_Vector tmp1 = ark_mem->ark_sdata;
  N_Vector tmp2 = ark_mem->ark_ycur;

  /* Compute updated solution and error estimate based on whether
     a non-identity mass matrix is present */
  if (ark_mem->ark_mass_matrix) {   /* M != I */

    /* setup mass matrix, using y, tmp1, tmp2 as temporaries */
    if (ark_mem->ark_MassSetupNonNull) {
      ier = ark_mem->ark_msetup(ark_mem, y, tmp1, tmp2);
      if (ier != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
    }

    /* compute y RHS (store in y) */
    N_VConst(ZERO, y);
    for (j=0; j<ark_mem->ark_stages; j++) {
      hb = ark_mem->ark_h * ark_mem->ark_b[j];
      if (!ark_mem->ark_implicit)
	N_VLinearSum(hb, ark_mem->ark_Fe[j], ONE, y, y);
      if (!ark_mem->ark_explicit)
	N_VLinearSum(hb, ark_mem->ark_Fi[j], ONE, y, y);
    }
    
    /* solve for y update (stored in y) */
    /* ier = ark_mem->ark_msolve(ark_mem, y, ark_mem->ark_ewt);  */
    ier = ark_mem->ark_msolve(ark_mem, y, ark_mem->ark_rwt); 
    ark_mem->ark_mass_solves++;
    if (ier < 0) {
      ark_mem->ark_nmassfails++;
      *dsm = 2.0;         /* indicate too much error, step with smaller step */
      N_VScale(ONE, ark_mem->ark_ynew, y);      /* place old solution into y */
      return(CONV_FAIL);
    }

    /* compute y = ynew + update */
    N_VLinearSum(ONE, ark_mem->ark_ynew, ONE, y, y);


    /* compute yerr RHS (store in yerr) */
    N_VConst(ZERO, yerr);
    for (j=0; j<ark_mem->ark_stages; j++) {
      hb = ark_mem->ark_h * (ark_mem->ark_b[j] - ark_mem->ark_b2[j]);
      if (!ark_mem->ark_implicit)
	N_VLinearSum(hb, ark_mem->ark_Fe[j], ONE, yerr, yerr);
      if (!ark_mem->ark_explicit)
	N_VLinearSum(hb, ark_mem->ark_Fi[j], ONE, yerr, yerr);
    }

    /* solve for yerr */
    /* ier = ark_mem->ark_msolve(ark_mem, yerr, ark_mem->ark_ewt);  */
    ier = ark_mem->ark_msolve(ark_mem, yerr, ark_mem->ark_rwt); 
    ark_mem->ark_mass_solves++;
    if (ier < 0) {
      ark_mem->ark_nmassfails++;
      *dsm = 2.0;         /* indicate too much error, step with smaller step */
      return(CONV_FAIL);
    }

  } else {                          /* M == I */

    /* Initialize solution to yn, error estimate to zero */
    N_VScale(ONE, ark_mem->ark_ynew, y);
    N_VConst(ZERO, yerr);

    /* Iterate over each stage, updating solution and error estimate */
    for (j=0; j<ark_mem->ark_stages; j++) {
      
      /* solution */
      hb = ark_mem->ark_h * ark_mem->ark_b[j];
      if (!ark_mem->ark_implicit)
	N_VLinearSum(hb, ark_mem->ark_Fe[j], ONE, y, y);
      if (!ark_mem->ark_explicit)
	N_VLinearSum(hb, ark_mem->ark_Fi[j], ONE, y, y);
      
      /* error */
      hb = ark_mem->ark_h * (ark_mem->ark_b[j] - ark_mem->ark_b2[j]);
      if (!ark_mem->ark_implicit)
	N_VLinearSum(hb, ark_mem->ark_Fe[j], ONE, yerr, yerr);
      if (!ark_mem->ark_explicit)
	N_VLinearSum(hb, ark_mem->ark_Fi[j], ONE, yerr, yerr);
      
    }

  }

  /* fill error norm and return */
  *dsm = N_VWrmsNorm(yerr, ark_mem->ark_ewt);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkDoErrorTest

 This routine performs the local error test for the ARK method. 
 The weighted local error norm dsm is passed in, and 
 the test dsm ?<= 1 is made.

 If the test passes, arkDoErrorTest returns ARK_SUCCESS. 

 If the test fails, we revert to the last successful solution 
 time, and:
   - if maxnef error test failures have occurred or if 
     SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
   - otherwise: update time step factor eta based on local error 
     estimate and reduce h.  Then set *nflagPtr to PREV_ERR_FAIL, 
     and return TRY_AGAIN. 
---------------------------------------------------------------*/
static int arkDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
			  realtype saved_t, int *nefPtr, realtype dsm)
{
  realtype ehist2, hhist2;
  int retval;

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */  
  if (dsm <= ONE) return(ARK_SUCCESS);
  
  /* Test failed; increment counters, set nflag, and restore tn */
  (*nefPtr)++;
  ark_mem->ark_netf++;
  *nflagPtr = PREV_ERR_FAIL;
  ark_mem->ark_tn = saved_t;

  /* At maxnef failures or |h| = hmin, return ARK_ERR_FAILURE */
  if ((SUNRabs(ark_mem->ark_h) <= ark_mem->ark_hmin*ONEPSM) ||
      (*nefPtr == ark_mem->ark_maxnef)) return(ARK_ERR_FAILURE);

  /* Set etamax=1 to prevent step size increase at end of this step */
  ark_mem->ark_etamax = ONE;

  /* Temporarily update error history array for recomputation of h */
  ehist2 = ark_mem->ark_hadapt_ehist[2];
  ark_mem->ark_hadapt_ehist[2] = ark_mem->ark_hadapt_ehist[1];
  ark_mem->ark_hadapt_ehist[1] = ark_mem->ark_hadapt_ehist[0];
  ark_mem->ark_hadapt_ehist[0] = dsm*ark_mem->ark_hadapt_bias;

  /* Temporarily update step history array for recomputation of h */
  hhist2 = ark_mem->ark_hadapt_hhist[2];
  ark_mem->ark_hadapt_hhist[2] = ark_mem->ark_hadapt_hhist[1];
  ark_mem->ark_hadapt_hhist[1] = ark_mem->ark_hadapt_hhist[0];
  ark_mem->ark_hadapt_hhist[0] = ark_mem->ark_h;

  /* Compute accuracy-based time step estimate (updated ark_eta) */
  retval = arkAdapt(ark_mem);
  if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);

  /* Revert error history array */
  ark_mem->ark_hadapt_ehist[0] = ark_mem->ark_hadapt_ehist[1];
  ark_mem->ark_hadapt_ehist[1] = ark_mem->ark_hadapt_ehist[2];
  ark_mem->ark_hadapt_ehist[2] = ehist2;

  /* Revert step history array */
  ark_mem->ark_hadapt_hhist[0] = ark_mem->ark_hadapt_hhist[1];
  ark_mem->ark_hadapt_hhist[1] = ark_mem->ark_hadapt_hhist[2];
  ark_mem->ark_hadapt_hhist[2] = hhist2;

  /* Enforce failure bounds on eta, update h, and return for retry of step */
  if (*nefPtr >= ark_mem->ark_small_nef) 
    ark_mem->ark_eta = SUNMIN(ark_mem->ark_eta, ark_mem->ark_etamxf);
  ark_mem->ark_h *= ark_mem->ark_eta;
  ark_mem->ark_next_h = ark_mem->ark_h;
  return(TRY_AGAIN);
}


/*---------------------------------------------------------------
 arkCompleteStep

 This routine performs various update operations when the step 
 solution has passed the local error test.  We update the current
 time, increment the step counter nst, record the values hold and
 tnew, reset the resized flag, update the error and time step 
 history arrays, update ycur to hold the current solution, and 
 update the yold, ynew, fold and fnew arrays.
---------------------------------------------------------------*/
static int arkCompleteStep(ARKodeMem ark_mem, realtype dsm)
{
  int retval;
  N_Vector tempvec;

  /* Set current time to the end of the step (in case the last 
     stage time does not coincide with the step solution time).  
     If tstop is enabled, it is possible for tn + h to be past 
     tstop by roundoff, and in that case, we reset tn (after 
     incrementing by h) to tstop. */
  ark_mem->ark_tn = ark_mem->ark_tnew + ark_mem->ark_h;
  if (ark_mem->ark_tstopset) {
    if ((ark_mem->ark_tn - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO) 
      ark_mem->ark_tn = ark_mem->ark_tstop;
  }

  /* update scalar quantities */
  ark_mem->ark_nst++;
  ark_mem->ark_hold = ark_mem->ark_h;
  ark_mem->ark_tnew = ark_mem->ark_tn;

  /* turn off flag regarding resized problem */
  ark_mem->ark_resized = FALSE;

  /* update error history array */
  ark_mem->ark_hadapt_ehist[2] = ark_mem->ark_hadapt_ehist[1];
  ark_mem->ark_hadapt_ehist[1] = ark_mem->ark_hadapt_ehist[0];
  ark_mem->ark_hadapt_ehist[0] = dsm*ark_mem->ark_hadapt_bias;

  /* update step history array */
  ark_mem->ark_hadapt_hhist[2] = ark_mem->ark_hadapt_hhist[1];
  ark_mem->ark_hadapt_hhist[1] = ark_mem->ark_hadapt_hhist[0];
  ark_mem->ark_hadapt_hhist[0] = ark_mem->ark_h;

  /* update ycur to current solution */
  N_VScale(ONE, ark_mem->ark_y, ark_mem->ark_ycur);
  
  /* swap yold and ynew arrays */
  tempvec = ark_mem->ark_yold;
  ark_mem->ark_yold = ark_mem->ark_ynew;
  ark_mem->ark_ynew = tempvec;

  /* swap fold and fnew arrays */
  tempvec = ark_mem->ark_fold;
  ark_mem->ark_fold = ark_mem->ark_fnew;
  ark_mem->ark_fnew = tempvec;

  /* update fnew array using explicit and implicit RHS:
     if c[s-1] = 1.0 use already-computed RHS values, 
     otherwise compute from scratch */
  N_VConst(ZERO, ark_mem->ark_fnew);
  if (SUNRabs(ark_mem->ark_c[ark_mem->ark_stages-1] - ONE) < TINY) {
    if (!ark_mem->ark_explicit) 
      N_VLinearSum(ONE, ark_mem->ark_Fi[ark_mem->ark_stages-1], 
		   ONE, ark_mem->ark_fnew, ark_mem->ark_fnew);
    if (!ark_mem->ark_implicit) 
      N_VLinearSum(ONE, ark_mem->ark_Fe[ark_mem->ark_stages-1], 
		   ONE, ark_mem->ark_fnew, ark_mem->ark_fnew);
  } else {
    /* NOTE: new y value is currently in yold, due to swap above */
    retval = arkFullRHS(ark_mem, ark_mem->ark_tn, ark_mem->ark_yold, 
			ark_mem->ark_ftemp, ark_mem->ark_fnew);
    if (retval != 0) return(ARK_RHSFUNC_FAIL);
  }

  /* if M!=I, update fnew with M^{-1}*fnew (note, mass matrix already current) */
  if (ark_mem->ark_mass_matrix) {   /* M != I */
    N_VScale(ark_mem->ark_h, ark_mem->ark_fnew, ark_mem->ark_fnew);      /* scale RHS */
    /* retval = ark_mem->ark_msolve(ark_mem, ark_mem->ark_fnew, ark_mem->ark_ewt);  */
    retval = ark_mem->ark_msolve(ark_mem, ark_mem->ark_fnew, ark_mem->ark_rwt); 
    N_VScale(ONE/ark_mem->ark_h, ark_mem->ark_fnew, ark_mem->ark_fnew);  /* scale result */
    ark_mem->ark_mass_solves++;
    if (retval != ARK_SUCCESS) {
      ark_mem->ark_nmassfails++;
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE", 
		      "arkCompleteStep", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }

  /* update ynew (cannot swap since ark_y is a pointer to user-supplied data) */
  N_VScale(ONE, ark_mem->ark_y, ark_mem->ark_ynew);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkPrepareNextStep

 This routine handles the setting of stepsize, hprime, for the 
 next step.  Along with hprime, it sets the ratio eta=hprime/h.
 It also updates other state variables related to a change of 
 step size.
---------------------------------------------------------------*/
 static int arkPrepareNextStep(ARKodeMem ark_mem)
{
  int retval;

  /* If etamax = 1, or if fixed time-stepping requesetd, defer 
     step size changes until next step */
  if ((ark_mem->ark_etamax == ONE) || ark_mem->ark_fixedstep){
    ark_mem->ark_hprime = ark_mem->ark_h;
    ark_mem->ark_eta = ONE;
    return(ARK_SUCCESS);
  }

  /* Adjust ark_eta in arkAdapt */
  retval = arkAdapt(ark_mem);
  
  /* set hprime value for next step size */
  ark_mem->ark_hprime = ark_mem->ark_h * ark_mem->ark_eta;

  if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkNls

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 It calls one of arkNlsAccelFP, arkLs or arkNlsNewton to do the 
 work.

 Upon a successful solve, the solution is held in ark_mem->ark_y.
---------------------------------------------------------------*/
static int arkNls(ARKodeMem ark_mem, int nflag)
{

  /* call the appropriate solver */

  /*   fixed point */
  if (ark_mem->ark_use_fp)
    return(arkNlsAccelFP(ark_mem, nflag));

  /*   linearly implicit (one Newton iteration) */
  if (ark_mem->ark_linear)
    return(arkLs(ark_mem, nflag));
  
  /*   Newton */
  return(arkNlsNewton(ark_mem, nflag));

}


/*---------------------------------------------------------------
 arkNlsResid

 This routine evaluates the negative nonlinear residual for the 
 additive Runge-Kutta method.  It assumes that any data from 
 previous time steps/stages is contained in ark_mem, and merely 
 combines this old data with the current ODE RHS vector, 
 fy = f(t,y), to compute the nonlinear residual r.

 At the ith stage, we compute 
    r = -M*zi + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j) 
                     + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = -M*zi + gamma*Fi(zi) + M*yn + data
 where
    zi is stored in y
    Fi(zi) is stored in fy
    (M*yn+data) is stored in ark_mem->ark_sdata,
 so we really just compute
    r = -M*y + gamma*fy + ark_mem->ark_sdata

 Possible return values:  ARK_SUCCESS  or  ARK_RHSFUNC_FAIL
---------------------------------------------------------------*/
static int arkNlsResid(ARKodeMem ark_mem, N_Vector y, 
		       N_Vector fy, N_Vector r)
{

  /* temporary variables */
  int retval;

  /* put M*y in r */
  if (ark_mem->ark_mass_matrix) {
    retval = ark_mem->ark_mtimes(y, r, ark_mem->ark_tn, 
				 ark_mem->ark_mtimes_data);
    ark_mem->ark_mass_mult++;
    if (retval != ARK_SUCCESS)  
      return (ARK_MASSMULT_FAIL);
  } else {
    N_VScale(ONE, y, r);
  }

  /* update with sdata */
  N_VLinearSum(-ONE, r, ONE, ark_mem->ark_sdata, r);

  /* update with gamma*fy */
  N_VLinearSum(ark_mem->ark_gamma, fy, ONE, r, r);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkNlsNewton

 This routine handles the Newton iteration for implicit portions 
 of the ARK and DIRK methods. It calls lsetup if indicated, 
 performs a modified Newton iteration (using a fixed 
 Jacobian / precond), and retries a failed attempt at Newton 
 iteration if that is indicated. 

 Upon entry, the predicted solution is held in ark_mem->ark_ycur; 
 this array is never changed throughout this routine.  If an 
 initial attempt at solving the nonlinear system fail (e.g. due to 
 a stale Jacobian), this allows for new attempts that first revert 
 ark_mem->ark_y back to this initial guess before trying again. 

 Upon a successful solve, the solution is held in ark_mem->ark_y.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+  
   ARK_LSETUP_FAIL    |-> halt the integration 
   ARK_LSOLVE_FAIL   -+

   CONV_FAIL        -+
   RHSFUNC_RECVR    -+-> predict again or stop if too many
---------------------------------------------------------------*/
static int arkNlsNewton(ARKodeMem ark_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3, b;
  int convfail, retval, ier, m;
  booleantype callSetup;
  realtype del, delp, dcon;
  
  vtemp1 = ark_mem->ark_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = ark_mem->ark_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = ark_mem->ark_tempv; /* rename tempv as vtemp3 for readability */
  b      = ark_mem->ark_tempv; /* also rename tempv as b for readability */

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (ark_mem->ark_setupNonNull) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (ark_mem->ark_firststage) || (ark_mem->ark_msbp < 0) ||
      (ark_mem->ark_nst >= ark_mem->ark_nstlp + abs(ark_mem->ark_msbp)) || 
      (SUNRabs(ark_mem->ark_gamrat-ONE) > ark_mem->ark_dgmax);
  } else {  
    ark_mem->ark_crate = ONE;
    callSetup = FALSE;
  }
  
  /* Looping point for attempts at solution of the nonlinear system:
       Evaluate f at predicted y, store result in ark_mem->ark_ftemp.
       Call lsetup if indicated, setting statistics and gamma factors.
       Zero out the correction array (ark_mem->ark_acor).
       Copy the predicted y (ycur) into the output (ark_mem->ark_y).
       Performs the modified Newton iteration using the existing lsetup.
       Repeat process if a recoverable failure occurred (convergence
	  failure with stale Jacobian). */
  for(;;) {

    if (!ark_mem->ark_explicit) {
      retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_ycur, 
			       ark_mem->ark_ftemp, ark_mem->ark_user_data);
      ark_mem->ark_nfi++; 
      if (retval < 0) return(ARK_RHSFUNC_FAIL);
      if (retval > 0) return(RHSFUNC_RECVR);
    }
    
    if (callSetup) {

      /* Solver diagnostics reporting */
      if (ark_mem->ark_report)  fprintf(ark_mem->ark_diagfp, "  lsetup\n");

      ier = ark_mem->ark_lsetup(ark_mem, convfail, ark_mem->ark_ycur, 
				ark_mem->ark_ftemp, &ark_mem->ark_jcur, 
				vtemp1, vtemp2, vtemp3);
      ark_mem->ark_nsetups++;
      callSetup = FALSE;
      ark_mem->ark_firststage = FALSE;
      ark_mem->ark_gamrat = ark_mem->ark_crate = ONE; 
      ark_mem->ark_gammap = ark_mem->ark_gamma;
      ark_mem->ark_nstlp  = ark_mem->ark_nst;

      /* Return if lsetup failed */
      if (ier < 0) return(ARK_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, ark_mem->ark_acor);
    N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_y);


    /***********************************
     Do the modified Newton iteration:
     - Reuses a single call to ark_mem->ark_lsetup (called by the enclosing loop). 
     - Upon entry, the predicted solution is held in ark_mem->ark_ycur.
     - Here, ark_mem->ark_acor accumulates the difference between the predictor 
       and the final solution to the time step.
     - Upon a successful solve, the solution is held in ark_mem->ark_y. 
    */

    /* Initialize temporary variables for use in iteration */
    ark_mem->ark_mnewt = m = 0;
    del = delp = ZERO;

    /* Reset the local truncation error estimate */
    ark_mem->ark_eLTE = ark_mem->ark_nlscoef * RCONST(0.1);

    /* Looping point for Newton iteration */
    for(;;) {

      /* Evaluate the nonlinear system residual, put result into b */
      retval = arkNlsResid(ark_mem, ark_mem->ark_acor, 
			   ark_mem->ark_ftemp, b);
      if (retval != ARK_SUCCESS) {
	ier = ARK_RHSFUNC_FAIL;
	break;
      }

      /* Call the lsolve function */
      retval = ark_mem->ark_lsolve(ark_mem, b, ark_mem->ark_rwt, 
				   ark_mem->ark_y, ark_mem->ark_ftemp); 
      /* retval = ark_mem->ark_lsolve(ark_mem, b, ark_mem->ark_ewt,  */
      /* 				   ark_mem->ark_y, ark_mem->ark_ftemp);  */
      ark_mem->ark_nni++;
    
      if (retval < 0) {
	ier = ARK_LSOLVE_FAIL;
	break;
      }
    
      /* If lsolve had a recoverable failure and Jacobian data is
	 not current, signal to try the solution again */
      if (retval > 0) { 
	if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	  ier = TRY_AGAIN;
	else 
	  ier = CONV_FAIL;
	break;
      }

      /* Get WRMS norm of correction; add correction to acor and y */
      del = N_VWrmsNorm(b, ark_mem->ark_ewt);
      N_VLinearSum(ONE, ark_mem->ark_acor, ONE, b, ark_mem->ark_acor);
      N_VLinearSum(ONE, ark_mem->ark_ycur, ONE, ark_mem->ark_acor, ark_mem->ark_y);

      /* Test for convergence.  If m > 0, an estimate of the convergence
	 rate constant is stored in crate, and used in the test */
      if (m > 0) 
	ark_mem->ark_crate = SUNMAX(ark_mem->ark_crdown*ark_mem->ark_crate, del/delp);
      if (ark_mem->ark_crate < ONE)
      	ark_mem->ark_eLTE = ark_mem->ark_nlscoef * (ONE - ark_mem->ark_crate);
      else
      	ark_mem->ark_eLTE = ark_mem->ark_nlscoef * RCONST(0.1);
      dcon = del * SUNMIN(ONE, ark_mem->ark_crate) / ark_mem->ark_eLTE;

#ifdef DEBUG_OUTPUT
 printf("Newton iter %i,  del = %19.16g,  crate = %19.16g\n", m, del, ark_mem->ark_crate);
 printf("   eLTE = %19.16g,  dcon = %19.16g\n", ark_mem->ark_eLTE, dcon);
 printf("Newton correction:\n");
 N_VPrint_Serial(ark_mem->ark_acor);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->ark_report) 
	fprintf(ark_mem->ark_diagfp, "    newt  %i  %19.16g  %19.16g\n", m, del, dcon);
    
      if (dcon <= ONE) {
	ark_mem->ark_jcur = FALSE;
	ier = ARK_SUCCESS;
	break;
      }

      /* update Newton iteration counter */
      ark_mem->ark_mnewt = ++m;
    
      /* Stop at maxcor iterations or if iteration seems to be diverging.
	 If still not converged and Jacobian data is not current, signal 
	 to try the solution again */
      if ((m == ark_mem->ark_maxcor) || 
	  ((m >= 2) && (del > ark_mem->ark_rdiv*delp))) {
	if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	  ier = TRY_AGAIN;
	else
	  ier = CONV_FAIL;
	break;
      }
    
      /* Save norm of correction, evaluate fi, and loop again */
      delp = del;
      if (!ark_mem->ark_explicit) {
	retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_y, 
				 ark_mem->ark_ftemp, ark_mem->ark_user_data);
	ark_mem->ark_nfi++;
	if (retval < 0) {
	  ier = ARK_RHSFUNC_FAIL;
	  break;
	}
	if (retval > 0) {
	  if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	    ier = TRY_AGAIN;
	  else
	    ier = RHSFUNC_RECVR;
	  break;
	}
      }
    } 
    /* end modified Newton iteration */
    /*********************************/
    
    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=ARK_FAIL_BAD_J.  Otherwise return. */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = TRUE;
    convfail = ARK_FAIL_BAD_J;
  }
}


/*---------------------------------------------------------------
 arkNlsAccelFP

 This routine handles the Anderson-accelerated fixed-point 
 iteration for implicit portions of the ARK and DIRK methods. It 
 performs an accelerated fixed-point iteration (without need for 
 a Jacobian or preconditioner. 

 Upon entry, the predicted solution is held in ark_mem->ark_ycur; 
 this array is never changed throughout this routine.  If an 
 initial attempt at solving the nonlinear system fail, this 
 allows for new attempts that first revert ark_mem->ark_y back 
 to this initial guess before trying again. 

 Upon a successful solve, the solution is held in ark_mem->ark_y.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  ---> halt the integration 

   CONV_FAIL         -+
   RHSFUNC_RECVR     -+-> predict again or stop if too many
---------------------------------------------------------------*/
static int arkNlsAccelFP(ARKodeMem ark_mem, int nflag)
{
  /* local variables */
  int retval;
  realtype del, delp, dcon;

  /* local shortcut variables */
  realtype tn     = ark_mem->ark_tn;
  realtype *R     = ark_mem->ark_fp_R;
  realtype *gamma = ark_mem->ark_fp_gamma;
  long int maa    = ark_mem->ark_fp_m;
  void *udata     = ark_mem->ark_user_data;
  N_Vector ypred  = ark_mem->ark_acor;
  N_Vector y      = ark_mem->ark_y;
  N_Vector ycur   = ark_mem->ark_ycur;
  N_Vector ftemp  = ark_mem->ark_ftemp;
  N_Vector fval   = ark_mem->ark_fp_fval;
  N_Vector tempv  = ark_mem->ark_tempv;

  /* Initialize temporary variables for use in iteration */
  del = delp = ZERO;
  ark_mem->ark_crate = ONE;

  /* Initialize iteration counter */
  ark_mem->ark_mnewt = 0;

  /* Load prediction into y vector, ypred */
  N_VScale(ONE, ycur, y);
  N_VScale(ONE, ycur, ypred);

  /* Looping point for attempts at solution of the nonlinear system:
       Evaluate f at predicted y, store result in ark_mem->ark_ftemp.
       Performs the accelerated fixed-point iteration.
       Repeat process if a recoverable failure occurred (convergence
	  failure with stale Jacobian). */
  for(;;) {

    /* evaluate implicit ODE RHS function, store in ftemp */
    if (!ark_mem->ark_explicit) {
      retval = ark_mem->ark_fi(tn, y, ftemp, udata);
      ark_mem->ark_nfi++;
      if (retval < 0) return(ARK_RHSFUNC_FAIL);
      if (retval > 0) return(RHSFUNC_RECVR);
    }

    /* store difference between current guess and prediction in tempv */
    N_VLinearSum(ONE, y, -ONE, ypred, tempv);

    /* evaluate nonlinear residual, store in fval */
    retval = arkNlsResid(ark_mem, tempv, ftemp, fval);
    if (retval != ARK_SUCCESS) return(ARK_RHSFUNC_FAIL);

    /* convert nonlinear residual result to a fixed-point function result */
    N_VLinearSum(ONE, y, ONE, fval, fval);

    /* perform fixed point update */
    if (maa == 0) {
      /* plain fixed-point solver, copy residual into y */
      N_VScale(ONE, fval, y);
    } else {
      /* Anderson-accelerated solver */
      N_VScale(ONE, ycur, y);   /* update guess */
      retval = arkAndersonAcc(ark_mem, fval, tempv, y,  /* accelerate guess */
			      ycur, (int)(ark_mem->ark_mnewt), R, gamma);
    }
    ark_mem->ark_nni++;

    /* measure convergence. If ark_mem->ark_mnewt > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test */
    N_VLinearSum(ONE, y, -ONE, ycur, tempv);
    del = N_VWrmsNorm(tempv, ark_mem->ark_ewt);
    if (ark_mem->ark_mnewt > 0)
      ark_mem->ark_crate = SUNMAX(ark_mem->ark_crdown*ark_mem->ark_crate, del/delp);
    if (ark_mem->ark_crate < ONE)
      ark_mem->ark_eLTE = ark_mem->ark_nlscoef * (ONE - ark_mem->ark_crate);
    else
      ark_mem->ark_eLTE = ark_mem->ark_nlscoef * RCONST(0.1);
    dcon = del * SUNMIN(ONE, ark_mem->ark_crate) / ark_mem->ark_eLTE;

#ifdef DEBUG_OUTPUT
 printf("FP iter %i,  del = %19.16g,  crate = %19.16g\n", ark_mem->ark_mnewt, del, ark_mem->ark_crate);
 printf("   eLTE = %19.16g,  dcon = %19.16g\n", ark_mem->ark_eLTE, dcon);
 printf("Fixed-point correction:\n");
 N_VPrint_Serial(tempv);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->ark_report)
      fprintf(ark_mem->ark_diagfp, "    fp  %i  %19.16g  %19.16g\n", ark_mem->ark_mnewt, del, dcon);

    /* update iteration counter */
    ark_mem->ark_mnewt++;

    if (dcon <= ONE)  return(ARK_SUCCESS);

    /* Stop at maxcor iterations or if iteration seems to be diverging */
    if ((ark_mem->ark_mnewt == ark_mem->ark_maxcor) ||
	((ark_mem->ark_mnewt >= 2) && (del > ark_mem->ark_rdiv*delp))) 
      return(CONV_FAIL);

    /* Update current solution guess */
    /*** I DON'T LIKE THAT I WAIT TO UPDATE YCUR UNTIL THIS POINT; THIS MEANS THAT 
	 I CAN RETURN WITH A LESS-ACCURATE SOLUTION THAN THE ONE I HAVE IN HAND.  
	 MOREOVER, I ESTIMATE DCON USING THE EXPECTED ERROR IN THE NEW ITERATE, 
	 NOT THE COMPUTED ERROR IN THE PREVIOUS ITERATE, SO MY STOP CRITERIA ARE 
	 EVEN POTENTIALLY INCORRECT. ***/
    N_VScale(ONE, y, ycur);
    
    /* Save norm of correction and loop again */
    delp = del;

  }
}


/*---------------------------------------------------------------
 arkAndersonAcc

 This routine computes the Anderson-accelerated fixed point 
 iterate itself, as used by the nonlinear solver arkNlsAccelFP().  

 Upon entry, the predicted solution is held in xold;
 this array is never changed throughout this routine.  

 The result of the routine is held in x.  

 Possible return values:
   ARK_SUCCESS   ---> successful completion

---------------------------------------------------------------*/
static int arkAndersonAcc(ARKodeMem ark_mem, N_Vector gval, 
			  N_Vector fv, N_Vector x, N_Vector xold, 
			  int iter, realtype *R, realtype *gamma)
{
  /* local variables */
  long int i_pt, i, j, lAA;
  realtype alfa, a, b, temp, c, s;

  /* local shortcut variables */
  N_Vector vtemp2 = ark_mem->ark_y;       /* rename y as vtemp2 for readability */
  long int *ipt_map = ark_mem->ark_fp_imap;
  long int maa = ark_mem->ark_fp_m;
  N_Vector gold = ark_mem->ark_fp_gold;
  N_Vector fold = ark_mem->ark_fp_fold;
  N_Vector *df = ark_mem->ark_fp_df;
  N_Vector *dg = ark_mem->ark_fp_dg;
  N_Vector *Q = ark_mem->ark_fp_q;
  
  for (i=0; i<maa; i++)  ipt_map[i]=0;
  i_pt = iter-1 - ((iter-1)/maa)*maa;
  N_VLinearSum(ONE, gval, -ONE, xold, fv);
  if (iter > 0) {
    /* compute dg_new = gval - gval_old*/
    N_VLinearSum(ONE, gval, -ONE, gold, dg[i_pt]);
    /* compute df_new = fval - fval_old */
    N_VLinearSum(ONE, fv, -ONE, fold, df[i_pt]);
  }
    
  N_VScale(ONE, gval, gold);
  N_VScale(ONE, fv, fold);
  
  if (iter == 0) {
    N_VScale(ONE, gval, x);
  } else {
    if (iter == 1) {
      R[0] = SUNRsqrt(N_VDotProd(df[i_pt], df[i_pt])); 
      alfa = ONE/R[0];
      N_VScale(alfa, df[i_pt], Q[i_pt]);
      ipt_map[0] = 0;
    } else if (iter <= maa) {
      N_VScale(ONE, df[i_pt], vtemp2);
      for (j=0; j<(iter-1); j++) {
	ipt_map[j] = j;
	R[(iter-1)*maa+j] = N_VDotProd(Q[j], vtemp2);
	N_VLinearSum(ONE, vtemp2, -R[(iter-1)*maa+j], Q[j], vtemp2);
      }
      R[(iter-1)*maa+iter-1] = SUNRsqrt(N_VDotProd(vtemp2, vtemp2)); 
      if (R[(iter-1)*maa+iter-1] == ZERO) {
	N_VScale(ZERO, vtemp2, Q[i_pt]);
      } else {
	N_VScale((ONE/R[(iter-1)*maa+iter-1]), vtemp2, Q[i_pt]);
      }
      ipt_map[iter-1] = iter-1;
    } else {
      /* Delete left-most column vector from QR factorization */
      for (i=0; i < maa-1; i++) {
        a = R[(i+1)*maa + i];
        b = R[(i+1)*maa + i+1];
        temp = SUNRsqrt(a*a + b*b);
        c = a / temp;
        s = b / temp;
        R[(i+1)*maa + i] = temp;
        R[(i+1)*maa + i+1] = 0.0;      
	/* OK to re-use temp */
        if (i < maa-1) {
          for (j = i+2; j < maa; j++) {
            a = R[j*maa + i];
            b = R[j*maa + i+1];
            temp = c * a + s * b;
            R[j*maa + i+1] = -s*a + c*b;
            R[j*maa + i] = temp;
	  }
	}
        N_VLinearSum(c, Q[i], s, Q[i+1], vtemp2);
        N_VLinearSum(-s, Q[i], c, Q[i+1], Q[i+1]);
        N_VScale(ONE, vtemp2, Q[i]);
      }

      /* Shift R to the left by one. */
      for (i = 1; i < maa; i++) {
        for (j = 0; j < maa-1; j++) {
          R[(i-1)*maa + j] = R[i*maa + j];
        }
      }

      /* Add the new df vector */
      N_VScale(ONE,df[i_pt],vtemp2);
      for (j=0; j < (maa-1); j++) {
        R[(maa-1)*maa+j] = N_VDotProd(Q[j],vtemp2);
        N_VLinearSum(ONE,vtemp2,-R[(maa-1)*maa+j],Q[j],vtemp2);
      }
      R[(maa-1)*maa+maa-1] = SUNRsqrt(N_VDotProd(vtemp2,vtemp2));
      N_VScale((1/R[(maa-1)*maa+maa-1]),vtemp2,Q[maa-1]);

      /* Update the iteration map */
      j = 0;
      for (i=i_pt+1; i<maa; i++)
	ipt_map[j++] = i;
      for (i=0; i<(i_pt+1); i++)
	ipt_map[j++] = i;
    }

    /* Solve least squares problem and update solution */
    lAA = iter;
    if (maa < iter)  lAA = maa;
    N_VScale(ONE, gval, x);
    for (i=0; i<lAA; i++)
      gamma[i] = N_VDotProd(fv, Q[i]);
    for (i=lAA-1; i>-1; i--) {
      for (j=i+1; j<lAA; j++) 
	gamma[i] = gamma[i] - R[j*maa+i]*gamma[j]; 
      if (gamma[i] == ZERO) {
	gamma[i] = ZERO;
      } else {
	gamma[i] = gamma[i]/R[i*maa+i];
      }
      N_VLinearSum(ONE, x, -gamma[i], dg[ipt_map[i]], x);
    }
  }

  return 0;
}

/*---------------------------------------------------------------
 arkLs

 This routine attempts to solve the linear system associated
 with a single implicit step of the ARK method.  This should 
 only be called if the user has specified that the implicit 
 problem is linear.  In this routine, we assume that the problem 
 depends linearly on the solution.  Additionally, if the Jacobian
 is not time dependent we only call lsetup on changes to gamma; 
 otherwise we call lsetup at every call.  In all cases, we then 
 call the user-specified linear solver (with a tight tolerance) 
 to compute the time-evolved solution.  

 Upon entry, the predicted solution is held in ark_mem->ark_ycur.

 Upon a successful solve, the solution is held in ark_mem->ark_y.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+  
   ARK_LSETUP_FAIL    |-> halt the integration 
   ARK_LSOLVE_FAIL   -+

   RHSFUNC_RECVR     --> predict again or stop if too many
---------------------------------------------------------------*/
static int arkLs(ARKodeMem ark_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3, b;
  int convfail, retval, ier;
  booleantype callSetup;
  realtype del;
  
  vtemp1 = ark_mem->ark_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = ark_mem->ark_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = ark_mem->ark_tempv; /* rename tempv as vtemp3 for readability */
  b      = ark_mem->ark_tempv; /* also rename tempv as b for readability */

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = (nflag == FIRST_CALL) ? ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (ark_mem->ark_setupNonNull) {      
    callSetup = (ark_mem->ark_firststage) || 
      (ark_mem->ark_linear_timedep) || (ark_mem->ark_msbp < 0) ||
      (SUNRabs(ark_mem->ark_gamrat-ONE) > ark_mem->ark_dgmax);
  } else {  
    callSetup = FALSE;
  }
  
  /* update implicit RHS, store in ark_mem->ark_ftemp */
  if (!ark_mem->ark_explicit) {
    retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_ycur, 
			     ark_mem->ark_ftemp, ark_mem->ark_user_data);
    ark_mem->ark_nfi++; 
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);
  }
  
  /* update system matrix if necessary */
  if (callSetup) {

    /* Solver diagnostics reporting */
    if (ark_mem->ark_report)  fprintf(ark_mem->ark_diagfp, "  lsetup\n");

    ier = ark_mem->ark_lsetup(ark_mem, convfail, ark_mem->ark_ycur, 
			      ark_mem->ark_ftemp, &ark_mem->ark_jcur, 
			      vtemp1, vtemp2, vtemp3);
    ark_mem->ark_nsetups++;
    callSetup = FALSE;
    ark_mem->ark_firststage = FALSE;
    ark_mem->ark_gamrat = ark_mem->ark_crate = ONE; 
    ark_mem->ark_gammap = ark_mem->ark_gamma;
    ark_mem->ark_nstlp  = ark_mem->ark_nst;

    /* Return if lsetup failed */
    if (ier < 0) return(ARK_LSETUP_FAIL);
    if (ier > 0) return(CONV_FAIL);
  }

  /* Set acor to zero and load prediction into y vector */
  N_VConst(ZERO, ark_mem->ark_acor);
  N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_y);
  

  /* Do a single Newton iteration */

  /*   Initialize temporary variables for use in iteration */
  ark_mem->ark_mnewt = 0;
  del = ZERO;

  /*   Set the local truncation error estimate to force an "accurate" linear solve */
  ark_mem->ark_eLTE = ark_mem->ark_nlscoef * RCONST(0.1);

  /*   Evaluate the nonlinear system residual, put result into b */
  retval = arkNlsResid(ark_mem, ark_mem->ark_acor, ark_mem->ark_ftemp, b);
  if (retval != ARK_SUCCESS)  return (ARK_RHSFUNC_FAIL);

  /*   Call the lsolve function */
  retval = ark_mem->ark_lsolve(ark_mem, b, ark_mem->ark_rwt, 
			       ark_mem->ark_y, ark_mem->ark_ftemp); 
  ark_mem->ark_nni++;
  if (retval != 0)  return (ARK_LSOLVE_FAIL);
    
  /*   Get WRMS norm of correction; add correction to acor and y */
  del = N_VWrmsNorm(b, ark_mem->ark_ewt);
  N_VLinearSum(ONE, ark_mem->ark_acor, ONE, b, ark_mem->ark_acor);
  N_VLinearSum(ONE, ark_mem->ark_ycur, ONE, ark_mem->ark_acor, ark_mem->ark_y);

  /*   Solver diagnostics reporting */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "    newt  %i  %19.16g  %19.16g\n", 0, del, 0.0);

  /* clean up and return */ 
  ark_mem->ark_jcur = FALSE;
  return (ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkHandleNFlag

 This routine takes action on the return value nflag = *nflagPtr
 returned by arkNls, as follows:

 If arkNls succeeded in solving the nonlinear system, then
 arkHandleNFlag returns the constant SOLVE_SUCCESS, which tells 
 arkStep it is safe to continue with other stage solves, or to 
 perform the error test.

 If the nonlinear system was not solved successfully, then ncfn and
 ncf = *ncfPtr are incremented.

 If the solution of the nonlinear system failed due to an
 unrecoverable failure by setup, we return the value ARK_LSETUP_FAIL.
 
 If it failed due to an unrecoverable failure in solve, then we return
 the value ARK_LSOLVE_FAIL.

 If it failed due to an unrecoverable failure in rhs, then we return
 the value ARK_RHSFUNC_FAIL.

 Otherwise, a recoverable failure occurred when solving the 
 nonlinear system (arkNls returned nflag == CONV_FAIL or RHSFUNC_RECVR). 
 In this case, if using fixed time step sizes, or if ncf is now equal 
 to maxncf, or if |h| = hmin, then we return the value ARK_CONV_FAILURE 
 (if nflag=CONV_FAIL) or ARK_REPTD_RHSFUNC_ERR (if nflag=RHSFUNC_RECVR).
 If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
 PREDICT_AGAIN, telling arkStep to reattempt the step.
---------------------------------------------------------------*/
static int arkHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, 
			  realtype saved_t, int *ncfPtr)
{
  int nflag;
  
  nflag = *nflagPtr;
  
  if (nflag == ARK_SUCCESS) return(SOLVE_SUCCESS);

  /* The nonlinear soln. failed; increment ncfn */
  ark_mem->ark_ncfn++;
  
  /* If fixed time stepping, then return with convergence failure */
  if (ark_mem->ark_fixedstep)    return(ARK_CONV_FAILURE);

  /* Restore the previous step time, tn */
  ark_mem->ark_tn = saved_t;
  
  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == ARK_LSETUP_FAIL)  return(ARK_LSETUP_FAIL);
  if (nflag == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
  if (nflag == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);
  
  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  ark_mem->ark_etamax = ONE;

  /* If we had maxncf failures or |h| = hmin, 
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((SUNRabs(ark_mem->ark_h) <= ark_mem->ark_hmin*ONEPSM) ||
      (*ncfPtr == ark_mem->ark_maxncf)) {
    if (nflag == CONV_FAIL)     return(ARK_CONV_FAILURE);
    if (nflag == RHSFUNC_RECVR) return(ARK_REPTD_RHSFUNC_ERR);    
  }

  /* Reduce step size; return to reattempt the step */
  ark_mem->ark_eta = SUNMAX(ark_mem->ark_etacf,
			 ark_mem->ark_hmin / SUNRabs(ark_mem->ark_h));
  ark_mem->ark_h *= ark_mem->ark_eta;
  ark_mem->ark_next_h = ark_mem->ark_h;
  *nflagPtr = PREV_CONV_FAIL;


  return(PREDICT_AGAIN);
}


/*---------------------------------------------------------------
 arkHandleFailure

 This routine prints error messages for all cases of failure by
 arkHin and arkStep. It returns to ARKode the value that ARKode is 
 to return to the user.
---------------------------------------------------------------*/
static int arkHandleFailure(ARKodeMem ark_mem, int flag)
{

  /* Depending on flag, print error message and return error flag */
  switch (flag) {
  case ARK_ERR_FAILURE: 
    arkProcessError(ark_mem, ARK_ERR_FAILURE, "ARKODE", "ARKode", 
		    MSGARK_ERR_FAILS, ark_mem->ark_tn, ark_mem->ark_h);
    break;
  case ARK_CONV_FAILURE:
    arkProcessError(ark_mem, ARK_CONV_FAILURE, "ARKODE", "ARKode", 
		    MSGARK_CONV_FAILS, ark_mem->ark_tn, ark_mem->ark_h);
    break;
  case ARK_LSETUP_FAIL:
    arkProcessError(ark_mem, ARK_LSETUP_FAIL, "ARKODE", "ARKode", 
		    MSGARK_SETUP_FAILED, ark_mem->ark_tn);
    break;
  case ARK_LSOLVE_FAIL:
    arkProcessError(ark_mem, ARK_LSOLVE_FAIL, "ARKODE", "ARKode", 
		    MSGARK_SOLVE_FAILED, ark_mem->ark_tn);
    break;
  case ARK_RHSFUNC_FAIL:
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "ARKode", 
		    MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    arkProcessError(ark_mem, ARK_UNREC_RHSFUNC_ERR, "ARKODE", "ARKode", 
		    MSGARK_RHSFUNC_UNREC, ark_mem->ark_tn);
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    arkProcessError(ark_mem, ARK_REPTD_RHSFUNC_ERR, "ARKODE", "ARKode", 
		    MSGARK_RHSFUNC_REPTD, ark_mem->ark_tn);
    break;
  case ARK_RTFUNC_FAIL:    
    arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKode", 
		    MSGARK_RTFUNC_FAILED, ark_mem->ark_tn);
    break;
  case ARK_TOO_CLOSE:
    arkProcessError(ark_mem, ARK_TOO_CLOSE, "ARKODE", "ARKode", 
		    MSGARK_TOO_CLOSE);
    break;
  case ARK_MASSSOLVE_FAIL:
    arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE", "ARKode", 
		    MSGARK_MASSSOLVE_FAIL);
    break;
  default:
    return(ARK_SUCCESS);   
  }

  return(flag);
}


/*---------------------------------------------------------------
 arkDenseEval

 This routine evaluates the dense output formula for the method,
 using a cubic Hermite interpolating formula with the data 
 {yn,yp,fn,fp}, where yn,yp are the stored solution values and 
 fn,fp are the sotred RHS function values at the endpoints of the
 last successful time step [tn,tp].  If greater polynomial order 
 than 3 is requested (and the method order allows it, then we 
 can bootstrap up to a 5th-order accurate interpolant.  For lower 
 order interpolants than cubic, we use {yn,yp,fp} for 
 quadratic, {yn,yp} for linear, and {0.5*(yn+yp)} for constant.
 
 Derivatives have reduced accuracy than the dense output function
 itself, losing one order per derivative.  We will provide 
 derivatives up to d = min(5,q).

 The input 'tau' specifies the time at which to return derivative
 information, the formula is 
    t = tn + tau*h,
 where h = tp-tn, i.e. values 0<tau<1 provide interpolation, 
 other values result in extrapolation.
---------------------------------------------------------------*/
static int arkDenseEval(ARKodeMem ark_mem, realtype tau, 
			int d, int order, N_Vector yout)
{
  /* local variables */
  int q;
  realtype h, a0, a1, a2, a3, tau2, tau3;
  h = ark_mem->ark_hold;
  tau2 = tau*tau;
  tau3 = tau*tau2;

  /* determine polynomial order q */
  q = SUNMIN(order, ark_mem->ark_dense_q);   /* respect Set routine  */
  q = SUNMIN(q, ark_mem->ark_q);             /* respect method order */
  q = SUNMAX(q, 0);                          /* respect lower bound  */

  /* check that d is possible */
  if ((d > SUNMIN(5,q)) || (d < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "arkDenseEval", "Requested illegal derivative.");
    return (ARK_ILL_INPUT);
  }

  /* build polynomial based on order */
  switch (q) {

  case(0):    /* constant interpolant, yout = 0.5*(yn+yp) */
    N_VLinearSum(HALF, ark_mem->ark_yold, HALF, ark_mem->ark_ynew, yout);
    break;

  case(1):    /* linear interpolant */
    if (d == 0) {
      a0 = -tau;
      a1 = ONE+tau;
    } else {  /* d=1 */
      a0 = -ONE/h;
      a1 =  ONE/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    break;

  case(2):    /* quadratic interpolant */
    if (d == 0) {
      a0 = tau2;
      a1 = ONE -tau2;
      a2 = h*(tau2 + tau);
    } else if (d == 1) {
      a0 = TWO*tau/h;
      a1 = -TWO*tau/h;
      a2 = (ONE + TWO*tau);
    } else {  /* d == 2 */
      a0 = TWO/h/h;
      a1 = -TWO/h/h;
      a2 = TWO/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    N_VLinearSum(a2, ark_mem->ark_fnew, ONE, yout, yout);
    break;

  case(3):    /* cubic interpolant */
    if (d == 0) {
      a0 = THREE*tau2 + TWO*tau3;
      a1 = ONE - THREE*tau2 - TWO*tau3;
      a2 = h*(tau2 + tau3);
      a3 = h*(tau + TWO*tau2 + tau3);
    } else if (d == 1) {
      a0 = (SIX*tau + SIX*tau2)/h;
      a1 = -(SIX*tau + SIX*tau2)/h;
      a2 = TWO*tau + THREE*tau2;
      a3 = ONE + FOUR*tau + THREE*tau2;
    } else if (d == 2) {
      a0 = (SIX + TWELVE*tau)/h/h;
      a1 = -(SIX + TWELVE*tau)/h/h;
      a2 = (TWO + SIX*tau)/h;
      a3 = (FOUR + SIX*tau)/h;
    } else {  /* d == 3 */
      a0 = TWELVE/h/h/h;
      a1 = -TWELVE/h/h/h;
      a2 = SIX/h/h;
      a3 = SIX/h/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    N_VLinearSum(a2, ark_mem->ark_fold, ONE, yout, yout);
    N_VLinearSum(a3, ark_mem->ark_fnew, ONE, yout, yout);
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkDenseEval", 
		    "Illegal polynomial order.");
    return (ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkFullRHS

 This routine computes the full ODE RHS [fe(t,y) + fi(t,y)].  If
 a non-identity mass matrix is defined for the problem, it fills
 the return array with M^{-1}*[fe(t,y) + fi(t,y)]

 Inputs (unchanged): ark_mem, t, y
 Output: f
 Temporary: tmp
---------------------------------------------------------------*/
static int arkFullRHS(ARKodeMem ark_mem, realtype t, 
		      N_Vector y, N_Vector tmp, N_Vector f)
{
  int retval;

  /* if the problem involves a non-identity mass matrix and setup is
     required, do so here (use f, tmp and ark_sdata as a temporaries) */
  if (ark_mem->ark_mass_matrix && ark_mem->ark_MassSetupNonNull) {
    retval = ark_mem->ark_msetup(ark_mem, f, tmp, ark_mem->ark_sdata);
    if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
  }


  /* explicit problem -- only call fe */
  if (ark_mem->ark_explicit) {

    retval = ark_mem->ark_fe(t, y, f, ark_mem->ark_user_data);
    ark_mem->ark_nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "arkFullRHS", 
		      MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }


  /* implicit problem -- only call fi */
  } else if (ark_mem->ark_implicit) {

    retval = ark_mem->ark_fi(t, y, f, ark_mem->ark_user_data);
    ark_mem->ark_nfi++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", 
		      "arkFullRHS", MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }


  /* imex problem */
  } else {
    
    /* explicit portion */
    retval = ark_mem->ark_fe(t, y, f, ark_mem->ark_user_data);
    ark_mem->ark_nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "arkFullRHS", 
		      MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }

    /* implicit portion */
    retval = ark_mem->ark_fi(t, y, tmp, ark_mem->ark_user_data);
    ark_mem->ark_nfi++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", 
		      "arkFullRHS", MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }
    N_VLinearSum(ONE, tmp, ONE, f, f);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkEwtSetSS

 This routine sets ewt as decribed above in the case tol_type = ARK_SS.
 It tests for non-positive components before inverting. arkEwtSetSS
 returns 0 if ewt is successfully set to a positive vector
 and -1 otherwise. In the latter case, ewt is considered undefined.
---------------------------------------------------------------*/
static int arkEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, ark_mem->ark_tempv);
  N_VScale(ark_mem->ark_reltol, ark_mem->ark_tempv, ark_mem->ark_tempv);
  N_VAddConst(ark_mem->ark_tempv, ark_mem->ark_Sabstol, ark_mem->ark_tempv);
  if (N_VMin(ark_mem->ark_tempv) <= ZERO) return(-1);
  N_VInv(ark_mem->ark_tempv, weight);
  return(0);
}


/*---------------------------------------------------------------
 arkEwtSetSV

 This routine sets ewt as decribed above in the case tol_type = ARK_SV.
 It tests for non-positive components before inverting. arkEwtSetSV
 returns 0 if ewt is successfully set to a positive vector
 and -1 otherwise. In the latter case, ewt is considered undefined.
---------------------------------------------------------------*/
static int arkEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, ark_mem->ark_tempv);
  N_VLinearSum(ark_mem->ark_reltol, ark_mem->ark_tempv, ONE, 
	       ark_mem->ark_Vabstol, ark_mem->ark_tempv);
  if (N_VMin(ark_mem->ark_tempv) <= ZERO) return(-1);
  N_VInv(ark_mem->ark_tempv, weight);
  return(0);
}


/*---------------------------------------------------------------
 arkRwtSetSS

 This routine sets rwt as decribed above in the case tol_type = ARK_SS.
 It tests for non-positive components before inverting. arkRwtSetSS
 returns 0 if rwt is successfully set to a positive vector
 and -1 otherwise. In the latter case, rwt is considered undefined.
---------------------------------------------------------------*/
static int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->ark_tempv);
  N_VScale(ark_mem->ark_reltol, ark_mem->ark_tempv, ark_mem->ark_tempv);
  N_VAddConst(ark_mem->ark_tempv, ark_mem->ark_SRabstol, ark_mem->ark_tempv);
  if (N_VMin(ark_mem->ark_tempv) <= ZERO) return(-1);
  N_VInv(ark_mem->ark_tempv, weight);
  return(0);
}


/*---------------------------------------------------------------
 arkRwtSetSV

 This routine sets rwt as decribed above in the case tol_type = ARK_SV.
 It tests for non-positive components before inverting. arkRwtSetSV
 returns 0 if rwt is successfully set to a positive vector
 and -1 otherwise. In the latter case, rwt is considered undefined.
---------------------------------------------------------------*/
static int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->ark_tempv);
  N_VLinearSum(ark_mem->ark_reltol, ark_mem->ark_tempv, ONE, 
	       ark_mem->ark_VRabstol, ark_mem->ark_tempv);
  if (N_VMin(ark_mem->ark_tempv) <= ZERO) return(-1);
  N_VInv(ark_mem->ark_tempv, weight);
  return(0);
}


/*---------------------------------------------------------------
 arkAdapt is the time step adaptivity wrapper function.  This 
 should be called after arkCompleteStep, as it depends on the 
 values updated by that routine.  It computes and sets the value 
 of ark_eta inside ark_mem.
---------------------------------------------------------------*/
static int arkAdapt(ARKodeMem ark_mem)
{
  int ier;
  realtype h_acc, h_cfl, safety, int_dir;
  safety = ark_mem->ark_hadapt_safety;

  /* Call algorithm-specific error adaptivity method */
  switch (ark_mem->ark_hadapt_imethod) {
  case(0):    /* PID controller */
    ier = arkAdaptPID(ark_mem, &h_acc);
    break;
  case(1):    /* PI controller */
    ier = arkAdaptPI(ark_mem, &h_acc);
    break;
  case(2):    /* I controller */
    ier = arkAdaptI(ark_mem, &h_acc);
    break;
  case(3):    /* explicit Gustafsson controller */
    ier = arkAdaptExpGus(ark_mem, &h_acc);
    break;
  case(4):    /* implicit Gustafsson controller */
    ier = arkAdaptImpGus(ark_mem, &h_acc);
    break;
  case(5):    /* imex Gustafsson controller */
    ier = arkAdaptImExGus(ark_mem, &h_acc);
    break;
  case(-1):   /* user-supplied controller */
    ier = ark_mem->ark_hadapt(ark_mem->ark_ycur,
			      ark_mem->ark_tn, 
			      ark_mem->ark_hadapt_hhist[0], 
			      ark_mem->ark_hadapt_hhist[1], 
			      ark_mem->ark_hadapt_hhist[2], 
			      ark_mem->ark_hadapt_ehist[0],
			      ark_mem->ark_hadapt_ehist[1],
			      ark_mem->ark_hadapt_ehist[2],
			      ark_mem->ark_q, ark_mem->ark_p, 
			      &h_acc, ark_mem->ark_hadapt_data);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkAdapt", 
		    "Illegal adapt_imethod.");
    return (ARK_ILL_INPUT);
  }
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkAdapt", 
		    "Error in accuracy-based adaptivity function.");
    return (ARK_ILL_INPUT);
  }

  /* determine direction of integration */
  int_dir = ark_mem->ark_h / SUNRabs(ark_mem->ark_h);

  /* Call explicit stability function */
  ier = ark_mem->ark_expstab(ark_mem->ark_ycur, ark_mem->ark_tn,
			     &h_cfl, ark_mem->ark_estab_data);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkAdapt", 
		    "Error in explicit stability function.");
    return (ARK_ILL_INPUT);
  }
  if (h_cfl <= 0.0)  h_cfl = RCONST(1.0e30) * SUNRabs(ark_mem->ark_h);

  /* Solver diagnostics reporting */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "  adapt  %19.16g  %19.16g  %19.16g  %19.16g  %19.16g  %19.16g  %19.16g  %19.16g  ",
	    ark_mem->ark_hadapt_ehist[0], ark_mem->ark_hadapt_ehist[1], 
	    ark_mem->ark_hadapt_ehist[2], ark_mem->ark_hadapt_hhist[0], 
	    ark_mem->ark_hadapt_hhist[1], ark_mem->ark_hadapt_hhist[2], h_acc, h_cfl);

  /* enforce safety factors */
  h_acc *= safety;
  h_cfl *= ark_mem->ark_hadapt_cfl * int_dir;

  /* enforce maximum bound on time step growth */
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(ark_mem->ark_etamax*ark_mem->ark_h));

  /* enforce minimum bound time step reduction */
  h_acc = int_dir * SUNMAX(SUNRabs(h_acc), SUNRabs(ETAMIN*ark_mem->ark_h));

  /* Solver diagnostics reporting */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "%19.16g  %19.16g  ", h_acc, h_cfl);

  /* increment the relevant step counter, set desired step */
  if (SUNRabs(h_acc) < SUNRabs(h_cfl))
    ark_mem->ark_nst_acc++;
  else
    ark_mem->ark_nst_exp++;
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(h_cfl));

  /* enforce adaptivity bounds to retain Jacobian/preconditioner accuracy */
  if ( (SUNRabs(h_acc) > SUNRabs(ark_mem->ark_h*ark_mem->ark_hadapt_lbound*ONEMSM)) &&
       (SUNRabs(h_acc) < SUNRabs(ark_mem->ark_h*ark_mem->ark_hadapt_ubound*ONEPSM)) )
    h_acc = ark_mem->ark_h;

  /* set basic value of ark_eta */
  ark_mem->ark_eta = h_acc / ark_mem->ark_h;

  /* enforce minimum time step size */
  ark_mem->ark_eta = SUNMAX(ark_mem->ark_eta,
			 ark_mem->ark_hmin / SUNRabs(ark_mem->ark_h));

  /* enforce maximum time step size */
  ark_mem->ark_eta /= SUNMAX(ONE, SUNRabs(ark_mem->ark_h) *
			  ark_mem->ark_hmax_inv*ark_mem->ark_eta);

  /* Solver diagnostics reporting */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "%19.16g\n", ark_mem->ark_eta);

  return(ier);
}


/*---------------------------------------------------------------
 arkAdaptPID implements a PID time step control algorithm.
---------------------------------------------------------------*/
static int arkAdaptPID(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k, k1, k2, k3, e1, e2, e3, hcur, h_acc;

  /* set usable time-step adaptivity parameters */
  k = (ark_mem->ark_hadapt_pq) ? ark_mem->ark_q : ark_mem->ark_p;
  k1 = -ark_mem->ark_hadapt_k1 / k;
  k2 =  ark_mem->ark_hadapt_k2 / k;
  k3 = -ark_mem->ark_hadapt_k3 / k;
  e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
  e2 = SUNMAX(ark_mem->ark_hadapt_ehist[1], TINY);
  e3 = SUNMAX(ark_mem->ark_hadapt_ehist[2], TINY);
  hcur = ark_mem->ark_h;
  
  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2) * SUNRpowerR(e3,k3);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkAdaptPI implements a PI time step control algorithm.
---------------------------------------------------------------*/
static int arkAdaptPI(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k, k1, k2, e1, e2, hcur, h_acc;

  /* set usable time-step adaptivity parameters */
  k = (ark_mem->ark_hadapt_pq) ? ark_mem->ark_q : ark_mem->ark_p;
  k1 = -ark_mem->ark_hadapt_k1 / k;
  k2 =  ark_mem->ark_hadapt_k2 / k;
  e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
  e2 = SUNMAX(ark_mem->ark_hadapt_ehist[1], TINY);
  hcur = ark_mem->ark_h;
  
  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkAdaptI implements an I time step control algorithm.
---------------------------------------------------------------*/
static int arkAdaptI(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k, k1, e1, hcur, h_acc;

  /* set usable time-step adaptivity parameters */
  k = (ark_mem->ark_hadapt_pq) ? ark_mem->ark_q : ark_mem->ark_p;
  k1 = -ark_mem->ark_hadapt_k1 / k;
  e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
  hcur = ark_mem->ark_h;
  
  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * SUNRpowerR(e1,k1);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkAdaptExpGus implements the explicit Gustafsson time step
 control algorithm.
---------------------------------------------------------------*/
static int arkAdaptExpGus(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k, k1, k2, e1, e2, hcur, h_acc;
  k = (ark_mem->ark_hadapt_pq) ? ark_mem->ark_q : ark_mem->ark_p;

  /* modified method for first step */
  if (ark_mem->ark_nst < 2) {

    k1 = -ONE / k;
    hcur = ark_mem->ark_h;
    e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -ark_mem->ark_hadapt_k1 / k;
    k2 = -ark_mem->ark_hadapt_k2 / k;
    e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
    e2 = e1 / SUNMAX(ark_mem->ark_hadapt_ehist[1], TINY);
    hcur = ark_mem->ark_h;
    h_acc = hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkAdaptImpGus implements the implicit Gustafsson time step 
 control algorithm.
---------------------------------------------------------------*/
static int arkAdaptImpGus(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k, k1, k2, e1, e2, hcur, hrat, h_acc;
  k = (ark_mem->ark_hadapt_pq) ? ark_mem->ark_q : ark_mem->ark_p;

  /* modified method for first step */
  if (ark_mem->ark_nst < 2) {

    k1 = -ONE / k;
    hcur = ark_mem->ark_h;
    e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -ark_mem->ark_hadapt_k1 / k;
    k2 = -ark_mem->ark_hadapt_k2 / k;
    e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
    e2 = e1 / SUNMAX(ark_mem->ark_hadapt_ehist[1], TINY);
    hcur = ark_mem->ark_hadapt_hhist[0];
    hrat = hcur / ark_mem->ark_hadapt_hhist[1];
    h_acc = hcur * hrat * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkAdaptImExGus implements a combination implicit/explicit 
 Gustafsson time step control algorithm.
---------------------------------------------------------------*/
static int arkAdaptImExGus(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k, k1, k2, k3, e1, e2, hcur, hrat, h_acc;
  k = (ark_mem->ark_hadapt_pq) ? ark_mem->ark_q : ark_mem->ark_p;

  /* modified method for first step */
  if (ark_mem->ark_nst < 2) {

    k1 = -ONE / k;
    hcur = ark_mem->ark_h;
    e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -ark_mem->ark_hadapt_k1 / k;
    k2 = -ark_mem->ark_hadapt_k2 / k;
    k3 = -ark_mem->ark_hadapt_k3 / k;
    e1 = SUNMAX(ark_mem->ark_hadapt_ehist[0], TINY);
    e2 = e1 / SUNMAX(ark_mem->ark_hadapt_ehist[1], TINY);
    hcur = ark_mem->ark_hadapt_hhist[0];
    hrat = hcur / ark_mem->ark_hadapt_hhist[1];
    /* implicit estimate */
    h_acc = hcur * hrat * SUNRpowerR(e1,k3) * SUNRpowerR(e2,k3);
    /* explicit estimate */
    h_acc = SUNMIN(h_acc, hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2));

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkExpStab is the default explicit stability estimation function
---------------------------------------------------------------*/
int arkExpStab(N_Vector y, realtype t, realtype *hstab, void *data)
{
  /* explicit stability not used by default, 
     set to zero to disable */
  *hstab = RCONST(0.0);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkProcessError is a high level error handling function
 - if ark_mem==NULL it prints the error message to stderr
 - otherwise, it sets-up and calls the error handling function 
   pointed to by ark_ehfun 
---------------------------------------------------------------*/
void arkProcessError(ARKodeMem ark_mem, int error_code, 
		     const char *module, const char *fname, 
		     const char *msgfmt, ...)
{
  va_list ap;
  char msg[256];

  /* Initialize the argument pointer variable 
     (msgfmt is the last required argument to arkProcessError) */
  va_start(ap, msgfmt);

  /* Compose the message */
  vsprintf(msg, msgfmt, ap);

  if (ark_mem == NULL) {    /* We write to stderr */

#ifndef NO_FPRINTF_OUTPUT
    fprintf(stderr, "\n[%s ERROR]  %s\n  ", module, fname);
    fprintf(stderr, "%s\n\n", msg);
#endif

  } else {                 /* We can call ehfun */
    ark_mem->ark_ehfun(error_code, module, fname, msg, 
		       ark_mem->ark_eh_data);
  }

  /* Finalize argument processing */
  va_end(ap);

  return;
}


/*---------------------------------------------------------------
 arkRootCheck1

 This routine completes the initialization of rootfinding memory
 information, and checks whether g has a zero both at and very near
 the initial point of the IVP.

 This routine returns an int equal to:
  ARK_RTFUNC_FAIL < 0  if the g function failed, or
  ARK_SUCCESS     = 0  otherwise.
---------------------------------------------------------------*/
static int arkRootCheck1(ARKodeMem ark_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  for (i = 0; i < ark_mem->ark_nrtfn; i++) 
    ark_mem->ark_iroots[i] = 0;
  ark_mem->ark_tlo = ark_mem->ark_tn;
  ark_mem->ark_ttol = (SUNRabs(ark_mem->ark_tn) +
		       SUNRabs(ark_mem->ark_h))*ark_mem->ark_uround*HUND;

  /* Evaluate g at initial t and check for zero values. */
  retval = ark_mem->ark_gfun(ark_mem->ark_tlo, ark_mem->ark_ycur,
			     ark_mem->ark_glo, ark_mem->ark_user_data);
  ark_mem->ark_nge = 1;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (SUNRabs(ark_mem->ark_glo[i]) == ZERO) {
      zroot = TRUE;
      ark_mem->ark_gactive[i] = FALSE;
    }
  }
  if (!zroot) return(ARK_SUCCESS);

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  hratio = SUNMAX(ark_mem->ark_ttol/SUNRabs(ark_mem->ark_h), TENTH);
  smallh = hratio*ark_mem->ark_h;
  tplus = ark_mem->ark_tlo + smallh;
  N_VLinearSum(ONE, ark_mem->ark_ycur, smallh,
	       ark_mem->ark_fold, ark_mem->ark_y);
  retval = ark_mem->ark_gfun(tplus, ark_mem->ark_y, ark_mem->ark_ghi, 
			     ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* We check now only the components of g which were exactly 0.0 at t0
   * to see if we can 'activate' them. */
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i] && SUNRabs(ark_mem->ark_ghi[i]) != ZERO) {
      ark_mem->ark_gactive[i] = TRUE;
      ark_mem->ark_glo[i] = ark_mem->ark_ghi[i];
    }
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkRootCheck2

 This routine checks for exact zeros of g at the last root found,
 if the last return was a root.  It then checks for a close pair of
 zeros (an error condition), and for a new root at a nearby point.
 The array glo = g(tlo) at the left endpoint of the search interval
 is adjusted if necessary to assure that all g_i are nonzero
 there, before returning to do a root search in the interval.

 On entry, tlo = tretlast is the last value of tret returned by
 ARKode.  This may be the previous tn, the previous tout value, or
 the last root location.

 This routine returns an int equal to:
      ARK_RTFUNC_FAIL < 0 if the g function failed, or
      CLOSERT         = 3 if a close pair of zeros was found, or
      RTFOUND         = 1 if a new zero of g was found near tlo, or
      ARK_SUCCESS     = 0 otherwise.
---------------------------------------------------------------*/
static int arkRootCheck2(ARKodeMem ark_mem)
{
  int i, retval;
  /* realtype smallh, hratio, tplus; */
  realtype smallh, tplus;
  booleantype zroot;

  /* return if no roots in previous step */
  if (ark_mem->ark_irfnd == 0) return(ARK_SUCCESS);

  /* Set ark_y = y(tlo) */
  (void) ARKodeGetDky(ark_mem, ark_mem->ark_tlo, 0, ark_mem->ark_y);

  /* Evaluate root-finding function: glo = g(tlo, y(tlo)) */
  retval = ark_mem->ark_gfun(ark_mem->ark_tlo, ark_mem->ark_y, 
			     ark_mem->ark_glo, ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* reset root-finding flags (overall, and for specific eqns) */
  zroot = FALSE;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) 
    ark_mem->ark_iroots[i] = 0;

  /* for all active roots, check if glo_i == 0 to mark roots found */
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i]) continue;
    if (SUNRabs(ark_mem->ark_glo[i]) == ZERO) {
      zroot = TRUE;
      ark_mem->ark_iroots[i] = 1;
    }
  }
  if (!zroot) return(ARK_SUCCESS);  /* return if no roots */

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  /*     set time tolerance */
  ark_mem->ark_ttol = (SUNRabs(ark_mem->ark_tn) +
		       SUNRabs(ark_mem->ark_h))*ark_mem->ark_uround*HUND;
  /*     set tplus = tlo + smallh */
  smallh = (ark_mem->ark_h > ZERO) ? ark_mem->ark_ttol : -ark_mem->ark_ttol;
  tplus = ark_mem->ark_tlo + smallh;
  /*     update ark_y with small explicit Euler step (if tplus is past tn) */
  if ( (tplus - ark_mem->ark_tn)*ark_mem->ark_h >= ZERO ) {
    /* hratio = smallh/ark_mem->ark_h; */
    N_VLinearSum(ONE, ark_mem->ark_y, smallh, 
		 ark_mem->ark_fold, ark_mem->ark_y);
  } else {
    /*   set ark_y = y(tplus) via interpolation */
    (void) ARKodeGetDky(ark_mem, tplus, 0, ark_mem->ark_y);
  }
  /*     set ghi = g(tplus,y(tplus)) */
  retval = ark_mem->ark_gfun(tplus, ark_mem->ark_y, ark_mem->ark_ghi, 
			     ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* Check for close roots (error return), for a new zero at tlo+smallh,
  and for a g_i that changed from zero to nonzero. */
  zroot = FALSE;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i]) continue;
    if (SUNRabs(ark_mem->ark_ghi[i]) == ZERO) {
      if (ark_mem->ark_iroots[i] == 1) return(CLOSERT);
      zroot = TRUE;
      ark_mem->ark_iroots[i] = 1;
    } else {
      if (ark_mem->ark_iroots[i] == 1) 
	ark_mem->ark_glo[i] = ark_mem->ark_ghi[i];
    }
  }
  if (zroot) return(RTFOUND);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkRootCheck3

 This routine interfaces to arkRootfind to look for a root of g
 between tlo and either tn or tout, whichever comes first.
 Only roots beyond tlo in the direction of integration are sought.

 This routine returns an int equal to:
      ARK_RTFUNC_FAIL < 0 if the g function failed, or
      RTFOUND         = 1 if a root of g was found, or
      ARK_SUCCESS     = 0 otherwise.
---------------------------------------------------------------*/
static int arkRootCheck3(ARKodeMem ark_mem)
{
  int i, retval, ier;

  /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
  if (ark_mem->ark_taskc == ARK_ONE_STEP) {
    ark_mem->ark_thi = ark_mem->ark_tn;
    N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_y);
  }
  if (ark_mem->ark_taskc == ARK_NORMAL) {
    if ( (ark_mem->ark_toutc - ark_mem->ark_tn)*ark_mem->ark_h >= ZERO) {
      ark_mem->ark_thi = ark_mem->ark_tn; 
      N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_y);
    } else {
      ark_mem->ark_thi = ark_mem->ark_toutc;
      (void) ARKodeGetDky(ark_mem, ark_mem->ark_thi, 0, ark_mem->ark_y);
    }
  }

  /* Set ark_mem->ark_ghi = g(thi) and call arkRootfind to search (tlo,thi) for roots. */
  retval = ark_mem->ark_gfun(ark_mem->ark_thi, ark_mem->ark_y, 
			     ark_mem->ark_ghi, ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  ark_mem->ark_ttol = (SUNRabs(ark_mem->ark_tn) +
		       SUNRabs(ark_mem->ark_h))*ark_mem->ark_uround*HUND;
  ier = arkRootfind(ark_mem);
  if (ier == ARK_RTFUNC_FAIL) return(ARK_RTFUNC_FAIL);
  for(i=0; i<ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i] && ark_mem->ark_grout[i] != ZERO) 
      ark_mem->ark_gactive[i] = TRUE;
  }
  ark_mem->ark_tlo = ark_mem->ark_trout;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) 
    ark_mem->ark_glo[i] = ark_mem->ark_grout[i];

  /* If no root found, return ARK_SUCCESS. */  
  if (ier == ARK_SUCCESS) return(ARK_SUCCESS);

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void) ARKodeGetDky(ark_mem, ark_mem->ark_trout, 0, ark_mem->ark_y);
  return(RTFOUND);

}


/*---------------------------------------------------------------
 arkRootfind

 This routine solves for a root of g(t) between tlo and thi, if
 one exists.  Only roots of odd multiplicity (i.e. with a change
 of sign in one of the g_i), or exact zeros, are found.
 Here the sign of tlo - thi is arbitrary, but if multiple roots
 are found, the one closest to tlo is returned.

 The method used is the Illinois algorithm, a modified secant method.
 Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 Defined Output Points for Solutions of ODEs, Sandia National
 Laboratory Report SAND80-0180, February 1980.

 This routine uses the following parameters for communication:

 nrtfn    = number of functions g_i, or number of components of
            the vector-valued function g(t).  Input only.

 gfun     = user-defined function for g(t).  Its form is
            (void) gfun(t, y, gt, user_data)

 rootdir  = in array specifying the direction of zero-crossings.
            If rootdir[i] > 0, search for roots of g_i only if
            g_i is increasing; if rootdir[i] < 0, search for
            roots of g_i only if g_i is decreasing; otherwise
            always search for roots of g_i.

 gactive  = array specifying whether a component of g should
            or should not be monitored. gactive[i] is initially
            set to TRUE for all i=0,...,nrtfn-1, but it may be
            reset to FALSE if at the first step g[i] is 0.0
            both at the I.C. and at a small perturbation of them.
            gactive[i] is then set back on TRUE only after the 
            corresponding g function moves away from 0.0.

 nge      = cumulative counter for gfun calls.

 ttol     = a convergence tolerance for trout.  Input only.
            When a root at trout is found, it is located only to
            within a tolerance of ttol.  Typically, ttol should
            be set to a value on the order of
               100 * UROUND * max (SUNRabs(tlo), SUNRabs(thi))
            where UROUND is the unit roundoff of the machine.

 tlo, thi = endpoints of the interval in which roots are sought.
            On input, and must be distinct, but tlo - thi may
            be of either sign.  The direction of integration is
            assumed to be from tlo to thi.  On return, tlo and thi
            are the endpoints of the final relevant interval.

 glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
            and g(thi) respectively.  Input and output.  On input,
            none of the glo[i] should be zero.

 trout    = root location, if a root was found, or thi if not.
            Output only.  If a root was found other than an exact
            zero of g, trout is the endpoint thi of the final
            interval bracketing the root, with size at most ttol.

 grout    = array of length nrtfn containing g(trout) on return.

 iroots   = int array of length nrtfn with root information.
            Output only.  If a root was found, iroots indicates
            which components g_i have a root at trout.  For
            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
            and g_i is increasing, iroots[i] = -1 if g_i has a
            root and g_i is decreasing, and iroots[i] = 0 if g_i
            has no roots or g_i varies in the direction opposite
            to that indicated by rootdir[i].

 This routine returns an int equal to:
      ARK_RTFUNC_FAIL < 0 if the g function failed, or
      RTFOUND         = 1 if a root of g was found, or
      ARK_SUCCESS     = 0 otherwise.
---------------------------------------------------------------*/
static int arkRootfind(ARKodeMem ark_mem)
{
  realtype alpha, tmid, gfrac, maxfrac, fracint, fracsub;
  int i, retval, imax, side, sideprev;
  booleantype zroot, sgnchg;

  imax = 0;

  /* First check for change in sign in ghi or for a zero in ghi. */
  maxfrac = ZERO;
  zroot = FALSE;
  sgnchg = FALSE;
  for (i = 0;  i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i]) continue;
    if (SUNRabs(ark_mem->ark_ghi[i]) == ZERO) {
      if (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) {
        zroot = TRUE;
      }
    } else {
      if ( (ark_mem->ark_glo[i]*ark_mem->ark_ghi[i] < ZERO) && 
	   (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) {
        gfrac = SUNRabs(ark_mem->ark_ghi[i]/(ark_mem->ark_ghi[i] - ark_mem->ark_glo[i]));
        if (gfrac > maxfrac) {
          sgnchg = TRUE;
          maxfrac = gfrac;
          imax = i;
        }
      }
    }
  }

  /* If no sign change was found, reset trout and grout.  Then return
     ARK_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */ 
  if (!sgnchg) {
    ark_mem->ark_trout = ark_mem->ark_thi;
    for (i = 0; i < ark_mem->ark_nrtfn; i++) 
      ark_mem->ark_grout[i] = ark_mem->ark_ghi[i];
    if (!zroot) return(ARK_SUCCESS);
    for (i = 0; i < ark_mem->ark_nrtfn; i++) {
      ark_mem->ark_iroots[i] = 0;
      if (!ark_mem->ark_gactive[i]) continue;
      if (SUNRabs(ark_mem->ark_ghi[i]) == ZERO)
	ark_mem->ark_iroots[i] = ark_mem->ark_glo[i] > 0 ? -1:1;
    }
    return(RTFOUND);
  }

  /* Initialize alpha to avoid compiler warning */
  alpha = ONE;

  /* A sign change was found.  Loop to locate nearest root. */
  side = 0;  sideprev = -1;
  for(;;) {                                    /* Looping point */

    /* If interval size is already less than tolerance ttol, break. */
    if (SUNRabs(ark_mem->ark_thi - ark_mem->ark_tlo) <= ark_mem->ark_ttol) break;
    
    /* Set weight alpha.
       On the first two passes, set alpha = 1.  Thereafter, reset alpha
       according to the side (low vs high) of the subinterval in which
       the sign change was found in the previous two passes.
       If the sides were opposite, set alpha = 1.
       If the sides were the same, then double alpha (if high side),
       or halve alpha (if low side).
       The next guess tmid is the secant method value if alpha = 1, but
       is closer to tlo if alpha < 1, and closer to thi if alpha > 1.    */
    if (sideprev == side) {
      alpha = (side == 2) ? alpha*TWO : alpha*HALF;
    } else {
      alpha = ONE;
    }

    /* Set next root approximation tmid and get g(tmid).
       If tmid is too close to tlo or thi, adjust it inward,
       by a fractional distance that is between 0.1 and 0.5.  */
    tmid = ark_mem->ark_thi - (ark_mem->ark_thi - ark_mem->ark_tlo) *
      ark_mem->ark_ghi[imax]/(ark_mem->ark_ghi[imax] - alpha*ark_mem->ark_glo[imax]);
    if (SUNRabs(tmid - ark_mem->ark_tlo) < HALF*ark_mem->ark_ttol) {
      fracint = SUNRabs(ark_mem->ark_thi - ark_mem->ark_tlo)/ark_mem->ark_ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = ark_mem->ark_tlo + fracsub*(ark_mem->ark_thi - ark_mem->ark_tlo);
    }
    if (SUNRabs(ark_mem->ark_thi - tmid) < HALF*ark_mem->ark_ttol) {
      fracint = SUNRabs(ark_mem->ark_thi - ark_mem->ark_tlo)/ark_mem->ark_ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = ark_mem->ark_thi - fracsub*(ark_mem->ark_thi - ark_mem->ark_tlo);
    }

    (void) ARKodeGetDky(ark_mem, tmid, 0, ark_mem->ark_y);
    retval = ark_mem->ark_gfun(tmid, ark_mem->ark_y, ark_mem->ark_grout, 
			       ark_mem->ark_user_data);
    ark_mem->ark_nge++;
    if (retval != 0) return(ARK_RTFUNC_FAIL);

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */  
    maxfrac = ZERO;
    zroot = FALSE;
    sgnchg = FALSE;
    sideprev = side;
    for (i = 0;  i < ark_mem->ark_nrtfn; i++) {
      if (!ark_mem->ark_gactive[i]) continue;
      if (SUNRabs(ark_mem->ark_grout[i]) == ZERO) {
        if (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) {
          zroot = TRUE;
        }
      } else {
        if ( (ark_mem->ark_glo[i]*ark_mem->ark_grout[i] < ZERO) && 
	     (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) {
          gfrac = SUNRabs(ark_mem->ark_grout[i]/(ark_mem->ark_grout[i] - ark_mem->ark_glo[i]));
          if (gfrac > maxfrac) {
            sgnchg = TRUE;
            maxfrac = gfrac;
            imax = i;
          }
        }
      }
    }
    if (sgnchg) {
      /* Sign change found in (tlo,tmid); replace thi with tmid. */
      ark_mem->ark_thi = tmid;
      for (i = 0; i < ark_mem->ark_nrtfn; i++) 
	ark_mem->ark_ghi[i] = ark_mem->ark_grout[i];
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (SUNRabs(ark_mem->ark_thi - ark_mem->ark_tlo) <= ark_mem->ark_ttol) break;
      continue;  /* Return to looping point. */
    }

    if (zroot) {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      ark_mem->ark_thi = tmid;
      for (i = 0; i < ark_mem->ark_nrtfn; i++) 
	ark_mem->ark_ghi[i] = ark_mem->ark_grout[i];
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    ark_mem->ark_tlo = tmid;
    for (i = 0; i < ark_mem->ark_nrtfn; i++) 
      ark_mem->ark_glo[i] = ark_mem->ark_grout[i];
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (SUNRabs(ark_mem->ark_thi - ark_mem->ark_tlo) <= ark_mem->ark_ttol)
      break;

  } /* End of root-search loop */

  /* Reset trout and grout, set iroots, and return RTFOUND. */
  ark_mem->ark_trout = ark_mem->ark_thi;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    ark_mem->ark_grout[i] = ark_mem->ark_ghi[i];
    ark_mem->ark_iroots[i] = 0;
    if (!ark_mem->ark_gactive[i]) continue;
    if ( (SUNRabs(ark_mem->ark_ghi[i]) == ZERO) &&
	 (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) 
      ark_mem->ark_iroots[i] = ark_mem->ark_glo[i] > 0 ? -1:1;
    if ( (ark_mem->ark_glo[i]*ark_mem->ark_ghi[i] < ZERO) && 
	 (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) 
      ark_mem->ark_iroots[i] = ark_mem->ark_glo[i] > 0 ? -1:1;
  }
  return(RTFOUND);
}


/*===============================================================
   EOF
===============================================================*/
