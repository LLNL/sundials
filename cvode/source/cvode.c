/*
 * -----------------------------------------------------------------
 * $Revision: 1.42 $
 * $Date: 2005-04-05 01:59:46 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Dan Shumaker @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the main CVODE integrator.
 * It is independent of the CVODE linear solver in use.
 * -----------------------------------------------------------------
 */

/*
 * Header files
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "sundialsmath.h"
#include "sundialstypes.h"

/* 
 * Macro: loop 
*/

#define loop for(;;)

/*
 * CVODE Private Constants
 */

#define ZERO   RCONST(0.0)   /* real 0.0   */
#define TINY   RCONST(1.0e-10) /* small number */
#define TENTH  RCONST(0.1)   /* real 0.1   */
#define FOURTH RCONST(0.25)  /* real 0.25  */
#define HALF   RCONST(0.5)   /* real 0.5   */
#define ONE    RCONST(1.0)   /* real 1.0   */
#define TWO    RCONST(2.0)   /* real 2.0   */
#define THREE  RCONST(3.0)   /* real 3.0   */
#define FOUR   RCONST(4.0)   /* real 4.0   */
#define FIVE   RCONST(5.0)   /* real 5.0   */
#define TWELVE RCONST(12.0)  /* real 12.0  */
#define HUN    RCONST(100.0) /* real 100.0 */

/*
 * Routine-Specific Constants
 */

/* CVodeGetDky */

#define FUZZ_FACTOR RCONST(100.0)

/* CVHin */

#define HLB_FACTOR RCONST(100.0)
#define HUB_FACTOR RCONST(0.1)
#define H_BIAS     HALF
#define MAX_ITERS  4

/* CVCreate */

#define CORTES RCONST(0.1)

/* CVStep return values */

#define SUCCESS_STEP      0
#define REP_ERR_FAIL     -1
#define REP_CONV_FAIL    -2
#define SETUP_FAILED     -3
#define SOLVE_FAILED     -4

/* CVStep control constants */

#define PREDICT_AGAIN    -5
#define DO_ERROR_TEST     1

/* CVStep */

#define THRESH RCONST(1.5)
#define ETAMX1 RCONST(10000.0) 
#define ETAMX2 RCONST(10.0)
#define ETAMX3 RCONST(10.0)
#define ETAMXF RCONST(0.2)
#define ETAMIN RCONST(0.1)
#define ETACF  RCONST(0.25)
#define ADDON  RCONST(0.000001)
#define BIAS1  RCONST(6.0)
#define BIAS2  RCONST(6.0)
#define BIAS3  RCONST(10.0)
#define ONEPSM RCONST(1.000001)

#define SMALL_NST    10   /* nst > SMALL_NST => use ETAMX3          */
#define MXNCF        10   /* max no. of convergence failures during */
                          /* one step try                           */
#define MXNEF         7   /* max no. of error test failures during  */
                          /* one step try                           */
#define MXNEF1        3   /* max no. of error test failures before  */
                          /* forcing a reduction of order           */
#define SMALL_NEF     2   /* if an error failure occurs and         */
                          /* SMALL_NEF <= nef <= MXNEF1, then       */
                          /* reset eta =  MIN(eta, ETAMXF)          */
#define LONG_WAIT    10   /* number of steps to wait before         */
                          /* considering an order change when       */
                          /* q==1 and MXNEF1 error test failures    */
                          /* have occurred                          */

/* CVnls return values */

#define SOLVED            0
#define CONV_FAIL        -1 
#define SETUP_FAIL_UNREC -2
#define SOLVE_FAIL_UNREC -3

/* CVnls input flags */

#define FIRST_CALL      0
#define PREV_CONV_FAIL -1
#define PREV_ERR_FAIL  -2

/* CVnls other constants */

#define NLS_MAXCOR 3       /* maximum no. of corrector iterations for
                              the nonlinear solver                     */

#define CRDOWN RCONST(0.3) /* constant used in the estimation of the   */
                           /* convergence rate (crate) of the          */
                           /* iterates for the nonlinear equation      */
#define DGMAX  RCONST(0.3) /* iter == CV_NEWTON, |gamma/gammap-1| > DGMAX */
                           /* => call lsetup                           */

#define RDIV      TWO  /* declare divergence if ratio del/delp > RDIV  */
#define MSBP       20  /* max no. of steps between lsetup calls        */

#define TRY_AGAIN  99  /* control constant for CVnlsNewton - should be */
                       /* distinct from CVnls return values            */

/* CVRcheck* return values */

#define INITROOT -1
#define CLOSERT  -2
#define RTFOUND   1

/*
 * Private Helper Functions Prototypes
 */

static booleantype CVCheckNvector(N_Vector tmpl);

static int CVInitialSetup(CVodeMem cv_mem);

static booleantype CVAllocVectors(CVodeMem cv_mem, N_Vector tmpl, int tol);
static void CVFreeVectors(CVodeMem cv_mem);

static int CVEwtSetSS(N_Vector ycur, N_Vector weight, CVodeMem cv_mem);
static int CVEwtSetSV(N_Vector ycur, N_Vector weight, CVodeMem cv_mem);

static booleantype CVHin(CVodeMem cv_mem, realtype tout);
static realtype CVUpperBoundH0(CVodeMem cv_mem, realtype tdist);
static realtype CVYddNorm(CVodeMem cv_mem, realtype hg);

static int CVStep(CVodeMem cv_mem);

static int CVsldet(CVodeMem cv_mem);

static void CVAdjustParams(CVodeMem cv_mem);
static void CVAdjustOrder(CVodeMem cv_mem, int deltaq);
static void CVAdjustAdams(CVodeMem cv_mem, int deltaq);
static void CVAdjustBDF(CVodeMem cv_mem, int deltaq);
static void CVIncreaseBDF(CVodeMem cv_mem);
static void CVDecreaseBDF(CVodeMem cv_mem);

static void CVRescale(CVodeMem cv_mem);

static void CVPredict(CVodeMem cv_mem);

static void CVSet(CVodeMem cv_mem);
static void CVSetAdams(CVodeMem cv_mem);
static realtype CVAdamsStart(CVodeMem cv_mem, realtype m[]);
static void CVAdamsFinish(CVodeMem cv_mem, realtype m[], realtype M[],
                          realtype hsum);
static realtype CVAltSum(int iend, realtype a[], int k);
static void CVSetBDF(CVodeMem cv_mem);
static void CVSetTqBDF(CVodeMem cv_mem, realtype hsum, realtype alpha0,
                     realtype alpha0_hat, realtype xi_inv, realtype xistar_inv);

static int CVnls(CVodeMem cv_mem, int nflag);
static int CVnlsFunctional(CVodeMem cv_mem);
static int CVnlsNewton(CVodeMem cv_mem, int nflag);
static int CVNewtonIteration(CVodeMem cv_mem);

static int  CVHandleNFlag(CVodeMem cv_mem, int *nflagPtr, realtype saved_t,
                          int *ncfPtr);

static void CVRestore(CVodeMem cv_mem, realtype saved_t);

static booleantype CVDoErrorTest(CVodeMem cv_mem, int *nflagPtr, int *kflagPtr,
                               realtype saved_t, int *nefPtr, realtype *dsmPtr);

static void CVCompleteStep(CVodeMem cv_mem);

static void CVPrepareNextStep(CVodeMem cv_mem, realtype dsm);
static void CVSetEta(CVodeMem cv_mem);
static realtype CVComputeEtaqm1(CVodeMem cv_mem);
static realtype CVComputeEtaqp1(CVodeMem cv_mem);
static void CVChooseEta(CVodeMem cv_mem);
static void CVBDFStab(CVodeMem cv_mem);

static int  CVHandleFailure(CVodeMem cv_mem,int kflag);

static int CVRcheck1(CVodeMem cv_mem);
static int CVRcheck2(CVodeMem cv_mem);
static int CVRcheck3(CVodeMem cv_mem);
static int CVRootfind(CVodeMem cv_mem);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/* 
 * CVodeCreate
 *
 * CVodeCreate creates an internal memory block for a problem to 
 * be solved by CVODE.
 * If successful, CVodeCreate returns a pointer to the problem memory. 
 * This pointer should be passed to CVodeMalloc.  
 * If an initialization error occurs, CVodeCreate prints an error 
 * message to standard err and returns NULL. 
 */

void *CVodeCreate(int lmm, int iter)
{
  int maxord;
  CVodeMem cv_mem;

  /* Test inputs */

  if ((lmm != CV_ADAMS) && (lmm != CV_BDF)) {
    fprintf(stderr, MSGCV_BAD_LMM);
    return(NULL);
  }
  
  if ((iter != CV_FUNCTIONAL) && (iter != CV_NEWTON)) {
    fprintf(stderr, MSGCV_BAD_ITER);
    return(NULL);
  }

  cv_mem = (CVodeMem) malloc(sizeof(struct CVodeMemRec));
  if (cv_mem == NULL) {
    fprintf(stderr, MSGCV_CVMEM_FAIL);
    return(NULL);
  }

  maxord = (lmm == CV_ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

  /* copy input parameters into cv_mem */
  cv_mem->cv_lmm  = lmm;
  cv_mem->cv_iter = iter;

  /* Set uround */
  cv_mem->cv_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  cv_mem->cv_f        = NULL;
  cv_mem->cv_f_data   = NULL;
  cv_mem->cv_efun     = NULL;
  cv_mem->cv_e_data   = NULL;
  cv_mem->cv_errfp    = stderr;
  cv_mem->cv_qmax     = maxord;
  cv_mem->cv_mxstep   = MXSTEP_DEFAULT;
  cv_mem->cv_mxhnil   = MXHNIL_DEFAULT;
  cv_mem->cv_sldeton  = FALSE;
  cv_mem->cv_hin      = ZERO;
  cv_mem->cv_hmin     = HMIN_DEFAULT;
  cv_mem->cv_hmax_inv = HMAX_INV_DEFAULT;
  cv_mem->cv_tstopset = FALSE;
  cv_mem->cv_maxcor   = NLS_MAXCOR;
  cv_mem->cv_maxnef   = MXNEF;
  cv_mem->cv_maxncf   = MXNCF;
  cv_mem->cv_nlscoef  = CORTES;

  /* CVodeMalloc not done yet */
  cv_mem->cv_VabstolMallocDone = FALSE;
  cv_mem->cv_MallocDone = FALSE;

  /* Return pointer to CVODE memory block */
  return((void *)cv_mem);
}

#define iter (cv_mem->cv_iter)  
#define lmm (cv_mem->cv_lmm) 
#define errfp (cv_mem->cv_errfp)

/*------------------     CVodeMalloc     --------------------------*/
/* 
   CVodeMalloc allocates and initializes memory for a problem. All 
   problem inputs are checked for errors. If any error occurs during 
   initialization, it is reported to the file whose file pointer is 
   errfp and an error flag is returned. Otherwise, it returns CV_SUCCESS
*/
/*-----------------------------------------------------------------*/

int CVodeMalloc(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0, 
                int itol, realtype reltol, void *abstol)
{
  CVodeMem cv_mem;
  booleantype nvectorOK, allocOK, neg_abstol, ewtsetOK;
  long int lrw1, liw1;
  int i,k;
  
  /* Check cvode_mem */
  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_CVM_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check for legal input parameters */

  if (y0==NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_Y0_NULL);
    return(CV_ILL_INPUT);
  }

  if ((itol != CV_SS) && (itol != CV_SV) && (itol != CV_WF) ) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_ITOL);
    return(CV_ILL_INPUT);
  }

  if (f == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_F_NULL);
    return(CV_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */

  nvectorOK = CVCheckNvector(y0);
  if(!nvectorOK) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_NVECTOR);
    return(CV_ILL_INPUT);
  }

  /* Test tolerances */

  if (abstol == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_ABSTOL_NULL);
    return(CV_ILL_INPUT);
  }

  if (itol != CV_WF) {

    if (reltol < ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_RELTOL);
      return(CV_ILL_INPUT);
    }
   
    if (itol == CV_SS)
      neg_abstol = (*((realtype *)abstol) < ZERO);
    else 
      neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);

    if (neg_abstol) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_ABSTOL);
      return(CV_ILL_INPUT);
    }

  }

  /* Set space requirements for one N_Vector */

  if (y0->ops->nvspace != NULL) {
    N_VSpace(y0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  cv_mem->cv_lrw1 = lrw1;
  cv_mem->cv_liw1 = liw1;

  /* Allocate the vectors (using y0 as a template) */

  allocOK = CVAllocVectors(cv_mem, y0, itol);
  if (!allocOK) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }
 
  /* Copy tolerances into memory */

  cv_mem->cv_itol   = itol;
  cv_mem->cv_reltol = reltol;      

  if (itol == CV_WF)
    cv_mem->cv_efun = (CVEwtFn)abstol;
  else {
    cv_mem->cv_efun = CVEwtSet;
    if (itol == CV_SS)
      cv_mem->cv_Sabstol = *((realtype *)abstol);
    else 
      N_VScale(ONE, (N_Vector)abstol, cv_mem->cv_Vabstol);
  }

  /* All error checking is complete at this point */

  /* Copy the input parameters into CVODE state */

  cv_mem->cv_f  = f;
  cv_mem->cv_tn = t0;

  /* Set step parameters */

  cv_mem->cv_q      = 1;
  cv_mem->cv_L      = 2;
  cv_mem->cv_qwait  = cv_mem->cv_L;
  cv_mem->cv_etamax = ETAMX1;

  cv_mem->cv_qu    = 0;
  cv_mem->cv_hu    = ZERO;
  cv_mem->cv_tolsf = ONE;

  /* Set the linear solver addresses to NULL.
     (We check != NULL later, in CVode, if using CV_NEWTON.) */

  cv_mem->cv_linit  = NULL;
  cv_mem->cv_lsetup = NULL;
  cv_mem->cv_lsolve = NULL;
  cv_mem->cv_lfree  = NULL;
  cv_mem->cv_lmem   = NULL;

  /* Initialize zn[0] in the history array */

  N_VScale(ONE, y0, cv_mem->cv_zn[0]);

  /* Initialize all the counters */

  cv_mem->cv_nst     = 0;
  cv_mem->cv_nfe     = 0;
  cv_mem->cv_ncfn    = 0;
  cv_mem->cv_netf    = 0;
  cv_mem->cv_nni     = 0;
  cv_mem->cv_nsetups = 0;
  cv_mem->cv_nhnil   = 0;
  cv_mem->cv_nstlp   = 0;
  cv_mem->cv_nscon   = 0;
  cv_mem->cv_nge     = 0;

  /* Initialize root finding variables */

  cv_mem->cv_glo    = NULL;
  cv_mem->cv_ghi    = NULL;
  cv_mem->cv_groot  = NULL;
  cv_mem->cv_iroots = NULL;
  cv_mem->cv_gfun   = NULL;
  cv_mem->cv_g_data = NULL;
  cv_mem->cv_nrtfn  = 0;

  /* Initialize Stablilty Limit Detection data */
  /* NOTE: We do this even if stab lim det was not
     turned on yet. This way, the user can turn it
     on at any time */

  cv_mem->cv_nor = 0;
  for (i = 1; i <= 5; i++)
    for (k = 1; k <= 3; k++) 
      cv_mem->cv_ssdat[i-1][k-1] = ZERO;
  
  /* Problem has been successfully initialized */

  cv_mem->cv_MallocDone = TRUE;
  return(CV_SUCCESS);
}

/*------------------     CVodeReInit     --------------------------*/
/*
  CVodeReInit re-initializes CVODE's memory for a problem, assuming
  it has already been allocated in a prior CVodeMalloc call.
  All problem specification inputs are checked for errors.
  If any error occurs during initialization, it is reported to the
  file whose file pointer is errfp.
  The return value is CV_SUCCESS = 0 if no errors occurred, or
  a negative value otherwise.
*/
/*-----------------------------------------------------------------*/

int CVodeReInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0, 
                int itol, realtype reltol, void *abstol)
{
  CVodeMem cv_mem;
  booleantype neg_abstol, ewtsetOK;
  int i,k;
 
  /* Check cvode_mem */

  if (cvode_mem==NULL) {
    fprintf(stderr, MSGCV_CVM_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if cvode_mem was allocated */

  if (cv_mem->cv_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_CVREI_NO_MALLOC);
    return(CV_NO_MALLOC);
  }

  /* Check for legal input parameters */

  if (y0 == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_Y0_NULL);
    return(CV_ILL_INPUT);
  }
  
  if ((itol != CV_SS) && (itol != CV_SV) && (itol != CV_WF)) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_ITOL);
    return(CV_ILL_INPUT);
  }

  if (f == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_F_NULL);
    return(CV_ILL_INPUT);
  }

  /* Test tolerances */

  if (abstol == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_ABSTOL_NULL);
    return(CV_ILL_INPUT);
  }

  if (itol != CV_WF) {

    if (reltol < ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_RELTOL);
      return(CV_ILL_INPUT);
    }
    
    if (itol == CV_SS) {
      neg_abstol = (*((realtype *)abstol) < ZERO);
    } else {
      neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);
    }

    if (neg_abstol) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_ABSTOL);
      return(CV_ILL_INPUT);
    }

  }

  /* Copy tolerances into memory */

  if ( (itol != CV_SV) && (cv_mem->cv_VabstolMallocDone) ) {
    N_VDestroy(cv_mem->cv_Vabstol);
    cv_mem->cv_VabstolMallocDone = FALSE;
  }

  if ( (itol == CV_SV) && !(cv_mem->cv_VabstolMallocDone) ) {
    cv_mem->cv_Vabstol = N_VClone(y0);
    cv_mem->cv_VabstolMallocDone = TRUE;
  }

  cv_mem->cv_itol   = itol;
  cv_mem->cv_reltol = reltol;    

  if (itol == CV_WF)
    cv_mem->cv_efun = (CVEwtFn)abstol;
  else {
    cv_mem->cv_efun = CVEwtSet;
    if (itol == CV_SS)
      cv_mem->cv_Sabstol = *((realtype *)abstol);
    else
      N_VScale(ONE, (N_Vector)abstol, cv_mem->cv_Vabstol);
  }

  /* All error checking is complete at this point */
  
  /* Copy the input parameters into CVODE state */

  cv_mem->cv_f = f;
  cv_mem->cv_tn = t0;

  /* Set step parameters */

  cv_mem->cv_q      = 1;
  cv_mem->cv_L      = 2;
  cv_mem->cv_qwait  = cv_mem->cv_L;
  cv_mem->cv_etamax = ETAMX1;

  cv_mem->cv_qu    = 0;
  cv_mem->cv_hu    = ZERO;
  cv_mem->cv_tolsf = ONE;

  /* Initialize zn[0] in the history array */

  N_VScale(ONE, y0, cv_mem->cv_zn[0]);
 
  /* Initialize all the counters */

  cv_mem->cv_nst     = 0;
  cv_mem->cv_nfe     = 0;
  cv_mem->cv_ncfn    = 0;
  cv_mem->cv_netf    = 0;
  cv_mem->cv_nni     = 0;
  cv_mem->cv_nsetups = 0;
  cv_mem->cv_nhnil   = 0;
  cv_mem->cv_nstlp   = 0;
  cv_mem->cv_nscon   = 0;
  cv_mem->cv_nge     = 0;

  /* Initialize Stablilty Limit Detection data */

  cv_mem->cv_nor = 0;
  for (i = 1; i <= 5; i++)
    for (k = 1; k <= 3; k++) 
      cv_mem->cv_ssdat[i-1][k-1] = ZERO;
  
  /* Problem has been successfully re-initialized */

  return(CV_SUCCESS);
}

#define gfun   (cv_mem->cv_gfun)
#define glo    (cv_mem->cv_glo)
#define ghi    (cv_mem->cv_ghi)
#define groot  (cv_mem->cv_groot)
#define iroots (cv_mem->cv_iroots)

/*------------------     CVodeRootInit     ------------------------*/
/*
  CVodeRootInit initializes a rootfinding problem to be solved
  during the integration of the ODE system.  It loads the root
  function pointer and the number of root functions, and allocates
  workspace memory.  The return value is CV_SUCCESS = 0 if no errors
  occurred, or a negative value otherwise.
*/
/*-----------------------------------------------------------------*/

int CVodeRootInit(void *cvode_mem, CVRootFn g, int nrtfn)
{
  CVodeMem cv_mem;
  int nrt;

  /* Check cvode_mem pointer */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGCV_ROOT_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  nrt = (nrtfn < 0) ? 0 : nrtfn;

  /* If rerunning CVodeRootInit() with a different number of root
     functions (changing number of gfun components), then free
     currently held memory resources */
  if ((nrt != cv_mem->cv_nrtfn) && (cv_mem->cv_nrtfn > 0)) {
    free(glo);
    free(ghi);
    free(groot);
    free(iroots);

    /* Linux version of free() routine doesn't set pointer to NULL */
    glo = ghi = groot = NULL;
    iroots = NULL;
  }

  /* If CVodeRootInit() was called with nrtfn == 0, then set cv_nrtfn to
     zero and cv_gfun to NULL before returning */
  if (nrt == 0) {
    cv_mem->cv_nrtfn = nrt;
    gfun = NULL;
    return(CV_SUCCESS);
  }

  /* If rerunning CVodeRootInit() with the same number of root functions
     (not changing number of gfun components), then check if the root
     function argument has changed */
  /* If g != NULL then return as currently reserved memory resources
     will suffice */
  if (nrt == cv_mem->cv_nrtfn) {
    if (g != gfun) {
      if (g == NULL) {
	free(glo);
	free(ghi);
	free(groot);
	free(iroots);
	fprintf(errfp, MSGCV_ROOT_FUNC_NULL);
	return(CV_RTFUNC_NULL);
      }
      else {
	gfun = g;
	return(CV_SUCCESS);
      }
    }
    else return(CV_SUCCESS);
  }

  /* Set variable values in CVode memory block */
  cv_mem->cv_nrtfn = nrt;
  if (g == NULL) {
    fprintf(errfp, MSGCV_ROOT_FUNC_NULL);
    return(CV_RTFUNC_NULL);
  }
  else gfun = g;

  /* Allocate necessary memory and return */
  glo = (realtype *) malloc(nrt*sizeof(realtype));
  if (glo == NULL) {
    fprintf(stderr, MSGCV_ROOT_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  ghi = (realtype *) malloc(nrt*sizeof(realtype));
  if (ghi == NULL) {
    free(glo);
    fprintf(stderr, MSGCV_ROOT_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  groot = (realtype *) malloc(nrt*sizeof(realtype));
  if (groot == NULL) {
    free(glo); free(ghi);
    fprintf(stderr, MSGCV_ROOT_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  iroots = (int *) malloc(nrt*sizeof(int));
  if (iroots == NULL) {
    free(glo); free(ghi); free(groot);
    fprintf(stderr, MSGCV_ROOT_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  return(CV_SUCCESS);
}

/* 
 * =================================================================
 * Readibility Constants
 * =================================================================
 */

#define f              (cv_mem->cv_f)      
#define f_data         (cv_mem->cv_f_data) 
#define efun           (cv_mem->cv_efun)
#define e_data         (cv_mem->cv_e_data) 
#define g_data         (cv_mem->cv_g_data) 
#define qmax           (cv_mem->cv_qmax) 
#define mxstep         (cv_mem->cv_mxstep)
#define mxhnil         (cv_mem->cv_mxhnil)
#define sldeton        (cv_mem->cv_sldeton)
#define hin            (cv_mem->cv_hin)
#define hmin           (cv_mem->cv_hmin)
#define hmax_inv       (cv_mem->cv_hmax_inv)
#define tstop          (cv_mem->cv_tstop)
#define tstopset       (cv_mem->cv_tstopset)
#define maxnef         (cv_mem->cv_maxnef)
#define maxncf         (cv_mem->cv_maxncf)
#define maxcor         (cv_mem->cv_maxcor)
#define nlscoef        (cv_mem->cv_nlscoef)
#define itol           (cv_mem->cv_itol)         
#define reltol         (cv_mem->cv_reltol)       
#define Sabstol        (cv_mem->cv_Sabstol)     
#define Vabstol        (cv_mem->cv_Vabstol)     

#define uround         (cv_mem->cv_uround)  
#define zn             (cv_mem->cv_zn) 
#define ewt            (cv_mem->cv_ewt)  
#define y              (cv_mem->cv_y)
#define acor           (cv_mem->cv_acor)
#define tempv          (cv_mem->cv_tempv)
#define ftemp          (cv_mem->cv_ftemp) 
#define q              (cv_mem->cv_q)
#define qprime         (cv_mem->cv_qprime)
#define next_q         (cv_mem->cv_next_q)
#define qwait          (cv_mem->cv_qwait)
#define L              (cv_mem->cv_L)
#define h              (cv_mem->cv_h)
#define hprime         (cv_mem->cv_hprime)
#define next_h         (cv_mem->cv_next_h)
#define eta            (cv_mem->cv_eta) 
#define etaqm1         (cv_mem->cv_etaqm1) 
#define etaq           (cv_mem->cv_etaq) 
#define etaqp1         (cv_mem->cv_etaqp1) 
#define nscon          (cv_mem->cv_nscon)
#define hscale         (cv_mem->cv_hscale)
#define tn             (cv_mem->cv_tn)
#define tau            (cv_mem->cv_tau)
#define tq             (cv_mem->cv_tq)
#define l              (cv_mem->cv_l)
#define rl1            (cv_mem->cv_rl1)
#define gamma          (cv_mem->cv_gamma) 
#define gammap         (cv_mem->cv_gammap) 
#define gamrat         (cv_mem->cv_gamrat)
#define crate          (cv_mem->cv_crate)
#define acnrm          (cv_mem->cv_acnrm)
#define mnewt          (cv_mem->cv_mnewt)
#define etamax         (cv_mem->cv_etamax)
#define nst            (cv_mem->cv_nst)
#define nfe            (cv_mem->cv_nfe)
#define ncfn           (cv_mem->cv_ncfn)
#define netf           (cv_mem->cv_netf)
#define nni            (cv_mem->cv_nni)
#define nsetups        (cv_mem->cv_nsetups)
#define nhnil          (cv_mem->cv_nhnil)
#define lrw1           (cv_mem->cv_lrw1)
#define liw1           (cv_mem->cv_liw1)
#define lrw            (cv_mem->cv_lrw)
#define liw            (cv_mem->cv_liw)
#define linit          (cv_mem->cv_linit)
#define lsetup         (cv_mem->cv_lsetup)
#define lsolve         (cv_mem->cv_lsolve) 
#define lfree          (cv_mem->cv_lfree) 
#define lmem           (cv_mem->cv_lmem) 
#define qu             (cv_mem->cv_qu)          
#define nstlp          (cv_mem->cv_nstlp)  
#define h0u            (cv_mem->cv_h0u)
#define hu             (cv_mem->cv_hu)         
#define saved_tq5      (cv_mem->cv_saved_tq5)  
#define jcur           (cv_mem->cv_jcur)         
#define tolsf          (cv_mem->cv_tolsf)      
#define setupNonNull   (cv_mem->cv_setupNonNull) 
#define nor            (cv_mem->cv_nor)
#define ssdat          (cv_mem->cv_ssdat)
#define nrtfn          (cv_mem->cv_nrtfn)
#define tlo            (cv_mem->cv_tlo)
#define thi            (cv_mem->cv_thi)
#define tretlast       (cv_mem->cv_tretlast)
#define toutc          (cv_mem->cv_toutc)
#define troot          (cv_mem->cv_troot)
#define ttol           (cv_mem->cv_ttol)
#define taskc          (cv_mem->cv_taskc)
#define irfnd          (cv_mem->cv_irfnd)
#define nge            (cv_mem->cv_nge)

/********************* CVode ****************************************

 This routine is the main driver of the CVODE package. 

 It integrates over a time interval defined by the user, by calling
 CVStep to do internal time steps.

 The first time that CVode is called for a successfully initialized
 problem, it computes a tentative initial step size h.

 CVode supports four modes, specified by itask: CV_NORMAL, CV_ONE_STEP,
 CV_NORMAL_TSTOP, and CV_ONE_STEP_TSTOP.
 In the CV_NORMAL mode, the solver steps until it reaches or passes tout
 and then interpolates to obtain y(tout).
 In the CV_ONE_STEP mode, it takes one internal step and returns.
 CV_NORMAL_TSTOP and CV_ONE_STEP_TSTOP are similar to CV_NORMAL and CV_ONE_STEP,
 respectively, but the integration never proceeds past tstop (which
 must have been defined through a call to CVodeSetStopTime).

********************************************************************/

int CVode(void *cvode_mem, realtype tout, N_Vector yout, 
          realtype *tret, int itask)
{
  CVodeMem cv_mem;
  long int nstloc;
  int kflag, istate, ier, task, irfndp;
  realtype troundoff, rh;
  booleantype istop, hOK;
  int ewtsetOK;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGCV_CVODE_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if cvode_mem was allocated */
  if (cv_mem->cv_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_CVODE_NO_MALLOC);
    return(CV_NO_MALLOC);
  }
  
  /* Check for yout != NULL */
  if ((y = yout) == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_YOUT_NULL);       
    return(CV_ILL_INPUT);
  }
  
  /* Check for tret != NULL */
  if (tret == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_TRET_NULL);
    return(CV_ILL_INPUT);
  }

  /* Check for valid itask */
  if ((itask != CV_NORMAL)       && 
      (itask != CV_ONE_STEP)     &&
      (itask != CV_NORMAL_TSTOP) &&
      (itask != CV_ONE_STEP_TSTOP) ) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_ITASK);
    return(CV_ILL_INPUT);
  }

  /* Split itask into task and istop */
  if ((itask == CV_NORMAL_TSTOP) || (itask == CV_ONE_STEP_TSTOP)) {
    if ( tstopset == FALSE ) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_NO_TSTOP);
      return(CV_ILL_INPUT);
    }
    istop = TRUE;
  } else {
    istop = FALSE;
  }
  if ((itask == CV_NORMAL) || (itask == CV_NORMAL_TSTOP)) {
    task = CV_NORMAL; toutc = tout;
  } else {
    task = CV_ONE_STEP;
  }
  taskc = task;

  /* Begin first call block */ 

  if (nst == 0) {

    ier = CVInitialSetup(cv_mem);
    if (ier!= CV_SUCCESS) return (ier);
    
    /* Call f at (t0,y0), set zn[1] = y'(t0), 
       set initial h (from H0 or CVHin), and scale zn[1] by h.
       Also check for zeros of root function g at and near t0.    */
    
    f(tn, zn[0], zn[1], f_data); 
    nfe++;

    h = hin;
    if ( (h != ZERO) && ((tout-tn)*h < ZERO) ) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_H0);
      return(CV_ILL_INPUT);
    }
    if (h == ZERO) {
      hOK = CVHin(cv_mem, tout);
      if (!hOK) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_TOO_CLOSE);
        return(CV_ILL_INPUT);
      }
    }
    rh = ABS(h)*hmax_inv;
    if (rh > ONE) h /= rh;
    if (ABS(h) < hmin) h *= hmin/ABS(h);

    /* Check for approach to tstop */

    if (istop) {
      if ( (tstop - tn)*h < ZERO ) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_TSTOP, tn);
        return(CV_ILL_INPUT);
      }
      if ( (tn + h - tstop)*h > ZERO ) 
        h = tstop - tn;
    }

    hscale = h; 
    h0u    = h;
    hprime = h;

    N_VScale(h, zn[1], zn[1]);

    if (nrtfn > 0) {
      ier = CVRcheck1(cv_mem);
      if (ier != CV_SUCCESS) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_INIT_ROOT);
        return(CV_ILL_INPUT);
      }
    }

  } /* end of first call block */

  /* At following steps, perform stop tests */

  if (nst > 0) {

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If task = CV_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */

    if (nrtfn > 0) {

      irfndp = irfnd;
      
      ier = CVRcheck2(cv_mem);

      if (ier == CLOSERT) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_CLOSE_ROOTS, tlo);
        return(CV_ILL_INPUT);
      }

      if (ier == RTFOUND) {
        tretlast = *tret = tlo;
        return(CV_ROOT_RETURN);
      }

      if (tn != tretlast) {       /* Check remaining interval for roots */
        ier = CVRcheck3(cv_mem);
        if (ier == CV_SUCCESS) {     /* no root found */
          irfnd = 0;
          if ((irfndp == 1) && (task == CV_ONE_STEP)) {
            tretlast = *tret = tn;
            N_VScale(ONE, zn[0], yout);
            return(CV_SUCCESS);
          }
        }
        if (ier == RTFOUND) {  /* a new root was found */
          irfnd = 1;
          tretlast = *tret = tlo;
          return(CV_ROOT_RETURN);
        }
      }

    } /* end of root stop check */

    /* Test for tn past tstop */
    if ( istop && ((tstop - tn)*h < ZERO) ) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_TSTOP, tn);
      return(CV_ILL_INPUT);
    }

    /* In CV_NORMAL mode, test if tout was reached */
    if ( (task == CV_NORMAL) && ((tn-tout)*h >= ZERO) ) {
      tretlast = *tret = tout;
      ier =  CVodeGetDky(cv_mem, tout, 0, yout);
      if (ier != CV_SUCCESS) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_TOUT, tout);
        return(CV_ILL_INPUT);
      }
      return(CV_SUCCESS);
    }

    /* In CV_ONE_STEP mode, test if tn was returned */
    if (task == CV_ONE_STEP && tretlast != tn) {
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      return(CV_SUCCESS);
    }

    /* Test for tn at tstop or near tstop */
    if ( istop ) {

      troundoff = FUZZ_FACTOR*uround*(ABS(tn) + ABS(h));
      if ( ABS(tn - tstop) <= troundoff) {
        ier =  CVodeGetDky(cv_mem, tstop, 0, yout);
        if (ier != CV_SUCCESS) {
          if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_TSTOP, tn);
          return(CV_ILL_INPUT);
        }
        tretlast = *tret = tstop;
        return(CV_TSTOP_RETURN);
      }
      
      if ( (tn + hprime - tstop)*h > ZERO ) {
        hprime = tstop - tn;
        eta = hprime/h;
      }

    } /* end of istop tests block */
    
  } /* end stopping tests block */  

  /* Looping point for internal steps */

  nstloc = 0;
  loop {
   
    next_h = h;
    next_q = q;
    
    /* Reset and check ewt */

    if (nst > 0) {
      ewtsetOK = efun(zn[0], ewt, e_data);
      if (ewtsetOK != 0) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_EWT_NOW_BAD, tn);
        istate = CV_ILL_INPUT;
        tretlast = *tret = tn;
        N_VScale(ONE, zn[0], yout);
        break;
      }
    }
    
    /* Check for too many steps */
    
    if (nstloc >= mxstep) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_MAX_STEPS, tn);
      istate = CV_TOO_MUCH_WORK;
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      break;
    }

    /* Check for too much accuracy requested */

    if ((tolsf = uround * N_VWrmsNorm(zn[0], ewt)) > ONE) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_TOO_MUCH_ACC, tn);
      istate = CV_TOO_MUCH_ACC;
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      tolsf *= TWO;
      break;
    } else {
      tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */

    if (tn + h == tn) {
      nhnil++;
      if (nhnil <= mxhnil) if(errfp!=NULL) fprintf(errfp, MSGCV_HNIL, tn, h);
      if (nhnil == mxhnil) if(errfp!=NULL) fprintf(errfp, MSGCV_HNIL_DONE);
    }

    /* Call CVStep to take a step */

    kflag = CVStep(cv_mem);

    /* Process failed step cases, and exit loop */
   
    if (kflag != SUCCESS_STEP) {
      istate = CVHandleFailure(cv_mem, kflag);
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      break;
    }
    
    nstloc++;

    /* Check for root in last step taken. */

    if (nrtfn > 0) {

      ier = CVRcheck3(cv_mem);

      if (ier == RTFOUND) {  /* A new root was found */
        irfnd = 1;
        istate = CV_ROOT_RETURN;
        tretlast = *tret = tlo;
        break;
      }
    }

    /* Check if tn is at tstop or near tstop */

    if ( istop ) {

      troundoff = FUZZ_FACTOR*uround*(ABS(tn) + ABS(h));
      if ( ABS(tn - tstop) <= troundoff) {
        (void) CVodeGetDky(cv_mem, tstop, 0, yout);
        tretlast = *tret = tstop;
        istate = CV_TSTOP_RETURN;
        break;
      }

      if ( (tn + hprime - tstop)*h > ZERO ) {
        hprime = tstop - tn;
        eta = hprime/h;
      }

    }

    /* Check if in one-step mode, and if so copy y and exit loop */
    
    if (task == CV_ONE_STEP) {
      istate = CV_SUCCESS;
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      next_q = qprime;
      next_h = hprime;
      break;
    }

    /* Check if tout reached, and if so interpolate and exit loop */

    if ((tn-tout)*h >= ZERO) {
      istate = CV_SUCCESS;
      tretlast = *tret = tout;
      (void) CVodeGetDky(cv_mem, tout, 0, yout);
      next_q = qprime;
      next_h = hprime;
      break;
    }
  }

  return(istate);
}

/*------------------    CVodeGetDky      --------------------------*/
/*
  This routine computes the k-th derivative of the interpolating
  polynomial at the time t and stores the result in the vector dky.
  The formula is:
          q 
   dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] , 
         j=k 
  where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
  zn[j] is the j-th column of the Nordsieck history array.

  This function is called by CVode with k = 0 and t = tout, but
  may also be called directly by the user.
*/
/*-----------------------------------------------------------------*/

int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  CVodeMem cv_mem;
  
  /* Check all inputs for legality */
 
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGCV_DKY_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (dky == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_DKY);
    return(CV_BAD_DKY);
  }

  if ((k < 0) || (k > q)) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_K);
    return(CV_BAD_K);
  }
  
  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * uround * (ABS(tn) + ABS(hu));
  if (hu < ZERO) tfuzz = -tfuzz;
  tp = tn - hu - tfuzz;
  tn1 = tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_T, t, tn-hu, tn);
    return(CV_BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */

  s = (t - tn) / h;
  for (j=q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == q) {
      N_VScale(c, zn[q], dky);
    } else {
      N_VLinearSum(c, zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(CV_SUCCESS);
  r = RPowerI(h,-k);
  N_VScale(r, dky, dky);
  return(CV_SUCCESS);
}

/********************* CVodeFree **********************************

 This routine frees the problem memory allocated by CVodeMalloc.
 Such memory includes all the vectors allocated by CVAllocVectors,
 and the memory lmem for the linear solver (deallocated by a call
 to lfree).

*******************************************************************/

void CVodeFree(void *cvode_mem)
{
  CVodeMem cv_mem;

  cv_mem = (CVodeMem) cvode_mem;
  
  if (cvode_mem == NULL) return;

  CVFreeVectors(cv_mem);

  if (iter == CV_NEWTON && lfree != NULL) lfree(cv_mem);

  if (nrtfn > 0) {
    free(glo);
    free(ghi);
    free(groot);
    free(iroots); 
  }

  free(cv_mem);
}

/* 
 * =================================================================
 *  Private Functions Implementation
 * =================================================================
 */

/****************** CVCheckNvector ***********************************
 This routine checks if all required vector operations are present.
 If any of them is missing it returns FALSE.
**********************************************************************/

static booleantype CVCheckNvector(N_Vector tmpl)
{
  if((tmpl->ops->nvclone     == NULL) ||
     (tmpl->ops->nvdestroy   == NULL) ||
     (tmpl->ops->nvlinearsum == NULL) ||
     (tmpl->ops->nvconst     == NULL) ||
     (tmpl->ops->nvprod      == NULL) ||
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

/****************** CVAllocVectors ***********************************

 This routine allocates the CVODE vectors ewt, acor, tempv, ftemp, and
 zn[0], ..., zn[qmax]. If tol=CV_SV, it also allocates space for Vabstol.
 If all memory allocations are successful, CVAllocVectors returns TRUE. 
 Otherwise all allocated memory is freed and CVAllocVectors returns FALSE.
 This routine also sets the optional outputs lrw and liw, which are
 (respectively) the lengths of the real and integer work spaces
 allocated here.

**********************************************************************/

static booleantype CVAllocVectors(CVodeMem cv_mem, N_Vector tmpl, int tol)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */
  
  ewt = N_VClone(tmpl);
  if (ewt == NULL) return(FALSE);
  acor = N_VClone(tmpl);
  if (acor == NULL) {
    N_VDestroy(ewt);
    return(FALSE);
  }
  tempv = N_VClone(tmpl);
  if (tempv == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(acor);
    return(FALSE);
  }
  ftemp = N_VClone(tmpl);
  if (ftemp == NULL) {
    N_VDestroy(tempv);
    N_VDestroy(ewt);
    N_VDestroy(acor);
    return(FALSE);
  }

  /* Allocate zn[0] ... zn[qmax] */

  for (j=0; j <= qmax; j++) {
    zn[j] = N_VClone(tmpl);
    if (zn[j] == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(acor);
      N_VDestroy(tempv);
      N_VDestroy(ftemp);
      for (i=0; i < j; i++) N_VDestroy(zn[i]);
      return(FALSE);
    }
  }

  if (tol == CV_SV) {
    Vabstol = N_VClone(tmpl);
    if (Vabstol == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(acor);
      N_VDestroy(tempv);
      N_VDestroy(ftemp);
      for (i=0; i <= qmax; i++) N_VDestroy(zn[i]);
      return(FALSE);
    }
    cv_mem->cv_VabstolMallocDone = TRUE;
  }

  /* Set solver workspace lengths  */

  lrw = (qmax + 5)*lrw1;
  liw = (qmax + 5)*liw1;

  return(TRUE);
}

/***************** CVFreeVectors *********************************
  
 This routine frees the CVODE vectors allocated in CVAllocVectors.

******************************************************************/

static void CVFreeVectors(CVodeMem cv_mem)
{
  int j;
  
  N_VDestroy(ewt);
  N_VDestroy(acor);
  N_VDestroy(tempv);
  N_VDestroy(ftemp);
  for(j=0; j <= qmax; j++) N_VDestroy(zn[j]);
  if (cv_mem->cv_VabstolMallocDone) N_VDestroy(Vabstol);
}

/*-----------------------------------------------------------------*/

/*  
 * CVInitialSetup
 *
 * This routine performs input consistency checks at the first step.
 * If needed, it also checks the linear solver module and calls the
 * linear solver initialization routine.
 */

static int CVInitialSetup(CVodeMem cv_mem)
{
  int ier;
  int ewtsetOK;

  /* Solver initial setup */

  if (itol != CV_WF)
    e_data = (void *)cv_mem;

  ewtsetOK = efun(zn[0], ewt, e_data);
  if (ewtsetOK != 0) {
    if(errfp!=NULL) fprintf(errfp, MSGCV_BAD_EWT);
    return(CV_ILL_INPUT);
  }
  
  /* Check if lsolve function exists (if needed)
     and call linit function (if it exists) */

  if (iter == CV_NEWTON) {
    if (lsolve == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSGCV_LSOLVE_NULL);
      return (CV_ILL_INPUT);
    }
    if (linit != NULL) {
      ier = linit(cv_mem);
      if (ier != 0) {
        if(errfp!=NULL) fprintf(errfp, MSGCV_LINIT_FAIL);
        return (CV_ILL_INPUT);
      }
    }
  }
    
  return(CV_SUCCESS);
}


/******************* CVHin ***************************************

 This routine computes a tentative initial step size h0. 
 If tout is too close to tn (= t0), then CVHin returns FALSE and
 h remains uninitialized. Otherwise, CVHin sets h to the chosen 
 value h0 and returns TRUE.

 The algorithm used seeks to find h0 as a solution of
       (WRMS norm of (h0^2 ydd / 2)) = 1, 
 where ydd = estimated second derivative of y.

*****************************************************************/

static booleantype CVHin(CVodeMem cv_mem, realtype tout)
{
  int sign, count;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hnew, hrat, h0, yddnrm;

  /* Test for tout too close to tn */
  
  if ((tdiff = tout-tn) == ZERO) return(FALSE);
  
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = ABS(tdiff);
  tround = uround * MAX(ABS(tn), ABS(tout));
  if (tdist < TWO*tround) return(FALSE);
  
  /* Set lower and upper bounds on h0, and take geometric mean 
     Exit with this value if the bounds cross each other       */

  hlb = HLB_FACTOR * tround;
  hub = CVUpperBoundH0(cv_mem, tdist);
  hg  = RSqrt(hlb*hub);
  if (hub < hlb) {
    if (sign == -1) hg = -hg;
    h = hg;
    return(TRUE);
  }
  
  /* Loop up to MAX_ITERS times to find h0.
     Stop if new and previous values differ by a factor < 2.
     Stop if hnew/hg > 2 after one iteration, as this probably means
     that the ydd value is bad because of cancellation error.        */

  count = 0;
  loop {
    hgs = hg*sign;
    yddnrm = CVYddNorm(cv_mem, hgs);
    hnew =  (yddnrm*hub*hub > TWO) ? RSqrt(TWO/yddnrm) : RSqrt(hg*hub);
    count++;
    if (count >= MAX_ITERS) break;
    hrat = hnew/hg;
    if ((hrat > HALF) && (hrat < TWO)) break;
    if ((count >= 2) && (hrat > TWO)) {
      hnew = hg;
      break;
    }
    hg = hnew;
  }
  
  /* Apply bounds, bias factor, and attach sign */

  h0 = H_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  h = h0;
  return(TRUE);
}

/******************** CVUpperBoundH0 ******************************

 This routine sets an upper bound on abs(h0) based on
 tdist = abs(tout - t0) and the values of y[i]/y'[i].

******************************************************************/

static realtype CVUpperBoundH0(CVodeMem cv_mem, realtype tdist)
{
  realtype hub_inv, hub;
  N_Vector temp1, temp2;

  temp1 = tempv;
  temp2 = acor;

  /* 
   * Bound based on |y0|/|y0'| -- allow at most an increase of
   * HUB_FACTOR in y0 (based on a forward Euler step). The weight 
   * factor is used as a safeguard against zero components in y0. 
   */

  N_VAbs(zn[0], temp2);
  efun(zn[0], temp1, e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(zn[1], temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /*
   * bound based on tdist -- allow at most a step of magnitude
   * HUB_FACTOR * tdist
   */

  hub = HUB_FACTOR*tdist;

  /* Use the smaler of the two */

  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}

/****************** CVYddNorm *************************************

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.

******************************************************************/

static realtype CVYddNorm(CVodeMem cv_mem, realtype hg)
{
  realtype yddnrm;
  
  N_VLinearSum(hg, zn[1], ONE, zn[0], y);
  f(tn+hg, y, tempv, f_data);
  nfe++;
  N_VLinearSum(ONE, tempv, -ONE, zn[1], tempv);
  N_VScale(ONE/hg, tempv, tempv);

  yddnrm = N_VWrmsNorm(tempv, ewt);
  return(yddnrm);
}

/********************* CVStep **************************************
 
 This routine performs one internal cvode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
  * preliminary adjustments if a new step size was chosen;
  * prediction of the Nordsieck history array zn at tn + h;
  * setting of multistep method coefficients and test quantities;
  * solution of the nonlinear system;
  * testing the local error;
  * updating zn and other state data if successful;
  * resetting stepsize and order for the next step.
  * if SLDET is on, check for stability, reduce order if necessary.
 On a failure in the nonlinear system solution or error test, the
 step may be reattempted, depending on the nature of the failure.

********************************************************************/

static int CVStep(CVodeMem cv_mem)
{
  realtype saved_t, dsm;
  int ncf, nef, nflag;
  booleantype passed;

  int kflag;
  
  saved_t = tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;

  
  if ((nst > 0) && (hprime != h)) CVAdjustParams(cv_mem);
  
  /* Looping point for attempts to take a step */
  loop {  
    CVPredict(cv_mem);  
    CVSet(cv_mem);

    nflag = CVnls(cv_mem, nflag);
    kflag = CVHandleNFlag(cv_mem, &nflag, saved_t, &ncf);
    if (kflag == PREDICT_AGAIN) continue;
    if (kflag != DO_ERROR_TEST) return(kflag);
    /* Return if nonlinear solve failed and recovery not possible. */

    passed = CVDoErrorTest(cv_mem, &nflag, &kflag, saved_t, &nef, &dsm);
    /* Return if error test failed and recovery not possible. */
    if ((!passed) && (kflag == REP_ERR_FAIL)) return(kflag);
    if (passed) break;
    /* Retry step if error test failed, nflag == PREV_ERR_FAIL */
  }

  /* Nonlinear system solve and error test were both successful.
     Update data, and consider change of step and/or order.       */


  CVCompleteStep(cv_mem); 
  CVPrepareNextStep(cv_mem, dsm); 

  /* If Stablilty Limit Detection is turned on, call stability limit
     detection routine for possible order reduction. */

  if (sldeton) CVBDFStab(cv_mem);

  etamax = (nst <= SMALL_NST) ? ETAMX2 : ETAMX3;

  /*  Finally, we rescale the acor array to be the 
      estimated local error vector. */

  N_VScale(ONE/tq[2], acor, acor);
  return(SUCCESS_STEP);
      
}


/********************* CVAdjustParams ********************************

 This routine is called when a change in step size was decided upon,
 and it handles the required adjustments to the history array zn.
 If there is to be a change in order, we call CVAdjustOrder and reset
 q, L = q+1, and qwait.  Then in any case, we call CVRescale, which
 resets h and rescales the Nordsieck array.

**********************************************************************/

static void CVAdjustParams(CVodeMem cv_mem)
{
  if (qprime != q) {
    CVAdjustOrder(cv_mem, qprime-q);
    q = qprime;
    L = q+1;
    qwait = L;
  }
  CVRescale(cv_mem);
}

/********************* CVAdjustOrder *****************************

  This routine is a high level routine which handles an order
  change by an amount deltaq (= +1 or -1). If a decrease in order
  is requested and q==2, then the routine returns immediately.
  Otherwise CVAdjustAdams or CVAdjustBDF is called to handle the
  order change (depending on the value of lmm).

******************************************************************/

static void CVAdjustOrder(CVodeMem cv_mem, int deltaq)
{
  if ((q==2) && (deltaq != 1)) return;
  
  switch(lmm){
    case CV_ADAMS: CVAdjustAdams(cv_mem, deltaq);
                break;
    case CV_BDF:   CVAdjustBDF(cv_mem, deltaq);
                break;
  }
}

/*************** CVAdjustAdams ***********************************

 This routine adjusts the history array on a change of order q by
 deltaq, in the case that lmm == CV_ADAMS.

*****************************************************************/

static void CVAdjustAdams(CVodeMem cv_mem, int deltaq)
{
  int i, j;
  realtype xi, hsum;

  /* On an order increase, set new column of zn to zero and return */
  
  if (deltaq==1) {
    N_VConst(ZERO, zn[L]);
    return;
  }

  /* On an order decrease, each zn[j] is adjusted by a multiple
     of zn[q].  The coefficients in the adjustment are the 
     coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
     integrated, where xi_j = [t_n - t_(n-j)]/h.               */

  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[1] = ONE;
  hsum = ZERO;
  for (j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum / hscale;
    for (i=j+1; i >= 1; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for (j=1; j <= q-2; j++) l[j+1] = q * (l[j] / (j+1));
  
  for (j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);
}

/********************** CVAdjustBDF *******************************

 This is a high level routine which handles adjustments to the
 history array on a change of order by deltaq in the case that 
 lmm == CV_BDF.  CVAdjustBDF calls CVIncreaseBDF if deltaq = +1 and 
 CVDecreaseBDF if deltaq = -1 to do the actual work.

******************************************************************/

static void CVAdjustBDF(CVodeMem cv_mem, int deltaq)
{
  switch(deltaq) {
    case 1 : CVIncreaseBDF(cv_mem);
             return;
    case -1: CVDecreaseBDF(cv_mem);
             return;
  }
}

/******************** CVIncreaseBDF **********************************

 This routine adjusts the history array on an increase in the 
 order q in the case that lmm == CV_BDF.  
 A new column zn[q+1] is set equal to a multiple of the saved 
 vector (= acor) in zn[qmax].  Then each zn[j] is adjusted by
 a multiple of zn[q+1].  The coefficients in the adjustment are the 
 coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_(q-1)),
 where xi_j = [t_n - t_(n-j)]/h.

*********************************************************************/

static void CVIncreaseBDF(CVodeMem cv_mem)
{
  realtype alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = hscale;
  if (q > 1) {
    for (j=1; j < q; j++) {
      hsum += tau[j+1];
      xi = hsum / hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--) l[i] = l[i]*xiold + l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;
  N_VScale(A1, zn[qmax], zn[L]);
  for (j=2; j <= q; j++) {
    N_VLinearSum(l[j], zn[L], ONE, zn[j], zn[j]);
  }  
}

/********************* CVDecreaseBDF ******************************

 This routine adjusts the history array on a decrease in the 
 order q in the case that lmm == CV_BDF.  
 Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 in the adjustment are the coefficients of the polynomial
 x*x*(x+xi_1)*...*(x+xi_(q-2)), where xi_j = [t_n - t_(n-j)]/h.

******************************************************************/

static void CVDecreaseBDF(CVodeMem cv_mem)
{
  realtype hsum, xi;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = ONE;
  hsum = ZERO;
  for(j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum /hscale;
    for (i=j+2; i >= 2; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for(j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);
}

/**************** CVRescale ***********************************

  This routine rescales the Nordsieck array by multiplying the
  jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
  h is rescaled by eta, and hscale is reset to h.

***************************************************************/

static void CVRescale(CVodeMem cv_mem)
{
  int j;
  realtype factor;
  
  factor = eta;
  for (j=1; j <= q; j++) {
    N_VScale(factor, zn[j], zn[j]);
    factor *= eta;
  }
  h = hscale * eta;
  hscale = h;
  nscon = 0;
}

/********************* CVPredict *************************************

 This routine advances tn by the tentative step size h, and computes
 the predicted array z_n(0), which is overwritten on zn.  The
 prediction of zn is done by repeated additions.

*********************************************************************/

static void CVPredict(CVodeMem cv_mem)
{
  int j, k;
  
  tn += h;
  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--) 
      N_VLinearSum(ONE, zn[j-1], ONE, zn[j], zn[j-1]); 
}

/************************** CVSet *********************************

 This routine is a high level routine which calls CVSetAdams or
 CVSetBDF to set the polynomial l, the test quantity array tq, 
 and the related variables  rl1, gamma, and gamrat.

******************************************************************/

static void CVSet(CVodeMem cv_mem)
{
  switch(lmm) {
    case CV_ADAMS: CVSetAdams(cv_mem);
                break;
    case CV_BDF  : CVSetBDF(cv_mem);
                break;
  }
  rl1 = ONE / l[1];
  gamma = h * rl1;
  if (nst == 0) gammap = gamma;
  gamrat = (nst > 0) ? gamma / gammap : ONE;  /* protect x / x != 1.0 */
}

/******************** CVSetAdams *********************************

 This routine handles the computation of l and tq for the
 case lmm == CV_ADAMS.

 The components of the array l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                          q-1
 (d/dx) Lambda(x) = c * PRODUCT (1 + x / xi_i) , where
                          i=1
 Lambda(-1) = 0, Lambda(0) = 1, and c is a normalization factor.
 Here xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.

*****************************************************************/

static void CVSetAdams(CVodeMem cv_mem)
{
  realtype m[L_MAX], M[3], hsum;
  
  if (q == 1) {
    l[0] = l[1] = tq[1] = tq[5] = ONE;
    tq[2] = TWO;
    tq[3] = TWELVE;
    tq[4] = nlscoef * tq[2];       /* = 0.1 * tq[2] */
    return;
  }
  
  hsum = CVAdamsStart(cv_mem, m);
  
  M[0] = CVAltSum(q-1, m, 1);
  M[1] = CVAltSum(q-1, m, 2);
  
  CVAdamsFinish(cv_mem, m, M, hsum);
}

/****************** CVAdamsStart ********************************

 This routine generates in m[] the coefficients of the product
 polynomial needed for the Adams l and tq coefficients for q > 1.
  
******************************************************************/

static realtype CVAdamsStart(CVodeMem cv_mem, realtype m[])
{
  realtype hsum, xi_inv, sum;
  int i, j;
  
  hsum = h;
  m[0] = ONE;
  for (i=1; i <= q; i++) m[i] = ZERO;
  for (j=1; j < q; j++) {
    if ((j==q-1) && (qwait == 1)) {
      sum = CVAltSum(q-2, m, 2);
      tq[1] = m[q-2] / (q * sum);
    }
    xi_inv = h / hsum;
    for (i=j; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    hsum += tau[j];
    /* The m[i] are coefficients of product(1 to j) (1 + x/xi_i) */
  }
  return(hsum);
}

/****************** CVAdamsFinish  *******************************

 This routine completes the calculation of the Adams l and tq.

******************************************************************/

static void CVAdamsFinish(CVodeMem cv_mem, realtype m[], realtype M[], 
                          realtype hsum)
{
  int i;
  realtype M0_inv, xi, xi_inv;
  
  M0_inv = ONE / M[0];
  
  l[0] = ONE;
  for (i=1; i <= q; i++) l[i] = M0_inv * (m[i-1] / i);
  xi = hsum / h;
  xi_inv = ONE / xi;
  
  tq[2] = xi * M[0] / M[1];
  tq[5] = xi / l[q];

  if (qwait == 1) {
    for (i=q; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    M[2] = CVAltSum(q, m, 2);
    tq[3] = L * M[0] / M[2];
  }

  tq[4] = nlscoef * tq[2];
}

/****************** CVAltSum **************************************
  
 CVAltSum returns the value of the alternating sum
   sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 If iend < 0 then CVAltSum returns 0.
 This operation is needed to compute the integral, from -1 to 0,
 of a polynomial x^(k-1) M(x) given the coefficients of M(x).

******************************************************************/

static realtype CVAltSum(int iend, realtype a[], int k)
{
  int i, sign;
  realtype sum;
  
  if (iend < 0) return(ZERO);
  
  sum = ZERO;
  sign = 1;
  for (i=0; i <= iend; i++) {
    sum += sign * (a[i] / (i+k));
    sign = -sign;
  }
  return(sum);
}

/***************** CVSetBDF **************************************

 This routine computes the coefficients l and tq in the case
 lmm == CV_BDF.  CVSetBDF calls CVSetTqBDF to set the test
 quantity array tq. 

 The components of the array l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                                 q-1
 Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
                                 i=1
 xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.


*****************************************************************/

static void CVSetBDF(CVodeMem cv_mem)
{
  realtype alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;
  
  l[0] = l[1] = xi_inv = xistar_inv = ONE;
  for (i=2; i <= q; i++) l[i] = ZERO;
  alpha0 = alpha0_hat = -ONE;
  hsum = h;
  if (q > 1) {
    for (j=2; j < q; j++) {
      hsum += tau[j-1];
      xi_inv = h / hsum;
      alpha0 -= ONE / j;
      for(i=j; i >= 1; i--) l[i] += l[i-1]*xi_inv;
      /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }
    
    /* j = q */
    alpha0 -= ONE / q;
    xistar_inv = -l[1] - alpha0;
    hsum += tau[q-1];
    xi_inv = h / hsum;
    alpha0_hat = -l[1] - xi_inv;
    for (i=q; i >= 1; i--) l[i] += l[i-1]*xistar_inv;
  }

  CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/****************** CVSetTqBDF ************************************

 This routine sets the test quantity array tq when lmm == CV_BDF.

******************************************************************/

static void CVSetTqBDF(CVodeMem cv_mem, realtype hsum, realtype alpha0,
                      realtype alpha0_hat, realtype xi_inv, realtype xistar_inv)
{
  realtype A1, A2, A3, A4, A5, A6;
  realtype C, CPrime, CPrimePrime;
  
  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE + q * A1;
  tq[2] = ABS(alpha0 * (A2 / A1));
  tq[5] = ABS((A2) / (l[q] * xi_inv/xistar_inv));
  if (qwait == 1) {
    C = xistar_inv / l[q];
    A3 = alpha0 + ONE / q;
    A4 = alpha0_hat + xi_inv;
    CPrime = A3 / (ONE - A4 + A3);
    tq[1] = ABS(CPrime / C);
    hsum += tau[q];
    xi_inv = h / hsum;
    A5 = alpha0 - (ONE / (q+1));
    A6 = alpha0_hat - xi_inv;
    CPrimePrime = A2 / (ONE - A6 + A5);
    tq[3] = ABS(CPrimePrime * xi_inv * (q+2) * A5);
  }
  tq[4] = nlscoef * tq[2];
}

/****************** CVnls *****************************************

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 Depending on iter, it calls CVnlsFunctional or CVnlsNewton
 to do the work.

******************************************************************/

static int CVnls(CVodeMem cv_mem, int nflag)
{
  int flag = SOLVED;

  switch(iter) {
  case CV_FUNCTIONAL: 
    flag = CVnlsFunctional(cv_mem);
    break;
  case CV_NEWTON:
    flag = CVnlsNewton(cv_mem, nflag);
    break;
  }

  return(flag);
}

/***************** CVnlsFunctional ********************************

 This routine attempts to solve the nonlinear system using 
 functional iteration (no matrices involved).

******************************************************************/

static int CVnlsFunctional(CVodeMem cv_mem)
{
  int m;
  realtype del, delp, dcon;

  /* Initialize counter and evaluate f at predicted y */
  
  crate = ONE;
  m = 0;
  f(tn, zn[0], tempv, f_data);
  nfe++;
  N_VConst(ZERO, acor);

  /* Loop until convergence; accumulate corrections in acor */

  loop {
    /* Correct y directly from the last f value */
    N_VLinearSum(h, tempv, -ONE, zn[1], tempv);
    N_VScale(rl1, tempv, tempv);
    N_VLinearSum(ONE, zn[0], ONE, tempv, y);
    /* Get WRMS norm of current correction to use in convergence test */
    N_VLinearSum(ONE, tempv, -ONE, acor, acor);
    del = N_VWrmsNorm(acor, ewt);
    N_VScale(ONE, tempv, acor);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) crate = MAX(CRDOWN * crate, del / delp);
    dcon = del * MIN(ONE, crate) / tq[4];
    if (dcon <= ONE) {
      acnrm = (m == 0) ? del : N_VWrmsNorm(acor, ewt);
      return(SOLVED);  /* Convergence achieved */
    }

    /* Stop at maxcor iterations or if iter. seems to be diverging */
    m++;
    if ((m==maxcor) || ((m >= 2) && (del > RDIV * delp)))
      return(CONV_FAIL);
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    f(tn, y, tempv, f_data);
    nfe++;
  }
}

/*********************** CVnlsNewton **********************************

 This routine handles the Newton iteration. It calls lsetup if 
 indicated, calls CVNewtonIteration to perform the iteration, and 
 retries a failed attempt at Newton iteration if that is indicated.
 See return values at top of this file.

**********************************************************************/

static int CVnlsNewton(CVodeMem cv_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, ier;
  booleantype callSetup;
  
  vtemp1 = acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = y;     /* rename y as vtemp2 for readability     */
  vtemp3 = tempv; /* rename tempv as vtemp3 for readability */
  
  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    CV_NO_FAILURES : CV_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (setupNonNull) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (nst == 0) || (nst >= nstlp + MSBP) || (ABS(gamrat-ONE) > DGMAX);
  } else {  
    crate = ONE;
    callSetup = FALSE;
  }
  
  /* Looping point for the solution of the nonlinear system.
     Evaluate f at the predicted y, call lsetup if indicated, and
     call CVNewtonIteration for the Newton iteration itself.      */
  
  loop {

    f(tn, zn[0], ftemp, f_data);
    nfe++; 
    
    if (callSetup) {
      ier = lsetup(cv_mem, convfail, zn[0], ftemp, &jcur, 
                   vtemp1, vtemp2, vtemp3);
      nsetups++;
      callSetup = FALSE;
      gamrat = crate = ONE; 
      gammap = gamma;
      nstlp = nst;
      /* Return if lsetup failed */
      if (ier < 0) return(SETUP_FAIL_UNREC);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, acor);
    N_VScale(ONE, zn[0], y);

    /* Do the Newton iteration */
    ier = CVNewtonIteration(cv_mem);

    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=CV_FAIL_BAD_J.  Otherwise return.                 */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = TRUE;
    convfail = CV_FAIL_BAD_J;
  }
}

/********************** CVNewtonIteration ****************************

 This routine performs the Newton iteration. If the iteration succeeds,
 it returns the value SOLVED. If not, it may signal the CVnlsNewton 
 routine to call lsetup again and reattempt the iteration, by
 returning the value TRY_AGAIN. (In this case, CVnlsNewton must set 
 convfail to CV_FAIL_BAD_J before calling setup again). 
 Otherwise, this routine returns one of the appropriate values 
 SOLVE_FAIL_UNREC or CONV_FAIL back to CVnlsNewton.

*********************************************************************/

static int CVNewtonIteration(CVodeMem cv_mem)
{
  int m, ret;
  realtype del, delp, dcon;
  N_Vector b;
  
  
  mnewt = m = 0;

  /* Looping point for Newton iteration */
  loop {

    /* Evaluate the residual of the nonlinear system*/
    N_VLinearSum(rl1, zn[1], ONE, acor, tempv);
    N_VLinearSum(gamma, ftemp, -ONE, tempv, tempv);

    /* Call the lsolve function */
    b = tempv;
    ret = lsolve(cv_mem, b, ewt, y, ftemp); 
    nni++;
    
    if (ret < 0) return(SOLVE_FAIL_UNREC);
    
    /* If lsolve had a recoverable failure and Jacobian data is
       not current, signal to try the solution again            */
    if (ret > 0) { 
      if ((!jcur) && (setupNonNull)) return(TRY_AGAIN);
      return(CONV_FAIL);
    }

    /* Get WRMS norm of correction; add correction to acor and y */
    del = N_VWrmsNorm(b, ewt);
    N_VLinearSum(ONE, acor, ONE, b, acor);
    N_VLinearSum(ONE, zn[0], ONE, acor, y);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) {
      crate = MAX(CRDOWN * crate, del/delp);
    }
    dcon = del * MIN(ONE, crate) / tq[4];
    
    if (dcon <= ONE) {
      acnrm = (m==0) ? del : N_VWrmsNorm(acor, ewt);
      jcur = FALSE;
      return(SOLVED); /* Nonlinear system was solved successfully */
    }
    
    mnewt = ++m;
    
    /* Stop at maxcor iterations or if iter. seems to be diverging.
       If still not converged and Jacobian data is not current, 
       signal to try the solution again                            */
    if ((m == maxcor) || ((m >= 2) && (del > RDIV*delp))) {
      if ((!jcur) && (setupNonNull)) return(TRY_AGAIN);
      return(CONV_FAIL);
    }
    
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    f(tn, y, ftemp, f_data);
    nfe++;
  }
}

/********************** CVHandleNFlag *******************************

 This routine takes action on the return value nflag = *nflagPtr
 returned by CVnls, as follows:
 
 If CVnls succeeded in solving the nonlinear system, then
 CVHandleNFlag returns the constant DO_ERROR_TEST, which tells CVStep
 to perform the error test.

 If the nonlinear system was not solved successfully, then ncfn and
 ncf = *ncfPtr are incremented and Nordsieck array zn is restored.

 If the solution of the nonlinear system failed due to an
 unrecoverable failure by setup, we return the value SETUP_FAILED.

 If it failed due to an unrecoverable failure in solve, then we return
 the value SOLVE_FAILED.

 Otherwise, a recoverable failure occurred when solving the 
 nonlinear system (CVnls returned nflag == CONV_FAIL). 
   In this case, we return the value REP_CONV_FAIL if ncf is now
   equal to maxncf or |h| = hmin. 
   If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
   PREDICT_AGAIN, telling CVStep to reattempt the step.

*********************************************************************/

static int CVHandleNFlag(CVodeMem cv_mem, int *nflagPtr, realtype saved_t,
                         int *ncfPtr)
{
  int nflag;
  
  nflag = *nflagPtr;
  
  if (nflag == SOLVED) return(DO_ERROR_TEST);

  /* The nonlinear soln. failed; increment ncfn and restore zn */
  ncfn++;
  CVRestore(cv_mem, saved_t);
  
  /* Return if lsetup or lsolve failed unrecoverably */
  if (nflag == SETUP_FAIL_UNREC) return(SETUP_FAILED);
  if (nflag == SOLVE_FAIL_UNREC) return(SOLVE_FAILED);
  
  /* At this point, nflag == CONV_FAIL; increment ncf */
  
  (*ncfPtr)++;
  etamax = ONE;
  /* If we had maxncf failures or |h| = hmin, return REP_CONV_FAIL */
  if ((ABS(h) <= hmin*ONEPSM) || (*ncfPtr == maxncf))
    return(REP_CONV_FAIL);

  /* Reduce step size; return to reattempt the step */
  eta = MAX(ETACF, hmin / ABS(h));
  *nflagPtr = PREV_CONV_FAIL;
  CVRescale(cv_mem);
  return(PREDICT_AGAIN);
}

/********************** CVRestore ************************************

 This routine restores the value of tn to saved_t and undoes the
 prediction.  After execution of CVRestore, the Nordsieck array zn has
 the same values as before the call to CVPredict.

********************************************************************/

static void CVRestore(CVodeMem cv_mem, realtype saved_t)
{
  int j, k;
  
  tn = saved_t;
  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--)
      N_VLinearSum(ONE, zn[j-1], -ONE, zn[j], zn[j-1]);
}

/******************* CVDoErrorTest ********************************

 This routine performs the local error test. 
 The weighted local error norm dsm is loaded into *dsmPtr, and 
 the test dsm ?<= 1 is made.

 If the test passes, CVDoErrorTest returns TRUE. 

 If the test fails, we undo the step just taken (call CVRestore), 
 set *nflagPtr to PREV_ERR_FAIL, and return FALSE. 

 If maxnef error test failures have occurred or if ABS(h) = hmin,
 we set *kflagPtr = REP_ERR_FAIL. (Otherwise *kflagPtr has the
 value last returned by CVHandleNflag.)

 If more than MXNEF1 error test failures have occurred, an order
 reduction is forced.

******************************************************************/

static booleantype CVDoErrorTest(CVodeMem cv_mem, int *nflagPtr, int *kflagPtr,
                                realtype saved_t, int *nefPtr, realtype *dsmPtr)
{
  realtype dsm;
  
  dsm = acnrm / tq[2];

  /* If est. local error norm dsm passes test, return TRUE */  
  *dsmPtr = dsm; 
  if (dsm <= ONE) return(TRUE);
  
  /* Test failed; increment counters, set nflag, and restore zn array */
  (*nefPtr)++;
  netf++;
  *nflagPtr = PREV_ERR_FAIL;
  CVRestore(cv_mem, saved_t);

  /* At maxnef failures or |h| = hmin, return with kflag = REP_ERR_FAIL */
  if ((ABS(h) <= hmin*ONEPSM) || (*nefPtr == maxnef)) {
    *kflagPtr = REP_ERR_FAIL;
    return(FALSE);
  }

  /* Set etamax = 1 to prevent step size increase at end of this step */
  etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  if (*nefPtr <= MXNEF1) {
    eta = ONE / (RPowerR(BIAS2*dsm,ONE/L) + ADDON);
    eta = MAX(ETAMIN, MAX(eta, hmin / ABS(h)));
    if (*nefPtr >= SMALL_NEF) eta = MIN(eta, ETAMXF);
    CVRescale(cv_mem);
    return(FALSE);
  }
  
  /* After MXNEF1 failures, force an order reduction and retry step */
  if (q > 1) {
    eta = MAX(ETAMIN, hmin / ABS(h));
    CVAdjustOrder(cv_mem,-1);
    L = q;
    q--;
    qwait = L;
    CVRescale(cv_mem);
    return(FALSE);
  }

  /* If already at order 1, restart: reload zn from scratch */
  eta = MAX(ETAMIN, hmin / ABS(h));
  h *= eta;
  hscale = h;
  qwait = LONG_WAIT;
  nscon = 0;
  f(tn, zn[0], tempv, f_data);
  nfe++;
  N_VScale(h, tempv, zn[1]);
  return(FALSE);
}

/*************** CVCompleteStep **********************************

 This routine performs various update operations when the solution
 to the nonlinear system has passed the local error test. 
 We increment the step counter nst, record the values hu and qu,
 update the tau array, and apply the corrections to the zn array.
 The tau[i] are the last q values of h, with tau[1] the most recent.
 The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 we save acor and tq[5] for a possible order increase.

******************************************************************/

static void CVCompleteStep(CVodeMem cv_mem)
{
  int i, j;
  
  nst++;
  nscon++;
  hu = h;
  qu = q;

  for (i=q; i >= 2; i--)  tau[i] = tau[i-1];
  if ((q==1) && (nst > 1)) tau[2] = tau[1];
  tau[1] = h;

  for (j=0; j <= q; j++) 
    N_VLinearSum(l[j], acor, ONE, zn[j], zn[j]);
  qwait--;
  if ((qwait == 1) && (q != qmax)) {
    N_VScale(ONE, acor, zn[qmax]);
    saved_tq5 = tq[5];
  }
}

/************* CVPrepareNextStep **********************************

 This routine handles the setting of stepsize and order for the
 next step -- hprime and qprime.  Along with hprime, it sets the
 ratio eta = hprime/h.  It also updates other state variables 
 related to a change of step size or order. 

******************************************************************/

 static void CVPrepareNextStep(CVodeMem cv_mem, realtype dsm)
{
  /* If etamax = 1, defer step size or order changes */
  if (etamax == ONE) {
    qwait = MAX(qwait, 2);
    qprime = q;
    hprime = h;
    eta = ONE;
    return;
  }

  /* etaq is the ratio of new to old h at the current order */  
  etaq = ONE /(RPowerR(BIAS2*dsm,ONE/L) + ADDON);
  
  /* If no order change, adjust eta and acor in CVSetEta and return */
  if (qwait != 0) {
    eta = etaq;
    qprime = q;
    CVSetEta(cv_mem);
    return;
  }
  
  /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are 
     the ratios of new to old h at orders q-1 and q+1, respectively.
     CVChooseEta selects the largest; CVSetEta adjusts eta and acor */
  qwait = 2;
  etaqm1 = CVComputeEtaqm1(cv_mem);
  etaqp1 = CVComputeEtaqp1(cv_mem);  
  CVChooseEta(cv_mem); 
  CVSetEta(cv_mem);
}

/***************** CVSetEta ***************************************

 This routine adjusts the value of eta according to the various
 heuristic limits and the optional input hmax.  It also resets
 etamax to be the estimated local error vector.

*******************************************************************/

static void CVSetEta(CVodeMem cv_mem)
{

  /* If eta below the threshhold THRESH, reject a change of step size */
  if (eta < THRESH) {
    eta = ONE;
    hprime = h;
  } else {
    /* Limit eta by etamax and hmax, then set hprime */
    eta = MIN(eta, etamax);
    eta /= MAX(ONE, ABS(h)*hmax_inv*eta);
    hprime = h * eta;
    if (qprime < q) nscon = 0;
  }
  
  /* Reset etamax for the next step size change, and scale acor */
}

/*************** CVComputeEtaqm1 **********************************

 This routine computes and returns the value of etaqm1 for a
 possible decrease in order by 1.

******************************************************************/

static realtype CVComputeEtaqm1(CVodeMem cv_mem)
{
  realtype ddn;
  
  etaqm1 = ZERO;
  if (q > 1) {
    ddn = N_VWrmsNorm(zn[q], ewt) / tq[1];
    etaqm1 = ONE/(RPowerR(BIAS1*ddn, ONE/q) + ADDON);
  }
  return(etaqm1);
}

/*************** CVComputeEtaqp1 **********************************

 This routine computes and returns the value of etaqp1 for a
 possible increase in order by 1.

******************************************************************/

static realtype CVComputeEtaqp1(CVodeMem cv_mem)
{
  realtype dup, cquot;
  
  etaqp1 = ZERO;
  if (q != qmax) {
    cquot = (tq[5] / saved_tq5) * RPowerI(h/tau[2], L);
    N_VLinearSum(-cquot, zn[qmax], ONE, acor, tempv);
    dup = N_VWrmsNorm(tempv, ewt) /tq[3];
    etaqp1 = ONE / (RPowerR(BIAS3*dup, ONE/(L+1)) + ADDON);
  }
  return(etaqp1);
}

/******************* CVChooseEta **********************************

 Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 q - 1, q, or q + 1, respectively), this routine chooses the 
 maximum eta value, sets eta to that value, and sets qprime to the
 corresponding value of q.  If there is a tie, the preference
 order is to (1) keep the same order, then (2) decrease the order,
 and finally (3) increase the order.  If the maximum eta value
 is below the threshhold THRESH, the order is kept unchanged and
 eta is set to 1.

******************************************************************/

static void CVChooseEta(CVodeMem cv_mem)
{
  realtype etam;
  
  etam = MAX(etaqm1, MAX(etaq, etaqp1));
  
  if (etam < THRESH) {
    eta = ONE;
    qprime = q;
    return;
  }

  if (etam == etaq) {
    eta = etaq;
    qprime = q;
  } else if (etam == etaqm1) {
    eta = etaqm1;
    qprime = q - 1;
  } else {
    eta = etaqp1;
    qprime = q + 1;
    if (lmm == CV_BDF) N_VScale(ONE, acor, zn[qmax]);
  }
}

/****************** CVHandleFailure ******************************

 This routine prints error messages for all cases of failure by
 CVStep. It returns to CVode the value that CVode is to return to
 the user.

*****************************************************************/

static int CVHandleFailure(CVodeMem cv_mem, int kflag)
{

  /* Set vector of  absolute weighted local errors */
  N_VProd(acor, ewt, tempv);
  N_VAbs(tempv, tempv);

  /* Depending on kflag, print error message and return error flag */
  switch (kflag) {
    case REP_ERR_FAIL: 
      if(errfp!=NULL) fprintf(errfp, MSGCV_ERR_FAILS, tn, h);
      return(CV_ERR_FAILURE);
    case REP_CONV_FAIL:
      if(errfp!=NULL) fprintf(errfp, MSGCV_CONV_FAILS, tn, h);
      return(CV_CONV_FAILURE);
    case SETUP_FAILED:
      if(errfp!=NULL) fprintf(errfp, MSGCV_SETUP_FAILED, tn);
      return(CV_LSETUP_FAIL);
    case SOLVE_FAILED:
      if(errfp!=NULL) fprintf(errfp, MSGCV_SOLVE_FAILED, tn);
      return(CV_LSOLVE_FAIL);
  }
  return(0);
}

/****************** CVBDFStab ***********************************
 This routine handles the BDF Stability Limit Detection Algorithm
 STALD.  It is called if lmm = CV_BDF and the SLDET option is on.
 If the order is 3 or more, the required norm data is saved.
 If a decision to reduce order has not already been made, and
 enough data has been saved, CVsldet is called.  If it signals
 a stability limit violation, the order is reduced, and the step
 size is reset accordingly.

*****************************************************************/

void CVBDFStab(CVodeMem cv_mem)
{
  int i,k, ldflag, factorial;
  realtype sq, sqm1, sqm2;
      
  /* If order is 3 or greater, then save scaled derivative data,
     push old data down in i, then add current values to top.    */

  if (q >= 3) {
    for (k = 1; k <= 3; k++)
      { for (i = 5; i >= 2; i--) ssdat[i][k] = ssdat[i-1][k]; }
    factorial = 1;
    for (i = 1; i <= q-1; i++) factorial *= i;
    sq = factorial*q*(q+1)*acnrm/tq[5];
    sqm1 = factorial*q*N_VWrmsNorm(zn[q], ewt);
    sqm2 = factorial*N_VWrmsNorm(zn[q-1], ewt);
    ssdat[1][1] = sqm2*sqm2;
    ssdat[1][2] = sqm1*sqm1;
    ssdat[1][3] = sq*sq;
  }  

  if (qprime >= q) {

    /* If order is 3 or greater, and enough ssdat has been saved,
       nscon >= q+5, then call stability limit detection routine.  */

    if ( (q >= 3) && (nscon >= q+5) ) {
      ldflag = CVsldet(cv_mem);
      if (ldflag > 3) {
        /* A stability limit violation is indicated by
           a return flag of 4, 5, or 6.
           Reduce new order.                     */
        qprime = q-1;
        eta = etaqm1; 
        eta = MIN(eta,etamax);
        eta = eta/MAX(ONE,ABS(h)*hmax_inv*eta);
        hprime = h*eta;
        nor = nor + 1;
      }
    }
  }
  else {
    /* Otherwise, let order increase happen, and 
       reset stability limit counter, nscon.     */
    nscon = 0;
  }
}

/********************* CVsldet ************************************
  This routine detects stability limitation using stored scaled 
  derivatives data. CVsldet returns the magnitude of the
  dominate characteristic root, rr. The presents of a stability
  limit is indicated by rr > "something a little less then 1.0",  
  and a positive kflag. This routine should only be called if
  order is greater than or equal to 3, and data has been collected
  for 5 time steps. 
 
  Returned values:
     kflag = 1 -> Found stable characteristic root, normal matrix case
     kflag = 2 -> Found stable characteristic root, quartic solution
     kflag = 3 -> Found stable characteristic root, quartic solution,
                  with Newton correction
     kflag = 4 -> Found stability violation, normal matrix case
     kflag = 5 -> Found stability violation, quartic solution
     kflag = 6 -> Found stability violation, quartic solution,
                  with Newton correction

     kflag < 0 -> No stability limitation, 
                  or could not compute limitation.

     kflag = -1 -> Min/max ratio of ssdat too small.
     kflag = -2 -> For normal matrix case, vmax > vrrt2*vrrt2
     kflag = -3 -> For normal matrix case, The three ratios
                   are inconsistent.
     kflag = -4 -> Small coefficient prevents elimination of quartics.  
     kflag = -5 -> R value from quartics not consistent.
     kflag = -6 -> No corrected root passes test on qk values
     kflag = -7 -> Trouble solving for sigsq.
     kflag = -8 -> Trouble solving for B, or R via B.
     kflag = -9 -> R via sigsq[k] disagrees with R from data.

********************************************************************/

static int CVsldet(CVodeMem cv_mem)
{
  int i, k, j, it, kmin, kflag = 0;
  realtype rat[5][4], rav[4], qkr[4], sigsq[4], smax[4], ssmax[4];
  realtype drr[4], rrc[4],sqmx[4], qjk[4][4], vrat[5], qc[6][4], qco[6][4];
  realtype rr, rrcut, vrrtol, vrrt2, sqtol, rrtol;
  realtype smink, smaxk, sumrat, sumrsq, vmin, vmax, drrmax, adrr;
  realtype small, tem, sqmax, saqk, qp, s, sqmaxk, saqj, sqmin;
  realtype rsa, rsb, rsc, rsd, rse, rd1a, rd1b, rd1c, rd1d;
  realtype rd2a, rd2b, rd2c, rd3a, rd3b, cest1, corr1; 
  realtype ratp, ratm, qfac1, qfac2, bb, rrb;

  /* The following are cutoffs and tolerances used by this routine */

  rrcut  = RCONST(0.98);
  vrrtol = RCONST(1.0e-4);
  vrrt2  = RCONST(5.0e-4);
  sqtol  = RCONST(1.0e-3);
  rrtol  = RCONST(1.0e-2);
  
  rr = ZERO;
  
  /*  Index k corresponds to the degree of the interpolating polynomial. */
  /*      k = 1 -> q-1          */
  /*      k = 2 -> q            */
  /*      k = 3 -> q+1          */
  
  /*  Index i is a backward-in-time index, i = 1 -> current time, */
  /*      i = 2 -> previous step, etc    */
  
  /* get maxima, minima, and variances, and form quartic coefficients  */
  
  for (k=1; k<=3; k++) {
    smink = ssdat[1][k];
    smaxk = ZERO;
    
    for (i=1; i<=5; i++) {
      smink = MIN(smink,ssdat[i][k]);
      smaxk = MAX(smaxk,ssdat[i][k]);
    }
    
    if (smink < TINY*smaxk) {
      kflag = -1;  
      return(kflag);
    }
    smax[k] = smaxk;
    ssmax[k] = smaxk*smaxk;
    
    sumrat = ZERO;
    sumrsq = ZERO;
    for (i=1; i<=4; i++) {
      rat[i][k] = ssdat[i][k]/ssdat[i+1][k];
      sumrat = sumrat + rat[i][k];
      sumrsq = sumrsq + rat[i][k]*rat[i][k];
    } 
    rav[k] = FOURTH*sumrat;
    vrat[k] = ABS(FOURTH*sumrsq - rav[k]*rav[k]);
    
    qc[5][k] = ssdat[1][k]*ssdat[3][k] - ssdat[2][k]*ssdat[2][k];
    qc[4][k] = ssdat[2][k]*ssdat[3][k] - ssdat[1][k]*ssdat[4][k];
    qc[3][k] = ZERO;
    qc[2][k] = ssdat[2][k]*ssdat[5][k] - ssdat[3][k]*ssdat[4][k];
    qc[1][k] = ssdat[4][k]*ssdat[4][k] - ssdat[3][k]*ssdat[5][k];
    
    for (i=1; i<=5; i++) {
      qco[i][k] = qc[i][k];
    }
  }                            /* End of k loop */
  
  /* Isolate normal or nearly-normal matrix case. Three quartic will
     have common or nearly-common roots in this case. 
     Return a kflag = 1 if this procedure works. If three root 
     differ more than vrrt2, return error kflag = -3.    */
  
  vmin = MIN(vrat[1],MIN(vrat[2],vrat[3]));
  vmax = MAX(vrat[1],MAX(vrat[2],vrat[3]));
  
  if(vmin < vrrtol*vrrtol) {
    if (vmax > vrrt2*vrrt2) {
      kflag = -2;  
      return(kflag);
    } else {
      rr = (rav[1] + rav[2] + rav[3])/THREE;
      
      drrmax = ZERO;
      for(k = 1;k<=3;k++) {
        adrr = ABS(rav[k] - rr);
        drrmax = MAX(drrmax, adrr);
      }
      if (drrmax > vrrt2) {
        kflag = -3;    
      }
      
      kflag = 1;

      /*  can compute charactistic root, drop to next section   */
      
    }
  } else {

    /* use the quartics to get rr. */
    
    if (ABS(qco[1][1]) < TINY*ssmax[1]) {
      small = qco[1][1];
      kflag = -4;    
      return(kflag);
    }
    
    tem = qco[1][2]/qco[1][1];
    for(i=2; i<=5; i++) {
      qco[i][2] = qco[i][2] - tem*qco[i][1];
    }

    qco[1][2] = ZERO;
    tem = qco[1][3]/qco[1][1];
    for(i=2; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][1];
    }
    qco[1][3] = ZERO;
    
    if (ABS(qco[2][2]) < TINY*ssmax[2]) {
      small = qco[2][2];
      kflag = -4;    
      return(kflag);
    }
    
    tem = qco[2][3]/qco[2][2];
    for(i=3; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][2];
    }
    
    if (ABS(qco[4][3]) < TINY*ssmax[3]) {
      small = qco[4][3];
      kflag = -4;    
      return(kflag);
    }
    
    rr = -qco[5][3]/qco[4][3];
    
    if (rr < TINY || rr > HUN) {
      kflag = -5;   
      return(kflag);
    }
    
    for(k=1; k<=3; k++) {
      qkr[k] = qc[5][k] + rr*(qc[4][k] + rr*rr*(qc[2][k] + rr*qc[1][k]));
    }  
    
    sqmax = ZERO;
    for(k=1; k<=3; k++) {
      saqk = ABS(qkr[k])/ssmax[k];
      if (saqk > sqmax) sqmax = saqk;
    } 
    
    if (sqmax < sqtol) {
      kflag = 2;
      
      /*  can compute charactistic root, drop to "given rr,etc"   */
      
    } else {

      /* do Newton corrections to improve rr.  */
      
      for(it=1; it<=3; it++) {
        for(k=1; k<=3; k++) {
          qp = qc[4][k] + rr*rr*(THREE*qc[2][k] + rr*FOUR*qc[1][k]);
          drr[k] = ZERO;
          if (ABS(qp) > TINY*ssmax[k]) drr[k] = -qkr[k]/qp;
          rrc[k] = rr + drr[k];
        } 
        
        for(k=1; k<=3; k++) {
          s = rrc[k];
          sqmaxk = ZERO;
          for(j=1; j<=3; j++) {
            qjk[j][k] = qc[5][j] + s*(qc[4][j] + 
                                      s*s*(qc[2][j] + s*qc[1][j]));
            saqj = ABS(qjk[j][k])/ssmax[j];
            if (saqj > sqmaxk) sqmaxk = saqj;
          } 
          sqmx[k] = sqmaxk;
        } 

        sqmin = sqmx[1]; kmin = 1;
        for(k=2; k<=3; k++) {
          if (sqmx[k] < sqmin) {
            kmin = k;
            sqmin = sqmx[k];
          }
        } 
        rr = rrc[kmin];
        
        if (sqmin < sqtol) {
          kflag = 3;
          /*  can compute charactistic root   */
          /*  break out of Newton correction loop and drop to "given rr,etc" */ 
          break;
        } else {
          for(j=1; j<=3; j++) {
            qkr[j] = qjk[j][kmin];
          }
        }     
      }          /*  end of Newton correction loop  */ 
      
      if (sqmin > sqtol) {
        kflag = -6;
        return(kflag);
      }
    }     /*  end of if (sqmax < sqtol) else   */
  }      /*  end of if(vmin < vrrtol*vrrtol) else, quartics to get rr. */
  
  /* given rr, find sigsq[k] and verify rr.  */
  /* All positive kflag drop to this section  */
  
  for(k=1; k<=3; k++) {
    rsa = ssdat[1][k];
    rsb = ssdat[2][k]*rr;
    rsc = ssdat[3][k]*rr*rr;
    rsd = ssdat[4][k]*rr*rr*rr;
    rse = ssdat[5][k]*rr*rr*rr*rr;
    rd1a = rsa - rsb;
    rd1b = rsb - rsc;
    rd1c = rsc - rsd;
    rd1d = rsd - rse;
    rd2a = rd1a - rd1b;
    rd2b = rd1b - rd1c;
    rd2c = rd1c - rd1d;
    rd3a = rd2a - rd2b;
    rd3b = rd2b - rd2c;
    
    if (ABS(rd1b) < TINY*smax[k]) {
      kflag = -7;
      return(kflag);
    }
    
    cest1 = -rd3a/rd1b;
    if (cest1 < TINY || cest1 > FOUR) {
      kflag = -7;
      return(kflag);
    }
    corr1 = (rd2b/cest1)/(rr*rr);
    sigsq[k] = ssdat[3][k] + corr1;
  }
  
  if (sigsq[2] < TINY) {
    kflag = -8;
    return(kflag);
  }
  
  ratp = sigsq[3]/sigsq[2];
  ratm = sigsq[1]/sigsq[2];
  qfac1 = FOURTH*(q*q - ONE);
  qfac2 = TWO/(q - ONE);
  bb = ratp*ratm - ONE - qfac1*ratp;
  tem = ONE - qfac2*bb;
  
  if (ABS(tem) < TINY) {
    kflag = -8;
    return(kflag);
  }
  
  rrb = ONE/tem;
  
  if (ABS(rrb - rr) > rrtol) {
    kflag = -9;
    return(kflag);
  }
  
  /* Check to see if rr is above cutoff rrcut  */
  if (rr > rrcut) {
    if (kflag == 1) kflag = 4;
    if (kflag == 2) kflag = 5;
    if (kflag == 3) kflag = 6;
  }
  
  /* All positive kflag returned at this point  */
  
  return(kflag);
  
}

/********************* CVRcheck1 ***********************************
 
 This routine completes the initialization of rootfinding memory
 information, and checks whether g has a zero both at and very near
 the initial point of the IVP.

 This routine returns an int equal to:
   INITROOT = -1 if a close pair of zeros was found, and
   CV_SUCCESS     =  0 otherwise.

********************************************************************/

static int CVRcheck1(CVodeMem cv_mem)
{
  int i;
  realtype smallh, hratio;
  booleantype zroot;

  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  tlo = tn;
  ttol = (ABS(tn) + ABS(h))*uround*HUN;

  /* Evaluate g at initial t and check for zero values. */
  gfun (tlo, zn[0], glo, g_data);
  nge = 1;
  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(glo[i]) == ZERO) zroot = TRUE;
  }
  if (!zroot) return(CV_SUCCESS);

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  smallh = (h > ZERO) ? ttol : -ttol;
  tlo += smallh;
  hratio = smallh/h;
  N_VLinearSum(ONE, zn[0], hratio, zn[1], y);
  gfun (tlo, y, glo, g_data);  nge++;
  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(glo[i]) == ZERO) {
      zroot = TRUE;
      iroots[i] = 1;
    }
  }
  if (zroot) return(INITROOT);
  return(CV_SUCCESS);

}

/********************* CVRcheck2 ***********************************
 
 This routine checks for exact zeros of g at the last root found,
 if the last return was a root.  It then checks for a close
 pair of zeros (an error condition), and for a new root at a
 nearby point.  The left endpoint (tlo) of the search interval
 is adjusted if necessary to assure that all g_i are nonzero
 there, before returning to do a root search in the interval.

 On entry, tlo = tretlast is the last value of tret returned by
 CVode.  This may be the previous tn, the previous tout value, or
 the last root location.

 This routine returns an int equal to:
      CLOSERT = -2 if a close pair of zeros was found,
      RTFOUND =  1 if a new zero of g was found near tlo, or
      CV_SUCCESS    =  0 otherwise.

********************************************************************/

static int CVRcheck2(CVodeMem cv_mem)
{
  int i;
  realtype smallh, hratio;
  booleantype zroot;

  if (irfnd == 0) return(CV_SUCCESS);

  (void) CVodeGetDky(cv_mem, tlo, 0, y);
  gfun (tlo, y, glo, g_data);  nge++;
  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(glo[i]) == ZERO) {
      zroot = TRUE;
      iroots[i] = 1;
    }
  }
  if (!zroot) return(CV_SUCCESS);

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  ttol = (ABS(tn) + ABS(h))*uround*HUN;
  smallh = (h > ZERO) ? ttol : -ttol;
  tlo += smallh;
  if ( (tlo - tn)*h >= ZERO) {
    hratio = smallh/h;
    N_VLinearSum(ONE, y, hratio, zn[1], y);
  } else {
    (void) CVodeGetDky(cv_mem, tlo, 0, y);
  }
  gfun (tlo, y, glo, g_data);  nge++;
  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(glo[i]) == ZERO) {
      if (iroots[i] == 1) return(CLOSERT);
      zroot = TRUE;
      iroots[i] = 1;
    }
  }
  if (zroot) return(RTFOUND);
  return(CV_SUCCESS);

}

/********************* CVRcheck3 ***********************************
 
 This routine interfaces to CVRootfind to look for a root of g
 between tlo and either tn or tout, whichever comes first.
 Only roots beyond tlo in the direction of integration are sought.

 This routine returns an int equal to:
      RTFOUND =  1 if a root of g was found, or
      CV_SUCCESS    =  0 otherwise.

********************************************************************/

static int CVRcheck3(CVodeMem cv_mem)
{
  int i, ier;

  /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
  if (taskc == CV_ONE_STEP) {
    thi = tn;
    N_VScale(ONE, zn[0], y);
  }
  if (taskc == CV_NORMAL) {
    if ( (toutc - tn)*h >= ZERO) {
      thi = tn; 
      N_VScale(ONE, zn[0], y);
    } else {
      thi = toutc;
      (void) CVodeGetDky(cv_mem, thi, 0, y);
    }
  }

  /* Set ghi = g(thi) and call CVRootfind to search (tlo,thi) for roots. */
  gfun (thi, y, ghi, g_data);  nge++;
  ttol = (ABS(tn) + ABS(h))*uround*HUN;
  ier = CVRootfind(cv_mem);
  tlo = troot;
  for (i = 0; i < nrtfn; i++) glo[i] = groot[i];

  /* If no root found, return CV_SUCCESS. */  
  if (ier == CV_SUCCESS) return(CV_SUCCESS);

  /* If a root was found, interpolate to get y(troot) and return.  */
  (void) CVodeGetDky(cv_mem, troot, 0, y);
  return(RTFOUND);

}

/********************* CVRootfind **********************************
 
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
            (void) gfun(t, y, gt, g_data)

 nge      = cumulative counter for gfun calls.

 ttol     = a convergence tolerance for troot.  Input only.
            When a root at troot is found, it is located only to
            within a tolerance of ttol.  Typically, ttol should
            be set to a value on the order of
               100 * UROUND * max (ABS(tlo), ABS(thi))
            where UROUND is the unit roundoff of the machine.

 tlo, thi = endpoints of the interval in which roots are sought.
            On input, and must be distinct, but tlo - thi may
            be of either sign.  The direction of integration is
            assumed to be from tlo to thi.  On return, tlo and thi
            are the endpoints of the final relevant interval.

 glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
            and g(thi) respectively.  Input and output.  On input,
            none of the glo[i] should be zero.

 troot    = root location, if a root was found, or thi if not.
            Output only.  If a root was found other than an exact
            zero of g, troot is the endpoint thi of the final
            interval bracketing the root, with size at most ttol.

 groot    = array of length nrtfn containing g(troot) on return.

 iroots   = int array of length nrtfn with root information.
            Output only.  If a root was found, iroots indicates
            which components g_i have a root at troot.  For
            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
            and iroots[i] = 0 otherwise.

 This routine returns an int equal to:
      RTFOUND = 1 if a root of g was found, or
      CV_SUCCESS = 0 otherwise.

********************************************************************/

static int CVRootfind(CVodeMem cv_mem)
{
  realtype alpha, tmid, gfrac, maxfrac, fracint, fracsub;
  int i, imax, side, sideprev;
  booleantype zroot, sgnchg;

  imax = 0;

  /* First check for change in sign in ghi or for a zero in ghi. */
  maxfrac = ZERO;
  zroot = FALSE;
  sgnchg = FALSE;
  for (i = 0;  i < nrtfn; i++) {
    if (ABS(ghi[i]) == ZERO) {
      zroot = TRUE;
    } else {
      if (glo[i]*ghi[i] < ZERO) {
        gfrac = ABS(ghi[i]/(ghi[i] - glo[i]));
        if (gfrac > maxfrac) {
          sgnchg = TRUE;
          maxfrac = gfrac;
          imax = i;
        }
      }
    }
  }

  /* If no sign change was found, reset troot and groot.  Then return
     CV_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */ 
  if (!sgnchg) {
    troot = thi;
    for (i = 0; i < nrtfn; i++) groot[i] = ghi[i];
    if (!zroot) return(CV_SUCCESS);
    for (i = 0; i < nrtfn; i++) {
      iroots[i] = 0;
      if (ABS(ghi[i]) == ZERO) iroots[i] = 1;
    }
    return(RTFOUND);
  }

  /* A sign change was found.  Loop to locate nearest root. */

  side = 0;  sideprev = -1;
  loop {                                    /* Looping point */

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
    tmid = thi - (thi - tlo)*ghi[imax]/(ghi[imax] - alpha*glo[imax]);
    if (ABS(tmid - tlo) < HALF*ttol) {
      fracint = ABS(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = tlo + fracsub*(thi - tlo);
    }
    if (ABS(thi - tmid) < HALF*ttol) {
      fracint = ABS(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = thi - fracsub*(thi - tlo);
    }

    (void) CVodeGetDky(cv_mem, tmid, 0, y);
    gfun (tmid, y, groot, g_data);  nge++;

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */  
    maxfrac = ZERO;
    zroot = FALSE;
    sgnchg = FALSE;
    sideprev = side;
    for (i = 0;  i < nrtfn; i++) {
      if (ABS(groot[i]) == ZERO) {
        zroot = TRUE;
      } else {
        if (glo[i]*groot[i] < ZERO) {
          gfrac = ABS(groot[i]/(groot[i] - glo[i]));
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
      thi = tmid;
      for (i = 0; i < nrtfn; i++) ghi[i] = groot[i];
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (ABS(thi - tlo) <= ttol) break;
    continue;  /* Return to looping point. */
    }

    if (zroot) {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      thi = tmid;
      for (i = 0; i < nrtfn; i++) ghi[i] = groot[i];
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    tlo = tmid;
    for (i = 0; i < nrtfn; i++) glo[i] = groot[i];
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (ABS(thi - tlo) <= ttol) break;

  } /* End of root-search loop */

  /* Reset troot and groot, set iroots, and return RTFOUND. */
  troot = thi;
  for (i = 0; i < nrtfn; i++) {
    groot[i] = ghi[i];
    iroots[i] = 0;
    if (ABS(ghi[i]) == ZERO) iroots[i] = 1;
    if (glo[i]*ghi[i] < ZERO) iroots[i] = 1;
  }
  return(RTFOUND);
}

/*=================================================================*/
/*      Internal EWT function                                      */
/*=================================================================*/

/*
 * CVEwtSet
 *
 * This routine is responsible for setting the error weight vector ewt,
 * according to tol_type, as follows:
 *
 * (1) ewt[i] = 1 / (reltol * ABS(ycur[i]) + *abstol), i=0,...,neq-1
 *     if tol_type = CV_SS
 * (2) ewt[i] = 1 / (reltol * ABS(ycur[i]) + abstol[i]), i=0,...,neq-1
 *     if tol_type = CV_SV
 *
 * CVEwtSet returns 0 if ewt is successfully set as above to a
 * positive vector and -1 otherwise. In the latter case, ewt is
 * considered undefined.
 *
 * All the real work is done in the routines CVEwtSetSS, CVEwtSetSV.
 */

int CVEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  CVodeMem cv_mem;
  int flag = 0;

  /* data points to cv_mem here */

  cv_mem = (CVodeMem) data;

  switch(itol) {
  case CV_SS: 
    flag = CVEwtSetSS(ycur, weight, cv_mem);
    break;
  case CV_SV: 
    flag = CVEwtSetSV(ycur, weight, cv_mem);
    break;
  }
  
  return(flag);
}

/*
 * CVEwtSetSS
 *
 * This routine sets ewt as decribed above in the case tol_type = CV_SS.
 * It tests for non-positive components before inverting. CVEwtSetSS
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */

static int CVEwtSetSS(N_Vector ycur, N_Vector weight, CVodeMem cv_mem)
{
  N_VAbs(ycur, tempv);
  N_VScale(reltol, tempv, tempv);
  N_VAddConst(tempv, Sabstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(-1);
  N_VInv(tempv, weight);
  return(0);
}

/*
 * CVEwtSetSV
 *
 * This routine sets ewt as decribed above in the case tol_type = CV_SV.
 * It tests for non-positive components before inverting. CVEwtSetSV
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */

static int CVEwtSetSV(N_Vector ycur, N_Vector weight, CVodeMem cv_mem)
{
  N_VAbs(ycur, tempv);
  N_VLinearSum(reltol, tempv, ONE, Vabstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(-1);
  N_VInv(tempv, weight);
  return(0);
}

