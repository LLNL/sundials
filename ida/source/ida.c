/*
 * -----------------------------------------------------------------
 * $Revision: 1.44 $
 * $Date: 2005-04-19 21:13:57 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the main IDA solver.
 * It is independent of the linear solver in use.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "sundialsmath.h"

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define ZERO       RCONST(0.0)    /* real 0.0    */
#define HALF       RCONST(0.5)    /* real 0.5    */
#define QUARTER    RCONST(0.25)   /* real 0.25   */
#define TWOTHIRDS  RCONST(0.667)  /* real 2/3    */
#define ONE        RCONST(1.0)    /* real 1.0    */
#define ONEPT5     RCONST(1.5)    /* real 1.5    */
#define TWO        RCONST(2.0)    /* real 2.0    */
#define TEN        RCONST(10.0)   /* real 10.0   */
#define TWELVE     RCONST(12.0)   /* real 12.0   */
#define TWENTY     RCONST(20.0)   /* real 20.0   */
#define HUNDRED    RCONST(100.0)  /* real 100.0  */
#define PT9        RCONST(0.9)    /* real 0.9    */
#define PT99       RCONST(0.99)   /* real 0.99   */
#define PT1        RCONST(0.1)    /* real 0.1    */
#define PT01       RCONST(0.01)   /* real 0.01   */
#define PT001      RCONST(0.001)  /* real 0.001  */
#define PT0001     RCONST(0.0001) /* real 0.0001 */

/*
 * -----------------------------------------------------------------
 * default constants
 * -----------------------------------------------------------------
 */

#define MXNCF           10  /* max number of convergence failures allowed */
#define MXNEF           10  /* max number of error test failures allowed  */
#define MAXNH            5  /* max. number of h tries in IC calc. */
#define MAXNJ            4  /* max. number of J tries in IC calc. */
#define MAXNI           10  /* max. Newton iterations in IC calc. */
#define EPCON RCONST(0.33)  /* Newton convergence test constant */

/*
 * -----------------------------------------------------------------
 * routine-specific constants
 * -----------------------------------------------------------------
 */

/* IDAStep control constants */

#define PREDICT_AGAIN    20

/* IDANewtonIter constants */

#define MAXIT    4
#define RATEMAX  RCONST(0.9)
#define XRATE    RCONST(0.25)        

/* Return values for lower level routines used by IDASolve */

#define IDA_RES_RECVR    +1
#define IDA_LSETUP_RECVR +2
#define IDA_LSOLVE_RECVR +3

#define IDA_NCONV_RECVR  +4
#define IDA_CONSTR_RECVR +5
#define CONTINUE_STEPS   +99

/* IDACompleteStep constants */

#define UNSET    -1
#define LOWER     1 
#define RAISE     2 
#define MAINTAIN  3

/* IDATestError constants */

#define ERROR_TEST_FAIL +7

/* Macro: loop */

#define loop for(;;)

/*
 * -----------------------------------------------------------------
 * private helper function prototypes
 * -----------------------------------------------------------------
 */

static booleantype IDACheckNvector(N_Vector tmpl);

static booleantype IDAAllocVectors(IDAMem IDA_mem, N_Vector tmpl, int tol);
static void IDAFreeVectors(IDAMem IDA_mem);

realtype IDAWrmsNorm(IDAMem IDA_mem, N_Vector x, N_Vector w, 
                     booleantype mask);
int IDAInitialSetup(IDAMem IDA_mem);

static int IDAEwtSetSS(N_Vector ycur, N_Vector weight, IDAMem IDA_mem);
static int IDAEwtSetSV(N_Vector ycur, N_Vector weight, IDAMem IDA_mem);

static int IDAStopTest1(IDAMem IDA_mem, realtype tout,realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask);
static int IDAStopTest2(IDAMem IDA_mem, realtype tout, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask);
static int IDAHandleFailure(IDAMem IDA_mem, int sflag);

static int IDAStep(IDAMem IDA_mem);
static void IDASetCoeffs(IDAMem IDA_mem, realtype *ck);
static int IDAnls(IDAMem IDA_mem);
static int IDAPredict(IDAMem IDA_mem);
static int IDANewtonIter(IDAMem IDA_mem);
static int IDATestError(IDAMem IDA_mem, realtype *ck, realtype *est,
                        realtype *terk, realtype *terkm1, realtype *erkm1);
static int IDAHandleNFlag(IDAMem IDA_mem, int nflag, realtype saved_t,
                          int *ncfPtr, int *nefPtr, realtype *est);
static int IDACompleteStep(IDAMem IDA_mem, realtype *est, 
                           realtype *terk, realtype *terkm1, realtype *erkm1);

/*
 * -----------------------------------------------------------------
 * user-callable functions
 * -----------------------------------------------------------------
 */

/* 
 * -----------------------------------------------------------------
 * IDACreate
 * -----------------------------------------------------------------
 * IDACreate creates an internal memory block for a problem to 
 * be solved by IDA.
 * If successful, IDACreate returns a pointer to the problem memory. 
 * This pointer should be passed to IDAMalloc.  
 * If an initialization error occurs, IDACreate prints an error 
 * message to standard err and returns NULL. 
 * -----------------------------------------------------------------
*/

void *IDACreate(void)
{
  IDAMem IDA_mem;

  IDA_mem = (IDAMem) malloc(sizeof(struct IDAMemRec));
  if (IDA_mem == NULL) {
    fprintf(stderr, MSG_MEM_FAIL);
    return (NULL);
  }

  /* Set unit roundoff in IDA_mem */
  IDA_mem->ida_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  IDA_mem->ida_res         = NULL;
  IDA_mem->ida_rdata       = NULL;
  IDA_mem->ida_efun        = NULL;
  IDA_mem->ida_edata       = NULL;
  IDA_mem->ida_errfp       = stderr;
  IDA_mem->ida_maxord      = MAXORD_DEFAULT;
  IDA_mem->ida_mxstep      = MXSTEP_DEFAULT;
  IDA_mem->ida_hmax_inv    = HMAX_INV_DEFAULT;
  IDA_mem->ida_hin         = ZERO;
  IDA_mem->ida_epcon       = EPCON;
  IDA_mem->ida_maxnef      = MXNEF;
  IDA_mem->ida_maxncf      = MXNCF;
  IDA_mem->ida_maxcor      = MAXIT;
  IDA_mem->ida_suppressalg = FALSE;
  IDA_mem->ida_id          = NULL;
  IDA_mem->ida_constraints = NULL;
  IDA_mem->ida_constraintsSet = FALSE;
  IDA_mem->ida_tstopset    = FALSE;

  /* Set default values for IC optional inputs */
  IDA_mem->ida_epiccon = PT01 * EPCON;
  IDA_mem->ida_maxnh   = MAXNH;
  IDA_mem->ida_maxnj   = MAXNJ;
  IDA_mem->ida_maxnit  = MAXNI;
  IDA_mem->ida_lsoff   = FALSE;
  IDA_mem->ida_steptol = RPowerR(IDA_mem->ida_uround, TWOTHIRDS);

  /* Initialize lrw and liw */
  IDA_mem->ida_lrw = 25 + 5*MXORDP1;
  IDA_mem->ida_liw = 38;

  /* No mallocs have been done yet */
  IDA_mem->ida_VatolMallocDone = FALSE;
  IDA_mem->ida_constraintsMallocDone = FALSE;
  IDA_mem->ida_idMallocDone = FALSE;
  IDA_mem->ida_MallocDone = FALSE;

  /* Return pointer to IDA memory block */
  return((void *)IDA_mem);
}

/*-----------------------------------------------------------------*/

#define errfp (IDA_mem->ida_errfp)
#define lrw   (IDA_mem->ida_lrw)
#define liw   (IDA_mem->ida_liw)

/*-----------------------------------------------------------------*/

/*
 * -----------------------------------------------------------------
 * IDAMalloc
 * -----------------------------------------------------------------
 * IDAMalloc allocates and initializes memory for a problem. All
 * problem specification inputs are checked for errors. If any
 * error occurs during initialization, it is reported to the file
 * whose file pointer is errfp and an error flag is returned. 
 * -----------------------------------------------------------------
 */

int IDAMalloc(void *ida_mem, IDAResFn res,
              realtype t0, N_Vector yy0, N_Vector yp0, 
              int itol, realtype rtol, void *atol)
{
  IDAMem IDA_mem;
  booleantype nvectorOK, allocOK, neg_atol;
  long int lrw1, liw1;

  /* Check ida_mem */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAM_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  
  /* Check for legal input parameters */
  
  if (yy0 == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_Y0_NULL); 
    return(IDA_ILL_INPUT); 
  }
  
  if (yp0 == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_YP0_NULL); 
    return(IDA_ILL_INPUT); 
  }

  if ((itol != IDA_SS) && (itol != IDA_SV) && (itol != IDA_WF)) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_ITOL);
    return(IDA_ILL_INPUT);
  }

  if (res == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_RES_NULL); 
    return(IDA_ILL_INPUT); 
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = IDACheckNvector(yy0);
  if(!nvectorOK) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_NVECTOR);
    return(IDA_ILL_INPUT);
  }

  /* Test tolerances */

  if (itol != IDA_WF) {

    if (atol == NULL) { 
      if(errfp!=NULL) fprintf(errfp, MSG_ATOL_NULL); 
      return(IDA_ILL_INPUT); 
    }

    if (rtol < ZERO) { 
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_RTOL); 
      return(IDA_ILL_INPUT); 
    }
   
    if (itol == IDA_SS) { 
      neg_atol = (*((realtype *)atol) < ZERO); 
    } else { 
      neg_atol = (N_VMin((N_Vector)atol) < ZERO); 
    }

    if (neg_atol) { 
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_ATOL); 
      return(IDA_ILL_INPUT); 
    }

  }

  /* Set space requirements for one N_Vector */
  if (yy0->ops->nvspace != NULL) {
    N_VSpace(yy0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  IDA_mem->ida_lrw1 = lrw1;
  IDA_mem->ida_liw1 = liw1;

  /* Allocate the vectors (using yy0 as a template) */
  allocOK = IDAAllocVectors(IDA_mem, yy0, itol);
  if (!allocOK) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }
 
  /* All error checking is complete at this point */

  /* Copy the input parameters into IDA memory block */

  IDA_mem->ida_res = res;
  IDA_mem->ida_tn  = t0;

  /* Copy tolerances into memory */

  IDA_mem->ida_itol = itol;
  IDA_mem->ida_rtol = rtol;      

  if (itol == IDA_SS)
    IDA_mem->ida_Satol = *((realtype *)atol);
  else if (itol == IDA_SV) 
    N_VScale(ONE, (N_Vector)atol, IDA_mem->ida_Vatol);

  /* Set the linear solver addresses to NULL */
  IDA_mem->ida_linit  = NULL;
  IDA_mem->ida_lsetup = NULL;
  IDA_mem->ida_lsolve = NULL;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = NULL;
  IDA_mem->ida_lmem   = NULL;

  /* Initialize the phi array */
  N_VScale(ONE, yy0, IDA_mem->ida_phi[0]);  
  N_VScale(ONE, yp0, IDA_mem->ida_phi[1]);  
 
  /* Initialize all the counters and other optional output values */
  IDA_mem->ida_nst     = 0;
  IDA_mem->ida_nre     = 0;
  IDA_mem->ida_ncfn    = 0;
  IDA_mem->ida_netf    = 0;
  IDA_mem->ida_nni     = 0;
  IDA_mem->ida_nsetups = 0;
  
  IDA_mem->ida_kused = 0;
  IDA_mem->ida_hused = ZERO;
  IDA_mem->ida_tolsf = ONE;

  /* Initial setup not done yet */
  IDA_mem->ida_SetupDone = FALSE;

  /* Problem memory has been successfully allocated */
  IDA_mem->ida_MallocDone = TRUE;
  return(IDA_SUCCESS);
}

#define lrw1 (IDA_mem->ida_lrw1)
#define liw1 (IDA_mem->ida_liw1)

/*
 * -----------------------------------------------------------------
 * IDAReInit
 * -----------------------------------------------------------------
 * IDAReInit re-initializes IDA's memory for a problem, assuming
 * it has already beeen allocated in a prior IDAMalloc call.
 * All problem specification inputs are checked for errors.
 * The problem size Neq is assumed to be unchaged since the call
 * to IDAMalloc, and the maximum order maxord must not be larger.
 * If any error occurs during reinitialization, it is reported to
 * the file whose file pointer is errfp.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 * -----------------------------------------------------------------
 */

int IDAReInit(void *ida_mem, IDAResFn res,
              realtype t0, N_Vector yy0, N_Vector yp0,
              int itol, realtype rtol, void *atol)
{
  IDAMem IDA_mem;
  booleantype neg_atol;

  /* Check for legal input parameters */
  
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAM_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if problem was malloc'ed */
  
  if (IDA_mem->ida_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSG_REI_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check for legal input parameters */
  
  if (yy0 == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_Y0_NULL); 
    return(IDA_ILL_INPUT); 
  }
  
  if (yp0 == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_YP0_NULL); 
    return(IDA_ILL_INPUT); 
  }

  if ((itol != IDA_SS) && (itol != IDA_SV) && (itol != IDA_WF)) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_ITOL);
    return(IDA_ILL_INPUT);
  }

  if (res == NULL) { 
    if(errfp!=NULL) fprintf(errfp, MSG_RES_NULL); 
    return(IDA_ILL_INPUT); 
  }

  /* Test tolerances */

  if (itol != IDA_WF) {

    if (atol == NULL) { 
      if(errfp!=NULL) fprintf(errfp, MSG_ATOL_NULL); 
      return(IDA_ILL_INPUT); 
    }

    if (rtol < ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_RTOL); 
      return(IDA_ILL_INPUT); 
    }
   
    if (itol == IDA_SS) { 
      neg_atol = (*((realtype *)atol) < ZERO); 
    } else { 
      neg_atol = (N_VMin((N_Vector)atol) < ZERO); 
    }
    if (neg_atol) { 
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_ATOL); 
      return(IDA_ILL_INPUT); 
    }

  }

  /* Copy the input parameters into IDA memory block */
  IDA_mem->ida_res = res;
  IDA_mem->ida_tn  = t0;

  if ( (itol != IDA_SV) && (IDA_mem->ida_VatolMallocDone) ) {
    N_VDestroy(IDA_mem->ida_Vatol);
    lrw -= lrw1;
    liw -= liw1;
    IDA_mem->ida_VatolMallocDone = FALSE;
  }

  if ( (itol == IDA_SV) && !(IDA_mem->ida_VatolMallocDone) ) {
    IDA_mem->ida_Vatol = N_VClone(yy0);
    lrw += lrw1;
    liw += liw1;
    IDA_mem->ida_VatolMallocDone = TRUE;
  }

  IDA_mem->ida_itol = itol;
  IDA_mem->ida_rtol = rtol;      
  if (itol == IDA_SS)
    IDA_mem->ida_Satol = *((realtype *)atol);
  else if (itol == IDA_SV)
    N_VScale(ONE, (N_Vector)atol, IDA_mem->ida_Vatol);

  /* Initialize the phi array */
  N_VScale(ONE, yy0, IDA_mem->ida_phi[0]);  
  N_VScale(ONE, yp0, IDA_mem->ida_phi[1]);  
 
  /* Initialize all the counters and other optional output values */
 
  IDA_mem->ida_nst     = 0;
  IDA_mem->ida_nre     = 0;
  IDA_mem->ida_ncfn    = 0;
  IDA_mem->ida_netf    = 0;
  IDA_mem->ida_nni     = 0;
  IDA_mem->ida_nsetups = 0;
  
  IDA_mem->ida_kused = 0;
  IDA_mem->ida_hused = ZERO;
  IDA_mem->ida_tolsf = ONE;

  /* Initial setup not done yet */
  IDA_mem->ida_SetupDone = FALSE;
      
  /* Problem has been successfully re-initialized */

  return(IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * readability constants
 * -----------------------------------------------------------------
 */

#define res         (IDA_mem->ida_res)
#define y0          (IDA_mem->ida_y0)
#define yp0         (IDA_mem->ida_yp0)

#define itol        (IDA_mem->ida_itol)
#define rtol        (IDA_mem->ida_rtol)
#define Satol       (IDA_mem->ida_Satol)
#define Vatol       (IDA_mem->ida_Vatol)
#define efun        (IDA_mem->ida_efun)
#define edata       (IDA_mem->ida_edata)

#define rdata       (IDA_mem->ida_rdata)
#define maxord      (IDA_mem->ida_maxord)
#define mxstep      (IDA_mem->ida_mxstep)
#define hin         (IDA_mem->ida_hin)
#define hmax_inv    (IDA_mem->ida_hmax_inv)
#define tstop       (IDA_mem->ida_tstop)
#define tstopset    (IDA_mem->ida_tstopset)
#define epcon       (IDA_mem->ida_epcon)
#define maxnef      (IDA_mem->ida_maxnef)
#define maxncf      (IDA_mem->ida_maxncf)
#define maxcor      (IDA_mem->ida_maxcor)
#define suppressalg (IDA_mem->ida_suppressalg)
#define id          (IDA_mem->ida_id)
#define constraints (IDA_mem->ida_constraints)

#define epiccon     (IDA_mem->ida_epiccon)
#define maxnh       (IDA_mem->ida_maxnh)
#define maxnj       (IDA_mem->ida_maxnj)
#define maxnit      (IDA_mem->ida_maxnit)
#define lsoff       (IDA_mem->ida_lsoff)
#define steptol     (IDA_mem->ida_steptol)

#define uround      (IDA_mem->ida_uround)  
#define phi         (IDA_mem->ida_phi) 
#define ewt         (IDA_mem->ida_ewt)  
#define yy          (IDA_mem->ida_yy)
#define yp          (IDA_mem->ida_yp)
#define delta       (IDA_mem->ida_delta)
#define mm          (IDA_mem->ida_mm)
#define ee          (IDA_mem->ida_ee)
#define savres      (IDA_mem->ida_savres)
#define tempv1      (IDA_mem->ida_tempv1)
#define tempv2      (IDA_mem->ida_tempv2) 
#define kk          (IDA_mem->ida_kk)
#define hh          (IDA_mem->ida_hh)
#define h0u         (IDA_mem->ida_h0u)
#define tn          (IDA_mem->ida_tn)
#define tretp       (IDA_mem->ida_tretp)
#define cj          (IDA_mem->ida_cj)
#define cjold       (IDA_mem->ida_cjold)
#define cjratio     (IDA_mem->ida_cjratio)
#define cjlast      (IDA_mem->ida_cjlast)
#define nbacktr     (IDA_mem->ida_nbacktr)
#define nst         (IDA_mem->ida_nst)
#define nre         (IDA_mem->ida_nre)
#define ncfn        (IDA_mem->ida_ncfn)
#define netf        (IDA_mem->ida_netf)
#define nni         (IDA_mem->ida_nni)
#define nsetups     (IDA_mem->ida_nsetups)
#define ns          (IDA_mem->ida_ns)
#define linit       (IDA_mem->ida_linit)
#define lsetup      (IDA_mem->ida_lsetup)
#define lsolve      (IDA_mem->ida_lsolve) 
#define lperf       (IDA_mem->ida_lperf)
#define lfree       (IDA_mem->ida_lfree) 
#define lmem        (IDA_mem->ida_lmem) 
#define knew        (IDA_mem->ida_knew)
#define kused       (IDA_mem->ida_kused)          
#define hused       (IDA_mem->ida_hused)         
#define tolsf       (IDA_mem->ida_tolsf)      
#define phase       (IDA_mem->ida_phase)
#define epsNewt     (IDA_mem->ida_epsNewt)
#define toldel      (IDA_mem->ida_toldel)
#define ss          (IDA_mem->ida_ss)
#define rr          (IDA_mem->ida_rr)
#define psi         (IDA_mem->ida_psi)
#define alpha       (IDA_mem->ida_alpha)
#define beta        (IDA_mem->ida_beta)
#define sigma       (IDA_mem->ida_sigma)
#define gamma       (IDA_mem->ida_gamma)
#define setupNonNull (IDA_mem->ida_setupNonNull) 
#define constraintsSet (IDA_mem->ida_constraintsSet)

/*
 * -----------------------------------------------------------------
 * IDASolve
 * -----------------------------------------------------------------
 * This routine is the main driver of the IDA package. 
 *
 * It integrates over an independent variable interval defined by the user, 
 * by calling IDAStep to take internal independent variable steps.
 *
 * The first time that IDASolve is called for a successfully initialized
 * problem, it computes a tentative initial step size.
 *
 * IDASolve supports four modes, specified by itask:
 * IDA_NORMAL,  IDA_ONE_STEP,  IDA_NORMAL_TSTOP,  and  IDA_ONE_STEP_TSTOP.
 * In the IDA_NORMAL and IDA_NORMAL_TSTOP modes, the solver steps until it 
 * passes tout and then interpolates to obtain y(tout) and yp(tout).
 * In the IDA_ONE_STEP and IDA_ONE_STEP_TSTOP modes, it takes one internal step
 * and returns.  In the IDA_NORMAL_TSTOP and IDA_ONE_STEP_TSTOP modes, it also
 * takes steps so as to reach tstop exactly and never to go past it.
 *
 * IDASolve returns integer values corresponding to success and failure as below:
 *
 * successful returns: 
 *
 * IDA_SUCCESS        
 * IDA_TSTOP_RETURN   
 *
 * failed returns:
 *
 * IDA_ILL_INPUT
 * IDA_TOO_MUCH_WORK
 * IDA_MEM_NULL
 * IDA_TOO_MUCH_ACC
 * IDA_CONV_FAIL
 * IDA_LSETUP_FAIL
 * IDA_LSOLVE_FAIL    
 * IDA_CONSTR_FAIL
 * IDA_ERR_FAIL   
 * IDA_REP_RES_ERR
 * IDA_RES_FAIL
 * -----------------------------------------------------------------
 */

int IDASolve(void *ida_mem, realtype tout, realtype *tret,
             N_Vector yret, N_Vector ypret, int itask)
{
  long int nstloc;
  int sflag, istate, ier;
  realtype tdist, troundoff, ypnorm, rh, nrm;
  booleantype istop;
  int ewtsetOK;
  IDAMem IDA_mem;

  /* Check for legal inputs in all cases. */

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDA_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if problem was malloc'ed */
  
  if (IDA_mem->ida_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check for legal arguments */

  if (yret == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_YRET_NULL);       
    return(IDA_ILL_INPUT);
  }
  yy = yret;  

  if (ypret == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_YPRET_NULL);       
    return(IDA_ILL_INPUT);
  }
  yp = ypret;
  
  if (tret == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_TRET_NULL);
    return(IDA_ILL_INPUT);
  }
  *tret = tretp = tn; /* Set tret now in case of illegal-input return. */

  if ((itask < IDA_NORMAL) || (itask > IDA_ONE_STEP_TSTOP)) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_ITASK);
    return(IDA_ILL_INPUT);
  }
  
  if ( (itask == IDA_NORMAL_TSTOP) || (itask == IDA_ONE_STEP_TSTOP) ) {
    if ( tstopset == FALSE ) {
      if(errfp!=NULL) fprintf(errfp, MSG_NO_TSTOP);
      return(IDA_ILL_INPUT);
    }
    istop = TRUE;
  } else {
      istop = FALSE;
  }


  if (nst == 0) {       /* THIS IS THE FIRST CALL */

    /* Check inputs to the IDA for correctness and consistency */

    if (IDA_mem->ida_SetupDone == FALSE) {
      ier = IDAInitialSetup(IDA_mem);
      if(ier != IDA_SUCCESS) return(IDA_ILL_INPUT);
      IDA_mem->ida_SetupDone = TRUE;
    }

    /* On the first call, check for tout - tn too small,
       set initial hh,
       check for approach to tstop, and scale phi[1] by hh. */

    tdist = ABS(tout - tn);
    troundoff = TWO*uround*(ABS(tn) + ABS(tout));    
    if (tdist < troundoff) {
      if(errfp!=NULL) fprintf(errfp, MSG_TOO_CLOSE);
      return(IDA_ILL_INPUT);
    }

    hh = hin;
    if ( (hh != ZERO) && ((tout-tn)*hh < ZERO) ) {
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_HINIT);
      return(IDA_ILL_INPUT);
    }

    if (hh == ZERO) {
      hh = PT001*tdist;
      ypnorm = IDAWrmsNorm(IDA_mem, phi[1], ewt, suppressalg);
      if (ypnorm > HALF/hh) hh = HALF/ypnorm;
      if(tout < tn) hh = -hh;
    }

    rh = ABS(hh)*hmax_inv;
    if (rh > ONE) hh /= rh;

    if(istop) {
      if ( (tstop - tn)*hh < ZERO) {
        if(errfp!=NULL) fprintf(errfp, MSG_BAD_TSTOP, tn);
        return(IDA_ILL_INPUT);
      }
      if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
    }

    h0u = hh;

    N_VScale(hh, phi[1], phi[1]);
    kk = 0; kused = 0;  /* set in case of an error return before a step */

    /* Set the convergence test constants epsNewt and toldel */

    epsNewt = epcon;
    toldel = PT0001 * epsNewt;

  } /* end of first-call block. */

  /* Call lperf function and set nstloc for later performance testing. */

  if (lperf != NULL) lperf(IDA_mem, 0);
  nstloc = 0;

  /* If not the first call, check for stop conditions. */

  if (nst > 0) {
    istate = IDAStopTest1(IDA_mem, tout, tret, yret, ypret, itask);
    if (istate != CONTINUE_STEPS) return(istate);
  }

  /* Looping point for internal steps. */

  loop {
   
    /* Check for too many steps taken. */
    
    if (nstloc >= mxstep) {
      if(errfp!=NULL) fprintf(errfp, MSG_MAX_STEPS, tn);
      istate = IDA_TOO_MUCH_WORK;
      *tret = tretp = tn;
      break; /* Here yy=yret and yp=ypret already have the current solution. */
    }

    /* Call lperf to generate warnings of poor performance. */

    if (lperf != NULL) lperf(IDA_mem, 1);

    /* Reset and check ewt (if not first call). */

    if (nst > 0) {
      ewtsetOK = efun(phi[0], ewt, edata);
      if (ewtsetOK != 0) {
	if(errfp!=NULL) {
          if (itol == IDA_WF) fprintf(errfp, MSG_EWT_NOW_FAIL, tn);
          else fprintf(errfp, MSG_EWT_NOW_BAD, tn);
	}
        istate = IDA_ILL_INPUT;
        ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
        *tret = tretp = tn;
        break;
      }
    }
    
    /* Check for too much accuracy requested. */
    
    nrm = IDAWrmsNorm(IDA_mem, phi[0], ewt, suppressalg);
    tolsf = uround * nrm;
    if (tolsf > ONE) {
      tolsf *= TEN;
      if(errfp!=NULL) fprintf(errfp, MSG_TOO_MUCH_ACC, tn);
      istate = IDA_TOO_MUCH_ACC;
      *tret = tretp = tn;
      if (nst > 0) ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      break;
    }

    /* Call IDAStep to take a step. */

    sflag = IDAStep(IDA_mem);

    /* Process all failed-step cases, and exit loop. */
   
    if (sflag != IDA_SUCCESS) {
      istate = IDAHandleFailure(IDA_mem, sflag);
      *tret = tretp = tn;
      ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      break;
    }
    
    nstloc++;

    /* After successful step, check for stop conditions; continue or break. */

    istate = IDAStopTest2(IDA_mem, tout, tret, yret, ypret, itask);
    if (istate != CONTINUE_STEPS) break;

  } /* End of step loop */

  return(istate);    
}

/* 
 * -----------------------------------------------------------------
 * IDAGetSolution
 * -----------------------------------------------------------------
 * This routine evaluates y(t) and y'(t) as the value and derivative of 
 * the interpolating polynomial at the independent variable t, and stores
 * the results in the vectors yret and ypret.  It uses the current
 * independent variable value, tn, and the method order last used, kused.
 * This function is called by IDASolve with t = tout, t = tn, or t = tstop.
 * 
 * If kused = 0 (no step has been taken), or if t = tn, then the order used
 * here is taken to be 1, giving yret = phi[0], ypret = phi[1]/psi[0].
 * 
 * The return values are:
 *   IDA_SUCCESS  if t is legal, or
 *   IDA_BAD_T    if t is not within the interval of the last step taken.
 * -----------------------------------------------------------------
 */

int IDAGetSolution(void *ida_mem, realtype t, N_Vector yret, N_Vector ypret)
{
  IDAMem IDA_mem;
  realtype tfuzz, tp, delt, c, d, gam;
  int j, kord;

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem) ida_mem; 
  /* Check t for legality.  Here tn - hused is t_{n-1}. */
 
  tfuzz = HUNDRED * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  if ( (t - tp)*hh < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAG_BAD_T, t, tn-hused, tn);
    return(IDA_BAD_T);
  }

  /* Initialize yret = phi[0], ypret = 0, and kord = (kused or 1). */

  N_VScale (ONE, phi[0], yret);
  N_VConst (ZERO, ypret);
  kord = kused; 
  if (kused == 0) kord = 1;

 /* Accumulate multiples of columns phi[j] into yret and ypret. */

  delt = t - tn;
  c = ONE; d = ZERO;
  gam = delt/psi[0];
  for (j=1; j <= kord; j++) {
    d = d*gam + c/psi[j-1];
    c = c*gam;
    gam = (delt + psi[j-1])/psi[j];
    N_VLinearSum(ONE,  yret, c, phi[j],  yret);
    N_VLinearSum(ONE, ypret, d, phi[j], ypret);
  }
  return(IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDAFree
 * -----------------------------------------------------------------
 * This routine frees the problem memory allocated by IDAMalloc
 * Such memory includes all the vectors allocated by IDAAllocVectors,
 * and the memory lmem for the linear solver (deallocated by a call
 * to lfree).
 * -----------------------------------------------------------------
 */

void IDAFree(void *ida_mem)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL) return;

  IDA_mem = (IDAMem) ida_mem;
  
  IDAFreeVectors(IDA_mem);
  if (lfree != NULL) lfree(IDA_mem);
  free(IDA_mem);
}

/*
 * -----------------------------------------------------------------
 * private helper functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * IDACheckNvector
 * -----------------------------------------------------------------
 * This routine checks if all required vector operations are present.
 * If any of them is missing it returns FALSE.
 * -----------------------------------------------------------------
 */

static booleantype IDACheckNvector(N_Vector tmpl)
{
  if((tmpl->ops->nvclone        == NULL) ||
     (tmpl->ops->nvdestroy      == NULL) ||
     (tmpl->ops->nvlinearsum    == NULL) ||
     (tmpl->ops->nvconst        == NULL) ||
     (tmpl->ops->nvprod         == NULL) ||
     (tmpl->ops->nvscale        == NULL) ||
     (tmpl->ops->nvabs          == NULL) ||
     (tmpl->ops->nvinv          == NULL) ||
     (tmpl->ops->nvaddconst     == NULL) ||
     (tmpl->ops->nvwrmsnorm     == NULL) ||
     (tmpl->ops->nvmin          == NULL))
    return(FALSE);
  else
    return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * IDAAllocVectors
 * -----------------------------------------------------------------
 * This routine allocates the IDA vectors ewt, tempv1, tempv2, and
 * phi[0], ..., phi[maxord]. If tol=IDA_SV, it also allocates space 
 * for Vatol.
 * If all memory allocations are successful, IDAAllocVectors returns 
 * TRUE. Otherwise all allocated memory is freed and IDAAllocVectors 
 * returns FALSE.
 * This routine also sets the optional outputs lrw and liw, which are
 * (respectively) the lengths of the real and integer work spaces
 * allocated here.
 * -----------------------------------------------------------------
 */

static booleantype IDAAllocVectors(IDAMem IDA_mem, N_Vector tmpl, int tol)
{
  int i, j, maxcol;

  /* Allocate ewt, ee, delta, tempv1, tempv2 */
  
  ewt = N_VClone(tmpl);
  if (ewt == NULL) return(FALSE);

  ee = N_VClone(tmpl);
  if (ee == NULL) {
    N_VDestroy(ewt);
    return(FALSE);
  }
  delta = N_VClone(tmpl);
  if (delta == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(ee);
    return(FALSE);
  }
  tempv1 = N_VClone(tmpl);
  if (tempv1 == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(ee);
    N_VDestroy(delta);
    return(FALSE);
  }
  tempv2= N_VClone(tmpl);
  if (tempv2 == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(ee);
    N_VDestroy(delta);
    N_VDestroy(tempv1);
    return(FALSE);
  }

  savres = tempv1;

  /* Allocate phi[0] ... phi[maxord].  Make sure phi[2] and phi[3] are
  allocated (for use as temporary vectors), regardless of maxord.       */

  maxcol = MAX(maxord,3);
  for (j=0; j <= maxcol; j++) {
    phi[j] = N_VClone(tmpl);
    if (phi[j] == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(ee);
      N_VDestroy(delta);
      N_VDestroy(tempv1);
      N_VDestroy(tempv2);
      for (i=0; i < j; i++) N_VDestroy(phi[i]);
      return(FALSE);
    }
  }

  /* Update solver workspace lengths  */
  lrw += (maxcol + 6)*lrw1;
  liw += (maxcol + 6)*liw1;

  if (tol == IDA_SV) {
    Vatol = N_VClone(tmpl);
    if (Vatol == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(ee);
      N_VDestroy(delta);
      N_VDestroy(tempv1);
      N_VDestroy(tempv2);
      for (i=0; i <= maxcol; i++) N_VDestroy(phi[i]);
    }
    lrw += lrw1;
    liw += liw1;
    IDA_mem->ida_VatolMallocDone = TRUE;
  }

  return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * IDAfreeVectors
 * -----------------------------------------------------------------
 * This routine frees the IDA vectors allocated for IDA.
 * -----------------------------------------------------------------
 */

static void IDAFreeVectors(IDAMem IDA_mem)
{
  int j, maxcol;
  
  N_VDestroy(ewt);
  N_VDestroy(ee);
  N_VDestroy(delta);
  N_VDestroy(tempv1);
  N_VDestroy(tempv2);
  maxcol = MAX(maxord,3);
  for(j=0; j <= maxcol; j++) N_VDestroy(phi[j]);

  lrw -= (maxcol + 6)*lrw1;
  liw -= (maxcol + 6)*liw1;

  if (IDA_mem->ida_VatolMallocDone) {
    N_VDestroy(Vatol);
    lrw -= lrw1;
    liw -= liw1;
  }

  if (IDA_mem->ida_constraintsMallocDone) {
    N_VDestroy(constraints);
    lrw -= lrw1;
    liw -= liw1;
  }

  if (IDA_mem->ida_idMallocDone) {
    N_VDestroy(id);
    lrw -= lrw1;
    liw -= liw1;
  }

}

/*
 * -----------------------------------------------------------------
 * IDAInitialSetup
 * -----------------------------------------------------------------
 * This routine is called by IDASolve once at the first step. It performs
 * all checks on optional inputs and inputs to IDAMalloc/IDAReInit that
 * could not be done before.
 *
 * If no merror is encountered, IDAInitialSetup returns IDA_SUCCESS. Otherwise,
 * it returns an error flag and prints a message to errfp.
 * -----------------------------------------------------------------
 */

int IDAInitialSetup(IDAMem IDA_mem)
{
  booleantype conOK;
  int ewtsetOK;
  int ier;
  
  /* Test for more vector operations, depending on options */

  if (suppressalg)
    if (id->ops->nvwrmsnormmask == NULL) {
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_NVECTOR);
      return(IDA_ILL_INPUT);
  }

  /* Test id vector for legality */
  
  if(suppressalg && (id==NULL)){ 
    if(errfp!=NULL) fprintf(errfp, MSG_MISSING_ID); 
    return(IDA_ILL_INPUT); 
  }

  /* Load ewt */

  if (itol != IDA_WF) {
    efun = IDAEwtSet;
    edata = (void *)IDA_mem;
  } else {
    if (efun == NULL) {
      if (errfp != NULL) fprintf(errfp, MSG_NO_EFUN);
      return(IDA_ILL_INPUT);
    }
  }

  ewtsetOK = efun(phi[0], ewt, edata);
  if (ewtsetOK != 0) {
    if(errfp!=NULL) {
      if (itol == IDA_WF) fprintf(errfp, MSG_FAIL_EWT);
      else fprintf(errfp, MSG_BAD_EWT);
    }
    return(IDA_ILL_INPUT);
  }

  /* Check to see if y0 satisfies constraints. */

  if (constraintsSet) {
    conOK = N_VConstrMask (constraints, phi[0], tempv2);
    if (!conOK) { 
      if(errfp!=NULL) fprintf(errfp, MSG_Y0_FAIL_CONSTR); 
      return(IDA_ILL_INPUT); 
    }
  }

  /* Check that lsolve exists and call linit function if it exists. */

  if (lsolve == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_LSOLVE_NULL);
    return(IDA_ILL_INPUT);
  }

  if (linit != NULL) {
    ier = linit(IDA_mem);
    if (ier < 0) {
      if(errfp!=NULL) fprintf(errfp, MSG_LINIT_FAIL);
      return(IDA_LINIT_FAIL);
    }
  }

  return(IDA_SUCCESS);

}

/*
 * -----------------------------------------------------------------
 * IDAStopTest1
 * -----------------------------------------------------------------
 * This routine tests for stop conditions before taking a step.
 * The tests depend on the value of itask.
 * The variable tretp is the previously returned value of tret.
 *
 * The return values are:
 * CONTINUE_STEPS       if no stop conditions were found
 * IDA_SUCCESS          for a normal return to the user
 * IDA_TSTOP_RETURN     for a tstop-reached return to the user
 * IDA_ILL_INPUT        for an illegal-input return to the user 
 *
 * In the tstop cases, this routine may adjust the stepsize hh to cause
 * the next step to reach tstop exactly.
 * -----------------------------------------------------------------
 */

static int IDAStopTest1(IDAMem IDA_mem, realtype tout, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask)
{

  int ier;
  realtype troundoff;

  switch (itask) {
    
  case IDA_NORMAL:  
    /* Test for tout = tretp, and for tn past tout. */
    if (tout == tretp) {
      *tret = tretp = tout;
      return(IDA_SUCCESS);
    }
    if ( (tn - tout)*hh >= ZERO) {
      ier = IDAGetSolution(IDA_mem, tout, yret, ypret);
      if (ier != IDA_SUCCESS) {
        if(errfp!=NULL) fprintf(errfp, MSG_BAD_TOUT, tout);
        return(IDA_ILL_INPUT);
      }
      *tret = tretp = tout;
      return(IDA_SUCCESS);
    }
    return(CONTINUE_STEPS);
    
  case IDA_ONE_STEP:
    /* Test for tn past tretp. */
    if ( (tn - tretp)*hh > ZERO) {
      ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      *tret = tretp = tn;
      return(IDA_SUCCESS);
    }
    return(CONTINUE_STEPS);
    
  case IDA_NORMAL_TSTOP:
    /* Test for tn past tstop, tn = tretp, tn past tout, tn near tstop. */
    if ( (tn - tstop)*hh > ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_TSTOP, tn);
      return(IDA_ILL_INPUT);
    }
    if (tout == tretp) {
      *tret = tretp = tout;
      return(IDA_SUCCESS);
    }
    if ( (tn - tout)*hh >= ZERO) {
      ier = IDAGetSolution(IDA_mem, tout, yret, ypret);
      if (ier != IDA_SUCCESS) {
        if(errfp!=NULL) fprintf(errfp, MSG_BAD_TOUT, tout);
        return(IDA_ILL_INPUT);
      }
      *tret = tretp = tout;
      return(IDA_SUCCESS);
    }
    troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
    if ( ABS(tn - tstop) <= troundoff) {
      ier = IDAGetSolution(IDA_mem, tstop, yret, ypret);
      if (ier != IDA_SUCCESS) {
	  if(errfp!=NULL) fprintf(errfp, MSG_BAD_TSTOP, tn);
        return(IDA_ILL_INPUT);
      }
      *tret = tretp = tstop;
      return(IDA_TSTOP_RETURN);
    }
    if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
    return(CONTINUE_STEPS);
    
  case IDA_ONE_STEP_TSTOP:
    /* Test for tn past tstop, tn past tretp, and tn near tstop. */
    if ( (tn - tstop)*hh > ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSG_BAD_TSTOP, tn);
      return(IDA_ILL_INPUT);
    }
    if ( (tn - tretp)*hh > ZERO) {
      ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      *tret = tretp = tn;
      return(IDA_SUCCESS);
    }
    troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
    if ( ABS(tn - tstop) <= troundoff) {
      ier = IDAGetSolution(IDA_mem, tstop, yret, ypret);
      if (ier != IDA_SUCCESS) {
        if(errfp!=NULL) fprintf(errfp, MSG_BAD_TSTOP, tn);
        return(IDA_ILL_INPUT);
      }
      *tret = tretp = tstop;
      return(IDA_TSTOP_RETURN);
    }
    if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
    return(CONTINUE_STEPS);
    
  }
  return(-99);
}

/*
 * -----------------------------------------------------------------
 * IDAStopTest2
 * -----------------------------------------------------------------
 * This routine tests for stop conditions after taking a step.
 * The tests depend on the value of itask.
 *
 * The return values are:
 *  CONTINUE_STEPS     if no stop conditions were found
 *  IDA_SUCCESS        for a normal return to the user
 *  IDA_TSTOP_RETURN   for a tstop-reached return to the user
 *
 * In the two cases with tstop, this routine may reset the stepsize hh
 * to cause the next step to reach tstop exactly.
 *
 * In the two cases with ONE_STEP mode, no interpolation to tn is needed
 * because yret and ypret already contain the current y and y' values.
 *
 * Note: No test is made for an error return from IDAGetSolution here,
 * because the same test was made prior to the step.
 * -----------------------------------------------------------------
 */

static int IDAStopTest2(IDAMem IDA_mem, realtype tout, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask)
{

  int ier;
  realtype troundoff;

  switch (itask) {

    case IDA_NORMAL:  
      /* Test for tn past tout. */
      if ( (tn - tout)*hh >= ZERO) {
        ier = IDAGetSolution(IDA_mem, tout, yret, ypret);
        *tret = tretp = tout;
        return(IDA_SUCCESS);
      }
      return(CONTINUE_STEPS);

    case IDA_ONE_STEP:
      *tret = tretp = tn;
      return(IDA_SUCCESS);

    case IDA_NORMAL_TSTOP:
      /* Test for tn at tstop, for tn past tout, and for tn near tstop. */
      troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
      if ( ABS(tn - tstop) <= troundoff) {
        ier = IDAGetSolution(IDA_mem, tstop, yret, ypret);
        *tret = tretp = tstop;
        return(IDA_TSTOP_RETURN);
      }
      if ( (tn - tout)*hh >= ZERO) {
        ier = IDAGetSolution(IDA_mem, tout, yret, ypret);
        *tret = tretp = tout;
        return(IDA_SUCCESS);
      }
      if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
      return(CONTINUE_STEPS);

    case IDA_ONE_STEP_TSTOP:
      /* Test for tn at tstop. */
      troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
      if ( ABS(tn - tstop) <= troundoff) {
        ier = IDAGetSolution(IDA_mem, tstop, yret, ypret);
        *tret = tretp = tstop;
        return(IDA_TSTOP_RETURN);
      }
      if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
      *tret = tretp = tn;
      return(IDA_SUCCESS);

  }
  return -99;
}

/*
 * -----------------------------------------------------------------
 * IDAHandleFailure
 * -----------------------------------------------------------------
 * This routine prints error messages for all cases of failure by
 * IDAStep.  It returns to IDASolve the value that it is to return to
 * the user.
 * -----------------------------------------------------------------
 */

static int IDAHandleFailure(IDAMem IDA_mem, int sflag)
{

  /* Depending on sflag, print error message and return error flag */
  switch (sflag) {

    case IDA_ERR_FAIL:  
      if(errfp!=NULL) fprintf(errfp, MSG_ERR_FAILS, tn, hh);
      return(IDA_ERR_FAIL);

    case IDA_CONV_FAIL:
      if(errfp!=NULL) fprintf(errfp, MSG_CONV_FAILS, tn, hh);
      return(IDA_CONV_FAIL);

    case IDA_LSETUP_FAIL:  
      if(errfp!=NULL) fprintf(errfp, MSG_SETUP_FAILED, tn);
      return(IDA_LSETUP_FAIL);

    case IDA_LSOLVE_FAIL: 
      if(errfp!=NULL) fprintf(errfp, MSG_SOLVE_FAILED, tn);
      return(IDA_LSOLVE_FAIL);

    case IDA_REP_RES_ERR:
      if(errfp!=NULL) fprintf(errfp, MSG_REP_RES_ERR, tn);
      return(IDA_REP_RES_ERR);

    case IDA_RES_FAIL: 
      if(errfp!=NULL) fprintf(errfp, MSG_RES_NONRECOV, tn);
      return(IDA_RES_FAIL);

    case IDA_CONSTR_FAIL: 
      if(errfp!=NULL) fprintf(errfp, MSG_FAILED_CONSTR, tn);
      return(IDA_CONSTR_FAIL);

  }

  return -99;

}

/*
 * -----------------------------------------------------------------
 * IDAStep
 * -----------------------------------------------------------------
 * This routine performs one internal IDA step, from tn to tn + hh.
 * It calls other routines to do all the work.
 *
 * It solves a system of differential/algebraic equations of the form
 *       F(t,y,y') = 0, for one step. In IDA, tt is used for t,
 * yy is used for y, and yp is used for y'. The function F is supplied as 'res'
 * by the user.
 *
 * The methods used are modified divided difference, fixed leading 
 * coefficient forms of backward differentiation formulas.
 * The code adjusts the stepsize and order to control the local error per step.
 *
 * The main operations done here are as follows:
 *  * initialize various quantities;
 *  * setting of multistep method coefficients;
 *  * solution of the nonlinear system for yy at t = tn + hh;
 *  * deciding on order reduction and testing the local error;
 *  * attempting to recover from failure in nonlinear solver or error test;
 *  * resetting stepsize and order for the next step.
 *  * updating phi and other state data if successful;
 *
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 *
 * Variables or arrays (all in the IDAMem structure) used in IDAStep are:
 *
 * tt -- Independent variable.
 * yy -- Solution vector at tt.
 * yp -- Derivative of solution vector after successful stelp.
 * res -- User-supplied function to evaluate the residual. See the 
 *        description given in file ida.h .
 * lsetup -- Routine to prepare for the linear solver call. It may either
 *        save or recalculate quantities used by lsolve. (Optional)
 * lsolve -- Routine to solve a linear system. A prior call to lsetup
 *        may be required. 
 * hh  -- Appropriate step size for next step.
 * ewt -- Vector of weights used in all convergence tests.
 * phi -- Array of divided differences used by IDAStep. This array is composed 
 *       of  (maxord+1) nvectors (each of size Neq). (maxord+1) is the maximum 
 *       order for the problem, maxord, plus 1.
 *
 *       Return values are:
 *       IDA_SUCCESS   IDA_RES_FAIL        LSETUP_ERROR_NONRECVR       
 *                     IDA_LSOLVE_FAIL   IDA_ERR_FAIL            
 *                     IDA_CONSTR_FAIL               IDA_CONV_FAIL          
 *                     IDA_REP_RES_ERR            
 * -----------------------------------------------------------------
 */

static int IDAStep(IDAMem IDA_mem)
{
  realtype saved_t, ck, est;
  realtype terk, terkm1, erkm1;
  int ncf, nef, nflag, kflag;
  
  saved_t = tn;
  ncf = nef = 0;

  if(nst == ZERO){
    kk = 1;
    kused = 0;
    hused = ZERO;
    psi[0] = hh;
    cj = ONE/hh;
    phase = 0;
    ns = 0;
  }
  
  /* Looping point for attempts to take a step */

  loop {  
    IDASetCoeffs(IDA_mem, &ck);
    kflag = IDA_SUCCESS;

    nflag = IDAnls(IDA_mem);

    if(nflag == IDA_SUCCESS) 
      nflag = IDATestError(IDA_mem, &ck, &est, &terk, &terkm1, &erkm1);

    if(nflag != IDA_SUCCESS) 
      kflag = IDAHandleNFlag(IDA_mem, nflag, saved_t, &ncf, &nef, &est);

    if (kflag == PREDICT_AGAIN) continue;
    else if(kflag == IDA_SUCCESS) break;
    else return(kflag);
  }

  /* Nonlinear system solve and error test were both successful;
     update data, and consider change of step and/or order       */

  IDACompleteStep(IDA_mem, &est, &terk, &terkm1, &erkm1);

  return(IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDASetCoeffs
 * -----------------------------------------------------------------
 *  This routine computes the coefficients relevant to the current step.
 *  The counter ns counts the number of consecutive steps taken at
 *  constant stepsize h and order k, up to a maximum of k + 2.
 *  Then the first ns components of beta will be one, and on a step  
 *  with ns = k + 2, the coefficients alpha, etc. need not be reset here.
 *  Also, IDACompleteStep prohibits an order increase until ns = k + 2.
 * -----------------------------------------------------------------
 */

static void IDASetCoeffs(IDAMem IDA_mem, realtype *ck)
{
  int i;
  realtype temp1, temp2, alpha0, alphas;

  /* Set coefficients for the current stepsize h */

  if(hh != hused || kk != kused) ns = 0;
  ns = MIN(ns+1,kused+2);
  if(kk+1 >= ns){
    beta[0] = ONE;
    alpha[0] = ONE;
    temp1 = hh;
    gamma[0] = ZERO;
    sigma[0] = ONE;
    for(i=1;i<=kk;i++){
      temp2 = psi[i-1];
      psi[i-1] = temp1;
      beta[i] = beta[i-1] * psi[i-1] / temp2;
      temp1 = temp2 + hh;
      alpha[i] = hh / temp1;
      sigma[i] = i * sigma[i-1] * alpha[i]; 
      gamma[i] = gamma[i-1] + alpha[i-1] / hh;
   }
    psi[kk] = temp1;
  }
  /* compute alphas, alpha0 */
  alphas = ZERO;
  alpha0 = ZERO;
  for(i=0;i<kk;i++){
    alphas = alphas - ONE/(i+1);
    alpha0 = alpha0 - alpha[i];
  }

  /* compute leading coefficient cj  */
  cjlast = cj;
  cj = -alphas/hh;
  
  /* compute variable stepsize error coefficient ck */

    *ck = ABS(alpha[kk] + alphas - alpha0);
    *ck = MAX(*ck, alpha[kk]);

 /* change phi to phi-star  */

  for(i=ns;i<=kk;i++) N_VScale(beta[i], phi[i], phi[i]);

  /* update independent variable */

  tn = tn + hh;
}


/*
 * -----------------------------------------------------------------
 * IDAnls
 * -----------------------------------------------------------------
 * This routine attempts to solve the nonlinear system using the linear
 * solver specified. NOTE: this routine uses N_Vector ee as the scratch
 * vector tempv3 passed to lsetup.
 *
 *  Possible return values:
 *
 *  IDA_SUCCESS
 *
 *  IDA_RES_RECVR       IDA_RES_FAIL
 *  IDA_LSETUP_RECVR    IDA_LSETUP_FAIL
 *  IDA_LSOLVE_RECVR    IDA_LSOLVE_FAIL
 *
 *  IDA_CONSTR_RECVR
 *  IDA_NCONV_RECVR
 * -----------------------------------------------------------------
 */

static int IDAnls(IDAMem IDA_mem)
{
  int retval, ier;
  booleantype constraintsPassed, callSetup, tryAgain;
  realtype temp1, temp2, vnorm;
  N_Vector tempv3;

  callSetup = FALSE;

  /* Initialize if the first time called */

  if(nst == 0){
    cjold = cj;
    ss = TWENTY;
    if(setupNonNull) callSetup = TRUE;
  }

  mm = tempv2;
  tempv3 = ee;


  /* Decide if lsetup is to be called */

  if(setupNonNull){
    cjratio = cj / cjold;
    temp1 = (ONE - XRATE) / (ONE + XRATE);
    temp2 = ONE/temp1;
    {if(cjratio < temp1 || cjratio > temp2) callSetup = TRUE;}
    {if(cj != cjlast) ss=HUNDRED;}
  }

  /* Begin the main loop. This loop is traversed at most twice. 
     The second pass only occurs when the first pass had a recoverable
     failure with old Jacobian data */
  loop{


    /* Compute predicted values for yy and yp, and compute residual there. */
    ier = IDAPredict(IDA_mem);

    retval = res(tn, yy, yp, delta, rdata);
    nre++;
    if(retval < 0) return(IDA_RES_FAIL);
    if(retval > 0) return(IDA_RES_RECVR);

    /* If indicated, call linear solver setup function and reset parameters. */
    if(callSetup){
      nsetups++;
      retval = lsetup(IDA_mem, yy, yp, delta, tempv1, tempv2, tempv3);
      cjold = cj;
      cjratio = ONE;
      ss = TWENTY;
      if (retval < 0) return(IDA_LSETUP_FAIL);
      if (retval > 0) return(IDA_LSETUP_RECVR);
    }

    /* Call the Newton iteration routine.  */

    retval = IDANewtonIter(IDA_mem);


    /* Retry the current step on recoverable failure with old Jacobian data. */

    tryAgain = (retval>0)&&(setupNonNull) &&(!callSetup);

    if(tryAgain){
      callSetup = TRUE;
      continue;
    }
    else break;

  }  /* end of loop */

  if(retval != IDA_SUCCESS) return(retval);


  /* If otherwise successful, check and enforce inequality constraints. */

  if(constraintsSet){  /* Check constraints and get mask vector mm, 
                          set where constraints failed */
    constraintsPassed = N_VConstrMask(constraints,yy,mm);
    if(constraintsPassed) return(IDA_SUCCESS);
    else {
      N_VCompare(ONEPT5, constraints, tempv1);  
      /* a , where a[i] =1. when |c[i]| = 2 ,  c the vector of constraints */
      N_VProd(tempv1, constraints, tempv1);       /* a * c */
      N_VDiv(tempv1, ewt, tempv1);                /* a * c * wt */
      N_VLinearSum(ONE, yy, -PT1, tempv1, tempv1);/* y - 0.1 * a * c * wt */
      N_VProd(tempv1, mm, tempv1);               /*  v = mm*(y-.1*a*c*wt) */
      vnorm = IDAWrmsNorm(IDA_mem, tempv1, ewt, FALSE); /*  ||v|| */
      
      /* If vector v of constraint corrections is small
         in norm, correct and accept this step */      
      if(vnorm <= epsNewt){  
        N_VLinearSum(ONE, ee, -ONE, tempv1, ee);  /* ee <- ee - v */
        return(IDA_SUCCESS);
      }
      else {
        /* Constraints not met -- reduce h by computing rr = h'/h */
        N_VLinearSum(ONE, phi[0], -ONE, yy, tempv1);
        N_VProd(mm, tempv1, tempv1);
        rr = PT9*N_VMinQuotient(phi[0], tempv1);
        rr = MAX(rr,PT1);
        return(IDA_CONSTR_RECVR);
      }
    }
  }
  return(IDA_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * IDAPredict
 * -----------------------------------------------------------------
 * This routine predicts the new values for vectors yy and yp.
 * -----------------------------------------------------------------
 */

static int IDAPredict(IDAMem IDA_mem)
{
  int j;

  N_VScale(ONE, phi[0], yy);
  N_VConst(ZERO, yp);
  
  for(j=1; j<=kk; j++) {
    N_VLinearSum(ONE,      phi[j], ONE, yy, yy);
    N_VLinearSum(gamma[j], phi[j], ONE, yp, yp);
  }

  return(IDA_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * IDANewtonIter
 * -----------------------------------------------------------------
 * This routine performs the Newton iteration.  
 * It assumes that delta contains the initial residual vector on entry.
 * If the iteration succeeds, it returns the value IDA_SUCCESS = 0.
 * If not, it returns either:
 *   a positive value (for a recoverable failure), namely one of:
 *     IDA_RES_RECVR
 *     IDA_LSOLVE_RECVR
 *     IDA_NCONV_RECVR
 * or
 *   a negative value (for a nonrecoverable failure), namely one of:
 *     IDA_RES_FAIL
 *     IDA_LSOLVE_FAIL
 *
 *     NOTE: This routine uses N_Vector savres, which is preset to tempv1.
 * -----------------------------------------------------------------
 */

static int IDANewtonIter(IDAMem IDA_mem)
{
  int mnewt, retval;
  realtype delnrm, oldnrm, rate;

  /* Initialize counter mnewt and cumulative correction vector ee. */
  mnewt = 0;
  N_VConst(ZERO, ee);

  /* Initialize oldnrm to avoid compiler warning message */
  oldnrm = ZERO;

  /* Looping point for Newton iteration.  Break out on any error. */
  loop {

    nni++;

    /* Save a copy of the residual vector in savres. */
    N_VScale(ONE, delta, savres);

    /* Call the lsolve function to get correction vector delta. */
    retval = lsolve(IDA_mem, delta, ewt, yy, yp, savres); 
    if(retval < 0) return(IDA_LSOLVE_FAIL);
    if(retval > 0) return(IDA_LSOLVE_RECVR);

    /* Apply delta to yy, yp, and ee, and get norm(delta). */
    N_VLinearSum(ONE, yy, -ONE, delta, yy);
    N_VLinearSum(ONE, ee, -ONE, delta, ee);
    N_VLinearSum(ONE, yp, -cj,  delta, yp);
    delnrm = IDAWrmsNorm(IDA_mem, delta, ewt, FALSE);

    /* Test for convergence, first directly, then with rate estimate. */

    if (mnewt == 0){ 
       oldnrm = delnrm;
       if (delnrm <= toldel) return(IDA_SUCCESS);
    }
    else {
      rate = RPowerR( delnrm/oldnrm, ONE/mnewt );
      if (rate > RATEMAX) return(IDA_NCONV_RECVR); 
      ss = rate/(ONE - rate);
    }

    if (ss*delnrm <= epsNewt) return(IDA_SUCCESS);

    /* Not yet converged.  Increment mnewt and test for max allowed. */
    mnewt++;
    if (mnewt >= maxcor) {retval = IDA_NCONV_RECVR; break;}

    /* Call res for new residual and check error flag from res. */
    retval = res(tn, yy, yp, delta, rdata);
    nre++;
    if(retval < 0) return(IDA_RES_FAIL);
    if(retval > 0) return(IDA_RES_RECVR);

    /* Loop for next iteration. */

  } /* end of Newton iteration loop */

  /* All error returns exit here. */
  return(retval);

}

/*
 * -----------------------------------------------------------------
 * IDATestError
 * -----------------------------------------------------------------
 * This routine estimates errors at orders k, k-1, k-2, decides whether 
 * or not to reduce order, and performs the local error test. 
 *
 * IDATestError returns either  IDA_SUCCESS   or ERROR_TEST_FAIL
 * -----------------------------------------------------------------
 */

static int IDATestError(IDAMem IDA_mem, realtype *ck, realtype *est,
                        realtype *terk, realtype *terkm1, realtype *erkm1)
{
  int retval;
  realtype enorm;
  realtype terkm2;
  realtype erk, erkm2;

  /* Compute error for order k. */
  enorm = IDAWrmsNorm(IDA_mem, ee, ewt, suppressalg);
  erk = sigma[kk] * enorm;
  *terk = (kk+1) * erk;
  *est = erk;
  knew = kk;

  /* Now compute the errors for orders k-1 and k-2, and decide whether to 
     reduce the order k to k-1 */
  
  if(kk > 1){
    N_VLinearSum(ONE, phi[kk], ONE, ee, delta);
    *erkm1 = sigma[kk-1] * IDAWrmsNorm(IDA_mem, delta, ewt, suppressalg);
    *terkm1 = kk * *erkm1;
    {
      if(kk > 2){
        N_VLinearSum(ONE, phi[kk-1], ONE, delta, delta);
        erkm2 = sigma[kk-2] * IDAWrmsNorm(IDA_mem, delta, ewt, suppressalg);
        terkm2 = (kk-1) * erkm2;
        if(MAX(*terkm1, terkm2) > *terk) goto evaltest;
      }
      
      else if(*terkm1 > (HALF * (*terk))) goto evaltest; /* executed for kk=2 only */
    }
    /* end of "kk>2" if/else block */
    
    knew = kk-1;
    *est = *erkm1;
    
  } /* end kk>1 if block */ 
  
  
 evaltest:
  retval = IDA_SUCCESS;
  
  if ((*ck * enorm) > ONE) retval = ERROR_TEST_FAIL;
  return(retval);
}


/*
 * -----------------------------------------------------------------
 * IDAHandleNFlag
 * -----------------------------------------------------------------
 * This routine handles failures indicated by the input variable nflag. 
 * Positive values indicate various recoverable failures while negative
 * values indicate nonrecoverable failures. This routine adjusts the
 * step size for recoverable failures. 
 *
 *  Possible nflag values (input):
 *
 *   --convergence failures--
 *   IDA_RES_RECVR              > 0
 *   IDA_LSOLVE_RECVR           > 0
 *   IDA_CONSTR_RECVR           > 0
 *   IDA_NCONV_RECVR            > 0
 *   IDA_RES_FAIL               < 0
 *   IDA_LSOLVE_FAIL            < 0
 *   IDA_LSETUP_FAIL            < 0
 *
 *   --error test failure--
 *   ERROR_TEST_FAIL            > 0
 *
 *  Possible kflag values (output):
 *
 *   --recoverable--
 *   PREDICT_AGAIN
 *
 *   --nonrecoverable--
 *   IDA_CONSTR_FAIL   
 *   IDA_REP_RES_ERR    
 *   IDA_ERR_FAIL  
 *   IDA_CONV_FAIL 
 *   IDA_RES_FAIL
 *   IDA_LSETUP_FAIL
 *   IDA_LSOLVE_FAIL
 * -----------------------------------------------------------------
 */

static int IDAHandleNFlag(IDAMem IDA_mem, int nflag, realtype saved_t,
                          int *ncfPtr, int *nefPtr, realtype *est)
{
  int j;
  int *ncf, *nef;
  
  ncf = ncfPtr; nef = nefPtr;
  phase = 1;
    
  /* restore tn, phi, and psi */
  tn = saved_t;
  for (j = ns; j <= kk; j++) N_VScale(ONE/beta[j], phi[j], phi[j]);
  for (j = 1; j <= kk; j++) psi[j-1] = psi[j] - hh;
  
  /* NLS FAILURE  */

  if (nflag != ERROR_TEST_FAIL) {

    (*ncf)++; ncfn++;

    /* Nonrecoverable failure */
    if (nflag < 0) return(nflag);

    /* If there were too many convergence failures */
    if (*ncf >= maxncf) {
      if (nflag == IDA_RES_RECVR)    return(IDA_REP_RES_ERR);
      if (nflag == IDA_CONSTR_RECVR) return(IDA_CONSTR_FAIL);
      return(IDA_CONV_FAIL);
    }
    
    /* Prepare to predict again */
    rr = QUARTER;
    hh *= rr;
    
  } 

  /* ERROR TEST FAILURE */

  else { 

    (*nef)++; netf++;
    
    /* If there were too many error test failures */
    if (*nef >= maxnef) return(IDA_ERR_FAIL);
    
    /* Prepare to predict again */
    if (*nef == 1){
      /* On first error test failure, keep current order or lower order 
         by one. Compute new stepsize based on differences of the solution. */
      kk = knew;
      
      rr = PT9 * RPowerR( TWO*(*est) + PT0001,(-ONE/(kk+1)) );
      rr = MAX(QUARTER, MIN(PT9,rr));
      hh *=rr;  /* adjust step size */
    } else if (*nef == 2){
      /* On second error test failure, use current order or decrease order 
         by one. Reduce stepsize by factor of 1/4. */
      kk = knew;
      rr = QUARTER;
      hh *= rr;
      
    } else if (*nef > 2){
      /* On third and subsequent error test failures, set order to 1 and
         reduce stepsize h by factor of 1/4. */
      kk = 1;
      rr = QUARTER;
      hh *= rr;
    }

  } /* end of nflag  if block */
      
  if (nst == 0){
    psi[0] = hh;
    N_VScale(rr, phi[1], phi[1]);
  }

  return(PREDICT_AGAIN);

}

/*
 * -----------------------------------------------------------------
 * IDACompleteStep
 * -----------------------------------------------------------------
 * This routine completes a successful step.  It increments nst,
 * saves the stepsize and order used, makes the final selection of
 * stepsize and order for the next step, and updates the phi array.
 * Its return value is IDA_SUCCESS= 0.
 * -----------------------------------------------------------------
 */

static int IDACompleteStep(IDAMem IDA_mem, realtype *est, 
                           realtype *terk, realtype *terkm1, realtype *erkm1)
{
  int j, kdiff, action;
  realtype terkp1, erkp1, temp, hnew;

  nst++;
  kdiff = kk - kused;
  kused = kk;
  hused = hh;

  if ( (knew == kk-1) || (kk == maxord) ) phase = 1;

  /* For the first few steps, until either a step fails, or the order is 
  reduced, or the order reaches its maximum, we raise the order and double 
  the stepsize. During these steps, phase = 0. Thereafter, phase = 1, and
  stepsize and order are set by the usual local error algorithm. 

  Note that, after the first step, the order is not increased, as not all
  of the neccessary information is available yet.
  */
  
  if (phase == 0) {

    if (nst > 1) {
      kk++;
      hnew = TWO * hh;
      hh = hnew;
    }

  } else {

    action = UNSET;
    
    /* Set action = LOWER/MAINTAIN/RAISE to specify order decision */
    
    if (knew == kk-1)   {action = LOWER; goto takeaction;}
    if (kk == maxord)   {action = MAINTAIN; goto takeaction;}
    if ( (kk+1 >= ns ) || (kdiff == 1)) {action = MAINTAIN;goto takeaction;}
    
    /* Estimate the error at order k+1, unless already decided to
       reduce order, or already using maximum order, or stepsize has not
       been constant, or order was just raised. */
    
    N_VLinearSum (ONE, ee, -ONE, phi[kk+1], delta);
    terkp1 = IDAWrmsNorm(IDA_mem, delta, ewt, suppressalg);
    erkp1= terkp1/(kk+2);
    
    /* Choose among orders k-1, k, k+1 using local truncation error norms. */
    
    if (kk == 1) {
      if (terkp1 >= HALF * (*terk)) {action = MAINTAIN; goto takeaction;}
      else                          {action = RAISE;    goto takeaction;}
    }
    else {
      if (*terkm1 <= MIN(*terk, terkp1)) {action = LOWER;    goto takeaction;}
      if (terkp1  >= *terk)             {action = MAINTAIN; goto takeaction;}
      action = RAISE;
      goto takeaction;
    }
    
  takeaction:
    
    /* On change of order, reset kk and the estimated error norm. */
    
    if (action == RAISE) { kk++; *est = erkp1;}
    else if (action == LOWER) { kk--; *est = *erkm1;}
    
    /* Compute rr = tentative ratio hnew/hh from error norm.
       Reduce hh if rr <= 1, double hh if rr >= 2, else leave hh as is.
       If hh is reduced, hnew/hh is restricted to be between .5 and .9. */
    
    hnew = hh;
    rr = RPowerR( (TWO * (*est) + PT0001) , (-ONE/(kk+1) ) );
    
    if (rr >= TWO) {
      hnew = TWO * hh;
      if( (temp = ABS(hnew)*hmax_inv) > ONE ) hnew /= temp;
    }
    else if (rr <= ONE ) { 
      rr = MAX(HALF, MIN(PT9,rr));
      hnew = hh * rr;
    }
    
    hh = hnew;
    
  } /* end of phase if block */
  
  /* Save ee for possible order increase on next step; update phi array. */
  
  if (kused < maxord) N_VScale(ONE, ee, phi[kused+1]);
  
  N_VLinearSum(ONE, ee, ONE, phi[kused], phi[kused]);
  for (j= kused-1; j>=0; j--)
    N_VLinearSum(ONE, phi[j], ONE, phi[j+1], phi[j]);

  return (IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDAWrmsNorm
 * -----------------------------------------------------------------
 *  Returns the WRMS norm of vector x with weights w.
 *  If mask = TRUE, the weight vector w is masked by id, i.e.,
 *      nrm = N_VWrmsNormMask(x,w,id);
 *  Otherwise,
 *      nrm = N_VWrmsNorm(x,w);
 * 
 * mask = FALSE       when the call is made from the nonlinear solver.
 * mask = suppressalg otherwise.
 * -----------------------------------------------------------------
*/

realtype IDAWrmsNorm(IDAMem IDA_mem, N_Vector x, N_Vector w, 
                     booleantype mask)
{
  realtype nrm;

  if (mask) nrm = N_VWrmsNormMask(x, w, id);
  else      nrm = N_VWrmsNorm(x, w);

  return(nrm);
}


/*=================================================================*/
/*      Internal EWT function                                      */
/*=================================================================*/


/*  
 * -----------------------------------------------------------------
 * IDAEwtSet
 * -----------------------------------------------------------------
 * This routine is responsible for loading the error weight vector
 * ewt, according to itol, as follows:
 * (1) ewt[i] = 1 / (rtol * ABS(ycur[i]) + atol), i=0,...,Neq-1
 *     if itol = IDA_SS
 * (2) ewt[i] = 1 / (rtol * ABS(ycur[i]) + atol[i]), i=0,...,Neq-1
 *     if itol = IDA_SV
 *
 *  IDAEwtSet returns 0 if ewt is successfully set as above to a
 *  positive vector and -1 otherwise. In the latter case, ewt is
 *  considered undefined.
 *
 * All the real work is done in the routines IDAEwtSetSS, IDAEwtSetSV.
 * -----------------------------------------------------------------
 */

int IDAEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  IDAMem IDA_mem;
  int flag = 0;

  /* data points to IDA_mem here */

  IDA_mem = (IDAMem) data;

  switch(itol) {
  case IDA_SS: 
    flag = IDAEwtSetSS(ycur, weight, IDA_mem); 
    break;
  case IDA_SV: 
    flag = IDAEwtSetSV(ycur, weight, IDA_mem); 
    break;
  }
  return(flag);
}

/*
 * -----------------------------------------------------------------
 * IDAEwtSetSS
 * -----------------------------------------------------------------
 * This routine sets ewt as decribed above in the case itol=IDA_SS.
 * It tests for non-positive components before inverting. IDAEwtSetSS
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 * -----------------------------------------------------------------
 */

static int IDAEwtSetSS(N_Vector ycur, N_Vector weight, IDAMem IDA_mem)
{
  N_VAbs(ycur, tempv1);
  N_VScale(rtol, tempv1, tempv1);
  N_VAddConst(tempv1, Satol, tempv1);
  if (N_VMin(tempv1) <= ZERO) return(-1);
  N_VInv(tempv1, ewt);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * IDAEwtSetSV
 * -----------------------------------------------------------------
 * This routine sets ewt as decribed above in the case itol=IDA_SV.
 * It tests for non-positive components before inverting. IDAEwtSetSV
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 * -----------------------------------------------------------------
 */

static int IDAEwtSetSV(N_Vector ycur, N_Vector weight, IDAMem IDA_mem)
{
  N_VAbs(ycur, tempv1);
  N_VLinearSum(rtol, tempv1, ONE, Vatol, tempv1);
  if (N_VMin(tempv1) <= ZERO) return(-1);
  N_VInv(tempv1, ewt);
  return(0);
}

