/******************************************************************
 * File          : ida.c                                          *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and        *
 *                 Radu Serban @ LLNL                             *
 * Version of    : 3 July 2002                                    *
 *----------------------------------------------------------------*
 * This is the implementation file for the main IDA solver.       *
 * It is independent of the linear solver in use.                 *
 *                                                                *
 ******************************************************************/


/************************************************************/
/******************* BEGIN Imports **************************/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "ida.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"

/************************************************************/
/******************** END Imports ***************************/
/************************************************************/


/***************************************************************/
/*********************** BEGIN Macros **************************/
/**************************************************************/

/* Macro: loop */

#define loop for(;;)

/***************************************************************/
/************************ END Macros ***************************/
/***************************************************************/



/************************************************************/
/************** BEGIN IDA Private Constants ***************/
/************************************************************/

#define ZERO       RCONST(0.0)    /* real 0.0    */
#define HALF       RCONST(0.5)    /* real 0.5    */
#define QUARTER    RCONST(0.25)   /* real 0.25   */
#define TWOTHIRDS  RCONST(0.667)  /* real 2/3 for default steptol */
#define ONE        RCONST(1.0)    /* real 1.0    */
#define ONEPT5     RCONST(1.5)    /* real 1.5    */
#define TWO        RCONST(2.0)    /* real 2.0    */
#define TWOPT5     RCONST(2.5)    /* real 2.5    */
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

/***************************************************************/
/************** BEGIN Default Constants ************************/
/***************************************************************/

#define MXSTEP_DEFAULT   500   /* mxstep default value   */
#define MAXORD_DEFAULT   5     /* maxord default value   */
#define MXNCF           10     /* max number of convergence failures allowed */
#define MXETF           10     /* max number of error test failures allowed  */
#define EPCON      RCONST(0.33)   /* Newton convergence test constant */
#define EPICFAC    PT01 /* convergence test factor in IC calc. */
#define MAXNH      5    /* max. number of h tries in IC calc. */
#define MAXNJ      4    /* max. number of J tries in IC calc. */
#define MAXNI      10   /* max. Newton iterations in IC calc. */

/***************************************************************/
/*************** END Default Constants *************************/
/***************************************************************/


/***************************************************************/
/************ BEGIN Routine-Specific Constants *****************/
/***************************************************************/

 /* IDASolve return values are defined by an enum statement in the 
    include file ida.h.  They are: 
 
       (successful returns:)
 
       NORMAL_RETURN
       INTERMEDIATE_RETURN
       TSTOP_RETURN
 
       (failed returns:)

       IDA_NO_MEM
       ILL_INPUT
       TOO_MUCH_WORK
       TOO_MUCH_ACC
       ERR_FAILURE
       CONV_FAILURE
       SETUP_FAILURE
       SOLVE_FAILURE
       CONSTR_FAILURE   
       REP_RES_REC_ERR
       RES_NONREC_ERR
                       */

/* IDAStep return values */

/* The first three IDAStep return values are defined in an enum statement
   given in file ida.h. They are listed below for continuity with the other
   return values. The next one is defined here as an alias for the value
   defined in ida.h (LSOLVE_ERROR_NONRECVR) */

/* SUCCESS                      0 */
/* RES_ERROR_NONRECVR          -1 */
/* LSETUP_ERROR_NONRECVR       -2 */

#define CONV_FAIL_LINR_NONRECVR  LSOLVE_ERROR_NONRECVR   /* = -3 */

#define REP_ERR_FAIL           -4 
#define CONSTR_FAIL            -5
#define REP_CONV_FAIL          -6
#define REP_RES_ERR            -7


/* IDAInterp return values: */

#define OKAY     0
#define BAD_T    -1


/* IDACalcIC control constants */

#define ICRATEMAX  RCONST(0.9)    /* max. Newton conv. rate */
#define ALPHALS    RCONST(0.0001) /* alpha in linesearch conv. test */


/* IDAStep control constants */

#define PREDICT_AGAIN    20


/* IDANewtonIter constants */

#define RATEMAX RCONST(0.9)
#define MAXIT    4


/* Return values from various low-level routines */

/* Return values for lower level routines used by IDACalcIC */

enum { IC_FAIL_RECOV = 1,  IC_CONSTR_FAILED = 2,  IC_LINESRCH_FAILED = 3,
       IC_CONV_FAIL =  4,  IC_SLOW_CONVRG =   5 };

/* The following two return values are set by an enumeration in ida.h.
   They are included here for clarity. The values are potentially
   encountered in:  'res', 'lsetup', 'lsolve', IDANewtonIter, IDAnls,
   IDAHandleNFlag, IDAStep, and IDASolve. The third value is defined
   as an alias for the value returned from lsolve (LSOLVE_ERROR_RECVR) */

/*RES_ERROR_RECVR              +1    */
/*LSETUP_ERROR_RECVR           +2    */

#define CONV_FAIL_LINR_RECVR         LSOLVE_ERROR_RECVR  /*  = +3  */

#define CONV_FAIL_NLINR_RECVR  +4 /* IDANewtonIter, IDAnls                   */
#define CONSTRAINT_FAIL_RECVR  +5 /* IDAnls                                  */
#define CONTINUE_STEPS        +99 /* IDASolve, IDAStopTest1, IDAStopTest2    */

#define UNSET -1    /* IDACompleteStep */
#define LOWER 1     /* IDACompleteStep */
#define RAISE 2     /* IDACompleteStep */
#define MAINTAIN 3  /* IDACompleteStep */

#define ERROR_TEST_FAIL       +7

#define XRATE                RCONST(0.25)        


/***************************************************************/
/*************** END Routine-Specific Constants  ***************/
/***************************************************************/


/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* IDAMalloc error messages */

#define IDAM               "IDAMalloc/IDAReInit-- "

#define MSG_Y0_NULL        IDAM "y0=NULL illegal.\n\n"
#define MSG_YP0_NULL       IDAM "yp0=NULL illegal.\n\n"

#define MSG_BAD_NEQ        IDAM "Neq=%ld < 1 illegal.\n\n"

#define MSG_BAD_ITOL_1     IDAM "itol=%d illegal.\n"
#define MSG_BAD_ITOL_2     "The legal values are SS=%d and SV=%d.\n\n"
#define MSG_BAD_ITOL       MSG_BAD_ITOL_1 MSG_BAD_ITOL_2

#define MSG_RES_NULL       IDAM "res=NULL illegal.\n\n"

#define MSG_RTOL_NULL      IDAM "rtol=NULL illegal.\n\n"

#define MSG_BAD_CONSTRAINTS  IDAM "illegal values in constraints vector.\n\n"
#define MSG_MISSING_ID     IDAM "id = NULL but suppressalg option on.\n\n"
#define MSG_BAD_ID         IDAM "illegal values in id vector.\n\n"
 
#define MSG_BAD_RTOL       IDAM "*rtol=%g < 0 illegal.\n\n"

#define MSG_ATOL_NULL      IDAM "atol=NULL illegal.\n\n"

#define MSG_BAD_ATOL       IDAM "Some atol component < 0.0 illegal.\n\n"

#define MSG_BAD_OPTIN_1    IDAM "optIn=%d illegal.\n"
#define MSG_BAD_OPTIN_2    "The legal values are FALSE=%d and TRUE=%d.\n\n"
#define MSG_BAD_OPTIN      MSG_BAD_OPTIN_1 MSG_BAD_OPTIN_2

#define MSG_BAD_OPT        IDAM "optIn=TRUE, but iopt=ropt=NULL.\n\n"
#define MSG_BAD_NCONFAC    IDAM "optIn=TRUE, but ropt[NCONFAC] negative\n\n"

#define MSG_BAD_HMAX       IDAM "Negative hmax: %g\n"

#define MSG_MEM_FAIL       IDAM "A memory request failed.\n\n"

#define MSG_BAD_EWT        IDAM "Some initial ewt component = 0.0 illegal.\n\n"

#define MSG_Y0_FAIL_CONSTR IDAM "y0 fails to satisfy constraints.\n\n"

#define MSG_REI_NO_MEM  "IDAReInit-- ida_mem = NULL illegal.\n\n"

#define MSG_REI_MAXORD1 "IDAReInit-- Illegal attempt to increase "
#define MSG_REI_MAXORD2 "maximum method order from %d to %d.\n\n"
#define MSG_REI_MAXORD  MSG_REI_MAXORD1 MSG_REI_MAXORD2 


/* IDACalcIC error messages */

#define IDAIC              "IDACalcIC-- "

#define NO_MEM             "IDA_mem = NULL illegal.\n\n"

#define MSG_IC_IDA_NO_MEM  IDAIC NO_MEM
 
#define MSG_BAD_ICOPT      IDAIC "icopt = %d is illegal.\n\n"

#define MSG_IC_MISSING_ID  IDAIC "id = NULL conflicts with icopt.\n\n"

#define MSG_IC_BAD_ID      IDAIC "id has illegal values.\n\n"

#define MSG_IC_TOO_CLOSE1  IDAIC "tout1 = %g too close to t0 = %g to attempt"
#define MSG_IC_TOO_CLOSE2  " initial condition calculation.\n\n"
#define MSG_IC_TOO_CLOSE   MSG_IC_TOO_CLOSE1 MSG_IC_TOO_CLOSE2

#define MSG_BAD_EPICFAC    IDAIC "epicfac < 0.0 illegal.\n\n"

#define MSG_BAD_MAXNH      IDAIC "maxnh < 0 illegal.\n\n"

#define MSG_BAD_MAXNJ      IDAIC "maxnj < 0 illegal.\n\n"

#define MSG_BAD_MAXNIT     IDAIC "maxnit < 0 illegal.\n\n"

#define MSG_BAD_LSOFF      IDAIC "lsoff = %d is illegal.\n\n"

#define MSG_BAD_STEPTOL    IDAIC "steptol < 0.0 illegal.\n\n"

#define MSG_IC_LINIT_FAIL  IDAIC "The linear solver's init routine failed.\n\n"

#define MSG_IC_BAD_EWT     IDAIC "Some ewt component = 0.0 illegal.\n\n"

#define MSG_IC_RES_NONR1   IDAIC "Non-recoverable error return from"
#define MSG_IC_RES_NONR2   " ResFn residual routine. \n\n"
#define MSG_IC_RES_NONREC  MSG_IC_RES_NONR1 MSG_IC_RES_NONR2

#define MSG_IC_RES_FAIL1   IDAIC "Recoverable error in first call to"
#define MSG_IC_RES_FAIL2   " ResFn residual routine. Cannot recover. \n\n"
#define MSG_IC_RES_FAIL    MSG_IC_RES_FAIL1 MSG_IC_RES_FAIL2

#define MSG_IC_SETUP_FL1   IDAIC "The linear solver setup routine"
#define MSG_IC_SETUP_FL2   " failed non-recoverably.\n\n"
#define MSG_IC_SETUP_FAIL  MSG_IC_SETUP_FL1 MSG_IC_SETUP_FL2

#define MSG_IC_SOLVE_FL1   IDAIC "The linear solver solve routine"
#define MSG_IC_SOLVE_FL2   " failed non-recoverably.\n\n"
#define MSG_IC_SOLVE_FAIL  MSG_IC_SOLVE_FL1 MSG_IC_SOLVE_FL2

#define MSG_IC_NO_RECOV1   IDAIC "The residual routine or the linear"
#define MSG_IC_NO_RECOV2   " setup or solve routine had a recoverable"
#define MSG_IC_NO_RECOV3   " error, but IDACalcIC was unable to recover.\n\n"
#define MSG_IC_NO_RECOVERY MSG_IC_NO_RECOV1 MSG_IC_NO_RECOV2 MSG_IC_NO_RECOV3

#define MSG_IC_FAIL_CON1   IDAIC "Unable to satisfy the inequality"
#define MSG_IC_FAIL_CON2   " constraints.\n\n"
#define MSG_IC_FAIL_CONSTR MSG_IC_FAIL_CON1 MSG_IC_FAIL_CON2

#define MSG_IC_FAILED_LS1  IDAIC "The Linesearch algorithm failed"
#define MSG_IC_FAILED_LS2  " with too small a step.\n\n"
#define MSG_IC_FAILED_LINS MSG_IC_FAILED_LS1 MSG_IC_FAILED_LS2

#define MSG_IC_CONV_FAIL1  IDAIC "Failed to get convergence in"
#define MSG_IC_CONV_FAIL2  " Newton/Linesearch algorithm.\n\n"
#define MSG_IC_CONV_FAILED MSG_IC_CONV_FAIL1 MSG_IC_CONV_FAIL2



/* IDASolve error messages */

#define IDAS               "IDASolve-- "

#define MSG_IDA_NO_MEM     IDAS NO_MEM
 
#define MSG_LINIT_NULL     IDAS "The linear solver's init routine is NULL.\n\n"

#define MSG_LSETUP_NULL    IDAS "The linear solver's setup routine is NULL.\n\n"

#define MSG_LSOLVE_NULL    IDAS "The linear solver's solve routine is NULL.\n\n"

#define MSG_LFREE_NULL     IDAS "The linear solver's free routine is NULL.\n\n"

#define MSG_LINIT_FAIL     IDAS "The linear solver's init routine failed.\n\n"

#define MSG_BAD_HINIT      IDAS "hinit=%g and tout-t0=%g inconsistent.\n\n"

#define MSG_BAD_TOUT_1     IDAS "Trouble interpolating at tout = %g.\n"
#define MSG_BAD_TOUT_2     "tout too far back in direction of integration.\n\n"
#define MSG_BAD_TOUT       MSG_BAD_TOUT_1 MSG_BAD_TOUT_2

#define MSG_BAD_TSTOP_1    IDAS "tstop = %g is behind  current t = %g \n"
#define MSG_BAD_TSTOP_2    "in the direction of integration.\n\n"
#define MSG_BAD_TSTOP      MSG_BAD_TSTOP_1 MSG_BAD_TSTOP_2


#define MSG_MAX_STEPS_1    IDAS "At t=%g, mxstep=%d steps taken on "
#define MSG_MAX_STEPS_2    "this call before\nreaching tout=%g.\n\n"
#define MSG_MAX_STEPS      MSG_MAX_STEPS_1 MSG_MAX_STEPS_2

#define MSG_EWT_NOW_BAD_1  IDAS "At t=%g, "
#define MSG_EWT_NOW_BAD_2  "some ewt component has become <= 0.0.\n\n"
#define MSG_EWT_NOW_BAD    MSG_EWT_NOW_BAD_1 MSG_EWT_NOW_BAD_2

#define MSG_TOO_MUCH_ACC   IDAS "At t=%g, too much accuracy requested.\n\n"

#define MSG_ERR_FAILS_1    IDAS "At t=%g and step size h=%g, the error test\n"
#define MSG_ERR_FAILS_2    "failed repeatedly or with |h| = hmin.\n\n"
#define MSG_ERR_FAILS      MSG_ERR_FAILS_1 MSG_ERR_FAILS_2

#define MSG_CONV_FAILS_1   IDAS "At t=%g and step size h=%g, the corrector\n"
#define MSG_CONV_FAILS_2   "convergence failed repeatedly.\n\n"
#define MSG_CONV_FAILS     MSG_CONV_FAILS_1 MSG_CONV_FAILS_2

#define MSG_SETUP_FAILED_1 IDAS "At t=%g, the linear solver setup routine "
#define MSG_SETUP_FAILED_2 "failed in an unrecoverable manner.\n\n"
#define MSG_SETUP_FAILED   MSG_SETUP_FAILED_1 MSG_SETUP_FAILED_2

#define MSG_SOLVE_FAILED_1 IDAS "At t=%g, the linear solver solve routine "
#define MSG_SOLVE_FAILED_2 "failed in an unrecoverable manner.\n\n"
#define MSG_SOLVE_FAILED   MSG_SOLVE_FAILED_1 MSG_SOLVE_FAILED_2

#define MSG_TOO_CLOSE_1    IDAS "tout=%g too close to t0=%g to start"
#define MSG_TOO_CLOSE_2    " integration.\n\n"
#define MSG_TOO_CLOSE      MSG_TOO_CLOSE_1 MSG_TOO_CLOSE_2

#define MSG_YRET_NULL      IDAS "yret=NULL illegal.\n\n"
#define MSG_YPRET_NULL     IDAS "ypret=NULL illegal.\n\n"
#define MSG_TRET_NULL      IDAS "tret=NULL illegal.\n\n"

#define MSG_BAD_ITASK      IDAS "itask=%d illegal.\n\n"

#define MSG_REP_RES_ERR1   IDAS "At t = %g, repeated recoverable error \n"
#define MSG_REP_RES_ERR2   "returns from ResFn residual function. \n\n"
#define MSG_REP_RES_ERR    MSG_REP_RES_ERR1 MSG_REP_RES_ERR2

#define MSG_RES_NONRECOV1  IDAS "At t = %g, nonrecoverable error \n"
#define MSG_RES_NONRECOV2  "return from ResFn residual function. \n\n"
#define MSG_RES_NONRECOV   MSG_RES_NONRECOV1 MSG_RES_NONRECOV2

#define MSG_FAILED_CONSTR1 IDAS "At t = %g, unable to satisfy \n"
#define MSG_FAILED_CONSTR2 "inequality constraints. \n\n"
#define MSG_FAILED_CONSTR  MSG_FAILED_CONSTR1 MSG_FAILED_CONSTR2

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/


/************************************************************/
/*************** END IDA Private Constants ****************/
/************************************************************/


/**************************************************************/
/********* BEGIN Private Helper Functions Prototypes **********/
/**************************************************************/

static booleantype IDAAllocVectors(IDAMem IDA_mem, integertype Neq,  M_Env machEnv);
static void IDAFreeVectors(IDAMem IDA_mem);

static int IDAnlsIC (IDAMem IDA_mem);
static int IDANewtonIC (IDAMem IDA_mem);
static int IDALineSrch (IDAMem IDA_mem, realtype *delnorm, realtype *fnorm);
static int IDAfnorm (IDAMem IDA_mem, realtype *fnorm);
static int IDANewyyp (IDAMem IDA_mem, realtype lambda);
static int IDANewy (IDAMem IDA_mem);
static int IDAICFailFlag (IDAMem IDA_mem, int retval);

static booleantype IDAEwtSet(IDAMem IDA_mem, N_Vector ycur);
static booleantype IDAEwtSet0(IDAMem IDA_mem, N_Vector ycur);
static booleantype IDAEwtSetSS(IDAMem IDA_mem, N_Vector ycur);
static booleantype IDAEwtSetSV(IDAMem IDA_mem, N_Vector ycur);

static int IDAStopTest1(IDAMem IDA_mem, realtype tout, realtype tstop, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask);
static int IDAStopTest2(IDAMem IDA_mem, realtype tout, realtype tstop, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask);
static int IDAInterp(IDAMem IDA_mem, realtype t, N_Vector yret, N_Vector ypret);
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


/**************************************************************/
/********** END Private Helper Functions Prototypes ***********/
/**************************************************************/


/**************************************************************/
/**************** BEGIN Readability Constants *****************/
/**************************************************************/

#define uround   (IDA_mem->ida_uround)  
#define phi      (IDA_mem->ida_phi) 
#define ewt      (IDA_mem->ida_ewt)  
#define yy       (IDA_mem->ida_yy)
#define yp       (IDA_mem->ida_yp)
#define delta    (IDA_mem->ida_delta)
#define mm       (IDA_mem->ida_mm)
#define ee       (IDA_mem->ida_ee)
#define savres   (IDA_mem->ida_savres)
#define tempv1   (IDA_mem->ida_tempv1)
#define tempv2   (IDA_mem->ida_tempv2) 
#define kk       (IDA_mem->ida_kk)
#define hh       (IDA_mem->ida_hh)
#define tn       (IDA_mem->ida_tn)
#define tretp    (IDA_mem->ida_tretp)
#define cj       (IDA_mem->ida_cj)
#define cjold    (IDA_mem->ida_cjold)
#define cjratio  (IDA_mem->ida_cjratio)
#define cjlast   (IDA_mem->ida_cjlast)
#define maxord   (IDA_mem->ida_maxord)
#define mxstep   (IDA_mem->ida_mxstep) 
#define hmax_inv (IDA_mem->ida_hmax_inv)
#define nbacktr  (IDA_mem->ida_nbacktr)
#define nst      (IDA_mem->ida_nst)
#define nre      (IDA_mem->ida_nre)
#define ncfn     (IDA_mem->ida_ncfn)
#define netf     (IDA_mem->ida_netf)
#define nni      (IDA_mem->ida_nni)
#define nsetups  (IDA_mem->ida_nsetups)
#define ns       (IDA_mem->ida_ns)
#define lrw      (IDA_mem->ida_lrw)
#define liw      (IDA_mem->ida_liw)
#define linit    (IDA_mem->ida_linit)
#define lsetup   (IDA_mem->ida_lsetup)
#define lsolve   (IDA_mem->ida_lsolve) 
#define lperf    (IDA_mem->ida_lperf)
#define lfree    (IDA_mem->ida_lfree) 
#define lmem     (IDA_mem->ida_lmem) 
#define linitOK  (IDA_mem->ida_linitOK)
#define knew     (IDA_mem->ida_knew)
#define kused    (IDA_mem->ida_kused)          
#define hused    (IDA_mem->ida_hused)         
#define tolsf    (IDA_mem->ida_tolsf)      
#define machenv  (IDA_mem->ida_machenv)
#define phase    (IDA_mem->ida_phase)
#define epsNewt  (IDA_mem->ida_epsNewt)
#define nconfac  (IDA_mem->ida_nconfac)
#define toldel   (IDA_mem->ida_toldel)
#define ss       (IDA_mem->ida_ss)
#define rr       (IDA_mem->ida_rr)
#define psi      (IDA_mem->ida_psi)
#define alpha    (IDA_mem->ida_alpha)
#define beta     (IDA_mem->ida_beta)
#define sigma    (IDA_mem->ida_sigma)
#define gamma    (IDA_mem->ida_gamma)
#define setupNonNull (IDA_mem->ida_setupNonNull) 
#define suppressalg (IDA_mem->ida_suppressalg) 
#define constraintsSet (IDA_mem->ida_constraintsSet)

/**************************************************************/
/***************** END Readability Constants ******************/
/**************************************************************/


/***************************************************************/
/************* BEGIN IDA Implementation ************************/
/***************************************************************/


/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/


/******************** IDAMalloc *******************************

 IDAMalloc allocates and initializes memory for a problem. All
 problem specification inputs are checked for errors. If any
 error occurs during initialization, it is reported to the file
 whose file pointer is errfp and NULL is returned. Otherwise, the
 pointer to successfully initialized problem memory is returned.
 
*****************************************************************/

void *IDAMalloc(integertype Neq, ResFn res, void *rdata, realtype t0,
      N_Vector y0, N_Vector yp0, int itol, realtype *rtol, void *atol, 
      N_Vector id, N_Vector constraints, FILE *errfp, booleantype optIn, 
      long int iopt[], realtype ropt[], M_Env machEnv)
{
  booleantype allocOK, ioptExists, roptExists, neg_atol, ewtsetOK, conOK;
  realtype temptest;
  IDAMem IDA_mem;
  N_Vector mskewt;
  FILE *fp;
  
  /* Check for legal input parameters */
  
  fp = (errfp == NULL) ? stdout : errfp;

  if (Neq <= 0) { fprintf(fp, MSG_BAD_NEQ, (long int) Neq); return(NULL); }

  if (y0==NULL) { fprintf(fp, MSG_Y0_NULL); return(NULL); }
  
  if (yp0==NULL) { fprintf(fp, MSG_YP0_NULL); return(NULL); }

  if ((itol != SS) && (itol != SV)) {
    fprintf(fp, MSG_BAD_ITOL, itol, SS, SV);
    return(NULL);
  }

  if (res == NULL) { fprintf(fp, MSG_RES_NULL); return(NULL); }

  if (rtol == NULL) { fprintf(fp, MSG_RTOL_NULL); return(NULL); }

  if (*rtol < ZERO) { fprintf(fp, MSG_BAD_RTOL, *rtol); return(NULL); }
   
  if (atol == NULL) { fprintf(fp, MSG_ATOL_NULL); return(NULL); }

  if (itol == SS) { neg_atol = (*((realtype *)atol) < ZERO); }
  else { neg_atol = (N_VMin((N_Vector)atol) < ZERO); }
  if (neg_atol) { fprintf(fp, MSG_BAD_ATOL); return(NULL); }

  if ((optIn != FALSE) && (optIn != TRUE)) {
    fprintf(fp, MSG_BAD_OPTIN, optIn, FALSE, TRUE);
    return(NULL);
  }

  if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
    fprintf(fp, MSG_BAD_OPT);
    return(NULL);
  } 

  ioptExists = (iopt != NULL);
  roptExists = (ropt != NULL);

  if(optIn && roptExists) {
    if (ropt[HMAX] < ZERO) {
      fprintf(fp, MSG_BAD_HMAX, ropt[HMAX]);
      return(NULL);
    }
    if(ropt[NCONFAC] < ZERO) { fprintf(fp, MSG_BAD_NCONFAC); return(NULL); }
  }


 /* Allocate memory for the IDA memory structure */

  IDA_mem = (IDAMem) malloc(sizeof(struct IDAMemRec));
  if (IDA_mem == NULL) { fprintf(fp, MSG_MEM_FAIL); return(NULL); }

  /* Set maxord and suppressalg prior to allocating vectors */

  maxord = MAXORD_DEFAULT;
  suppressalg = FALSE;
  if (optIn && ioptExists) {
    if(iopt[SUPPRESSALG] == ONE) suppressalg = TRUE;
    if (iopt[MAXORD] > 0)  maxord = MIN(maxord, iopt[MAXORD]);
  }
 
  /* Test id vector for legality */

  if(suppressalg && (id==NULL)){ fprintf(fp, MSG_MISSING_ID); return(NULL); }
  if(suppressalg) {
    temptest = N_VMin(id);
    if(temptest < ZERO){ fprintf(fp, MSG_BAD_ID); return(NULL); }
  }

  /* Allocate the vectors */

  allocOK = IDAAllocVectors(IDA_mem, Neq, machEnv);
  if (!allocOK) {
    fprintf(fp, MSG_MEM_FAIL);
    free(IDA_mem);
    return(NULL);
  }
 
  /* Set the mskewt vector, set relevant memory pointers, and load ewt */

  if (suppressalg) mskewt = id;
  else mskewt = ewt;

  IDA_mem->ida_itol = itol;
  IDA_mem->ida_rtol = rtol;      
  IDA_mem->ida_atol = atol;  
  IDA_mem->ida_mskewt = mskewt;
  IDA_mem->ida_id  = id;

  ewtsetOK = IDAEwtSet(IDA_mem, y0);
  if (!ewtsetOK) {
    fprintf(fp, MSG_BAD_EWT);
    IDAFreeVectors(IDA_mem);
    free(IDA_mem);
    return(NULL);
  }

  /*  Check the constraints pointer and vector */
  
  if (constraints == NULL) constraintsSet = FALSE;
  else {
    constraintsSet = TRUE;
    temptest = N_VMaxNorm(constraints);
    if(temptest > TWOPT5){ fprintf(fp, MSG_BAD_CONSTRAINTS); return(NULL); }

    else if(temptest < HALF) constraintsSet = FALSE; /* constraints empty */
  }

  /* Check to see if y0 satisfies constraints. */

  if (constraintsSet) {
    conOK = N_VConstrMask (constraints, y0, tempv2);
    if (!conOK) { fprintf(fp, MSG_Y0_FAIL_CONSTR); return(NULL); }
  }

  /* Copy the input parameters into IDA memory block
     (Corresponding readability constants are defined below) */

  IDA_mem->ida_Neq = Neq;
  IDA_mem->ida_res = res;
  IDA_mem->ida_constraints = constraints;
  IDA_mem->ida_rdata = rdata;
  IDA_mem->ida_iopt = iopt;
  IDA_mem->ida_ropt = ropt;
  IDA_mem->ida_errfp = fp;
  IDA_mem->ida_y0  = y0;
  IDA_mem->ida_yp0 = yp0;

  tn = t0;
  machenv = machEnv;

  /* Set unit roundoff uround */

  uround = UnitRoundoff();


  /* Set the linear solver addresses to NULL, linitOK to FALSE */

  linit  = NULL;
  lsetup = NULL;
  lsolve = NULL;
  lperf  = NULL;
  lfree  = NULL;
  lmem = NULL;
  linitOK = FALSE;

  /* Initialize the phi array */

  N_VScale(ONE, y0, phi[0]);  
  N_VScale(ONE, yp0, phi[1]);  
 
  /* Handle the remaining optional inputs */

  hmax_inv = ZERO;
  nconfac  = ONE;
  if (optIn && roptExists) {
    if (ropt[HMAX] > ZERO) hmax_inv = ONE/ropt[HMAX];
    if (ropt[NCONFAC] > ZERO) nconfac = ropt[NCONFAC];
  }

  mxstep = MXSTEP_DEFAULT;
  if (optIn && ioptExists) if (iopt[MXSTEP] > 0) mxstep = iopt[MXSTEP];

  if ((!optIn) && roptExists) ropt[HINIT] = ZERO;

  /* All error checking is complete at this point */

    
  /* Initialize all the counters and other optional output values */
 
  nst = nre = ncfn = netf = nni = nsetups  = 0;
  
  kused = 0;
  hused = ZERO;
  tolsf = ONE;


  /* Initialize optional output locations in iopt, ropt */

  if (ioptExists) {
    iopt[NST] = iopt[NRE] = iopt[NSETUPS] = iopt[NNI] = 0;
    iopt[NCFN] = iopt[NETF] = iopt[NBACKTR] = 0;
    iopt[KUSED] = 0;
    iopt[KNEXT] = 0;
    iopt[LENRW] = lrw;
    iopt[LENIW] = liw;
  }
  
  if (roptExists) {
    ropt[HUSED] = ZERO;
    ropt[HNEXT] = ZERO;
    ropt[TNOW]  = t0;
    ropt[TOLSF] = ONE;
  }
      
  /* Problem has been successfully initialized */

  return((void *)IDA_mem);
}

/******************** IDAReInit ********************************

 IDAReInit re-initializes IDA's memory for a problem, assuming
 it has already beeen allocated in a prior IDAMalloc call.
 All problem specification inputs are checked for errors.
 The problem size Neq is assumed to be unchaged since the call
 to IDAMalloc, and the maximum order maxord must not be larger.
 If any error occurs during reinitialization, it is reported to
 the file whose file pointer is errfp.
 The return value is SUCCESS = 0 if no errors occurred, or
 a negative value otherwise.
 
*****************************************************************/

int IDAReInit(void *ida_mem, ResFn res, void *rdata, realtype t0,
      N_Vector y0, N_Vector yp0, int itol, realtype *rtol, void *atol, 
      N_Vector id, N_Vector constraints, FILE *errfp, booleantype optIn, 
      long int iopt[], realtype ropt[], void *machEnv)
{
  booleantype ioptExists, roptExists, neg_atol, ewtsetOK, conOK;
  realtype temptest;
  IDAMem IDA_mem;
  N_Vector mskewt;
  FILE *fp;
  int oldmaxord;

  /* Check for legal input parameters */
  
  fp = (errfp == NULL) ? stdout : errfp;

  if (ida_mem == NULL) {
    fprintf(fp, MSG_REI_NO_MEM);
    return(IDAREI_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (y0 == NULL) { fprintf(fp, MSG_Y0_NULL); return(IDAREI_ILL_INPUT); }
  
  if (yp0 == NULL) { fprintf(fp, MSG_YP0_NULL); return(IDAREI_ILL_INPUT); }

  if ((itol != SS) && (itol != SV)) {
    fprintf(fp, MSG_BAD_ITOL, itol, SS, SV);
    return(IDAREI_ILL_INPUT);
  }

  if (res == NULL) { fprintf(fp, MSG_RES_NULL); return(IDAREI_ILL_INPUT); }

  if (rtol == NULL) { fprintf(fp, MSG_RTOL_NULL); return(IDAREI_ILL_INPUT); }

  if (*rtol < ZERO) { fprintf(fp, MSG_BAD_RTOL, *rtol); return(IDAREI_ILL_INPUT); }
   
  if (atol == NULL) { fprintf(fp, MSG_ATOL_NULL); return(IDAREI_ILL_INPUT); }

  if (itol == SS) { neg_atol = (*((realtype *)atol) < ZERO); }
  else { neg_atol = (N_VMin((N_Vector)atol) < ZERO); }
  if (neg_atol) { fprintf(fp, MSG_BAD_ATOL); return(IDAREI_ILL_INPUT); }

  if ((optIn != FALSE) && (optIn != TRUE)) {
    fprintf(fp, MSG_BAD_OPTIN, optIn, FALSE, TRUE);
    return(IDAREI_ILL_INPUT);
  }

  if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
    fprintf(fp, MSG_BAD_OPT);
    return(IDAREI_ILL_INPUT);
  } 

  ioptExists = (iopt != NULL);
  roptExists = (ropt != NULL);

  if(optIn && roptExists) {
    if (ropt[HMAX] < ZERO) {
      fprintf(fp, MSG_BAD_HMAX, ropt[HMAX]);
      return(IDAREI_ILL_INPUT);
    }
    if(ropt[NCONFAC] < ZERO) { fprintf(fp, MSG_BAD_NCONFAC); return(IDAREI_ILL_INPUT); }
  }

  /* Set maxord and suppressalg, and check new vs old maxord */

  oldmaxord = maxord;
  maxord = MAXORD_DEFAULT;
  suppressalg = FALSE;
  if (optIn && ioptExists) {
    if(iopt[SUPPRESSALG] == ONE) suppressalg = TRUE;
    if (iopt[MAXORD] > 0) maxord = MIN(maxord, iopt[MAXORD]);
  }
  if (maxord > oldmaxord) {
    fprintf(fp, MSG_REI_MAXORD, oldmaxord, maxord);
    return(IDAREI_ILL_INPUT);
  }
 
  /* Test id vector for legality */

  if(suppressalg && (id==NULL)){ fprintf(fp, MSG_MISSING_ID); return(IDAREI_ILL_INPUT); }
  if(suppressalg) {
    temptest = N_VMin(id);
    if(temptest < ZERO){ fprintf(fp, MSG_BAD_ID); return(IDAREI_ILL_INPUT); }
  }

  /* Set the mskewt vector, set relevant memory pointers, and load ewt */

  if (suppressalg) mskewt = id;
  else mskewt = ewt;

  IDA_mem->ida_itol = itol;
  IDA_mem->ida_rtol = rtol;      
  IDA_mem->ida_atol = atol;  
  IDA_mem->ida_mskewt = mskewt;
  IDA_mem->ida_id  = id;

  ewtsetOK = IDAEwtSet(IDA_mem, y0);
  if (!ewtsetOK) {
    fprintf(fp, MSG_BAD_EWT);
    return(IDAREI_ILL_INPUT);
  }

  /*  Check the constraints pointer and vector */
  
  if (constraints == NULL) constraintsSet = FALSE;
  else {
    constraintsSet = TRUE;
    temptest = N_VMaxNorm(constraints);
    if(temptest > TWOPT5){ fprintf(fp, MSG_BAD_CONSTRAINTS); return(IDAREI_ILL_INPUT); }

    else if(temptest < HALF) constraintsSet = FALSE; /* constraints empty */
  }

  /* Check to see if y0 satisfies constraints. */

  if (constraintsSet) {
    conOK = N_VConstrMask (constraints, y0, tempv2);
    if (!conOK) { fprintf(fp, MSG_Y0_FAIL_CONSTR); return(IDAREI_ILL_INPUT); }
  }

  /* Copy the input parameters into IDA memory block
     (Corresponding readability constants are defined below.) */

  IDA_mem->ida_res = res;
  IDA_mem->ida_constraints = constraints;
  IDA_mem->ida_rdata = rdata;
  IDA_mem->ida_iopt = iopt;
  IDA_mem->ida_ropt = ropt;
  IDA_mem->ida_errfp = fp;
  IDA_mem->ida_y0  = y0;
  IDA_mem->ida_yp0 = yp0;

  tn = t0;
  machenv = machEnv;

  /* Set unit roundoff uround */

  uround = UnitRoundoff();


  /* Set linitOK to FALSE */

  linitOK = FALSE;

  /* Initialize the phi array */

  N_VScale(ONE, y0, phi[0]);  
  N_VScale(ONE, yp0, phi[1]);  
 
  /* Handle the remaining optional inputs */

  hmax_inv = ZERO;
  nconfac  = ONE;
  if (optIn && roptExists) {
    if (ropt[HMAX] > ZERO) hmax_inv = ONE/ropt[HMAX];
    if (ropt[NCONFAC] > ZERO) nconfac = ropt[NCONFAC];
  }

  mxstep = MXSTEP_DEFAULT;
  if (optIn && ioptExists) if (iopt[MXSTEP] > 0) mxstep = iopt[MXSTEP];

  if ((!optIn) && roptExists) ropt[HINIT] = ZERO;

  /* All error checking is complete at this point */

    
  /* Initialize all the counters and other optional output values */
 
  nst = nre = ncfn = netf = nni = nsetups  = 0;
  
  kused = 0;
  hused = ZERO;
  tolsf = ONE;


  /* Initialize optional output locations in iopt, ropt */

  if (ioptExists) {
    iopt[NST] = iopt[NRE] = iopt[NSETUPS] = iopt[NNI] = 0;
    iopt[NCFN] = iopt[NETF] = iopt[NBACKTR] = 0;
    iopt[KUSED] = 0;
    iopt[KNEXT] = 0;
    iopt[LENRW] = lrw;
    iopt[LENIW] = liw;
  }
  
  if (roptExists) {
    ropt[HUSED] = ZERO;
    ropt[HNEXT] = ZERO;
    ropt[TNOW]  = t0;
    ropt[TOLSF] = ONE;
  }
      
  /* Problem has been successfully re-initialized */

  return(SUCCESS);
}


/**************************************************************/
/************** BEGIN More Readability Constants **************/
/**************************************************************/
#define y0       (IDA_mem->ida_y0)
#define yp0      (IDA_mem->ida_yp0)
#define res      (IDA_mem->ida_res)
#define rdata    (IDA_mem->ida_rdata)
#define itol     (IDA_mem->ida_itol)
#define rtol     (IDA_mem->ida_rtol)
#define atol     (IDA_mem->ida_atol)
#define iopt     (IDA_mem->ida_iopt)
#define ropt     (IDA_mem->ida_ropt)
#define errfp    (IDA_mem->ida_errfp)
#define id       (IDA_mem->ida_id)
#define mskewt   (IDA_mem->ida_mskewt)
#define constraints (IDA_mem->ida_constraints)
#define sysindex (IDA_mem->ida_sysindex)
#define tscale   (IDA_mem->ida_tscale)
/**************************************************************/
/*************** END More Readability Constants ***************/
/**************************************************************/


/******************** IDACalcIC *********************************

 IDACalcIC computes consistent initial conditions, given the 
 user's initial guess for unknown components of y0 and/or yp0.

 The return value is SUCCESS = 0 if no error occurred.

 The error return values (fully described in ida.h) are:
    IC_IDA_NO_MEM      ida_mem is NULL
    IC_ILL_INPUT       bad value for icopt, tout1, epicfac, maxnh,
                       maxnj, maxnit, lsoff, steptol, or id
    IC_LINIT_FAIL      the linear solver linit routine failed
    IC_BAD_EWT         zero value of some component of ewt
    RES_NONRECOV_ERR   res had a non-recoverable error
    IC_FIRST_RES_FAIL  res failed recoverably on the first call
    SETUP_FAILURE      lsetup had a non-recoverable error
    SOLVE_FAILURE      lsolve had a non-recoverable error
    IC_NO_RECOVERY     res, lsetup, or lsolve had a recoverable
                       error, but IDACalcIC could not recover
    IC_FAILED_CONSTR   the inequality constraints could not be met
    IC_FAILED_LINESRCH the linesearch failed (on steptol test)
    IC_CONV_FAILURE    the Newton iterations failed to converge

 
*****************************************************************/

int IDACalcIC (void *ida_mem, int icopt, realtype tout1, realtype epicfac, 
               int maxnh, int maxnj, int maxnit, int lsoff, realtype steptol)
{
  booleantype ewtsetOK;
  int nwt, nh, retval, maxnh1, mxnh, icret;
  realtype tdist, troundoff, minid, hic, ypnorm;
  IDAMem IDA_mem;

  /* Check legality of input arguments, and set IDA memory copies. */

  IDA_mem = (IDAMem) ida_mem;
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IC_IDA_NO_MEM);
    return(IC_IDA_NO_MEM);
  }

  if (icopt < CALC_YA_YDP_INIT || icopt > CALC_Y_INIT) {
    fprintf(errfp, MSG_BAD_ICOPT, icopt);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_icopt = icopt;

  if (icopt == CALC_YA_YDP_INIT && (id == NULL)) {
    fprintf(errfp, MSG_IC_MISSING_ID);
    return(IC_ILL_INPUT);
  }

  tdist = ABS(tout1 - tn);
  troundoff = TWO*uround*(ABS(tn) + ABS(tout1));    
  if (tdist < troundoff) {
    fprintf(errfp, MSG_IC_TOO_CLOSE, tout1, tn);
    return(IC_ILL_INPUT);
  }

  /* For use in the CALC_YA_YP_INIT case, set sysindex and tscale. */

  sysindex = 1;
  tscale = tdist;
  if (icopt == CALC_YA_YDP_INIT) {
    minid = N_VMin(id);
    if (minid < ZERO) {
      fprintf(errfp, MSG_IC_BAD_ID);
      return(IC_ILL_INPUT);
    }
    if (minid > HALF) sysindex = 0;
  }

  /* For optional inputs, check legality and apply defaults if zero. */

  if (epicfac < ZERO) {
    fprintf(errfp, MSG_BAD_EPICFAC);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_epsNewt = (epicfac == ZERO) ? EPICFAC*EPCON : epicfac*EPCON;

  if (icopt == CALC_YA_YDP_INIT) {
    if (maxnh < 0) {
      fprintf(errfp, MSG_BAD_MAXNH);
      return(IC_ILL_INPUT);
    }
    maxnh1 = (maxnh == 0) ? MAXNH : maxnh;
  }

  if (maxnj < 0) {
    fprintf(errfp, MSG_BAD_MAXNJ);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_maxnj = (maxnj == 0) ? MAXNJ : maxnj;

  if (maxnit < 0) {
    fprintf(errfp, MSG_BAD_MAXNIT);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_maxnit = (maxnit == 0) ? MAXNI : maxnit;

  if (lsoff < 0 || lsoff > 1) {
    fprintf(errfp, MSG_BAD_LSOFF, lsoff);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_lsoff = lsoff;


  if (steptol < ZERO) {
    fprintf(errfp, MSG_BAD_STEPTOL);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_steptol = (steptol == ZERO) ? 
                         RPowerR(uround,TWOTHIRDS) : steptol;

  /* Initializations: cjratio = 1 (for use in direct linear solvers); 
     set nbacktr = 0; call linit routine. */

  cjratio = ONE;
  nbacktr = 0;
  linitOK = (linit(IDA_mem) == LINIT_OK);
  if (!linitOK) {
    fprintf(errfp, MSG_IC_LINIT_FAIL);
    return(IC_LINIT_FAIL);
  }

  /* Set hic, hh, cj, and mxnh. */
  hic = PT001*tdist;
  ypnorm = N_VWrmsNorm (yp0, mskewt);
  if (ypnorm > HALF/hic) hic = HALF/ypnorm;
  if( tout1 < tn) hic = -hic;
  hh = hic;
  if (icopt == CALC_YA_YDP_INIT) {
    cj = ONE/hic;
    mxnh = maxnh1;
  }
  else {
    cj = ZERO;
    mxnh = 1;
  }

  /* If suppressalg is on, reset id to bit vector form. */
  if (suppressalg) N_VOneMask (id);

  /* Loop over nwt = number of evaluations of ewt vector. */

  for (nwt = 1; nwt <= 2; nwt++) {
 
    /* Loop over nh = number of h values. */
    for (nh = 1; nh <= mxnh; nh++) {

      /* Call the IC nonlinear solver function. */
      retval = IDAnlsIC(IDA_mem);

      /* Cut h and loop on recoverable CALC_YA_YDP_INIT failure; else break. */
      if (retval == SUCCESS) break;
      ncfn++;
      if (retval < 0) break;
      if (nh == mxnh) break;
      /* If looping to try again, reset y0 and yp0 if not converging. */
      if (retval != IC_SLOW_CONVRG) {
        N_VScale (ONE, phi[0], y0);
        N_VScale (ONE, phi[1], yp0);
      }
      hic *= PT1;
      cj = ONE/hic;
      hh = hic;
    }   /* End of nh loop */

    /* Break on failure; else reset ewt, save y0,yp0 in phi, and loop. */
    if (retval != SUCCESS) break;
    ewtsetOK = IDAEwtSet0(IDA_mem, y0);
    if (!ewtsetOK) { retval = IC_BAD_EWT; break; }
    N_VScale (ONE, y0,  phi[0]);
    N_VScale (ONE, yp0, phi[1]);

  }   /* End of nwt loop */


  /* If suppressalg is on, reset mskewt. */
  if (suppressalg) N_VProd (id, ewt, mskewt);

  /* Load the optional outputs. */
  if (icopt == CALC_YA_YDP_INIT && ropt != NULL) ropt[HUSED] = hic;
  if (iopt != NULL) {
    iopt[NRE] = nre;
    iopt[NNI] = nni;
    iopt[NCFN] = ncfn;
    iopt[NSETUPS] = nsetups;
    iopt[NBACKTR] = nbacktr;
  }

  /* On any failure, print message and return proper flag. */
  if (retval != SUCCESS) {
    icret = IDAICFailFlag(IDA_mem, retval);
    return(icret);
  }

  /* Otherwise return success flag. */
  return(SUCCESS);

}


/********************* IDASolve ***********************************

 This routine is the main driver of the IDA package. 

 It integrates over an independent variable interval defined by the user, 
 by calling IDAStep to take internal independent variable steps.

 The first time that IDASolve is called for a successfully initialized
 problem, it computes a tentative initial step size.

 IDASolve supports four modes, specified by itask:
    NORMAL,  ONE_STEP,  NORMAL_TSTOP,  and  ONE_STEP_TSTOP.
 In the NORMAL and NORMAL-TSTOP modes, the solver steps until it 
 passes tout and then interpolates to obtain y(tout) and yp(tout).
 In the ONE_STEP and ONE_STEP_TSTOP modes, it takes one internal step
 and returns.  In the NORMAL_TSTOP and ONE_STEP_TSTOP modes, it also
 takes steps so as to reach tstop exactly and never to go past it.

 IDASolve returns integer values corresponding to success and failure as below:

       (successful returns:)     (failed returns:)
 
       NORMAL_RETURN             ILL_INPUT               TOO_MUCH_WORK
       INTERMEDIATE_RETURN       IDA_NO_MEM              TOO_MUCH_ACC
       TSTOP_RETURN              CONV_FAILURE            SETUP_FAILURE
                                 SOLVE_FAILURE           CONSTR_FAILURE
                                 ERR_FAILURE             REP_RES_REC_ERR
                                 RES_NONREC_ERR

********************************************************************/

int IDASolve(void *ida_mem, realtype tout, realtype tstop, realtype *tret,
             N_Vector yret, N_Vector ypret, int itask)
{
  int nstloc, sflag, istate, ier;
  realtype tdist, troundoff, ypnorm, rh;
  booleantype ewtsetOK;
  IDAMem IDA_mem;

  /* Check for legal inputs in all cases. */

  IDA_mem = (IDAMem) ida_mem;
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDA_NO_MEM);
    return(IDA_NO_MEM);
  }
  
  if (yret == NULL) {
    fprintf(errfp, MSG_YRET_NULL);       
    return(ILL_INPUT);
  }
  yy = yret;  

  if (ypret == NULL) {
    fprintf(errfp, MSG_YPRET_NULL);       
    return(ILL_INPUT);
  }
  yp = ypret;
  
  if (tret == NULL) {
    fprintf(errfp, MSG_TRET_NULL);
    return(ILL_INPUT);
  }
  *tret = tretp = tn; /* Set tret now in case of illegal-input return. */

  if ((itask < NORMAL) || (itask > ONE_STEP_TSTOP)) {
    fprintf(errfp, MSG_BAD_ITASK, itask);
    return(ILL_INPUT);
  }

  /* On first call, check linear solver functions and call linit function. */
  
  if (nst == 0) {
    if (linit == NULL) {
      fprintf(errfp, MSG_LINIT_NULL);
      return(ILL_INPUT);
    }
    if (lsetup == NULL) {
      fprintf(errfp, MSG_LSETUP_NULL);
      return(ILL_INPUT);
    }
    if (lsolve == NULL) {
      fprintf(errfp, MSG_LSOLVE_NULL);
      return(ILL_INPUT);
    }
    if (lfree == NULL) {
      fprintf(errfp, MSG_LFREE_NULL);
      return(ILL_INPUT);
    }
    /* Call linit if not already successfully called by IDACalcIC */
    if (!linitOK) {
      linitOK = (linit(IDA_mem) == LINIT_OK);
      if (!linitOK) {
        fprintf(errfp, MSG_LINIT_FAIL);
        return(ILL_INPUT);
      }
    }

    /* On the first call, check for tout - tn too small,
       set initial hh (from HINIT or internally),
       check for approach to tstop, and scale phi[1] by hh. */

    tdist = ABS(tout - tn);
    troundoff = TWO*uround*(ABS(tn) + ABS(tout));    
    if (tdist < troundoff) {
      fprintf(errfp, MSG_TOO_CLOSE, tout, tn);
      return(ILL_INPUT);
    }

    hh = ZERO;
    if (ropt != NULL) hh = ropt[HINIT];
    if ( (hh != ZERO) && ((tout-tn)*hh < ZERO) ) {
      fprintf(errfp, MSG_BAD_HINIT, hh, tout-tn);
      return(ILL_INPUT);
    }

    if (hh == ZERO) {
      hh = PT001*tdist;
      ypnorm = N_VWrmsNorm(phi[1], mskewt);
      if (ypnorm > HALF/hh) hh = HALF/ypnorm;
      if(tout < tn) hh = -hh;
    }

    rh = ABS(hh)*hmax_inv;
    if (rh > ONE) hh /= rh;

    if ( (itask == NORMAL_TSTOP) || (itask == ONE_STEP_TSTOP) ) {
      if ( (tstop - tn)*hh < ZERO) {
        fprintf(errfp, MSG_BAD_TSTOP, tstop, tn);
        return(ILL_INPUT);
      }
      if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
    }

    N_VScale(hh, phi[1], phi[1]);
    kk = 0; kused = 0;  /* set in case of an error return before a step */

    /* Set the convergence test constants epsNewt and toldel */

    epsNewt = nconfac * EPCON;
    toldel = PT0001 * epsNewt;

  } /* end of first-call block. */

  /* Call lperf function and set nstloc for later performance testing. */

  if (lperf != NULL) lperf(IDA_mem, 0);
  nstloc = 0;

  /* If not the first call, check for stop conditions. */

  if (nst > 0) {
    istate = IDAStopTest1(IDA_mem, tout, tstop, tret, yret, ypret, itask);
    if (istate != CONTINUE_STEPS) return(istate);
  }

  /* Looping point for internal steps. */

  loop {
   
    /* Check for too many steps taken. */
    
    if (nstloc >= mxstep) {
      fprintf(errfp, MSG_MAX_STEPS, tn, mxstep, tout);
      istate = TOO_MUCH_WORK;
      *tret = tretp = tn;
      break; /* Here yy=yret and yp=ypret already have the current solution. */
    }

    /* Call lperf to generate warnings of poor performance. */

    if (lperf != NULL) lperf(IDA_mem, 1);

    /* Reset and check ewt and mskewt (if not first call). */

    if (nst > 0) {
      ewtsetOK = IDAEwtSet(IDA_mem, phi[0]);
      if (!ewtsetOK) {
        fprintf(errfp, MSG_EWT_NOW_BAD, tn);
        istate = ILL_INPUT;
        ier = IDAInterp(IDA_mem, tn, yret, ypret);
        *tret = tretp = tn;
        break;
      }
    }
    
    /* Check for too much accuracy requested. */
    
    if ((tolsf = uround * N_VWrmsNorm(phi[0], ewt)) > ONE) {
      tolsf *= TEN;
      fprintf(errfp, MSG_TOO_MUCH_ACC, tn);
      istate = TOO_MUCH_ACC;
      *tret = tretp = tn;
      if (nst > 0) ier = IDAInterp(IDA_mem, tn, yret, ypret);
      break;
    }

    /* Call IDAStep to take a step. */

    sflag = IDAStep(IDA_mem);

    /* Process all failed-step cases, and exit loop. */
   
    if (sflag != SUCCESS) {
      istate = IDAHandleFailure(IDA_mem, sflag);
      *tret = tretp = tn;
      ier = IDAInterp(IDA_mem, tn, yret, ypret);
      break;
    }
    
    nstloc++;

    /* After successful step, check for stop conditions; continue or break. */

    istate = IDAStopTest2(IDA_mem, tout, tstop, tret, yret, ypret, itask);
    if (istate != CONTINUE_STEPS) break;

  } /* End of step loop */

  /* Load optional outputs and return. */

  if (iopt != NULL) {
    iopt[NST] = nst;
    iopt[NRE] = nre;
    iopt[NSETUPS] = nsetups;
    iopt[NNI] = nni;
    iopt[NCFN] = ncfn;
    iopt[NETF] = netf;
    iopt[KUSED] = kused;
    iopt[KNEXT] = kk;
  }
  
  if (ropt != NULL) {
    ropt[HUSED] = hused;
    ropt[HNEXT] = hh;
    ropt[TNOW]  = tn;
    ropt[TOLSF] = tolsf;
  }
  
  return(istate);    
}



/********************* IDAFree **********************************

 This routine frees the problem memory allocated by IDAMalloc
 Such memory includes all the vectors allocated by IDAAllocVectors,
 and the memory lmem for the linear solver (deallocated by a call
 to lfree).

*******************************************************************/

void IDAFree(void *ida_mem)
{
  IDAMem IDA_mem;

  IDA_mem = (IDAMem) ida_mem;
  
  if (IDA_mem == NULL) return;

  IDAFreeVectors(IDA_mem);
  lfree(IDA_mem);
  free(IDA_mem);
}


/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/


/*******************************************************************/
/******** BEGIN Private Helper Functions Implementation ************/
/*******************************************************************/
 
/****************** IDAAllocVectors ***********************************

 This routine allocates the IDA vectors ewt, tempv1, tempv2, and
 phi[0], ..., phi[maxord]. The length of the vectors is the input
 parameter Neq and the maximum order (needed to allocate phi) is the
 input parameter maxord. If all memory allocations are successful,
 IDAAllocVectors returns TRUE. Otherwise all allocated memory is freed
 and IDAAllocVectors returns FALSE.
 This routine also sets the optional outputs lrw and liw, which are
 (respectively) the lengths of the real and integer work spaces
 allocated here.

**********************************************************************/

static booleantype IDAAllocVectors(IDAMem IDA_mem, integertype Neq, M_Env machEnv)
{
  int i, j, maxcol;

  /* Allocate ewt, ee, delta, tempv1, tempv2 */
  
  ewt = N_VNew(Neq, machEnv);
  if (ewt == NULL) return(FALSE);
  ee = N_VNew(Neq, machEnv);
  if (ee == NULL) {
    N_VFree(ewt);
    return(FALSE);
  }
  delta = N_VNew(Neq, machEnv);
  if (delta == NULL) {
    N_VFree(ewt);
    N_VFree(ee);
    return(FALSE);
  }
  tempv1 = N_VNew(Neq, machEnv);
  if (tempv1 == NULL) {
    N_VFree(ewt);
    N_VFree(ee);
    N_VFree(delta);
    return(FALSE);
  }
   tempv2= N_VNew(Neq, machEnv);
  if (tempv2 == NULL) {
    N_VFree(ewt);
    N_VFree(ee);
    N_VFree(delta);
    N_VFree(tempv1);
    return(FALSE);
  }

  savres = tempv1;

  /* Allocate phi[0] ... phi[maxord].  Make sure phi[2] and phi[3] are
  allocated (for use as temporary vectors), regardless of maxord.       */

  maxcol = MAX(maxord,3);
  for (j=0; j <= maxcol; j++) {
    phi[j] = N_VNew(Neq, machEnv);
    if (phi[j] == NULL) {
      N_VFree(ewt);
      N_VFree(ee);
      N_VFree(delta);
      N_VFree(tempv1);
      N_VFree(tempv2);
      for (i=0; i < j; i++) N_VFree(phi[i]);
      return(FALSE);
    }
  }

  /* Set solver workspace lengths  */

  lrw = (maxcol + 6)*Neq;
  liw = 0;

  return(TRUE);
}


#define Neq    (IDA_mem->ida_Neq)


/***************** IDAFreeVectors *********************************
  
 This routine frees the IDA vectors allocated in IDAAllocVectors.

******************************************************************/

static void IDAFreeVectors(IDAMem IDA_mem)
{
  int j, maxcol;
  
  N_VFree(ewt);
  N_VFree(ee);
  N_VFree(delta);
  N_VFree(tempv1);
  N_VFree(tempv2);
  maxcol = MAX(maxord,3);
  for(j=0; j <= maxcol; j++) N_VFree(phi[j]);
}


/**************************************************************/
/*** BEGIN More Readability Constants for IDACalcIC routines **/
/**************************************************************/
#define icopt     (IDA_mem->ida_icopt)
#define epsic     (IDA_mem->ida_epsic)
#define maxnj     (IDA_mem->ida_maxnj)
#define maxnit    (IDA_mem->ida_maxnit)
#define lsoff     (IDA_mem->ida_lsoff)
#define steptol   (IDA_mem->ida_steptol)
#define ynew      (IDA_mem->ida_ynew)
#define ypnew     (IDA_mem->ida_ypnew)
#define delnew    (IDA_mem->ida_delnew)
#define dtemp     (IDA_mem->ida_dtemp)
/**************************************************************/
/*** END More Readability Constants for IDACalcIC routines ****/
/**************************************************************/


/******************** IDAnlsIC **********************************

 IDAnlsIC solves a nonlinear system for consistent initial 
 conditions.  It calls IDANewtonIC to do most of the work.

 The return value is SUCCESS = 0 if no error occurred.
 The error return values (positive) considered recoverable are:
    IC_FAIL_RECOV      if res, lsetup, or lsolve failed recoverably
    IC_CONSTR_FAILED   if the constraints could not be met
    IC_LINESRCH_FAILED if the linesearch failed (on steptol test)
    IC_CONV_FAIL       if the Newton iterations failed to converge
    IC_SLOW_CONVRG     if the iterations are converging slowly
                       (failed the convergence test, but showed
                       norm reduction or convergence rate < 1)
 The error return values (negative) considered non-recoverable are:
    RES_NONRECOV_ERR   if res had a non-recoverable error
    IC_FIRST_RES_FAIL  if res failed recoverably on the first call
    SETUP_FAILURE      if lsetup had a non-recoverable error
    SOLVE_FAILURE      if lsolve had a non-recoverable error
 
*****************************************************************/

static int IDAnlsIC (IDAMem IDA_mem)
{
  int retval, nj;
  N_Vector tv1, tv2, tv3;

  tv1 = ee;
  tv2 = tempv2;
  tv3 = phi[2];

  retval = res (Neq, tn, y0, yp0, delta, rdata);
  nre++;
  if(retval < 0) return(RES_NONRECOV_ERR);
  if(retval > 0) return(IC_FIRST_RES_FAIL);

  N_VScale (ONE, delta, savres);

  /* Loop over nj = number of linear solve Jacobian setups. */

  for (nj = 1; nj <= maxnj; nj++) {

    /* If there is a setup routine, call it. */
    if (setupNonNull) {
      nsetups++;
      retval = lsetup (IDA_mem, y0, yp0, delta, tv1, tv2, tv3);
      if(retval < 0) return(SETUP_FAILURE);
      if(retval > 0) return(IC_FAIL_RECOV);
    }

    /* Call the Newton iteration routine, and return if successful.  */
    retval = IDANewtonIC(IDA_mem);
    if (retval == SUCCESS) return(SUCCESS);

    /* If converging slowly and lsetup is nontrivial, retry. */
    if (retval == IC_SLOW_CONVRG && setupNonNull) {
      N_VScale (ONE, savres, delta);
      continue;
    }

    else return(retval);

  }   /* End of nj loop */

  /* No convergence after maxnj tries; return failure flag. */
  return(retval);

}

/******************** IDANewtonIC ************************************

 IDANewtonIC performs the Newton iteration to solve for consistent
 initial conditions.  It calls IDALineSrch within each iteration.
 On return, savres contains the current residual vector.

 The return value is SUCCESS = 0 if no error occurred.
 The error return values (positive) considered recoverable are:
    IC_FAIL_RECOV      if res or lsolve failed recoverably
    IC_CONSTR_FAILED   if the constraints could not be met
    IC_LINESRCH_FAILED if the linesearch failed (on steptol test)
    IC_CONV_FAIL       if the Newton iterations failed to converge
    IC_SLOW_CONVRG     if the iterations appear to be converging slowly.
                       They failed the convergence test, but showed 
                       an overall norm reduction (by a factor of < 0.1)
                       or a convergence rate <= ICRATEMAX).
 The error return values (negative) considered non-recoverable are:
    RES_NONRECOV_ERR   if res had a non-recoverable error
    SOLVE_FAILURE      if lsolve had a non-recoverable error
 
    **********************************************************************/

static int IDANewtonIC (IDAMem IDA_mem)
{
  int retval, mnewt;
  realtype delnorm, fnorm, fnorm0, oldfnrm, rate;

  /* Set pointer for vector delnew */
  delnew = phi[2];

  /* Call the linear solve function to get the Newton step, delta. */
  retval = lsolve (IDA_mem, delta, y0, yp0, savres);
  if(retval < 0) return(SOLVE_FAILURE);
  if(retval > 0) return(IC_FAIL_RECOV);

  /* Compute the norm of the step; return now if this is small. */
  fnorm = N_VWrmsNorm (delta, ewt);
  if (sysindex == 0) fnorm *= tscale*abs(cj);
  if (fnorm <= epsNewt) return(SUCCESS);
  fnorm0 = fnorm;

  /* Newton iteration loop */

  for (mnewt = 0; mnewt < maxnit; mnewt++) {

    nni++;
    delnorm = fnorm;
    oldfnrm = fnorm;

    /* Call the Linesearch function and return if it failed. */
    retval = IDALineSrch(IDA_mem, &delnorm, &fnorm);
    if (retval != SUCCESS) return(retval);

    /* Set the observed convergence rate and test for convergence. */
    rate = fnorm/oldfnrm;
    if (fnorm <= epsNewt) return(SUCCESS);

    /* If not converged, copy new step vector, and loop. */
    N_VScale (ONE, delnew, delta);

  }   /* End of Newton iteration loop */

  /* Return either IC_SLOW_CONVRG or recoverable fail flag. */
  if (rate <= ICRATEMAX || fnorm < PT1*fnorm0) return(IC_SLOW_CONVRG);
  return(IC_CONV_FAIL);

}


/******************** IDALineSrch *******************************

 IDALineSrch performs the Linesearch algorithm with the 
 calculation of consistent initial conditions.

 On entry, y0 and yp0 are the current values of y and y', the 
 Newton step is delta, the current residual vector F is savres,
 delnorm is WRMS-norm(delta), and fnorm is the norm of the vector
 J-inverse F.

 On a successful return, y0, yp0, and savres have been updated, 
 delnew contains the current value of J-inverse F, and fnorm is
 WRMS-norm(delnew).
 
 The return value is SUCCESS = 0 if no error occurred.
 The error return values (positive) considered recoverable are:
    IC_FAIL_RECOV      if res or lsolve failed recoverably
    IC_CONSTR_FAILED   if the constraints could not be met
    IC_LINESRCH_FAILED if the linesearch failed (on steptol test)
 The error return values (negative) considered non-recoverable are:
    RES_NONRECOV_ERR   if res had a non-recoverable error
    SOLVE_FAILURE      if lsolve had a non-recoverable error
 
*****************************************************************/

static int IDALineSrch (IDAMem IDA_mem, realtype *delnorm, realtype *fnorm)
{
  booleantype conOK;
  int retval;
  realtype f1norm, fnormp, f1normp, ratio, lambda, minlam, slpi;
  N_Vector mc;

  /* Initialize work space pointers, f1norm, ratio.
     (Use of mc in constraint check does not conflict with ypnew.) */
  mc = ee;
  dtemp = phi[3];
  ynew = tempv2;
  ypnew = ee;
  f1norm = (*fnorm)*(*fnorm)*HALF;
  ratio = ONE;

  /* If there are constraints, check and reduce step if necessary. */
  if (constraintsSet) {

    /* Update y and check constraints. */
    IDANewy(IDA_mem);
    conOK = N_VConstrMask (constraints, ynew, mc);

    if (!conOK) {
      /* Not satisfied.  Compute scaled step to satisfy constraints. */
      N_VProd (mc, delta, dtemp);
      ratio = PT99*N_VMinQuotient (y0, dtemp);
      (*delnorm) *= ratio;
      if ((*delnorm) <= steptol) return(IC_CONSTR_FAILED);
      N_VScale (ratio, delta, delta);
    }

  } /* End of constraints check */

  slpi = -TWO*f1norm*ratio;
  minlam = steptol/(*delnorm);
  lambda = ONE;

  /* In CALC_Y_INIT case, set ypnew = yp0 (fixed) for linesearch. */
  if (icopt == CALC_Y_INIT) N_VScale (ONE, yp0, ypnew);

  /* Loop on linesearch variable lambda. */

  loop {

    /* Get new (y,y') = (ynew,ypnew) and norm of new function value. */
    IDANewyyp(IDA_mem, lambda);
    retval = IDAfnorm(IDA_mem, &fnormp);
    if (retval != SUCCESS) return(retval);

    /* If lsoff option is on, break out. */
    if (lsoff == 1) break;

    /* Do alpha-condition test. */
    f1normp = fnormp*fnormp*HALF;
    if (f1normp <= f1norm + ALPHALS*slpi*lambda) break;
    if (lambda < minlam) return(IC_LINESRCH_FAILED);
    lambda /= TWO;
    nbacktr++;

  }  /* End of breakout linesearch loop */

  /* Update y0, yp0, and fnorm, then return. */
  N_VScale (ONE, ynew,  y0);
  if (icopt == CALC_YA_YDP_INIT) N_VScale (ONE, ypnew, yp0);
  *fnorm = fnormp;
  return(SUCCESS);

}

/******************** IDAfnorm **********************************

 IDAfnorm computes the norm of the current function value, by
 evaluating the DAE residual function, calling the linear 
 system solver, and computing a WRMS-norm.
 
 On return, savres contains the current residual vector F, and
 delnew contains J-inverse F.

 The return value is SUCCESS = 0 if no error occurred, or
    IC_FAIL_RECOV    if res or lsolve failed recoverably, or
    RES_NONRECOV_ERR if res had a non-recoverable error, or
    SOLVE_FAILURE    if lsolve had a non-recoverable error.
 
*****************************************************************/

static int IDAfnorm (IDAMem IDA_mem, realtype *fnorm)
{

  int retval;

  /* Get residual vector F, return if failed, and save F in savres. */
  retval = res (Neq, tn, ynew, ypnew, delnew, rdata);
  nre++;
  if(retval < 0) return(RES_NONRECOV_ERR);
  if(retval > 0) return(IC_FAIL_RECOV);

  N_VScale (ONE, delnew, savres);

  /* Call the linear solve function to get J-inverse F; return if failed. */
  retval = lsolve (IDA_mem, delnew, ynew, ypnew, savres);
  if(retval < 0) return(SOLVE_FAILURE);
  if(retval > 0) return(IC_FAIL_RECOV);

  /* Compute the WRMS-norm; rescale if index = 0. */
  *fnorm = N_VWrmsNorm (delnew, ewt);
  if (sysindex == 0) (*fnorm) *= tscale*abs(cj);

  return(SUCCESS);

}

/******************** IDANewyyp *********************************

 IDANewyyp updates the vectors ynew and ypnew from y0 and yp0,
 using the current step vector lambda*delta, in a manner
 depending on icopt and the input id vector.
 
 The return value is always SUCCESS = 0.
 
*****************************************************************/

static int IDANewyyp (IDAMem IDA_mem, realtype lambda)
{
  
  /* CALC_YA_YDP_INIT case: ynew = y0  - lambda*delta    where id_i = 0
                           ypnew = yp0 - cj*lambda*delta where id_i = 1. */
  if (icopt == CALC_YA_YDP_INIT) {
    N_VProd (id, delta, dtemp);
    N_VLinearSum (ONE, yp0, -cj*lambda, dtemp, ypnew);
    N_VLinearSum (ONE, delta, -ONE, dtemp, dtemp);
    N_VLinearSum (ONE, y0, -lambda, dtemp, ynew);
    return(SUCCESS);
  }

  /* CALC_Y_INIT case: ynew = y0 - lambda*delta. (ypnew = yp0 preset.) */
  N_VLinearSum (ONE, y0, -lambda, delta, ynew);
  return(SUCCESS);

}


/******************** IDANewy ***********************************

 IDANewy updates the vector ynew from y0,
 using the current step vector delta, in a manner
 depending on icopt and the input id vector.
 
 The return value is always SUCCESS = 0.
 
*****************************************************************/

static int IDANewy (IDAMem IDA_mem)
{
  
  /* CALC_YA_YDP_INIT case: ynew = y0  - delta    where id_i = 0. */
  if (icopt == CALC_YA_YDP_INIT) {
    N_VProd (id, delta, dtemp);
    N_VLinearSum (ONE, delta, -ONE, dtemp, dtemp);
    N_VLinearSum (ONE, y0, -ONE, dtemp, ynew);
    return(SUCCESS);
  }

  /* CALC_Y_INIT case: ynew = y0 - delta. */
  N_VLinearSum (ONE, y0, -ONE, delta, ynew);
  return(SUCCESS);

}


/******************** IDAICFailFlag *****************************

 IDAICFailFlag prints a message and sets the IDACalcIC return
 value appropriate to the flag retval returned by IDAnlsIC.
 
*****************************************************************/

static int IDAICFailFlag (IDAMem IDA_mem, int retval)
{

  /* Depending on retval, print error message and return error flag. */
  switch (retval) {

    case RES_NONRECOV_ERR:  fprintf(errfp, MSG_IC_RES_NONREC);
                         return(RES_NONRECOV_ERR);

    case IC_FIRST_RES_FAIL:  fprintf(errfp, MSG_IC_RES_FAIL);
                         return(IC_FIRST_RES_FAIL);

    case SETUP_FAILURE:  fprintf(errfp, MSG_IC_SETUP_FAIL);
                         return(SETUP_FAILURE);

    case SOLVE_FAILURE:  fprintf(errfp, MSG_IC_SOLVE_FAIL);
                         return(SOLVE_FAILURE);

    case IC_FAIL_RECOV:  fprintf(errfp, MSG_IC_NO_RECOVERY);
                         return(IC_NO_RECOVERY);

    case IC_CONSTR_FAILED: fprintf(errfp, MSG_IC_FAIL_CONSTR);
                         return(IC_FAILED_CONSTR);

    case IC_LINESRCH_FAILED:  fprintf(errfp, MSG_IC_FAILED_LINS);
                         return(IC_FAILED_LINESRCH);

    case IC_CONV_FAIL:   fprintf(errfp, MSG_IC_CONV_FAILED);
                         return(IC_CONV_FAILURE);

    case IC_SLOW_CONVRG: fprintf(errfp, MSG_IC_CONV_FAILED);
                         return(IC_CONV_FAILURE);

    case IC_BAD_EWT:     fprintf(errfp, MSG_IC_BAD_EWT);
                         return(IC_BAD_EWT);

  }
  return -99;
}


/*********************** IDAEwtSet **************************************
  
 This routine is responsible for loading the error weight vector
 ewt, according to itol, as follows:
 (1) ewt[i] = 1 / (*rtol * ABS(ycur[i]) + *atol), i=0,...,Neq-1
     if itol = SS
 (2) ewt[i] = 1 / (*rtol * ABS(ycur[i]) + atol[i]), i=0,...,Neq-1
     if itol = SV

  It also loads the masked error weight vector mskewt if the 
  suppressalg option is on.

  IDAEwtSet returns TRUE if ewt is successfully set as above to a
  positive vector and FALSE otherwise. In the latter case, ewt is
  considered undefined after the FALSE return from IDAEwtSet.

  All the real work is done in the routines IDAEwtSetSS, IDAEwtSetSV,
  N_VOneMask, N_VProd.
 
***********************************************************************/

static booleantype IDAEwtSet(IDAMem IDA_mem, N_Vector ycur)
{
  booleantype ewtsetOK;

  switch(itol) {
  case SS: 
    ewtsetOK = IDAEwtSetSS(IDA_mem, ycur); 
    break;
  case SV: 
    ewtsetOK = IDAEwtSetSV(IDA_mem, ycur); 
    break;
  }
  if(!ewtsetOK) return(FALSE);
  if(suppressalg) {
    N_VOneMask(mskewt);  
    /* Note: mskewt is identical to id, so the vector id
       of the next line is what was just processed by N_VOneMask. It is called
       id below to emphasize that it now has been restored to its original
       settings (1.0 or 0.0) as received (id) in preparation to recreating 
       vector mskewt */
    N_VProd(ewt, id, mskewt);
   }
  return(TRUE);
}


/*********************** IDAEwtSet0 *************************************
  
 This routine is responsible for loading the error weight vector
 ewt, according to itol, as follows:
 (1) ewt[i] = 1 / (*rtol * ABS(ycur[i]) + *atol), i=0,...,Neq-1
     if itol = SS
 (2) ewt[i] = 1 / (*rtol * ABS(ycur[i]) + atol[i]), i=0,...,Neq-1
     if itol = SV

  IDAEwtSet0 returns TRUE if ewt is successfully set as above to a
  positive vector and FALSE otherwise. In the latter case, ewt is
  considered undefined after the FALSE return from IDAEwtSet0.

  All the real work is done in the routines IDAEwtSetSS, IDAEwtSetSV.
 
***********************************************************************/

static booleantype IDAEwtSet0(IDAMem IDA_mem, N_Vector ycur)
{
  booleantype ewtsetOK;

  switch(itol) {
  case SS: 
    ewtsetOK = IDAEwtSetSS(IDA_mem, ycur); 
    break;
  case SV: 
    ewtsetOK = IDAEwtSetSV(IDA_mem, ycur); 
    break;
  }
  return(ewtsetOK);
}


/*********************** IDAEwtSetSS *********************************

 This routine sets ewt as decribed above in the case itol=SS.
 It tests for non-positive components before inverting. IDAEwtSetSS
 returns TRUE if ewt is successfully set to a positive vector
 and FALSE otherwise. In the latter case, ewt is considered
 undefined after the FALSE return from IDAEwtSetSS.

********************************************************************/

static booleantype IDAEwtSetSS(IDAMem IDA_mem, N_Vector ycur)
{
  realtype rtoli, *atoli;
  
  rtoli = *rtol;
  atoli = (realtype *)atol;
  N_VAbs(ycur, tempv1);
  N_VScale(rtoli, tempv1, tempv1);
  N_VAddConst(tempv1, *atoli, tempv1);
  if (N_VMin(tempv1) <= ZERO) return(FALSE);
  N_VInv(tempv1, ewt);
  return(TRUE);
}

/*********************** IDAEwtSetSV *********************************

 This routine sets ewt as decribed above in the case itol=SV.
 It tests for non-positive components before inverting. IDAEwtSetSV
 returns TRUE if ewt is successfully set to a positive vector
 and FALSE otherwise. In the latter case, ewt is considered
 undefined after the FALSE return from IDAEwtSetSV.

********************************************************************/

static booleantype IDAEwtSetSV(IDAMem IDA_mem, N_Vector ycur)
{
  realtype rtoli;
  N_Vector atoli;
  
  rtoli = *rtol;
  atoli = (N_Vector)atol;
  N_VAbs(ycur, tempv1);
  N_VLinearSum(rtoli, tempv1, ONE, atoli, tempv1);
  if (N_VMin(tempv1) <= ZERO) return(FALSE);
  N_VInv(tempv1, ewt);
  return(TRUE);
}



/********************* IDAStopTest1 ********************************

 This routine tests for stop conditions before taking a step.
 The tests depend on the value of itask.
 The variable tretp is the previously returned value of tret.

 The return values are:
   CONTINUE_STEPS       if no stop conditions were found
   NORMAL_RETURN        for a normal return to the user
   INTERMEDIATE_RETURN  for an intermediate-output return to the user
   TSTOP_RETURN         for a tstop-reached return to the user
   ILL_INPUT            for an illegal-input return to the user 

 In the tstop cases, this routine may adjust the stepsize hh to cause
 the next step to reach tstop exactly.

********************************************************************/

static int IDAStopTest1(IDAMem IDA_mem, realtype tout, realtype tstop, realtype *tret, 
             N_Vector yret, N_Vector ypret, int itask)
{

  int ier;
  realtype troundoff;

  switch (itask) {
    
  case NORMAL:  
    /* Test for tout = tretp, and for tn past tout. */
    if (tout == tretp) {
      *tret = tretp = tout;
      return(NORMAL_RETURN);
    }
    if ( (tn - tout)*hh >= ZERO) {
      ier = IDAInterp(IDA_mem, tout, yret, ypret);
      if (ier != OKAY) {
        fprintf(errfp,MSG_BAD_TOUT, tout);
        return(ILL_INPUT);
      }
      *tret = tretp = tout;
      return(NORMAL_RETURN);
    }
    return(CONTINUE_STEPS);
    
  case ONE_STEP:
    /* Test for tn past tretp. */
    if ( (tn - tretp)*hh > ZERO) {
      ier = IDAInterp(IDA_mem, tn, yret, ypret);
      *tret = tretp = tn;
      return(INTERMEDIATE_RETURN);
    }
    return(CONTINUE_STEPS);
    
  case NORMAL_TSTOP:
    /* Test for tn past tstop, tn = tretp, tn past tout, tn near tstop. */
    if ( (tn - tstop)*hh > ZERO) {
      fprintf(errfp, MSG_BAD_TSTOP, tstop, tn);
      return(ILL_INPUT);
    }
    if (tout == tretp) {
      *tret = tretp = tout;
      return(NORMAL_RETURN);
    }
    if ( (tn - tout)*hh >= ZERO) {
      ier = IDAInterp(IDA_mem, tout, yret, ypret);
      if (ier != OKAY) {
        fprintf(errfp, MSG_BAD_TOUT, tout);
        return(ILL_INPUT);
      }
      *tret = tretp = tout;
      return(NORMAL_RETURN);
    }
    troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
    if ( ABS(tn - tstop) <= troundoff) {
      ier = IDAInterp(IDA_mem, tstop, yret, ypret);
      if (ier != OKAY) {
        fprintf(errfp, MSG_BAD_TSTOP, tstop, tn);
        return(ILL_INPUT);
      }
      *tret = tretp = tstop;
      return(TSTOP_RETURN);
    }
    if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
    return(CONTINUE_STEPS);
    
  case ONE_STEP_TSTOP:
    /* Test for tn past tstop, tn past tretp, and tn near tstop. */
    if ( (tn - tstop)*hh > ZERO) {
      fprintf(errfp, MSG_BAD_TSTOP, tstop, tn);
      return(ILL_INPUT);
    }
    if ( (tn - tretp)*hh > ZERO) {
      ier = IDAInterp(IDA_mem, tn, yret, ypret);
      *tret = tretp = tn;
      return(INTERMEDIATE_RETURN);
    }
    troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
    if ( ABS(tn - tstop) <= troundoff) {
      ier = IDAInterp(IDA_mem, tstop, yret, ypret);
      if (ier != OKAY) {
        fprintf(errfp, MSG_BAD_TSTOP, tstop, tn);
        return(ILL_INPUT);
      }
      *tret = tretp = tstop;
      return(TSTOP_RETURN);
    }
    if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
    return(CONTINUE_STEPS);
    
  }
  return -99;
}

/********************* IDAStopTest2 ********************************

 This routine tests for stop conditions after taking a step.
 The tests depend on the value of itask.

 The return values are:
   CONTINUE_STEPS       if no stop conditions were found
   NORMAL_RETURN        for a normal return to the user
   INTERMEDIATE_RETURN  for an intermediate-output return to the user
   TSTOP_RETURN         for a tstop-reached return to the user

 In the two cases with tstop, this routine may reset the stepsize hh
 to cause the next step to reach tstop exactly.

 In the two cases with ONE_STEP mode, no interpolation to tn is needed
 because yret and ypret already contain the current y and y' values.

 Note: No test is made for an error return from IDAInterp here,
 because the same test was made prior to the step.

********************************************************************/

static int IDAStopTest2(IDAMem IDA_mem, realtype tout, realtype tstop, realtype *tret, 
             N_Vector yret, N_Vector ypret, int itask)
{

  int ier;
  realtype troundoff;

  switch (itask) {

    case NORMAL:  
      /* Test for tn past tout. */
      if ( (tn - tout)*hh >= ZERO) {
        ier = IDAInterp(IDA_mem, tout, yret, ypret);
        *tret = tretp = tout;
        return(NORMAL_RETURN);
      }
      return(CONTINUE_STEPS);

    case ONE_STEP:
      *tret = tretp = tn;
      return(INTERMEDIATE_RETURN);

    case NORMAL_TSTOP:
      /* Test for tn at tstop, for tn past tout, and for tn near tstop. */
      troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
      if ( ABS(tn - tstop) <= troundoff) {
        ier = IDAInterp(IDA_mem, tstop, yret, ypret);
        *tret = tretp = tstop;
        return(TSTOP_RETURN);
      }
      if ( (tn - tout)*hh >= ZERO) {
        ier = IDAInterp(IDA_mem, tout, yret, ypret);
        *tret = tretp = tout;
        return(NORMAL_RETURN);
      }
      if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
      return(CONTINUE_STEPS);

    case ONE_STEP_TSTOP:
      /* Test for tn at tstop. */
      troundoff = HUNDRED*uround*(ABS(tn) + ABS(hh));
      if ( ABS(tn - tstop) <= troundoff) {
        ier = IDAInterp(IDA_mem, tstop, yret, ypret);
        *tret = tretp = tstop;
        return(TSTOP_RETURN);
      }
      if ( (tn + hh - tstop)*hh > ZERO) hh = tstop - tn;
      *tret = tretp = tn;
      return(INTERMEDIATE_RETURN);

  }
  return -99;
}


/*************** IDAInterp ********************************************

 This routine evaluates y(t) and y'(t) as the value and derivative of 
 the interpolating polynomial at the independent variable t, and stores
 the results in the vectors yret and ypret.  It uses the current
 independent variable value, tn, and the method order last used, kused.
 This function is called by IDASolve with t = tout, t = tn, or t = tstop.

 If kused = 0 (no step has been taken), or if t = tn, then the order used
 here is taken to be 1, giving yret = phi[0], ypret = phi[1]/psi[0].

 The return values are:
   OKAY  if t is legal, or
   BAD_T if t is not within the interval of the last step taken.

**********************************************************************/

static int IDAInterp(IDAMem IDA_mem, realtype t, N_Vector yret, N_Vector ypret)
{
  realtype tfuzz, tp, delt, c, d, gam;
  int j, kord;
  
  /* Check t for legality.  Here tn - hused is t_{n-1}. */
 
  tfuzz = HUNDRED * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  if ( (t - tp)*hh < ZERO) return(BAD_T);

  /* Initialize yret = phi[0], ypret = 0, and kord = (kused or 1). */

  N_VScale (ONE, phi[0], yret);
  N_VConst (ZERO, ypret);
  kord = kused; if (kused == 0 || t == tn) kord = 1;

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
  return(OKAY);
}


/****************** IDAHandleFailure ******************************

 This routine prints error messages for all cases of failure by
 IDAStep.  It returns to IDASolve the value that it is to return to
 the user.

*****************************************************************/

static int IDAHandleFailure(IDAMem IDA_mem, int sflag)
{

  /* Depending on sflag, print error message and return error flag */
  switch (sflag) {

    case REP_ERR_FAIL:  fprintf(errfp, MSG_ERR_FAILS, tn, hh);
                        return(ERR_FAILURE);

    case REP_CONV_FAIL: fprintf(errfp, MSG_CONV_FAILS, tn, hh);
                        return(CONV_FAILURE);

    case LSETUP_ERROR_NONRECVR:  fprintf(errfp, MSG_SETUP_FAILED, tn);
                        return(SETUP_FAILURE);

    case CONV_FAIL_LINR_NONRECVR:  fprintf(errfp, MSG_SOLVE_FAILED, tn);
                        return(SOLVE_FAILURE);

    case REP_RES_ERR:   fprintf(errfp, MSG_REP_RES_ERR, tn);
                        return(REP_RES_REC_ERR);

    case RES_ERROR_NONRECVR:  fprintf(errfp, MSG_RES_NONRECOV, tn);
                        return(RES_NONRECOV_ERR);

    case CONSTR_FAIL:   fprintf(errfp, MSG_FAILED_CONSTR, tn);
                        return(CONSTR_FAILURE);

  }
  return -99;
}

/********************* IDAStep **************************************
 
 This routine performs one internal IDA step, from tn to tn + hh.
 It calls other routines to do all the work.

 It solves a system of differential/algebraic equations of the form
       F(t,y,y') = 0, for one step. In IDA, tt is used for t,
 yy is used for y, and yp is used for y'. The function F is supplied as 'res'
 by the user.

 The methods used are modified divided difference, fixed leading 
 coefficient forms of backward differentiation formulas.
 The code adjusts the stepsize and order to control the local error per step.

 The main operations done here are as follows:
  * initialize various quantities;
  * setting of multistep method coefficients;
  * solution of the nonlinear system for yy at t = tn + hh;
  * deciding on order reduction and testing the local error;
  * attempting to recover from failure in nonlinear solver or error test;
  * resetting stepsize and order for the next step.
  * updating phi and other state data if successful;

 On a failure in the nonlinear system solution or error test, the
 step may be reattempted, depending on the nature of the failure.

 Variables or arrays (all in the IDAMem structure) used in IDAStep are:

 tt -- Independent variable.
 yy -- Solution vector at tt.
 yp -- Derivative of solution vector after successful stelp.
 Neq -- Number of equations to be integrated.
 res -- User-supplied function to evaluate the residual. See the 
        description given in file ida.h .
 lsetup -- Routine to prepare for the linear solver call. It may either
        save or recalculate quantities used by lsolve. (Optional)
 lsolve -- Routine to solve a linear system. A prior call to lsetup
        may be required. 
 hh  -- Appropriate step size for next step.
 ewt -- Vector of weights used in all convergence tests.
 mskewt -- Masked vector of weights used in error test.
 phi -- Array of divided differences used by IDAStep. This array is composed 
       of  (maxord+1) nvectors (each of size Neq). (maxord+1) is the maximum 
       order for the problem, maxord, plus 1.
 
       Return values are:
       SUCCESS       RES_ERROR_NONRECVR        LSETUP_ERROR_NONRECVR       
                     CONV_FAIL_LINR_NONRECVR   REP_ERR_FAIL            
                     CONSTR_FAIL               REP_CONV_FAIL          
                     REP_RES_ERR            

********************************************************************/

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
    kflag = SUCCESS;

    nflag = IDAnls(IDA_mem);

    if(nflag == SUCCESS) nflag=IDATestError(IDA_mem, &ck, &est,
                                             &terk, &terkm1, &erkm1);

    if(nflag != SUCCESS) kflag=IDAHandleNFlag(IDA_mem, nflag, saved_t, 
                                                &ncf, &nef, &est);
    if (kflag == PREDICT_AGAIN) continue;
    else if(kflag == SUCCESS) break;
    else return(kflag);
  }

  /* Nonlinear system solve and error test were both successful;
     update data, and consider change of step and/or order       */

  IDACompleteStep(IDA_mem, &est, &terk, &terkm1, &erkm1);

  return(SUCCESS);
}


/********************* IDASetCoeffs ********************************

  This routine computes the coefficients relevant to the current step.
  The counter ns counts the number of consecutive steps taken at
  constant stepsize h and order k, up to a maximum of k + 2.
  Then the first ns components of beta will be one, and on a step  
  with ns = k + 2, the coefficients alpha, etc. need not be reset here.
  Also, IDACompleteStep prohibits an order increase until ns = k + 2.

***********************************************************************/

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


/****************** IDAnls *****************************************

 This routine attempts to solve the nonlinear system using the linear
 solver specified. NOTE: this routine uses N_Vector ee as the scratch
 vector tempv3 passed to lsetup.

  Possible return values:

  SUCCESS

  RES_ERROR_RECVR       RES_ERROR_NONRECVR
  LSETUP_ERROR_RECVR    LSETUP_ERROR_NONRECVR
  CONV_FAIL_LINR_RECVR  CONV_FAIL_LINR_NONRECVR

  CONSTRAINT_FAIL_RECVR
  CONV_FAIL_NLINR_RECVR


**********************************************************************/

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

    retval = res(Neq, tn, yy, yp, delta, rdata);
    nre++;
    if(retval != SUCCESS) break;

    /* If indicated, call linear solver setup function and reset parameters. */
    if(callSetup){
      nsetups++;

      retval = lsetup(IDA_mem, yy, yp, delta, tempv1, tempv2, tempv3);

      cjold = cj;
      cjratio = ONE;
      ss = TWENTY;
      if(retval != SUCCESS) break;
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

  if(retval != SUCCESS) return(retval);


  /* If otherwise successful, check and enforce inequality constraints. */

  if(constraintsSet){  /* Check constraints and get mask vector mm, 
                          set where constraints failed */
    constraintsPassed = N_VConstrMask(constraints,yy,mm);
    if(constraintsPassed) return(SUCCESS);
    else {
      N_VCompare(ONEPT5, constraints, tempv1);  
      /* a , where a[i] =1. when |c[i]| = 2 ,  c the vector of constraints */
      N_VProd(tempv1, constraints, tempv1);       /* a * c */
      N_VDiv(tempv1, ewt, tempv1);                /* a * c * wt */
      N_VLinearSum(ONE, yy, -PT1, tempv1, tempv1);/* y - 0.1 * a * c * wt */
      N_VProd(tempv1, mm, tempv1);               /*  v = mm*(y-.1*a*c*wt) */
      vnorm = N_VWrmsNorm(tempv1, ewt);          /*  ||v|| */
      
      /* If vector v of constraint corrections is small
         in norm, correct and accept this step */      
      if(vnorm <= epsNewt){  
        N_VLinearSum(ONE, yy, -ONE, tempv1, yy);  /* y <- y - v */
        return(SUCCESS);
      }
      else {
        /* Constraints not met -- reduce h by computing rr = h'/h */
        N_VLinearSum(ONE, phi[0], -ONE, yy, tempv1);
        N_VProd(mm, tempv1, tempv1);
        rr = PT9*N_VMinQuotient(phi[0], tempv1);
        rr = MAX(rr,PT1);
        return(CONSTRAINT_FAIL_RECVR);
      }
    }
  }
  return(SUCCESS);
}


/***************** IDAPredict ********************************************

  This routine predicts the new values for vectors yy and yp.

*************************************************************************/

static int IDAPredict(IDAMem IDA_mem)
{
  int j;

  N_VScale(ONE, phi[0], yy);
  N_VConst(ZERO, yp);
  
  for(j=1; j<=kk; j++) {
    N_VLinearSum(ONE,      phi[j], ONE, yy, yy);
    N_VLinearSum(gamma[j], phi[j], ONE, yp, yp);
  }

  return(SUCCESS);
}


/********************** IDANewtonIter *********************************

 This routine performs the Newton iteration.  
 It assumes that delta contains the initial residual vector on entry.
 If the iteration succeeds, it returns the value SUCCESS = 0.
 If not, it returns either:
   a positive value (for a recoverable failure), namely one of:
     RES_ERROR_RECOVR
     CONV_FAIL_LINR_RECVR
     CONV_FAIL_NLINR_RECVR
 or
   a negative value (for a nonrecoverable failure), namely one of:
     RES_ERROR_NONRECOVR
     CONV_FAIL_LINR_NONRECVR

     NOTE: This routine uses N_Vector savres, which is preset to tempv1.

***********************************************************************/

static int IDANewtonIter(IDAMem IDA_mem)
{
  int mnewt, retval;
  realtype delnrm, oldnrm, rate;

  /* Initialize counter mnewt and cumulative correction vector ee. */
  mnewt = 0;
  N_VConst (ZERO, ee);

  /* Looping point for Newton iteration.  Break out on any error. */
  loop {

    /* Save a copy of the residual vector in savres. */
    N_VScale(ONE, delta, savres);

    /* Call the lsolve function to get correction vector delta. */
    retval = lsolve(IDA_mem, delta, yy, yp, savres); 
    nni++;
    if (retval != SUCCESS) break;

    /* Apply delta to yy, yp, and ee, and get norm(delta). */
    N_VLinearSum(ONE, yy, -ONE, delta, yy);
    N_VLinearSum(ONE, ee, -ONE, delta, ee);
    N_VLinearSum(ONE, yp, -cj,  delta, yp);
    delnrm = N_VWrmsNorm(delta, ewt);

    /* Test for convergence, first directly, then with rate estimate. */

    if (mnewt == 0){ 
       oldnrm = delnrm;
       if (delnrm <= toldel) return(SUCCESS);
    }
    else {
      rate = RPowerR( delnrm/oldnrm, ONE/mnewt );
      if (rate > RATEMAX) { retval = CONV_FAIL_NLINR_RECVR; break;}
      ss = rate/(ONE - rate);
    }

    if (ss*delnrm <= epsNewt) return(SUCCESS);

    /* Not yet converged.  Increment mnewt and test for max allowed. */
    mnewt++;
    if (mnewt >= MAXIT) {retval = CONV_FAIL_NLINR_RECVR; break;}

    /* Call res for new residual and check error flag from res. */
    retval = res(Neq, tn, yy, yp, delta, rdata);
    nre++;
    if (retval != SUCCESS) break;

    /* Loop for next iteration. */

  } /* end of Newton iteration loop */

  /* All error returns exit here. */
  return(retval);

}


/******************* IDATestError ********************************

 This routine estimates errors at orders k, k-1, k-2, decides whether 
 or not to reduce order, and performs the local error test. 

 IDATestError returns either  SUCCESS   or ERROR_TEST_FAIL

******************************************************************/

static int IDATestError(IDAMem IDA_mem, realtype *ck, realtype *est,
                         realtype *terk, realtype *terkm1, realtype *erkm1)
{
  int retval;
  realtype enorm;
  realtype terkm2;
  realtype erk, erkm2;

  /* Compute error for order k. */
  enorm = N_VWrmsNorm(ee, mskewt);
  erk = sigma[kk] * enorm;
  *terk = (kk+1) * erk;
  *est = erk;
  knew = kk;

  /* Now compute the errors for orders k-1 and k-2, and decide whether to 
     reduce the order k to k-1 */
  
  if(kk > 1){
    N_VLinearSum(ONE, phi[kk], ONE, ee, delta);
    *erkm1 = sigma[kk-1] * N_VWrmsNorm(delta, mskewt);
    *terkm1 = kk * *erkm1;
    {
      if(kk > 2){
        N_VLinearSum(ONE, phi[kk-1], ONE, delta, delta);
        erkm2 = sigma[kk-2] * N_VWrmsNorm(delta, mskewt);
        terkm2 = (kk-1) * erkm2;
        if(MAX(*terkm1, terkm2) > *terk) goto evaltest;
      }
      
      else if(*terkm1 > 0.5 * (*terk)) goto evaltest; /* executed for kk=2 only */
    }
    /* end of "kk>2" if/else block */
    
    knew = kk-1;
    *est = *erkm1;
    
  } /* end kk>1 if block */ 
  
  
 evaltest:
  retval = SUCCESS;
  
  if ((*ck * enorm) > ONE) retval = ERROR_TEST_FAIL;
  return(retval);
}


/********************** IDAHandleNFlag *******************************

This routine handles failures indicated by the input variable nflag. 
Positive values indicate various recoverable failures while negative
values indicate nonrecoverable failures. This routine adjusts the
step size for recoverable failures. The possible return values are:

(recoverable)
PREDICT_AGAIN

(nonrecoverable)
CONSTR_FAIL      REP_RES_ERR        RES_ERROR_NONRECVR      [from IDAnls] 
REP_ERR_FAIL     REP_CONV_FAIL      LSETUP_ERROR_NONRECVR   [from IDAnls]
                                    CONV_FAIL_LINR_NONRECVR [from IDAnls]
     
**********************************************************************/

static int IDAHandleNFlag(IDAMem IDA_mem, int nflag, realtype saved_t,
                          int *ncfPtr, int *nefPtr, realtype *est)
{
  int j, retval;
  int *ncf, *nef;
  
  ncf = ncfPtr; nef = nefPtr;
  phase = 1;
  
  tn = saved_t; /* restore tn */
  
  /* restore phi and psi */
  
  for (j = ns; j<=kk; j++) N_VScale(ONE/beta[j], phi[j], phi[j]);
  
  for (j = 1; j <= kk; j++) psi[j-1] = psi[j] - hh;
  
  loop{  /* begin 'breakout' loop */
    
    if (nflag < 0) {    /*  nonrecoverable failure cases */
      retval = nflag; ncfn++;
      break;
    }
    
    /* Only positive nflag values (apparently recoverable) will appear here*/
    
    else if (nflag != ERROR_TEST_FAIL) {   /*  misc. recoverable cases  */
      (*ncf)++; ncfn++;
      if (nflag != CONSTRAINT_FAIL_RECVR) rr = QUARTER;
      hh *= rr;
      if (*ncf < MXNCF){
        retval = PREDICT_AGAIN;
        break;
      }
      else if (nflag == RES_ERROR_RECVR) {
        retval= REP_RES_ERR;
        break;
      }
      else if (nflag == CONSTRAINT_FAIL_RECVR) {
        retval = CONSTR_FAIL;
        break;
      }
      else {
        retval = REP_CONV_FAIL;
        break;
      }
    }
    else {    /* error test failed case */
      (*nef)++; netf++;
      
      if (*nef == 1){
        /* On first error test failure, keep current order or lower order 
           by one. Compute new stepsize based on differences of the solution. */
        kk = knew;
        
        rr = PT9 * RPowerR( TWO*(*est) + PT0001,(-ONE/(kk+1)) );
        rr = MAX(QUARTER, MIN(PT9,rr));
        hh *=rr;  /* adjust step size */
        
        retval = PREDICT_AGAIN;
        break;
      }
      else if (*nef == 2){
        /* On second error test failure, use current order or decrease order 
           by one. Reduce stepsize by factor of 1/4. */
        kk = knew;
        rr = QUARTER;
        hh *= rr;
        
        retval = PREDICT_AGAIN;
        break;
      }
      else if (*nef > 2){
        /* On third and subsequent error test failures, set order to 1 and
           reduce stepsize h by factor of 1/4. */
        if (*nef<MXETF){
          kk = 1;
          rr = QUARTER;
          hh *= rr;
          retval = PREDICT_AGAIN;
          break;
        }
        else {
          retval = REP_ERR_FAIL;
          break;
        }
      }
    } /* end of nflag  if block */
    
  } /* end of 'breakout loop' */
  
  if (retval < 0) return(retval);
  
  if (retval == PREDICT_AGAIN) {
    if (nst == 0){
      psi[0] = hh;
      N_VScale(rr, phi[1], phi[1]);
    }
    return(PREDICT_AGAIN);
  }
  return -99;  
}


/********************* IDACompleteStep *********************************

 This routine completes a successful step.  It increments nst,
 saves the stepsize and order used, makes the final selection of
 stepsize and order for the next step, and updates the phi array.
 Its return value is SUCCESS= 0.

***********************************************************************/

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
  stepsize and order are set by the usual local error algorithm.         */
  
  if (phase  == 0) {
    kk++;
    hnew = TWO * hh;
    hh = hnew;
  }
  
  else {
    action = UNSET;
    
    /* Set action = LOWER/MAINTAIN/RAISE to specify order decision */
    
    if (knew == kk-1)   {action = LOWER; goto takeaction;}
    if (kk == maxord)   {action = MAINTAIN; goto takeaction;}
    if ( (kk+1 >= ns ) || (kdiff == 1)) {action = MAINTAIN;goto takeaction;}
    
    /* Estimate the error at order k+1, unless already decided to
       reduce order, or already using maximum order, or stepsize has not
       been constant, or order was just raised. */
    
    N_VLinearSum (ONE, ee, -ONE, phi[kk+1], delta);
    terkp1 = N_VWrmsNorm(delta, mskewt);
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
  
  return (SUCCESS);
}


/********************************************************************/
/********* END Private Helper Functions Implementation **************/
/********************************************************************/


/********************************************************************/
/************** END IDA Implementation ******************************/
/********************************************************************/
