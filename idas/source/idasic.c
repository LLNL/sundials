/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2004-06-18 21:37:59 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the implementation file for the IC calculation for IDAS.
 * It is independent of the linear solver in use.                  
 * -----------------------------------------------------------------
 */

/*=================================================================*/
/*BEGIN        Import Header Files                                 */
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "idas_impl.h"
#include "sundialsmath.h"

/*=================================================================*/
/*END          Import Header Files                                 */
/*=================================================================*/

/* Macro: loop */
#define loop for(;;)

/*=================================================================*/
/*BEGIN        IDAS Constants                                      */
/*=================================================================*/

/* Private Constants */

#define ZERO       RCONST(0.0)    /* real 0.0    */
#define HALF       RCONST(0.5)    /* real 0.5    */
#define ONE        RCONST(1.0)    /* real 1.0    */
#define TWO        RCONST(2.0)    /* real 2.0    */
#define PT99       RCONST(0.99)   /* real 0.99   */
#define PT1        RCONST(0.1)    /* real 0.1    */
#define PT001      RCONST(0.001)  /* real 0.001  */

/* IDACalcIC control constants */

#define ICRATEMAX  RCONST(0.9)    /* max. Newton conv. rate */
#define ALPHALS    RCONST(0.0001) /* alpha in linesearch conv. test */

/* Return values for lower level routines used by IDACalcIC */

enum { IC_FAIL_RECOV = 1,  IC_CONSTR_FAILED = 2,  IC_LINESRCH_FAILED = 3,
       IC_CONV_FAIL =  4,  IC_SLOW_CONVRG =   5 };

/*=================================================================*/
/*END          IDAS Constants                                      */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        IDAS Error Messages                                 */
/*=================================================================*/

/* IDACalcIC error messages */

#define IDAIC              "IDACalcIC-- "

#define MSG_IC_IDA_NO_MEM  IDAIC "IDA_mem = NULL illegal.\n\n"

#define MSG_IC_NO_MALLOC   IDAIC "Attempt to call before IDAMalloc. \n\n"
 
#define MSG_BAD_ICOPT      IDAIC "icopt = %d is illegal.\n\n"

#define MSG_IC_MISSING_ID  IDAIC "id = NULL conflicts with icopt.\n\n"

#define MSG_IC_BAD_ID      IDAIC "id has illegal values.\n\n"

#define MSG_IC_TOO_CLOSE1  IDAIC "tout1 = %g too close to t0 = %g to attempt"
#define MSG_IC_TOO_CLOSE2  " initial condition calculation.\n\n"
#define MSG_IC_TOO_CLOSE   MSG_IC_TOO_CLOSE1 MSG_IC_TOO_CLOSE2

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

/*=================================================================*/
/*END          IDAS Error Messages                                 */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        Private Helper Functions Prototypes                 */
/*=================================================================*/

extern int IDAInitialSetup(IDAMem IDA_mem);
extern booleantype IDAEwtSet(IDAMem IDA_mem, N_Vector ycur);
extern realtype IDAWrmsNorm(IDAMem IDA_mem, N_Vector x, N_Vector w, 
                            booleantype mask);

static int IDAnlsIC (IDAMem IDA_mem);
static int IDANewtonIC (IDAMem IDA_mem);
static int IDALineSrch (IDAMem IDA_mem, realtype *delnorm, realtype *fnorm);
static int IDAfnorm (IDAMem IDA_mem, realtype *fnorm);
static int IDANewyyp (IDAMem IDA_mem, realtype lambda);
static int IDANewy (IDAMem IDA_mem);
static int IDAICFailFlag (IDAMem IDA_mem, int retval);

/*=================================================================*/
/*END          Private Helper Functions Prototypes                 */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        Readibility Constants                               */
/*=================================================================*/

#define errfp    (IDA_mem->ida_errfp)
#define rdata    (IDA_mem->ida_rdata)
#define res      (IDA_mem->ida_res)
#define y0       (IDA_mem->ida_y0)
#define yp0      (IDA_mem->ida_yp0)
#define uround   (IDA_mem->ida_uround)  
#define phi      (IDA_mem->ida_phi) 
#define ewt      (IDA_mem->ida_ewt)  
#define delta    (IDA_mem->ida_delta)
#define ee       (IDA_mem->ida_ee)
#define savres   (IDA_mem->ida_savres)
#define tempv2   (IDA_mem->ida_tempv2) 
#define hh       (IDA_mem->ida_hh)
#define tn       (IDA_mem->ida_tn)
#define cj       (IDA_mem->ida_cj)
#define cjratio  (IDA_mem->ida_cjratio)
#define nbacktr  (IDA_mem->ida_nbacktr)
#define nre      (IDA_mem->ida_nre)
#define ncfn     (IDA_mem->ida_ncfn)
#define nni      (IDA_mem->ida_nni)
#define nsetups  (IDA_mem->ida_nsetups)
#define ns       (IDA_mem->ida_ns)
#define linit    (IDA_mem->ida_linit)
#define lsetup   (IDA_mem->ida_lsetup)
#define lsolve   (IDA_mem->ida_lsolve) 
#define lperf    (IDA_mem->ida_lperf)
#define lfree    (IDA_mem->ida_lfree) 
#define lmem     (IDA_mem->ida_lmem) 
#define linitOK  (IDA_mem->ida_linitOK)
#define hused    (IDA_mem->ida_hused)         
#define epsNewt  (IDA_mem->ida_epsNewt)
#define id       (IDA_mem->ida_id)
#define setupNonNull   (IDA_mem->ida_setupNonNull) 
#define suppressalg    (IDA_mem->ida_suppressalg)
#define constraints    (IDA_mem->ida_constraints)
#define constraintsSet (IDA_mem->ida_constraintsSet)

#define epiccon  (IDA_mem->ida_epiccon)
#define maxnh    (IDA_mem->ida_maxnh)
#define maxnj    (IDA_mem->ida_maxnj)
#define maxnit   (IDA_mem->ida_maxnit)
#define lsoff    (IDA_mem->ida_lsoff)
#define steptol  (IDA_mem->ida_steptol)

/*=================================================================*/
/*END          Readibility Constants                               */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        EXPORTED FUNCTIONS IMPLEMENTATION                   */
/*=================================================================*/

/*-------------------- IDACalcIC ----------------------------------*/
/*
 IDACalcIC computes consistent initial conditions, given the 
 user's initial guess for unknown components of y0 and/or yp0.

 The return value is SUCCESS = 0 if no error occurred.

 The error return values (fully described in ida.h) are:
    IC_IDA_NO_MEM      ida_mem is NULL
    IC_NO_MALLOC       ida_mem was not allocated
    IC_ILL_INPUT       bad value for icopt, tout1, or id
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
*/
/*-----------------------------------------------------------------*/

int IDACalcIC (void *ida_mem, int icopt, realtype tout1)
{
  booleantype ewtsetOK;
  int ier, nwt, nh, mxnh, icret, retval=0;
  realtype tdist, troundoff, minid, hic, ypnorm;
  IDAMem IDA_mem;

  /* Check if IDA memory exists */

  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IC_IDA_NO_MEM);
    return(IC_IDA_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if problem was malloc'ed */
  
  if (IDA_mem->ida_MallocDone == FALSE) {
    if(errfp!=NULL) fprintf(errfp, MSG_IC_NO_MALLOC);
    return(IC_NO_MALLOC);
  }

  /* Check inputs to IDA for correctness and consistency */

  ier = IDAInitialSetup(IDA_mem);
  if (ier != SUCCESS) return(IC_ILL_INPUT);
  IDA_mem->ida_SetupDone = TRUE;

  /* Check legality of input arguments, and set IDA memory copies. */

  if (icopt < CALC_YA_YDP_INIT || icopt > CALC_Y_INIT) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_ICOPT, icopt);
    return(IC_ILL_INPUT);
  }
  IDA_mem->ida_icopt = icopt;

  if (icopt == CALC_YA_YDP_INIT && (id == NULL)) {
    if(errfp!=NULL) fprintf(errfp, MSG_IC_MISSING_ID);
    return(IC_ILL_INPUT);
  }

  tdist = ABS(tout1 - tn);
  troundoff = TWO*uround*(ABS(tn) + ABS(tout1));    
  if (tdist < troundoff) {
    if(errfp!=NULL) fprintf(errfp, MSG_IC_TOO_CLOSE, tout1, tn);
    return(IC_ILL_INPUT);
  }

  /* For use in the CALC_YA_YP_INIT case, set sysindex and tscale. */

  IDA_mem->ida_sysindex = 1;
  IDA_mem->ida_tscale   = tdist;
  if (icopt == CALC_YA_YDP_INIT) {
    minid = N_VMin(id);
    if (minid < ZERO) {
      if(errfp!=NULL) fprintf(errfp, MSG_IC_BAD_ID);
      return(IC_ILL_INPUT);
    }
    if (minid > HALF) IDA_mem->ida_sysindex = 0;
  }

  /* Set the test constant in the Newton convergence test */

  IDA_mem->ida_epsNewt = epiccon;

  /* Initializations: cjratio = 1 (for use in direct linear solvers); 
     set nbacktr = 0; call linit routine. */

  cjratio = ONE;
  nbacktr = 0;
  linitOK = (linit(IDA_mem) == LINIT_OK);
  if (!linitOK) {
    if(errfp!=NULL) fprintf(errfp, MSG_IC_LINIT_FAIL);
    return(IC_LINIT_FAIL);
  }

  /* Set hic, hh, cj, and mxnh. */
  hic = PT001*tdist;
  ypnorm = IDAWrmsNorm(IDA_mem, yp0, ewt, suppressalg);
  if (ypnorm > HALF/hic) hic = HALF/ypnorm;
  if( tout1 < tn) hic = -hic;
  hh = hic;
  if (icopt == CALC_YA_YDP_INIT) {
    cj = ONE/hic;
    mxnh = maxnh;
  }
  else {
    cj = ZERO;
    mxnh = 1;
  }

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
    ewtsetOK = IDAEwtSet(IDA_mem, y0);
    if (!ewtsetOK) { retval = IC_BAD_EWT; break; }
    N_VScale (ONE, y0,  phi[0]);
    N_VScale (ONE, yp0, phi[1]);

  }   /* End of nwt loop */

  /* Load the optional outputs. */
  if (icopt == CALC_YA_YDP_INIT)   hused = hic;

  /* On any failure, print message and return proper flag. */
  if (retval != SUCCESS) {
    icret = IDAICFailFlag(IDA_mem, retval);
    return(icret);
  }

  /* Otherwise return success flag. */
  return(SUCCESS);

}

/*=================================================================*/
/*END          EXPORTED FUNCTIONS IMPLEMENTATION                   */
/*=================================================================*/

/*=================================================================*/
/*BEGIN        PRIVATE FUNCTIONS IMPLEMENTATION                    */
/*=================================================================*/

#define icopt    (IDA_mem->ida_icopt)
#define sysindex (IDA_mem->ida_sysindex)
#define tscale   (IDA_mem->ida_tscale)
#define ynew     (IDA_mem->ida_ynew)
#define ypnew    (IDA_mem->ida_ypnew)
#define delnew   (IDA_mem->ida_delnew)
#define dtemp    (IDA_mem->ida_dtemp)

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

static int IDAnlsIC(IDAMem IDA_mem)
{
  int retval, nj;
  N_Vector tv1, tv2, tv3;

  tv1 = ee;
  tv2 = tempv2;
  tv3 = phi[2];

  retval = res (tn, y0, yp0, delta, rdata);
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

static int IDANewtonIC(IDAMem IDA_mem)
{
  int retval, mnewt;
  realtype delnorm, fnorm, fnorm0, oldfnrm, rate;

  /* Set pointer for vector delnew */
  delnew = phi[2];

  /* Call the linear solve function to get the Newton step, delta. */
  retval = lsolve (IDA_mem, delta, ewt, y0, yp0, savres);
  if(retval < 0) return(SOLVE_FAILURE);
  if(retval > 0) return(IC_FAIL_RECOV);

  /* Compute the norm of the step; return now if this is small. */
  fnorm = IDAWrmsNorm(IDA_mem, delta, ewt, FALSE);
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

static int IDALineSrch(IDAMem IDA_mem, realtype *delnorm, realtype *fnorm)
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
    if (lsoff) break;

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

static int IDAfnorm(IDAMem IDA_mem, realtype *fnorm)
{

  int retval;

  /* Get residual vector F, return if failed, and save F in savres. */
  retval = res (tn, ynew, ypnew, delnew, rdata);
  nre++;
  if(retval < 0) return(RES_NONRECOV_ERR);
  if(retval > 0) return(IC_FAIL_RECOV);

  N_VScale (ONE, delnew, savres);

  /* Call the linear solve function to get J-inverse F; return if failed. */
  retval = lsolve (IDA_mem, delnew, ewt, ynew, ypnew, savres);
  if(retval < 0) return(SOLVE_FAILURE);
  if(retval > 0) return(IC_FAIL_RECOV);

  /* Compute the WRMS-norm; rescale if index = 0. */
  *fnorm = IDAWrmsNorm(IDA_mem, delnew, ewt, FALSE);
  if (sysindex == 0) (*fnorm) *= tscale*abs(cj);

  return(SUCCESS);

}

/******************** IDANewyyp *********************************

 IDANewyyp updates the vectors ynew and ypnew from y0 and yp0,
 using the current step vector lambda*delta, in a manner
 depending on icopt and the input id vector.
 
 The return value is always SUCCESS = 0.
 
*****************************************************************/

static int IDANewyyp(IDAMem IDA_mem, realtype lambda)
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

static int IDANewy(IDAMem IDA_mem)
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

static int IDAICFailFlag(IDAMem IDA_mem, int retval)
{

  /* Depending on retval, print error message and return error flag. */
  switch (retval) {

    case RES_NONRECOV_ERR:  if(errfp!=NULL) fprintf(errfp, MSG_IC_RES_NONREC);
                         return(RES_NONRECOV_ERR);

    case IC_FIRST_RES_FAIL:  if(errfp!=NULL) fprintf(errfp, MSG_IC_RES_FAIL);
                         return(IC_FIRST_RES_FAIL);

    case SETUP_FAILURE:  if(errfp!=NULL) fprintf(errfp, MSG_IC_SETUP_FAIL);
                         return(SETUP_FAILURE);

    case SOLVE_FAILURE:  if(errfp!=NULL) fprintf(errfp, MSG_IC_SOLVE_FAIL);
                         return(SOLVE_FAILURE);

    case IC_FAIL_RECOV:  if(errfp!=NULL) fprintf(errfp, MSG_IC_NO_RECOVERY);
                         return(IC_NO_RECOVERY);

    case IC_CONSTR_FAILED: if(errfp!=NULL) fprintf(errfp, MSG_IC_FAIL_CONSTR);
                         return(IC_FAILED_CONSTR);

    case IC_LINESRCH_FAILED:  if(errfp!=NULL) fprintf(errfp, MSG_IC_FAILED_LINS);
                         return(IC_FAILED_LINESRCH);

    case IC_CONV_FAIL:   if(errfp!=NULL) fprintf(errfp, MSG_IC_CONV_FAILED);
                         return(IC_CONV_FAILURE);

    case IC_SLOW_CONVRG: if(errfp!=NULL) fprintf(errfp, MSG_IC_CONV_FAILED);
                         return(IC_CONV_FAILURE);

    case IC_BAD_EWT:     if(errfp!=NULL) fprintf(errfp, MSG_IC_BAD_EWT);
                         return(IC_BAD_EWT);

  }
  return -99;
}

/*=================================================================*/
/*END          PRIVATE FUNCTIONS IMPLEMENTATION                    */
/*=================================================================*/
