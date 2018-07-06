/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * This the implementation file for the CVODE nonlinear solver interface.
 * ---------------------------------------------------------------------------*/

#include "cvode_impl.h"
#include "sundials/sundials_math.h"

/* nonlinear solver constants */
#define ZERO    RCONST(0.0)     /* real 0.0     */
#define ONE     RCONST(1.0)     /* real 1.0     */
#define TWO     RCONST(2.0)     /* real 2.0     */

#define CRDOWN  RCONST(0.3)
#define DGMAX   RCONST(0.3)

#define RDIV    TWO
#define MSBP    20

/* return values */
#define CONV_FAIL        +4 
#define TRY_AGAIN        +5

#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8

#define RHSFUNC_RECVR    +9

/* private functions */
static int cvNewtonIteration(CVodeMem cv_mem);
static int cvNlsRes(N_Vector y, N_Vector res, void* cvode_mem);

static N_Vector delta;


/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

int cvNlsInit(CVodeMem cvode_mem)
{
  delta = N_VClone(cvode_mem->cv_acor);

  return(CV_SUCCESS);
}


int cvNlsFree(CVodeMem cvode_mem)
{
  N_VDestroy(delta);

  return(CV_SUCCESS);
}


/*
 * cvNlsFunctional
 *
 * This routine attempts to solve the nonlinear system using 
 * functional iteration (no matrices involved).
 *
 * Possible return values are:
 *
 *   CV_SUCCESS      --->  continue with error test
 *
 *   CV_RHSFUNC_FAIL --->  halt the integration
 *
 *   CONV_FAIL       -+
 *   RHSFUNC_RECVR   -+->  predict again or stop if too many
 *
 */

int cvNlsFunctional(CVodeMem cv_mem)
{
  int retval, m;
  realtype del, delp, dcon;

  /* Initialize counter and evaluate f at predicted y */
  
  cv_mem->cv_crate = ONE;
  m = 0;

  retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_zn[0],
                        cv_mem->cv_tempv, cv_mem->cv_user_data);
  cv_mem->cv_nfe++;
  if (retval < 0) return(CV_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  N_VConst(ZERO, cv_mem->cv_acor);

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Loop until convergence; accumulate corrections in acor */

  for(;;) {

    cv_mem->cv_nni++;

    /* Correct y directly from the last f value */
    N_VLinearSum(cv_mem->cv_h, cv_mem->cv_tempv, -ONE,
                 cv_mem->cv_zn[1], cv_mem->cv_tempv);
    N_VScale(cv_mem->cv_rl1, cv_mem->cv_tempv, cv_mem->cv_tempv);
    N_VLinearSum(ONE, cv_mem->cv_zn[0], ONE, cv_mem->cv_tempv, cv_mem->cv_y);
    /* Get WRMS norm of current correction to use in convergence test */
    N_VLinearSum(ONE, cv_mem->cv_tempv, -ONE, cv_mem->cv_acor, cv_mem->cv_acor);
    del = N_VWrmsNorm(cv_mem->cv_acor, cv_mem->cv_ewt);
    N_VScale(ONE, cv_mem->cv_tempv, cv_mem->cv_acor);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) cv_mem->cv_crate = SUNMAX(CRDOWN * cv_mem->cv_crate, del / delp);
    dcon = del * SUNMIN(ONE, cv_mem->cv_crate) / cv_mem->cv_tq[4];
    if (dcon <= ONE) {
      cv_mem->cv_acnrm = (m == 0) ?
        del : N_VWrmsNorm(cv_mem->cv_acor, cv_mem->cv_ewt);
      return(CV_SUCCESS);  /* Convergence achieved */
    }

    /* Stop at maxcor iterations or if iter. seems to be diverging */
    m++;
    if ((m==cv_mem->cv_maxcor) || ((m >= 2) && (del > RDIV * delp)))
      return(CONV_FAIL);

    /* Save norm of correction, evaluate f, and loop again */
    delp = del;

    retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_y,
                          cv_mem->cv_tempv, cv_mem->cv_user_data);
    cv_mem->cv_nfe++;
    if (retval < 0) return(CV_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

  }
}


/*
 * cvNlsNewton
 *
 * This routine handles the Newton iteration. It calls lsetup if 
 * indicated, calls cvNewtonIteration to perform the iteration, and 
 * retries a failed attempt at Newton iteration if that is indicated.
 *
 * Possible return values:
 *
 *   CV_SUCCESS       ---> continue with error test
 *
 *   CV_RHSFUNC_FAIL  -+  
 *   CV_LSETUP_FAIL    |-> halt the integration 
 *   CV_LSOLVE_FAIL   -+
 *
 *   CONV_FAIL        -+
 *   RHSFUNC_RECVR    -+-> predict again or stop if too many
 *
 */

int cvNlsNewton(CVodeMem cv_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, retval, ier;
  booleantype callSetup;
  
  vtemp1 = cv_mem->cv_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = cv_mem->cv_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = cv_mem->cv_tempv; /* rename tempv as vtemp3 for readability */
  
  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    CV_NO_FAILURES : CV_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (cv_mem->cv_lsetup) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (cv_mem->cv_nst == 0) ||
      (cv_mem->cv_nst >= cv_mem->cv_nstlp + MSBP) ||
      (SUNRabs(cv_mem->cv_gamrat-ONE) > DGMAX);
  } else {  
    cv_mem->cv_crate = ONE;
    callSetup = SUNFALSE;
  }
  
  /* Looping point for the solution of the nonlinear system.
     Evaluate f at the predicted y, call lsetup if indicated, and
     call cvNewtonIteration for the Newton iteration itself.      */
  

  for(;;) {

    /* compute the residual */
    retval = cvNlsRes(cv_mem->cv_zn[0], delta, cv_mem);
    if (retval != CV_SUCCESS) return(retval);

    if (callSetup) {
      ier = cv_mem->cv_lsetup(cv_mem, convfail, cv_mem->cv_zn[0],
                              cv_mem->cv_ftemp, &(cv_mem->cv_jcur),
                              vtemp1, vtemp2, vtemp3);
      cv_mem->cv_nsetups++;
      callSetup = SUNFALSE;
      cv_mem->cv_gamrat = cv_mem->cv_crate = ONE; 
      cv_mem->cv_gammap = cv_mem->cv_gamma;
      cv_mem->cv_nstlp = cv_mem->cv_nst;
      /* Return if lsetup failed */
      if (ier < 0) return(CV_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, cv_mem->cv_acor);
    N_VScale(ONE, cv_mem->cv_zn[0], cv_mem->cv_y);

    /* Do the Newton iteration */
    ier = cvNewtonIteration(cv_mem);

    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=CV_FAIL_BAD_J.  Otherwise return.                 */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = SUNTRUE;
    convfail = CV_FAIL_BAD_J;
  }
}


/*
 * cvNewtonIteration
 *
 * This routine performs the Newton iteration. If the iteration succeeds,
 * it returns the value CV_SUCCESS. If not, it may signal the cvNlsNewton 
 * routine to call lsetup again and reattempt the iteration, by
 * returning the value TRY_AGAIN. (In this case, cvNlsNewton must set 
 * convfail to CV_FAIL_BAD_J before calling setup again). 
 * Otherwise, this routine returns one of the appropriate values 
 * CV_LSOLVE_FAIL, CV_RHSFUNC_FAIL, CONV_FAIL, or RHSFUNC_RECVR back 
 * to cvNlsNewton.
 */

static int cvNewtonIteration(CVodeMem cv_mem)
{
  int m, retval;
  realtype del, delp, dcon;

  cv_mem->cv_mnewt = m = 0;

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Looping point for Newton iteration */
  for(;;) {

    /* Call the lsolve function */
    retval = cv_mem->cv_lsolve(cv_mem, delta, cv_mem->cv_ewt,
                               cv_mem->cv_y, cv_mem->cv_ftemp); 
    cv_mem->cv_nni++;
    
    if (retval < 0) return(CV_LSOLVE_FAIL);
    
    /* If lsolve had a recoverable failure and Jacobian data is
       not current, signal to try the solution again            */
    if (retval > 0) { 
      if ((!cv_mem->cv_jcur) && (cv_mem->cv_lsetup))
        return(TRY_AGAIN);
      else
        return(CONV_FAIL);
    }

    /* Get WRMS norm of correction; add correction to acor and y */
    del = N_VWrmsNorm(delta, cv_mem->cv_ewt);
    N_VLinearSum(ONE, cv_mem->cv_acor, ONE, delta, cv_mem->cv_acor);
    N_VLinearSum(ONE, cv_mem->cv_zn[0], ONE, cv_mem->cv_acor, cv_mem->cv_y);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) {
      cv_mem->cv_crate = SUNMAX(CRDOWN * cv_mem->cv_crate, del/delp);
    }
    dcon = del * SUNMIN(ONE, cv_mem->cv_crate) / cv_mem->cv_tq[4];
    
    if (dcon <= ONE) {
      cv_mem->cv_acnrm = (m==0) ?
        del : N_VWrmsNorm(cv_mem->cv_acor, cv_mem->cv_ewt);
      cv_mem->cv_jcur = SUNFALSE;
      return(CV_SUCCESS); /* Nonlinear system was solved successfully */
    }
    
    cv_mem->cv_mnewt = ++m;
    
    /* Stop at maxcor iterations or if iter. seems to be diverging.
       If still not converged and Jacobian data is not current, 
       signal to try the solution again                            */
    if ((m == cv_mem->cv_maxcor) || ((m >= 2) && (del > RDIV*delp))) {
      if ((!cv_mem->cv_jcur) && (cv_mem->cv_lsetup))
        return(TRY_AGAIN);
      else
        return(CONV_FAIL);
    }
    
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;

    /* Evaluate the residual of the nonlinear system */
    retval = cvNlsRes(cv_mem->cv_y, delta, cv_mem);

    if (retval < 0) return(CV_RHSFUNC_FAIL);
    if (retval > 0) {
      if ((!cv_mem->cv_jcur) && (cv_mem->cv_lsetup))
        return(TRY_AGAIN);
      else
        return(RHSFUNC_RECVR);
    }

  } /* end loop */
}


static int cvNlsRes(N_Vector y, N_Vector res, void* cvode_mem)
{
  CVodeMem cv_mem;
  int retval;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "cvNlsRes", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  retval = cv_mem->cv_f(cv_mem->cv_tn, y, cv_mem->cv_ftemp,
                        cv_mem->cv_user_data);
  cv_mem->cv_nfe++;
  if (retval < 0) return(CV_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* compute the accumulated correction */
  N_VLinearSum(ONE, y, -ONE, cv_mem->cv_zn[0], res);

  /* evaluate the residual of the nonlinear system */
  N_VLinearSum(cv_mem->cv_rl1, cv_mem->cv_zn[1], ONE, res, res);
  N_VLinearSum(cv_mem->cv_gamma, cv_mem->cv_ftemp, -ONE, res, res);

  return(0);
}
