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

#define NLS_MAXCOR 3

/* return values */
#define CONV_FAIL        +4 
#define TRY_AGAIN        +5

#define RHSFUNC_RECVR    +9

/* NLS VALUES */

/* Recoverable */
#define SUN_NLS_CONTINUE   +1  /* not converged, keep iterating      */
#define SUN_NLS_CONV_RECVR +2  /* convergece failure, try to recover */

/* private functions */
static int cvNls_Res(N_Vector y, N_Vector res, void* cvode_mem);

static int cvNls_LSetup(N_Vector ycor, N_Vector res, booleantype jbad,
                        booleantype* jcur, void* cvode_mem);
static int cvNls_LSolve(N_Vector ycor, N_Vector delta, void* cvode_mem);
static int cvNls_ConvTest(SUNNonlinearSolver NLS, N_Vector ycor, N_Vector del,
                          realtype tol, N_Vector ewt, void* cvode_mem);

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

int CVodeSetNonlinearSolver(void *cvode_mem, SUNNonlinearSolver NLS)
{
  CVodeMem cv_mem;
  int retval;

  /* Return immediately if CVode memory is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "CVSetNonlinearSolver", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Return immediately if NLS memory is NULL */
  if (NLS == NULL) {
    cvProcessError(NULL, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                   "NLS must be non-NULL");
    return (CV_ILL_INPUT);
  }

  /* check for required nonlinear solver functions */
  if ( NLS->ops->gettype    == NULL ||
       NLS->ops->initialize == NULL ||
       NLS->ops->solve      == NULL ||
       NLS->ops->free       == NULL ||
       NLS->ops->setsysfn   == NULL ) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                   "NLS does not support required operations");
    return(CV_ILL_INPUT);
  }

  /* free any existing nonlinear solver */
  if (cv_mem->NLS)
    retval = SUNNonlinSolFree(cv_mem->NLS);

  /* set SUNNonlinearSolver pointer */
  cv_mem->NLS = NLS;

  /* set the nonlinear residual function */
  if (SUNNonlinSolGetType(NLS) == SUNNONLINEARSOLVER_ROOTFIND) {
    retval = SUNNonlinSolSetSysFn(cv_mem->NLS, cvNls_Res);
    if (retval != CV_SUCCESS) {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                     "Setting nonlinear system function failed");
      return(CV_ILL_INPUT);
    }
  } else if (SUNNonlinSolGetType(NLS) ==  SUNNONLINEARSOLVER_STATIONARY) {
    /* >>>>>>> NEED TO UPDATE FOR FIXED POINT <<<<<<< */
    retval = SUNNonlinSolSetSysFn(cv_mem->NLS, cvNls_Res);
    if (retval != CV_SUCCESS) {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                     "Setting nonlinear system function failed");
      return(CV_ILL_INPUT);
    }
  } else {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                   "Invalid nonlinear solver type");
    return(CV_ILL_INPUT);
  }



  /* set convergence test function */
  retval = SUNNonlinSolSetConvTestFn(cv_mem->NLS, cvNls_ConvTest);
  if (retval != CV_SUCCESS) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                   "Setting convergence test function failed");
    return(CV_ILL_INPUT);
  }

  /* set max allowed nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(cv_mem->NLS, NLS_MAXCOR);
  if (retval != CV_SUCCESS) {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVSetNonlinearSolver",
                   "Setting maximum number of nonlinear iterations failed");
    return(CV_ILL_INPUT);
  }

  return(CV_SUCCESS);
}


int cvNlsInit(CVodeMem cvode_mem)
{
  int retval;

  if (cvode_mem->cv_iter == CV_NEWTON) {
    retval = SUNNonlinSolSetLSetupFn(cvode_mem->NLS, cvNls_LSetup);
    retval = SUNNonlinSolSetLSolveFn(cvode_mem->NLS, cvNls_LSolve);
    retval = SUNNonlinSolInitialize(cvode_mem->NLS);
  }

  return(CV_SUCCESS);
}


int cvNlsFree(CVodeMem cvode_mem)
{
  return(CV_SUCCESS);
}


static int cvNls_LSetup(N_Vector ycor, N_Vector res, booleantype jbad,
                        booleantype* jcur, void* cvode_mem)
{
  CVodeMem cv_mem;
  int      retval;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "cvNls_LSetup", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* if the nonlinear solver marked the Jacobian as bad update convfail */
  if (jbad)
    cv_mem->convfail = CV_FAIL_BAD_J;

  /* setup the linear solver */
  retval = cv_mem->cv_lsetup(cv_mem, cv_mem->convfail, cv_mem->cv_y, cv_mem->cv_ftemp,
                             &(cv_mem->cv_jcur), cv_mem->cv_vtemp1, cv_mem->cv_vtemp2,
                             cv_mem->cv_vtemp3);
  cv_mem->cv_nsetups++;

  /* update Jacobian status */
  *jcur = cv_mem->cv_jcur;

  cv_mem->cv_gamrat = cv_mem->cv_crate = ONE;
  cv_mem->cv_gammap = cv_mem->cv_gamma;
  cv_mem->cv_nstlp  = cv_mem->cv_nst;

  if (retval < 0) return(CV_LSETUP_FAIL);
  if (retval > 0) return(SUN_NLS_CONV_RECVR);

  return(CV_SUCCESS);
}


static int cvNls_LSolve(N_Vector ycor, N_Vector delta, void* cvode_mem)
{
  CVodeMem cv_mem;
  int      retval;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "cvNls_LSolve", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  retval = cv_mem->cv_lsolve(cv_mem, delta, cv_mem->cv_ewt, cv_mem->cv_y, cv_mem->cv_ftemp);

  if (retval < 0) return(CV_LSOLVE_FAIL);
  if (retval > 0) return(SUN_NLS_CONV_RECVR);

  return(CV_SUCCESS);
}


static int cvNls_ConvTest(SUNNonlinearSolver NLS, N_Vector ycor, N_Vector delta,
                          realtype tol, N_Vector ewt, void* cvode_mem)
{
  CVodeMem cv_mem;
  int m, retval;
  realtype del;
  realtype dcon;
  static realtype delp;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "cvNlsConvTest", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* compute the norm of the correction */
  del = N_VWrmsNorm(delta, ewt);

  /* get the current nonlinear solver iteration count */
  retval = SUNNonlinSolGetCurIter(NLS, &m);
  if (retval != CV_SUCCESS) return(CV_MEM_NULL);

  /* Test for convergence. If m > 0, an estimate of the convergence
     rate constant is stored in crate, and used in the test.        */
  if (m > 0) {
    cv_mem->cv_crate = SUNMAX(CRDOWN * cv_mem->cv_crate, del/delp);
  }
  dcon = del * SUNMIN(ONE, cv_mem->cv_crate) / tol;

  if (dcon <= ONE) {
    cv_mem->cv_acnrm = (m==0) ? del : N_VWrmsNorm(ycor, cv_mem->cv_ewt);
    return(CV_SUCCESS); /* Nonlinear system was solved successfully */
  }

  /* check if the iteration seems to be diverging */
  if ((m >= 1) && (del > RDIV*delp)) return(SUN_NLS_CONV_RECVR);

  /* Save norm of correction and loop again */
  delp = del;

  /* Not yet converged */
  return(SUN_NLS_CONTINUE);
}


static int cvNls_Res(N_Vector ycor, N_Vector res, void* cvode_mem)
{
  CVodeMem cv_mem;
  int retval;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "cvNls_Res", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* update the state based on the current correction */
  N_VLinearSum(ONE, cv_mem->cv_zn[0], ONE, ycor, cv_mem->cv_y);

  /* evaluate the rhs function */
  retval = cv_mem->cv_f(cv_mem->cv_tn, cv_mem->cv_y, cv_mem->cv_ftemp,
                        cv_mem->cv_user_data);
  cv_mem->cv_nfe++;
  if (retval < 0) return(CV_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  N_VLinearSum(cv_mem->cv_rl1, cv_mem->cv_zn[1], ONE, ycor, res);
  N_VLinearSum(-cv_mem->cv_gamma, cv_mem->cv_ftemp, ONE, res, res);

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
