#include "stdlib.h"
#include "sunnls_newton.h"
#include "sundials/sundials_math.h"

#define ZERO   RCONST(0.0)    /* real 0.0    */
#define ONE    RCONST(1.0)    /* real 1.0    */

/* private function to perform Newton iteration */
static int NewtonIter(N_Vector yy_predict, N_Vector yp_predict,
                      N_Vector yy, N_Vector yp,
                      N_Vector ewt, realtype tol, void* mem);

/*
 * Proxy for NLS content
 */
static SUNNonlinSolSysFn Res;
static SUNNonlinSolLSetupFn LSetup;
static SUNNonlinSolLSolveFn LSolve;
static SUNNonlinSolCTestFn  CTest;

static N_Vector cor;
static N_Vector delta;

static long int nni;
static int maxiter;

static realtype cj; /* LMM scaling factor to update yp */

/*
 * Exported Functions
 */

int SUNNonlinSolInit_Newton(N_Vector tmpl)
{
  /* allocate solver-specific memory */

  /* correction vector */
  if (cor == NULL) {
    cor = N_VClone(tmpl);
    if (cor == NULL) {
      /* free nonlinear solver */
      /* set last flag         */
      return(SUN_NLS_MEM_FAIL);
    }
  }    

  /* iteration update vector */
  if (delta == NULL) {
    delta = N_VClone(tmpl);
    if (delta == NULL) {
      /* free nonlinear solver */
      /* set last flag         */
      return(SUN_NLS_MEM_FAIL);
    }
  }    

  /* initialize constants */
  nni     = 0;   /* number of nonlinear iterations         */
  maxiter = 3;   /* maximum number of nonlinear iterations */

  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetSysFn(SUNNonlinSolSysFn SysFn)
{
  Res = SysFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetLSetupFn(SUNNonlinSolLSetupFn LSetupFn)
{
  LSetup = LSetupFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetLSolveFn(SUNNonlinSolLSolveFn LSolveFn)
{
  LSolve = LSolveFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetConvTestFn(SUNNonlinSolCTestFn CTestFn)
{
  CTest = CTestFn;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolFree_Newton()
{
  /* free solver-specific memory */

  if (cor) {
    N_VDestroy(cor);
    cor = NULL;
  }

  if (delta) {
    N_VDestroy(delta);
    delta = NULL;
  }

  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetAlphaFactor_Newton(realtype alpha)
{
  cj = alpha;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolSetMaxNonlinIters_Newton(int max)
{
  maxiter = max;
  return(SUN_NLS_SUCCESS);
}

int SUNNonlinSolGetNumIters_Newton(long int *numiters)
{
  *numiters = nni;
  return(SUN_NLS_SUCCESS);
}

/*
 * This routine attempts to solve the nonlinear system using the linear
 * solver specified. NOTE: this routine uses N_Vector ee as the scratch
 * vector tempv3 passed to lsetup.
 *
 *  Possible return values:
 *
 *  SUN_NLS_SUCCESS
 *
 *  Residual recoverable or non-recoverable errors
 *  SUN_NLS_LSETUP_RECVR    SUN_NLS_LSETUP_FAIL
 *  SUN_NLS_LSOLVE_RECVR    SUN_NLS_LSOLVE_FAIL
 *
 *  SUN_NLS_CONSTR_RECVR
 *  SUN_NLS_NCONV_RECVR
 */

int SUNNonlinSolSolve_Newton(N_Vector yy_predict, N_Vector yp_predict,
                             N_Vector yy, N_Vector yp,
                             N_Vector ee,
                             N_Vector ewt, realtype tol,
                             booleantype callSetup, void* mem)
{
  int retval;
  booleantype tryAgain;

  /* Begin the main loop. This loop is traversed at most twice. 
     The second pass only occurs when the first pass had a recoverable
     failure with old Jacobian data */
  for(;;){

    /* Compute residual at predicted values */
    retval = Res(yy_predict, yp_predict, delta, mem);
    if (retval < 0) return(SUN_NLS_SYS_FAIL);
    if (retval > 0) return(SUN_NLS_SYS_RECVR);

    /* If indicated, call linear solver setup function */
    if (callSetup) {
      retval = LSetup(yy_predict, yp_predict, delta, mem);
      if (retval < 0) return(SUN_NLS_LSETUP_FAIL);
      if (retval > 0) return(SUN_NLS_LSETUP_RECVR);
    }

    /* Call the Newton iteration routine */
    retval = NewtonIter(yy_predict, yp_predict, yy, yp, ewt, tol, mem);

    /* Retry the current step on recoverable failure with old Jacobian data */
    tryAgain = (retval>0) && (!callSetup);

    if (tryAgain){
      callSetup = SUNTRUE;
      continue;
    } else {
      break;
    }

  }  /* end of setup loop */

  /* copy correction and solution to solver mem */
  N_VScale(ONE, cor, ee);

  if (retval != SUN_NLS_SUCCESS) return(retval);

  return(SUN_NLS_SUCCESS);
}


/*
 * NewtonIter
 *
 * This routine performs the Newton iteration.  
 *
 * It assumes that delta contains the initial residual vector on entry.
 *
 * If the iteration succeeds, it returns the value SUN_NLS_SUCCESS = 0.
 * If not, it returns either:
 *   a positive value (for a recoverable failure), namely one of:
 *     residual
 *     SUN_NLS_LSOLVE_RECVR
 *     SUN_NLS_NCONV_RECVR
 * or
 *   a negative value (for a nonrecoverable failure), namely one of:
 *     residual
 *     SUN_NLS_LSOLVE_FAIL
 */

static int NewtonIter(N_Vector yy_predict, N_Vector yp_predict,
                      N_Vector yy, N_Vector yp,
                      N_Vector ewt, realtype tol, void* mem)
{
  int mnewt, retval;
  realtype delnrm;

  /* local variables for fused vector operation */
  realtype cvals[3];
  N_Vector Xvecs[3];

  /* Initialize counter mnewt and cumulative correction vector cor. */
  mnewt = 0;
  N_VConst(ZERO, cor);

  /* Load predictions into y and yp vectors */
  N_VScale(ONE, yy_predict, yy);
  N_VScale(ONE, yp_predict, yp);

  /* Looping point for Newton iteration.  Break out on any error. */
  for(;;) {

    /* increment number of nonlinear solver iterations */
    nni++;

    /* Call the lsolve function to get correction vector delta. */
    retval = LSolve(yy, yp, delta, mem);
    if (retval < 0) return(SUN_NLS_LSOLVE_FAIL);
    if (retval > 0) return(SUN_NLS_LSOLVE_RECVR);

    /* Apply delta to yy, yp, and acro, and get norm(delta). */
    cvals[0] = -ONE;  Xvecs[0] = yy;
    cvals[1] = -ONE;  Xvecs[1] = cor;
    cvals[2] = -cj;   Xvecs[2] = yp;

    retval = N_VScaleAddMulti(3, cvals, delta, Xvecs, Xvecs);
    if (retval != SUN_NLS_SUCCESS) return(SUN_NLS_VECTOROP_ERR);

    /* compute the norm of the correction */
    delnrm = N_VWrmsNorm(delta, ewt);

    /* Test for convergence */
    retval = CTest(mnewt, delnrm, tol, mem);
    if (retval != SUN_NLS_CONV_CONTINUE) return(retval);

    /* Not yet converged. Increment mnewt and test for max allowed. */
    mnewt++;
    if (mnewt >= maxiter) return(SUN_NLS_CONV_RECVR);

    /* Call res for new residual and check error flag from res. */
    retval = Res(yy, yp, delta, mem);
    if (retval < 0) return(SUN_NLS_SYS_FAIL);
    if (retval > 0) return(SUN_NLS_SYS_RECVR);

  } /* end of Newton iteration loop */

  /* All error returns exit here. */
  return(retval);
}
