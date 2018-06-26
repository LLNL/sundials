#include "ida_nls.h"
#include "sunnls_newton.h"
#include "sundials/sundials_math.h"

/* nonlinear solver constants */
#define MAXIT   4              /* max number of nonlinear iterations */
#define RATEMAX RCONST(0.9)    /* max convergence rate               */
#define PT0001  RCONST(0.0001) /* real 0.0001                        */

/* private functions passed to nonlinear solver */
static int IDANls_Res(N_Vector yy, N_Vector yp, N_Vector res, void* ida_mem);
static int IDANls_LSetup(N_Vector yy, N_Vector yp, N_Vector res, void* ida_mem);
static int IDANls_LSolve(N_Vector yy, N_Vector yp, N_Vector delta, void* ida_mem);
static int IDANls_ConvergenceTest(int m, realtype delnrm, realtype tol, void* ida_mem);

/* attach nonlinear solver to IDA */
int IDASetNonlinearSolver(void *ida_mem)
{
  IDAMem IDA_mem;

  /* Return immediately if IDA memory is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", 
                    "IDASetNonlinearSolver", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* initialize nonlinear solver */
  SUNNonlinSolInit_Newton(IDA_mem->ida_yy);

  /* set the residual function */
  SUNNonlinSolSetSysFn(IDANls_Res);

  /* set the setup function */
  SUNNonlinSolSetLSetupFn(IDANls_LSetup);

  /* set the linear solve function */
  SUNNonlinSolSetLSolveFn(IDANls_LSolve);

  /* set convergence test function */
  SUNNonlinSolSetConvTestFn(IDANls_ConvergenceTest);

  /* set max allowed nonlinear iterations */
  SUNNonlinSolSetMaxNonlinIters_Newton(MAXIT);

  /* set initial convergence rate factor */
  SUNNonlinSolSetAlphaFactor_Newton(IDA_mem->ida_cj);

  return(IDA_SUCCESS);
}


static int IDANls_LSetup(N_Vector yy, N_Vector yp, N_Vector res, void* ida_mem)
{
  IDAMem   IDA_mem;
  N_Vector tempv3;
  int      retval;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDANlsRes", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  tempv3 = IDA_mem->ida_ee;

  IDA_mem->ida_nsetups++;
  retval = IDA_mem->ida_lsetup(IDA_mem, yy, yp, res,
                               IDA_mem->ida_tempv1, IDA_mem->ida_tempv2, tempv3);

  IDA_mem->ida_cjold = IDA_mem->ida_cj;
  IDA_mem->ida_cjratio = ONE;
  IDA_mem->ida_ss = TWENTY;

  return(retval);
}


static int IDANls_LSolve(N_Vector yy, N_Vector yp, N_Vector delta, void* ida_mem)
{
  IDAMem IDA_mem;
  int    retval;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDANlsRes", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  retval = IDA_mem->ida_lsolve(IDA_mem, delta, IDA_mem->ida_ewt, yy, yp,
                               IDA_mem->ida_savres);
  return(retval);
}


/* only called from within the nonlinear solver so can assume predictor
 * has already been called */
static int IDANls_Res(N_Vector yy, N_Vector yp, N_Vector res, void* ida_mem)
{
  IDAMem IDA_mem;
  int retval;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDANlsRes", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  retval = IDA_mem->ida_res(IDA_mem->ida_tn, yy, yp, res, IDA_mem->ida_user_data);

  /* increment the number of residual evaluations */
  IDA_mem->ida_nre++;

  /* save a copy of the residual vector in savres */
  N_VScale(ONE, res, IDA_mem->ida_savres);

  return(retval);
}


static int IDANls_ConvergenceTest(int m, realtype delnrm, realtype tol, void* ida_mem)
{
  IDAMem IDA_mem;
  realtype rate;
  static realtype oldnrm;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDANlsRes", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test for convergence, first directly, then with rate estimate. */
  if (m == 0){
    oldnrm = delnrm;
    if (delnrm <= PT0001 * IDA_mem->ida_toldel) return(SUN_NLS_SUCCESS);
  } else {
    rate = SUNRpowerR( delnrm/oldnrm, ONE/m );
    if (rate > RATEMAX) return(SUN_NLS_CONV_RECVR);
    IDA_mem->ida_ss = rate/(ONE - rate);
  }

  if (IDA_mem->ida_ss*delnrm <= tol) return(SUN_NLS_SUCCESS);

  /* Not yet converged */
  return(SUN_NLS_CONV_CONTINUE);
}
