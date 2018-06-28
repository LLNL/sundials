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
 * This the implementation file for the IDA nonlinear solver interface.
 * ---------------------------------------------------------------------------*/

#include "ida_impl.h"
#include "sunnls_newton.h" /* relace with sundial_nonlinearsolve.h */
#include "sundials/sundials_math.h"

/* nonlinear solver constants */
#define MAXIT   4              /* max number of nonlinear iterations */
#define RATEMAX RCONST(0.9)    /* max convergence rate               */
#define PT0001  RCONST(0.0001) /* real 0.0001                        */
#define ONE     RCONST(1.0)    /* real 1.0                           */
#define TWENTY  RCONST(20.0)   /* real 20.0                          */

/* private functions passed to nonlinear solver */
static int IDANls_Res(N_Vector yy, N_Vector res, void* ida_mem);
static int IDANls_LSetup(N_Vector yy, N_Vector res, void* ida_mem);
static int IDANls_LSolve(N_Vector yy, N_Vector delta, void* ida_mem);
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
  if (IDA_mem->ida_lsetup)
    SUNNonlinSolSetLSetupFn(IDANls_LSetup);

  /* set the linear solve function */
  SUNNonlinSolSetLSolveFn(IDANls_LSolve);

  /* set convergence test function */
  SUNNonlinSolSetConvTestFn(IDANls_ConvergenceTest);

  /* set max allowed nonlinear iterations */
  SUNNonlinSolSetMaxNonlinIters_Newton(MAXIT);

  return(IDA_SUCCESS);
}


static int IDANls_LSetup(N_Vector yy, N_Vector res, void* ida_mem)
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
  retval = IDA_mem->ida_lsetup(IDA_mem, yy, IDA_mem->ida_yp, res,
                               IDA_mem->ida_tempv1, IDA_mem->ida_tempv2, tempv3);

  IDA_mem->ida_cjold = IDA_mem->ida_cj;
  IDA_mem->ida_cjratio = ONE;
  IDA_mem->ida_ss = TWENTY;

  return(retval);
}


static int IDANls_LSolve(N_Vector yy, N_Vector delta, void* ida_mem)
{
  IDAMem IDA_mem;
  int    retval;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDANlsRes", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  retval = IDA_mem->ida_lsolve(IDA_mem, delta, IDA_mem->ida_ewt, yy, IDA_mem->ida_yp,
                               IDA_mem->ida_savres);
  return(retval);
}


static int IDANls_Res(N_Vector yy, N_Vector res, void* ida_mem)
{
  IDAMem IDA_mem;
  int retval;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDANlsRes", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* compute current yp */
  N_VLinearSum(IDA_mem->ida_cj, yy, ONE, IDA_mem->ida_ypbeta, IDA_mem->ida_yp);

  /* evaluate residual */
  retval = IDA_mem->ida_res(IDA_mem->ida_tn, yy, IDA_mem->ida_yp, res, IDA_mem->ida_user_data);

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
