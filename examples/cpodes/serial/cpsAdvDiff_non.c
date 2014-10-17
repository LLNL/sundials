/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *    ydot = A * y, where A is a banded lower triangular
 * matrix derived from 2-D advection PDE.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Problem Constants */

#define ATOL RCONST(1.0e-6)   /* 1.0e-6 */
#define RTOL RCONST(0.0)      /* 0.0 */

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define TWO    RCONST(2.0)
#define THIRTY RCONST(30.0)

#define MESHX      5
#define MESHY      5
#define NEQ        MESHX*MESHY
#define ALPH1      RCONST(1.0)
#define ALPH2      RCONST(1.0)
#define NOUT       5
#define ML         5
#define MU         0
#define T0         RCONST(0.0)
#define T1         RCONST(0.01)
#define TOUT_MULT  RCONST(10.0)

/* Private Helper Functions */
static realtype MaxError(N_Vector y, realtype t);
static void PrintFinalStats(void *cpode_mem);

/* Functions Called by the Solver */
static int res(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *f_data);
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);


int main()
{
  void *cpode_mem;
  N_Vector y, yp;
  realtype reltol=RTOL, abstol=ATOL, t, tout, erm, hu;
  int flag, iout, qu;

  y = NULL;
  yp = NULL;
  cpode_mem = NULL;

  y = N_VNew_Serial(NEQ);
  N_VConst(ZERO, y);
  NV_Ith_S(y,0) = ONE;

  yp = N_VNew_Serial(NEQ);

  cpode_mem = CPodeCreate(CP_ADAMS, CP_FUNCTIONAL);      
  
  flag = CPodeInitExpl(cpode_mem, f, T0, y);

  /*
  f(T0, y, yp, NULL);
  flag = CPodeInitimpl(cpode_mem, res, T0, y, yp);
  */
  
  flag = CPodeSStolerances(cpode_mem, reltol, abstol);
  
  printf("\n      t        max.err      qu     hu \n");
  for(iout=1, tout=T1; iout <= NOUT; iout++, tout*=TOUT_MULT) {
    flag = CPode(cpode_mem, tout, &t, y, yp, CP_NORMAL);
    if (flag != CP_SUCCESS) break;
    erm = MaxError(y, t);
    flag = CPodeGetLastOrder(cpode_mem, &qu);
    flag = CPodeGetLastStep(cpode_mem, &hu);
    printf("%10.3f  %12.4e   %2d   %12.4e\n", t, erm, qu, hu);
  }
  
  PrintFinalStats(cpode_mem);
  
  CPodeFree(&cpode_mem);
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yp);

  return 0;
}

static int res(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *user_data)
{
  f(t,y,res,user_data);
  N_VLinearSum(1.0, yp, -1.0, res, res);
  return(0);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  long int i, j, k;
  realtype d, *ydata, *dydata;
  
  ydata = NV_DATA_S(y);
  dydata = NV_DATA_S(ydot);

  /*
     Excluding boundaries, 

     ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
         i,j    i,j       i,j             i-1,j             i,j-1
  */

  for (j=0; j < MESHY; j++) {
    for (i=0; i < MESHX; i++) {
      k = i + j * MESHX;
      d = -TWO*ydata[k];
      if (i != 0) d += ALPH1 * ydata[k-1];
      if (j != 0) d += ALPH2 * ydata[k-MESHX];
      dydata[k] = d;
    }
  }

  return(0);
}

static realtype MaxError(N_Vector y, realtype t)
{
  long int i, j, k;
  realtype *ydata, er, ex=ZERO, yt, maxError=ZERO, ifact_inv, jfact_inv=ONE;
  
  if (t == ZERO) return(ZERO);

  ydata = NV_DATA_S(y);
  if (t <= THIRTY) ex = SUN_EXP(-TWO*t);
  
  for (j = 0; j < MESHY; j++) {
    ifact_inv = ONE;
    for (i = 0; i < MESHX; i++) {
      k = i + j * MESHX;
      yt = RPowerI(t,i+j) * ex * ifact_inv * jfact_inv;
      er = SUN_ABS(ydata[k] - yt);
      if (er > maxError) maxError = er;
      ifact_inv /= (i+1);
    }
    jfact_inv /= (j+1);
  }
  return(maxError);
}



static void PrintFinalStats(void *cpode_mem)
{
  long int nst, nfe, nni, ncfn, netf;
  realtype h0u;
  int flag;
  
  flag = CPodeGetActualInitStep(cpode_mem, &h0u);
  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  printf("\n Final statistics:\n\n");
  printf(" Number of steps                          = %4ld \n", nst);
  printf(" Number of f-s                            = %4ld \n", nfe);
  printf(" Number of nonlinear iterations           = %4ld \n", nni);
  printf(" Number of nonlinear convergence failures = %4ld \n", ncfn);
  printf(" Number of error test failures            = %4ld \n", netf);
  printf(" Initial step size                        = %g \n\n", h0u);

}

