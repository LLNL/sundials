/***********************************************************************
 * File       : irobx.c
 * Written by : Allan G. Taylor, Alan C. Hindmarsh, and Radu Serban
 * Version of : 23 July 2003
 *----------------------------------------------------------------------
 *
 * This simple example problem for IDAS, due to Robertson, is from
 * chemical kinetics, and consists of the following three equations:
 *
 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions
 * y1 = 1, y2 = y3 = 0.
 *
 * The problem is solved with IDAS using IDADENSE for the linear solver,
 * with a user-supplied Jacobian.
 * Output is printed at t = .4, 4, 40, ..., 4e10.
 *************************************************************************/

#include <stdio.h>
#include <math.h>
#include "sundialstypes.h"
#include "sundialsmath.h"
#include "nvector_serial.h" /* definition of type N_Vector and macros NV_DATA_S */
#include "idas.h"
#include "idasdense.h"

/* Macro to define dense matrix elements, indexed from 1. */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

int resrob(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *rdata);

int jacrob(integertype Neq, realtype tt, N_Vector yy, N_Vector yp,
           realtype cj, void *jdata, N_Vector resvec, DenseMat JJ,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

#define NOUT  12

int main()
{
  void *mem;
  NV_Spec nvSpec;
  N_Vector yy, yp, avtol;
  integertype SystemSize = 3;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, t1, tout, tret, hused;
  int iout, itol, itask;
  int retval, kused, nst, nni, nje, nre, nreD, netf, ncfn;

  /* Initialize serial machine environment */
  nvSpec = NV_SpecInit_Serial(SystemSize);

  /* Allocate N-vectors. */
  yy = N_VNew(nvSpec); 
  yp = N_VNew(nvSpec);
  avtol = N_VNew(nvSpec);

  /* Set various input parameters and initialize y and y'. */
  yval  = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);
  atval = NV_DATA_S(avtol);
  t0   = 0.0;
  t1   = 0.4;
  itol = SV;  /* scalar relative tolerance, vector absolute tolerance. */
  rtol = 1.e-4;
  atval[0] = 1.e-6;
  atval[1] = 1.e-10;
  atval[2] = 1.e-6;
  yval[0] = 1.0;
  yval[1] = 0.0;
  yval[2] = 0.0;
  ypval[0]  = -0.04;
  ypval[1]  = +0.04;
  ypval[2]  =   0.;  
  itask = NORMAL;

  /* Call IDACreate and IDAMalloc to initialize solution */
  mem = IDACreate();
  if (mem == NULL) { printf("IDACreate failed."); return(1); }
  retval = IDAMalloc(mem, resrob, t0, yy, yp, itol, &rtol, avtol, nvSpec);
  if (retval != SUCCESS) {printf("IDAMalloc failed.\n"); return(1); }

  /* Call IDADense and set up the linear solver package. */
  retval = IDADense(mem, SystemSize);
  if (retval != SUCCESS) {printf("IDADense failed.\n"); return(1); }
  retval = IDADenseSetJacFn(mem, jacrob);
  if (retval != SUCCESS) {printf("IDADenseSetJacFn failed.\n"); return(1); }

  printf("irobx: Robertson kinetics DAE serial example problem for IDAS \n");
  printf("       Three equation chemical kinetics problem. \n\n");
  printf("Linear solver: IDADENSE, with user-supplied Jacobian.\n");
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",yval[0], yval[1], yval[2]);
  printf("Constraints and id not used.\n");
  printf("...............................................\n");

  /* Loop over tout values and call IDASolve. */

  for (tout = t1, iout = 1; iout <= NOUT ; iout++, tout *= 10.0) {
    retval=IDASolve(mem, tout, &tret, yy, yp, itask);
    if (retval != SUCCESS) {printf("IDASolve returned %d\n",retval); return(1); }
    
    printf("\nt = %g,  y = (%g, %g, %g)\n", tret, yval[0],yval[1],yval[2]);
    
    IDAGetLastOrder(mem, &kused);
    IDAGetNumSteps(mem, &nst);
    IDAGetNumNonlinSolvIters(mem, &nni);
    IDAGetNumResEvals(mem, &nre);
    IDAGetLastStep(mem, &hused);
    IDADenseGetNumJacEvals(mem, &nje);
    IDADenseGetNumResEvals(mem, &nreD);

    printf("k = %d, nst = %d, nni = %d, nje = %d, nre = %d, nreD = %d  h = %g\n",
           kused, nst, nni, nje, nre, nreD, hused);

  } /* End of tout loop. */

  IDAGetNumErrTestFails(mem, &netf);
  IDAGetNumNonlinSolvConvFails(mem, &ncfn);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps =                %d\n",nst);
  printf("Number of residual evaluations = %d\n",nre+nreD);
  printf("Number of Jacobian evaluations = %d\n",nje);
  printf("Number of error test failures =  %d\n",netf);
  printf("Number of nonlinear (corrector) convergence failures = %d\n",ncfn);

  IDAFree(mem);
  N_VFree(yy);
  N_VFree(yp);
  N_VFree(avtol);
  
  NV_SpecFree_Serial(nvSpec);

  return(0);
  
} /* End of irobx main. */


/**************************************************************************/

/*  Define the system residual function resrob. */

int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *rdata)
{
  realtype *yval, *ypval, *rval;

  yval = NV_DATA_S(yy); ypval = NV_DATA_S(yp); rval = NV_DATA_S(rr);

  rval[0]  = -0.04*yval[0] + 1.e4*yval[1]*yval[2];
  rval[1]  = -rval[0]    - 3.e7*yval[1]*yval[1] - ypval[1];
  rval[0] -=  ypval[0];
  rval[2]  =  yval[0] + yval[1] + yval[2] - 1.0;

  return(SUCCESS);

} /* End of residual function resrob. */

/* Define the Jacobian function jacrob. */

int jacrob(integertype Neq, realtype tt, N_Vector yy, N_Vector yp,
           realtype cj, void *jdata, N_Vector resvec, DenseMat JJ,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{

  realtype *yval;
  
  yval = NV_DATA_S(yy);

  IJth(JJ,1,1) = -0.04 - cj;
  IJth(JJ,2,1) = +0.04;
  IJth(JJ,3,1) =  1.  ;
  IJth(JJ,1,2) =  1.e4*yval[2];
  IJth(JJ,2,2) = -1.e4*yval[2] - 6.e7*yval[1] - cj;
  IJth(JJ,3,2) =  1.;
  IJth(JJ,1,3) =  1.e4*yval[1];
  IJth(JJ,2,3) = -1.e4*yval[1];
  IJth(JJ,3,3) =  1.;

  return(SUCCESS);

} /* End of the Jacobian function jacrob. */
