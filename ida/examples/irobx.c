/***********************************************************************
 * File:       irobx.c   
 * Written by: Allan G. Taylor and Alan C. Hindmarsh
 * Version of: 2 March 2001
 *----------------------------------------------------------------------
 *
 * This simple example problem for IDA, due to Robertson, is from
 * chemical kinetics, and consists of the following three equations:
 *
 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions
 * y1 = 1, y2 = y3 = 0.
 *
 * The problem is solved with IDA using IDADENSE for the linear solver,
 * with a user-supplied Jacobian.
 * Output is printed at t = .4, 4, 40, ..., 4e10.
 *************************************************************************/

#include <stdio.h>
#include <math.h>
#include "llnltyps.h"
#include "llnlmath.h"
#include "nvector.h"
#include "ida.h"
#include "idadense.h"

/* Macro to define dense matrix elements, indexed from 1. */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

int resrob(integer Neq, real tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *rdata);

int jacrob(integer Neq, real tt, N_Vector yy, N_Vector yp,
           real cj, N_Vector constraints, ResFn res, void *rdata,
           void *jdata, N_Vector resvec, N_Vector ewt, real hh,
           real uround, DenseMat JJ, long int *nrePtr,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

#define NOUT  12

main()
{
  int itol, itask;
  integer SystemSize = 3, retval, i,iout;
  long int iopt[OPT_SIZE];
  boole optIn;
  real ropt[OPT_SIZE], rtol, *yval, *ypval, *atval;
  real t0, t1, tout, tret;
  N_Vector yy, yp, avtol;
  void *mem;


  /* Allocate N-vectors. */
  yy = N_VNew(SystemSize,NULL); 
  yp = N_VNew(SystemSize,NULL);
  avtol = N_VNew(SystemSize,NULL);

  /* Set various input parameters and initialize y and y'. */
  yval  = N_VDATA(yy);
  ypval = N_VDATA(yp);
  atval = N_VDATA(avtol);
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
  optIn = FALSE; 
  itask = NORMAL;

  /* Call IDAMalloc to set up problem memory.
  NULL arguments are: id, constraints, errfp, and machEnv, respectively. */

  mem = IDAMalloc(SystemSize, resrob, NULL, t0, yy, yp, itol, &rtol, avtol,
                  NULL , NULL , NULL, optIn, iopt, ropt,  NULL);
  mem = (IDAMem)mem;
  if (mem == NULL) {printf("IDAMalloc failed.\n"); return(1); }

  /* Call IDADense and set up the linear solver package. */

  retval = IDADense(mem, jacrob, NULL);
  if (retval != SUCCESS) {printf("IDADense failed.\n"); return(1); }

  printf("irobx: Robertson kinetics DAE serial example problem for IDA \n");
  printf("       Three equation chemical kinetics problem. \n\n");
  printf("Linear solver: IDADENSE, with user-supplied Jacobian.\n");
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
                                 rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",yval[0], yval[1], yval[2]);
  printf("Constraints and id not used.\n");
  printf("...............................................\n");

  /* Loop over tout values and call IDASolve. */

  for (tout = t1, iout = 1; iout <= NOUT ; iout++, tout *= 10.0) {
    retval=IDASolve(mem, tout, t0, &tret, yy, yp, itask);
    if (retval != SUCCESS) {printf("IDASolve returned %d\n",retval); return(1); }

    printf("\nt = %g,  y = (%g, %g, %g)\n", tret, yval[0],yval[1],yval[2]);
 
    printf("k = %d, nst = %d, nni = %d, nje = %d, nre = %d, h = %g\n",
           iopt[KUSED], iopt[NST],iopt[NNI], iopt[DENSE_NJE], 
           iopt[NRE], ropt[HUSED]);

  } /* End of tout loop. */

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps =                %d\n",iopt[NST]);
  printf("Number of residual evaluations = %d\n",iopt[NRE]);
  printf("Number of Jacobian evaluations = %d\n",iopt[DENSE_NJE]);
  printf("Number of error test failures =  %d\n",iopt[NETF]);
  printf("Number of nonlinear (corrector) convergence failures = %d\n",
         iopt[NCFN]);

  IDAFree(mem);
  N_VFree(yy);
  N_VFree(yp);
  N_VFree(avtol);

} /* End of irobx main. */


/**************************************************************************/

/*  Define the system residual function resrob. */

int resrob(integer Neq, real tres, N_Vector yy, N_Vector yp, N_Vector rr, 
           void *rdata)
{
  int i;
  real *yval, *ypval, *rval;

  yval = N_VDATA(yy); ypval = N_VDATA(yp); rval = N_VDATA(rr);

  rval[0]  = -0.04*yval[0] + 1.e4*yval[1]*yval[2];
  rval[1]  = -rval[0]    - 3.e7*yval[1]*yval[1] - ypval[1];
  rval[0] -=  ypval[0];
  rval[2]  =  yval[0] + yval[1] + yval[2] - 1.0;

  return(SUCCESS);

} /* End of residual function resrob. */

/* Define the Jacobian function jacrob. */

int jacrob(integer Neq, real tt, N_Vector yy, N_Vector yp,
                   real cj, N_Vector constraints, ResFn res, void *rdata,
                   void *jdata, N_Vector rr, N_Vector ewt, real hh,
                   real uround, DenseMat JJ, long int *nrePtr,
                   N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{

  real *yval;

  yval = N_VDATA(yy);

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
