/***********************************************************************
 * File:       sensirobx.c   
 * Written by: Steven L. Lee and Alan C. Hindmarsh
 * Version of: 1 September 2001
 *----------------------------------------------------------------------
 *
 * This simple example problem for SensIDA, due to Robertson, is from
 * chemical kinetics, and consists of the following three equations:
 *
 *      dy1/dt = -p1*y1 + p2*y2*y3
 *      dy2/dt =  p1*y1 - p2*y2*y3 - p3*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions
 * y1 = 1.0, y2 = y3 = 0.0, and parameters p1 = 0.04, p2 = 1.e4, p3 = 3.e7.
 * Note that this system of equations can be expressed in standard DAE form:
 *      F(t, y, y', p) = 0.
 * The time-derivatives for the sensitivity and/or state variables are 
 * initialized using IDACalcIC.
 *
 * The problem is solved with SensIDA using IDADENSE for the linear solver,
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
#include "sensida.h"

typedef struct {
  real *p;
} *UserData;

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
#define NY    3
#define NP    3
#define NS    2
#define ZERO  RCONST(0.0)

main()
{
  int itol, itask;
  integer retval, i, j, iout;
  long int iopt[OPT_SIZE];
  boole optIn;
  real ropt[OPT_SIZE], rtol, *yval, *ypval, *atval, *idval;
  real t0, t1, tout, tret;
  N_Vector yy, yp, avtol, id;
  N_Vector *yysub;
  void *mem;
  UserData data;
  int Ntotal;
  real *pbar;
  int *plist;
  real rhomax;
  machEnvType machEnv;
  real p1, p2, p3;

  Ntotal = (1+NS)*NY;
  machEnv = SVecInitSens(NY);

  /* Allocate N-vectors. */
  yy = N_VNew(Ntotal, machEnv); 
  yp = N_VNew(Ntotal, machEnv);
  avtol = N_VNew(Ntotal, machEnv);
  id = N_VNew(Ntotal, machEnv);

  /* Create pointers to subvectors */
  yysub = N_VSUB(yy);

  data = (UserData) malloc(sizeof *data);

  /* Store nominal parameter values in p */  
  data->p = (real *) malloc(NP * sizeof(real));
  data->p[0] = p1 = 0.04;
  data->p[1] = p2 = 1.e4;
  data->p[2] = p3 = 3.e7;
  
  /* Scaling factor for each sensitivity equation */
  pbar = (real *) malloc(NP * sizeof(real));
  pbar[0] = 0.04;
  pbar[1] = 1.e4;
  pbar[2] = 3.e7;

  /* Store ordering of parameters in plist */
  plist = (int *) malloc(NP * sizeof(int));
  plist[0] = 3;
  plist[1] = 1;
  plist[2] = 2;

  rhomax = 0.0;
  
  /* Set various input parameters and initialize y and y'. */
  yval  = N_VDATA(yy);
  ypval = N_VDATA(yp);
  atval = N_VDATA(avtol);
  idval = N_VDATA(id);
  t0   = 0.0;
  t1   = 0.4;
  itol = SV;  /* scalar relative tolerance, vector absolute tolerance. */
  rtol = 1.e-4;
  atval[0] = 1.e-6;
  atval[1] = 1.e-10;
  atval[2] = 1.e-6;
  yval[0] = 1.0;   /* initialize the state variables, y */
  yval[1] = 0.0;
  yval[2] = 0.0;
  ypval[0]  =  0.0; /* initial guess for y' */
  ypval[1]  =  0.0;
  ypval[2]  =  0.0;  

  /* Identify the differential and algebraic components for y and y' */
  idval[0] = 1.0;
  idval[1] = 1.0;
  idval[2] = 0.0;

  optIn = FALSE; 
  itask = NORMAL;

  /* Initialize sensitivity variables, s and s' */
  SensInitZero(yy, NS);
  SensInitZero(yp, NS);

  /* Set the absolute tolerances for s and s' */
  SensSetVecAtol(avtol, NS);

  /* Identify differential and algebraic components for s and s' */
  SensSetId(id, NS);

  /* Call SensIDAMalloc to set up problem memory.
  NULL arguments are: constraints and errfp, respectively. */

  mem = SensIDAMalloc(NY, NS, Ntotal, resrob, data, t0, yy, yp, 
		      itol, &rtol, avtol, id , NULL , NULL, optIn, 
		      iopt, ropt,  machEnv, data->p, pbar, plist, rhomax);
  mem = (IDAMem)mem;

  if (mem == NULL) {printf("SensIDAMalloc failed.\n"); return(1); }

  /* Call SensIDADense and set up the linear solver package. */
  retval = SensIDADense(mem, jacrob, NULL);

  if (retval != SUCCESS) {printf("SensIDADense failed.\n"); return(1); }

  /* Call IDACalcIC to calculate initial values for sensitivity derivatives */
  retval = IDACalcIC(mem, CALC_YA_YDP_INIT, t1, ZERO, 0, 0, 0, 0, ZERO);


  if (retval != SUCCESS) {
    printf("IDACalcIC failed. retval = %d\n", retval);
    return(1);
  }

  /* Print output heading */

  printf("sensrobx: Robertson kinetics DAE serial example problem for SensIDA \n");
  printf("          Three equation chemical kinetics problem. \n\n");
  printf("Linear solver: SENSIDADENSE, with user-supplied Jacobian.\n");
  printf("Tolerance parameters:  rtol = %g   atol = %g, %g, %g \n",
                                 rtol, atval[0], atval[1], atval[2]);
  printf("Initial conditions: y0 = (%g, %g, %g)\n", yval[0], yval[1], yval[2]);
  printf("Constraints not used.\n");
  printf("Derivatives initialized via IDACalcIC.\n");
  printf("Parameter values: p1 = %g, p2 = %g, p3 = %g\n", p1, p2, p3);
  printf("Number of sensitivity vectors: Ns  = %3d \n", NS);
  printf("s_k is the sensitivity of y with respect to parameter k\n"); 
  printf(".......................................................\n");

  /* Loop over tout values and call IDASolve. */

  for (tout = t1, iout = 1; iout <= NOUT ; iout++, tout *= 10.0) {
    retval=IDASolve(mem, tout, t0, &tret, yy, yp, itask);
    if (retval != SUCCESS) {printf("IDASolve returned %d\n",retval); return(1); }

    printf("\nt = %g\ny   = (%g, %g, %g)\n", tret, yval[0], yval[1], yval[2]);
    for (i = 1; i <= NS; i++) {
      j = (plist == NULL) ? i : plist[i-1];
      printf("s_%d = (%g, %g, %g)\n", j, yysub[i]->data[0]/pbar[j-1],
	     yysub[i]->data[1]/pbar[j-1], yysub[i]->data[2]/pbar[j-1]);
    }

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

  SensIDAFree(mem);
  N_VFree(yy);
  N_VFree(yp);
  N_VFree(avtol);
  N_VFree(id);

  if (plist != NULL) free(plist);
  free(pbar);
  free(data->p);
  free(data);

  SVecFreeSens(machEnv);

} /* End of sensrobx main. */


/**************************************************************************/

/*  Define the system residual function resrob. */

int resrob(integer Neq, real tres, N_Vector yy, N_Vector yp, N_Vector rr, 
           void *rdata)
{
  int i;
  real *yval, *ypval, *rval;
  UserData data;
  real p1, p2, p3;

  yval = N_VDATA(yy); ypval = N_VDATA(yp); rval = N_VDATA(rr);

  data = (UserData) rdata;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  rval[0]  = -p1*yval[0] + p2*yval[1]*yval[2];
  rval[1]  = -rval[0]    - p3*yval[1]*yval[1] - ypval[1];
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
  UserData data;
  real p1, p2, p3;

  yval = N_VDATA(yy);

  data = (UserData) rdata;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  IJth(JJ,1,1) = -p1 - cj;
  IJth(JJ,2,1) =  p1;
  IJth(JJ,3,1) =  1.  ;
  IJth(JJ,1,2) =  p2*yval[2];
  IJth(JJ,2,2) = -p2*yval[2] - 2.0*p3*yval[1] - cj;
  IJth(JJ,3,2) =  1.;
  IJth(JJ,1,3) =  p2*yval[1];
  IJth(JJ,2,3) = -p2*yval[1];
  IJth(JJ,3,3) =  1.;

  return(SUCCESS);

} /* End of the Jacobian function jacrob. */
