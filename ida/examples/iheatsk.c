/***********************************************************************
 * File       : iheatsk.c   
 * Written by : Allan G. Taylor and Alan C. Hindmarsh
 * Version of : 8 March 2002
 *----------------------------------------------------------------------
 * Modified by R. Serban to work with new serial nvector (8/3/2002)    
 *----------------------------------------------------------------------
 *
 * Example problem for IDA: 2D heat equation, serial, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver IDASpgmr.
 *
 * The DAE system solved is a spatial discretization of the PDE 
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square.  The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MxM grid.
 * The values of u at the interior points satisfy ODEs, and equations
 * u = 0 at the boundaries are appended, to form a DAE system of size
 * N = M^2.  Here M = 10.
 *
 * The system is solved with IDA using the Krylov linear solver IDASPGMR. 
 * The preconditioner uses the diagonal elements of the Jacobian only.
 * Routines for preconditioning, required by IDASPGMR, are supplied here.
 * The constraints u >= 0 are posed for all components.
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"
#include "llnlmath.h"
#include "nvector_serial.h" /* definitions of type N_Vector, macro NV_DATA_S   */
#include "ida.h"
#include "idaspgmr.h"
#include "iterativ.h"


typedef struct {  
  integer    neq;
  integer    mm; /* mm is the number of points in the grid. */
  real       dx;
  real       coeff;
  N_Vector   pp; /* pp is the vector of diagonal preconditioner elements. */
} *UserData;


static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res);

int heatres(integer Neq, real tres, N_Vector uu, N_Vector up,
            N_Vector resval, void *rdata);

int PrecondHeateq(integer Neq, real tt, N_Vector uu,
                  N_Vector up, N_Vector rr, real cj,
                  ResFn res, void *rdata, void *pdata,
                  N_Vector ewt, N_Vector constraints, real hh, 
                  real uround, long int *nrePtr,
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

int PSolveHeateq(integer Neq, real tt, N_Vector uu,
                 N_Vector up, N_Vector rr, real cj, ResFn res, void *rdata,
                 void *pdata, N_Vector ewt, real delta, N_Vector rvec,
                 N_Vector zvec, long int *nrePtr, N_Vector tempv);

#define NOUT  11
#define MGRID 10
#define NEQ   (MGRID*MGRID)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define FOUR  RCONST(4.0)

int main()
{
  M_Env machEnv;
  integer retval, iout, itol, itask;
  long int iopt[OPT_SIZE];
  boole optIn;
  real ropt[OPT_SIZE], rtol, atol, t0, t1, tout, tret, umax;
  N_Vector uu, up, constraints, id, res;
  UserData data;
  void *mem;

  /* Initialize serial machine environment */
  machEnv = M_EnvInit_Serial(NEQ);  

  /* Allocate N-vectors and the user data structure. */

  uu = N_VNew(NEQ, machEnv); 
  up = N_VNew(NEQ, machEnv);
  res = N_VNew(NEQ, machEnv);
  constraints = N_VNew(NEQ, machEnv);
  id = N_VNew(NEQ, machEnv);

  data = (UserData) malloc(sizeof *data);

  /* Assign parameters in the user data structure. */

  data->neq = NEQ;
  data->mm  = MGRID;
  data->dx = ONE/(MGRID-ONE);
  data->coeff = ONE/(data->dx * data->dx);
  data->pp = N_VNew(NEQ, machEnv);

  /* Initialize uu, up, id. */
  SetInitialProfile(data, uu, up, id, res);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  /* Assign various parameters. */
  t0   = ZERO;
  t1   = 0.01;
  itol = SS; /* Specify scalar relative and absolute tolerance. */
  rtol = ZERO;
  atol = 1.e-3; 
  optIn = FALSE; 
  itask = NORMAL;

  /* Call IDAMalloc to initialize solution.
     NULL arguments are errfp and machEnv, respectively. */

  mem = IDAMalloc(NEQ, heatres, data , t0, uu, up, itol, &rtol, &atol,
                  id , constraints , NULL, optIn, iopt, ropt,  machEnv);
  /*mem = (IDAMem)mem;*/
  if (mem == NULL) { printf("IDAMalloc failed."); return(1); }

  /* Call IDASpgmr to specify the linear solver. */

  retval = IDASpgmr(mem, PrecondHeateq, PSolveHeateq,  MODIFIED_GS, 
                    0, 0, 0., 0., data);
  if (retval != SUCCESS) {printf("IDASpgmr failed."); return(1); }
  
  /* Print output heading. */
  
  printf("iheatsk: Heat equation, serial example problem for IDA \n");
  printf("         Discretized heat equation on 2D unit square. \n");
  printf("         Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("         Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: IDASPGMR, preconditioner using diagonal elements. \n");

  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nli   nre     h       npe nps\n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");
  umax = N_VMaxNorm(uu);

  printf(" %5.2f %13.5e  %ld  %3ld  %3ld  %3ld  %4ld  %9.2e  %3ld %3ld\n",
         t0, umax, iopt[KUSED], iopt[NST], iopt[NNI], iopt[SPGMR_NLI], 
         iopt[NRE],  ropt[HUSED], iopt[SPGMR_NPE], iopt[SPGMR_NPS]);

  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1,iout = 1; iout <= NOUT ; iout++, tout *= TWO) {
    
    retval = IDASolve(mem, tout, t0, &tret, uu, up, itask);
    
    umax = N_VMaxNorm(uu);
    printf(" %5.2f %13.5e  %ld  %3ld  %3ld  %3ld  %4ld  %9.2e  %3ld %3ld\n",
           tret, umax, iopt[KUSED], iopt[NST], iopt[NNI], iopt[SPGMR_NLI], 
           iopt[NRE],  ropt[HUSED], iopt[SPGMR_NPE], iopt[SPGMR_NPS]);

    if (retval < 0) {printf("IDASolve returned %ld.\n",retval); return(1); }

    } /* End of tout loop. */

  /* Print remaining counters and free memory. */

  printf("\n netf = %ld,   ncfn = %ld,   ncfl = %ld \n", 
         iopt[NETF], iopt[NCFN], iopt[SPGMR_NCFL]);

  IDAFree(mem);
  N_VFree(uu);
  N_VFree(up);
  N_VFree(constraints);
  N_VFree(id);
  N_VFree(res);
  N_VFree(data->pp);
  free(data);
  M_EnvFree_Serial(machEnv);

  return(0);

} /* End of iheatsk main program. */


/*************************************************************************
 * SetInitialProfile: routine to initialize u, up, and id vectors.       */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res)
{
  integer mm, mm1, i, j, offset, loc;
  real xfact, yfact, *udata, *updata, *iddata;

  mm = data->mm;

  udata = NV_DATA_S(uu);
  updata = NV_DATA_S(up);
  iddata = NV_DATA_S(id);

  /* Initialize id to 1's. */
  N_VConst(ONE, id);

  /* Initialize uu on all grid points. */ 
  mm1 = mm - 1;
  for (j = 0; j < mm; j++) {
    yfact = data->dx * j;
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      xfact = data->dx * i;
      loc = offset + i;
      udata[loc] = 16. * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }
  
  /* Initialize up vector to 0. */
  N_VConst(ZERO, up);

  /* heatres sets res to negative of ODE RHS values at interior points. */
  heatres(data->neq, ZERO, uu, up, res, data);

  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Set up and id at boundary points to zero. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0; i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) {
        updata[loc] = ZERO; iddata[loc] = ZERO; }
    }
  }
  
  return(SUCCESS);
  
} /* End of SetInitialProfiles. */


/* Routines called by IDASPGMR. */

/*************************************************************************
 * heatres: heat equation system residual function (user-supplied)       *
 * This uses 5-point central differencing on the interior points, and    *
 * includes algebraic equations for the boundary values.                 *
 * So for each interior point, the residual component has the form       *
 *    res_i = u'_i - (central difference)_i                              *
 * while for each boundary point, it is res_i = u_i.                     */

int heatres(integer Neq, real tres, N_Vector uu, N_Vector up,
            N_Vector res, void *rdata)
{
  integer i, j, offset, loc, mm;
  real *uv, *upv, *resv, coeff;
  UserData data;
  
  uv = NV_DATA_S(uu); upv = NV_DATA_S(up); resv = NV_DATA_S(res);

  data = (UserData)rdata;
  
  coeff = data->coeff;
  mm    = data->mm;
  
  /* Initialize res to uu, to take care of boundary equations. */
  N_VScale(ONE, uu, res);
  
  /* Loop over interior points; set res = up - (central difference). */
  for (j = 1; j < MGRID-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      resv[loc]= upv[loc] - coeff * 
       (uv[loc-1] + uv[loc+1] + uv[loc-mm] + uv[loc+mm] - 4*uv[loc]);
    }
  }

  return(SUCCESS);

 } /* End of residual function heatres. */


/*******************************************************************
 * PrecondHeateq: setup for diagonal preconditioner for iheatsk.   *
 *                                                                 *
 * The optional user-supplied functions PrecondHeateq and          *
 * PSolveHeateq together must define the left preconditoner        *
 * matrix P approximating the system Jacobian matrix               *
 *                   J = dF/du + cj*dF/du'                         *
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   *
 * systems P z = r.   This is done in this case by keeping only    *
 * the diagonal elements of the J matrix above, storing them as    *
 * inverses in a vector pp, when computed in PrecondHeateq, for    *
 * subsequent use in PSolveHeateq.                                 *
 *                                                                 *
 * In this instance, only cj and data (user data structure, with   * 
 * pp etc.) are used from the PrecondHeateq argument list.         */
  
int PrecondHeateq(integer Neq, real tt, N_Vector uu,
                  N_Vector up, N_Vector rr, real cj,
                  ResFn res, void *rdata, void *pdata,
                  N_Vector ewt, N_Vector constraints, real hh, 
                  real uround, long int *nrePtr,
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  
  integer i, j, offset, loc, mm;
  real *ppv, pelinv;
  UserData data;
  
  data = (UserData) rdata;
  ppv = NV_DATA_S(data->pp);
  mm = data->mm;

  /* Initialize the entire vector to 1., then set the interior points to the
     correct value for preconditioning. */
  
  N_VConst(ONE,data->pp);
  
  /* Compute the inverse of the preconditioner diagonal elements. */
  pelinv = ONE/(cj + FOUR*data->coeff); 
  
  for (j = 1; j < mm-1; j++) {
    offset = mm * j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      ppv[loc] = pelinv;
    }
  }
  
  return(SUCCESS);
  
} /* End of PrecondHeateq. */


/******************************************************************
 * PSolveHeateq: solve preconditioner linear system.              *
 * This routine multiplies the input vector rvec by the vector pp *
 * containing the inverse diagonal Jacobian elements (previously  *
 * computed in PrecondHeateq), returning the result in zvec.      */

int PSolveHeateq(integer Neq, real tt, N_Vector uu,
                 N_Vector up, N_Vector rr, real cj, ResFn res, void *rdata,
                 void *pdata, N_Vector ewt, real delta, N_Vector rvec,
                 N_Vector zvec, long int *nrePtr, N_Vector tempv)
{
  UserData data;
  
  data = (UserData) rdata;
  
  N_VProd(data->pp, rvec, zvec);
  
  return(SUCCESS);

} /* End of PSolveHeateq. */
