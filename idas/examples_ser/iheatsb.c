/***********************************************************************
 * File       : iheatsb.c   
 * Written by : Allan G. Taylor, Alan C. Hindmarsh, and Radu Serban
 * Version of : 23 July 2003
 *----------------------------------------------------------------------
 *
 * Example problem for IDAS: 2D heat equation, serial, banded. 
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the band solver IDABand, and IDACalcIC.
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
 * The system is solved with IDAS using the banded linear system solver,
 * half-bandwidths equal to M, and default difference-quotient Jacobian.
 * For purposes of illustration, IDACalcIC is called to compute correct
 * values at the boundary, given incorrect values as input initial guesses.
 * The constraints u >= 0 are posed for all components.
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 * (Output at t = 0 is for IDACalcIC cost statistics only.)
 *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "sundialsmath.h"
#include "nvector_serial.h" /* definitions of type N_Vector and macro NV_DATA_S */
#include "idas.h"
#include "idasband.h"

typedef struct {
  integertype    mm;
  realtype       dx;
  realtype       coeff;
} *UserData;


static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res);

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval, void *rdata);


#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.1)

int main()
{
  void *mem;
  UserData data;
  NV_Spec nvSpec;
  N_Vector uu, up, constraints, id, res;
  int ier, iout, itol, itask;
  integertype mu, ml;
  realtype rtol, atol, t0, t1, tout, tret, umax, hused;
  int kused, nst, nni, nje, nre, nreB, netf, ncfn;

  /* Initialize serial vector specification */
  nvSpec = NV_SpecInit_Serial(NEQ);

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew(nvSpec); 
  up = N_VNew(nvSpec);
  res = N_VNew(nvSpec);
  constraints = N_VNew(nvSpec);
  id = N_VNew(nvSpec);

  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  data->mm = MGRID;
  data->dx = ONE/(MGRID - ONE);
  data->coeff = ONE/( (data->dx) * (data->dx) );

  /* Initialize uu, up, id. */
  SetInitialProfile(data, uu, up, id, res);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = 0.01;
  itol = SS;    /* Specify scalar relative and absolute tolerance. */
  rtol = ZERO;
  atol = 1.0e-3;
  itask = NORMAL;

  /* Call IDACreate and IDAMalloc to initialize solution */
  mem = IDACreate();
  if (mem == NULL) { printf("IDACreate failed."); return(1); }
  ier = IDASetRdata(mem, data);
  if (ier != SUCCESS) { printf("IDASetRdata failed. "); return(1); }
  ier = IDASetID(mem, id);
  if (ier != SUCCESS) { printf("IDASetID failed. "); return(1); }
  ier = IDASetConstraints(mem, constraints);
  if (ier != SUCCESS) { printf("IDASetConstraints failed. "); return(1); }
  ier = IDAMalloc(mem, heatres, t0, uu, up, itol, &rtol, &atol, nvSpec);
  if (ier != SUCCESS) { printf("IDAMalloc failed. "); return(1); }

  /* Call IDABand to specify the linear solver. */
  mu = MGRID; ml = MGRID;
  ier = IDABand(mem, NEQ, mu, ml);
  if (ier != SUCCESS) { printf("IDABand failed."); return(1); }
 
  /* Call IDACalcIC to correct the initial values. */
  
  ier = IDACalcIC(mem, CALC_YA_YDP_INIT, t1);
  if (ier != SUCCESS) { printf("IDACalcIC failed. "); return(1); }

  /* Print output heading. */

  printf("iheatsb: Heat equation, serial example problem for IDAS \n");
  printf("         Discretized heat equation on 2D unit square. \n");
  printf("         Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("         Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: IDABAND, banded direct solver \n");
  printf("       difference quotient Jacobian, half-bandwidths = %d \n",MGRID);
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);

  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time       umax     k  nst  nni  nje   nre   nreB     h      \n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");

  umax = N_VMaxNorm(uu);
  
  IDAGetLastOrder(mem, &kused);
  IDAGetNumSteps(mem, &nst);
  IDAGetNumNonlinSolvIters(mem, &nni);
  IDAGetNumResEvals(mem, &nre);
  IDAGetLastStep(mem, &hused);
  IDABandGetNumJacEvals(mem, &nje);
  IDABandGetNumResEvals(mem, &nreB);

  printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e \n",
         t0, umax, kused, nst, nni, nje, nre, nreB, hused);

  /* Loop over output times, call IDASolve, and print results. */
  
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {
    
    ier = IDASolve(mem, tout, &tret, uu, up, itask);
    
    umax = N_VMaxNorm(uu);
    IDAGetLastOrder(mem, &kused);
    IDAGetNumSteps(mem, &nst);
    IDAGetNumNonlinSolvIters(mem, &nni);
    IDAGetNumResEvals(mem, &nre);
    IDAGetLastStep(mem, &hused);
    IDABandGetNumJacEvals(mem, &nje);
    IDABandGetNumResEvals(mem, &nreB);

    printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e \n",
           tret, umax, kused, nst, nni, nje, nre, nreB, hused);
    
    if (ier < 0) { printf("IDASolve failed."); return(1); }

  } /* End of tout loop. */
  
  /* Print remaining counters and free memory. */
  IDAGetNumErrTestFails(mem, &netf);
  IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  printf("\n netf = %d,   ncfn = %d \n", netf, ncfn);

  IDAFree(mem);
  N_VFree(uu);
  N_VFree(up);
  N_VFree(constraints);
  N_VFree(id);
  N_VFree(res);
  free(data);
  NV_SpecFree_Serial(nvSpec);

  return(0);

} /* End of iheatsb main program. */


/*************************************************************************
 * SetInitialProfile: routine to initialize u, up, and id vectors.       */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res)
{
  realtype xfact, yfact, *udata, *updata, *iddata;
  integertype mm, mm1, i, j, offset, loc;
  
  mm = data->mm;
  mm1 = mm - 1;
  
  udata = NV_DATA_S(uu);
  updata = NV_DATA_S(up);
  iddata = NV_DATA_S(id);

  /* Initialize id to 1's. */
  N_VConst(ONE, id);

  /* Initialize uu on all grid points. */ 
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
  heatres(ZERO, uu, up, res, data);
  
  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Finally, set values of u, up, and id at boundary points. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) {
        udata[loc] = BVAL; updata[loc] = ZERO; iddata[loc] = ZERO; }
    }
  }
  
  return(SUCCESS);

} /* End of SetInitialProfiles */


/*************************************************************************
 * heatres: heat equation system residual function                       *
 * This uses 5-point central differencing on the interior points, and    *
 * includes algebraic equations for the boundary values.                 *
 * So for each interior point, the residual component has the form       *
 *    res_i = u'_i - (central difference)_i                              *
 * while for each boundary point, it is res_i = u_i.                     */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval, 
            void *rdata)
{
  integertype mm, i, j, offset, loc;
  realtype *uv, *upv, *resv, coeff;
  UserData data;
  
  uv = NV_DATA_S(uu); upv = NV_DATA_S(up); resv = NV_DATA_S(resval);

  data = (UserData)rdata;
  mm = data->mm;
  coeff = data->coeff;
  
  /* Initialize resval to uu, to take care of boundary equations. */
  N_VScale(ONE, uu, resval);
  
  /* Loop over interior points; set res = up - (central difference). */
  for (j = 1; j < mm-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      resv[loc] = upv[loc] - coeff * 
           (uv[loc-1] + uv[loc+1] + uv[loc-mm] + uv[loc+mm] - 4.0*uv[loc]);
    }
  }
  
  return(SUCCESS);

} /* End of residual function heatres. */
