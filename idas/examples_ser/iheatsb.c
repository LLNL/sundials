/***********************************************************************
 * File       : iheatsb.c   
 * Written by : Allan G. Taylor, Alan C. Hindmarsh, and Radu Serban
 * Version of : 11 February 2004
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
#include "nvector_serial.h"
#include "idas.h"
#include "idasband.h"

typedef struct {
  long int    mm;
  realtype       dx;
  realtype       coeff;
} *UserData;

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res);

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval, void *rdata);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.1)

int main(void)
{
  void *mem;
  UserData data;
  NV_Spec nvSpec;
  N_Vector uu, up, constraints, id, res;
  int ier, iout, itol, itask;
  long int mu, ml;
  realtype rtol, atol, t0, t1, tout, tret, umax, hused;
  long int nst, nni, nje, nre, nreB, netf, ncfn;
  int kused;

  mem = NULL;
  data = NULL;
  nvSpec = NULL;
  uu = up = constraints = id = res = NULL;

  /* Initialize serial vector specification */
  nvSpec = NV_SpecInit_Serial(NEQ);
  if(check_flag((void *)nvSpec, "NV_SpecInit", 0)) return(1);

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew(nvSpec);
  if(check_flag((void *)uu, "N_VNew", 0)) return(1);
  up = N_VNew(nvSpec);
  if(check_flag((void *)up, "N_VNew", 0)) return(1);
  res = N_VNew(nvSpec);
  if(check_flag((void *)res, "N_VNew", 0)) return(1);
  constraints = N_VNew(nvSpec);
  if(check_flag((void *)constraints, "N_VNew", 0)) return(1);
  id = N_VNew(nvSpec);
  if(check_flag((void *)id, "N_VNew", 0)) return(1);

  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2)) return(1);
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
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);
  ier = IDASetRdata(mem, data);
  if(check_flag(&ier, "IDASetRdata", 1)) return(1);
  ier = IDASetId(mem, id);
  if(check_flag(&ier, "IDASetId", 1)) return(1);
  ier = IDASetConstraints(mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1)) return(1);
  ier = IDAMalloc(mem, heatres, t0, uu, up, itol, &rtol, &atol, nvSpec);
  if(check_flag(&ier, "IDAMalloc", 1)) return(1);

  /* Call IDABand to specify the linear solver. */
  mu = MGRID; ml = MGRID;
  ier = IDABand(mem, NEQ, mu, ml);
  if(check_flag(&ier, "IDABand", 1)) return(1);
 
  /* Call IDACalcIC to correct the initial values. */
  
  ier = IDACalcIC(mem, CALC_YA_YDP_INIT, t1);
  if(check_flag(&ier, "IDACalcIC", 1)) return(1);

  /* Print output heading. */

  printf("iheatsb: Heat equation, serial example problem for IDA \n");
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
  
  ier = IDAGetLastOrder(mem, &kused);
  check_flag(&ier, "IDAGetLastOrder", 1);
  ier = IDAGetNumSteps(mem, &nst);
  check_flag(&ier, "IDAGetNumSteps", 1);
  ier = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&ier, "IDAGetNumNonlinSolvIters", 1);
  ier = IDAGetNumResEvals(mem, &nre);
  check_flag(&ier, "IDAGetNumResEvals", 1);
  ier = IDAGetLastStep(mem, &hused);
  check_flag(&ier, "IDAGetLastStep", 1);
  ier = IDABandGetNumJacEvals(mem, &nje);
  check_flag(&ier, "IDABandGetNumJacEvals", 1);
  ier = IDABandGetNumResEvals(mem, &nreB);
  check_flag(&ier, "IDABandGetNumResEvals", 1);

  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n",
         t0, umax, kused, nst, nni, nje, nre, nreB, hused);

  /* Loop over output times, call IDASolve, and print results. */
  
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {
    
    ier = IDASolve(mem, tout, &tret, uu, up, itask);
    if(check_flag(&ier, "IDASolve", 1)) return(1);

    umax = N_VMaxNorm(uu);
    ier = IDAGetLastOrder(mem, &kused);
    check_flag(&ier, "IDAGetLastOrder", 1);
    ier = IDAGetNumSteps(mem, &nst);
    check_flag(&ier, "IDAGetNumSteps", 1);
    ier = IDAGetNumNonlinSolvIters(mem, &nni);
    check_flag(&ier, "IDAGetNumNonlinSolvIters", 1);
    ier = IDAGetNumResEvals(mem, &nre);
    check_flag(&ier, "IDAGetNumResEvals", 1);
    ier = IDAGetLastStep(mem, &hused);
    check_flag(&ier, "IDAGetLastStep", 1);
    ier = IDABandGetNumJacEvals(mem, &nje);
    check_flag(&ier, "IDABandGetNumJacEvals", 1);
    ier = IDABandGetNumResEvals(mem, &nreB);
    check_flag(&ier, "IDABandGetNumResEvals", 1);

    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n",
           tret, umax, kused, nst, nni, nje, nre, nreB, hused);

  } /* End of tout loop. */
  
  /* Print remaining counters and free memory. */
  ier = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&ier, "IDAGetNumErrTestFails", 1);
  ier = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&ier, "IDAGetNumNonlinSolvConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);

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
  long int mm, mm1, i, j, offset, loc;
  
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
  long int mm, i, j, offset, loc;
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

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if flag == SUCCESS
     opt == 2 means function allocates memory so check if returned NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  /* Check if flag != SUCCESS */
  else if (opt == 1) {
    errflag = flagvalue;
    if (*errflag != SUCCESS) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  return(0);
}
