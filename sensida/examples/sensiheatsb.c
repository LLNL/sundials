/***********************************************************************
 * File:       sensiheatsb.c   
 * Written by: Steven L. Lee and Alan C. Hindmarsh
 * Version of: 20 March 2001
 *----------------------------------------------------------------------
 *
 * Example problem for SensIDA: 2D heat equation, serial, banded. 
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
 * The system is solved with IDA using the banded linear system solver,
 * half-bandwidths equal to M, and default difference-quotient Jacobian.
 * For purposes of illustration, IDACalcIC is called to compute correct
 * values at the boundary, given incorrect values as input initial guesses.
 * The constraints u >= 0 are posed for all components.
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 * (Output at t = 0 is for IDACalcIC cost statistics only.)
 *
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "llnltyps.h"
#include "llnlmath.h"
#include "nvector.h"
#include "ida.h"
#include "idaband.h"
#include "sensida.h"

typedef struct {
  real       *p;
  integer    ny;
  integer    mm;
  real       dx;
  real       coeff;
} *UserData;


static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res);

int heatres(integer Neq, real tres, N_Vector uu, N_Vector up,
            N_Vector resval, void *rdata);


#define NOUT  11
#define MGRID 10
#define NY    MGRID*MGRID
#define NP    2
#define NS    2
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.1)

main()
{
  integer retval, i, j, iout, itol, itask, mu, ml;
  long int iopt[OPT_SIZE];
  boole optIn;
  real ropt[OPT_SIZE], rtol, atol, t0, t1, tout, tret, umax;
  N_Vector uu, up, constraints, id, res;
  UserData data;
  void *mem;
  N_Vector *uusub, *upsub;
  N_Vector *constraintssub, *idsub, *ressub;
  integer Ntotal;
  real *pbar;
  real rhomax;
  int *plist;
  machEnvType machEnv;

  Ntotal = (1+NS)*NY;
  machEnv = SVecInitSens(NY);

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew(Ntotal, machEnv); 
  up = N_VNew(Ntotal, machEnv);
  res = N_VNew(Ntotal, machEnv);
  constraints = N_VNew(Ntotal, machEnv);
  id = N_VNew(Ntotal, machEnv);

  /* Create pointers to the subvectors */
  uusub = N_VSUB(uu);
  upsub = N_VSUB(up);
  constraintssub = N_VSUB(constraints);
  idsub = N_VSUB(id);
  ressub = N_VSUB(res);

  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  data->ny = NY;
  data->mm = MGRID;
  data->dx = ONE/(MGRID - ONE);
  data->coeff = ONE/( (data->dx) * (data->dx) );

  /* Store nominal parameter values in p */
  data->p = (real *) malloc(NP * sizeof(real));
  data->p[0] = 1.0;
  data->p[1] = 1.0;

  /* Scaling factor for each sensitivity equation */
  pbar = (real *) malloc(NP * sizeof(real));
  pbar[0] = 1.0;
  pbar[1] = 1.0;

  /* Store ordering of parameters in plist */
  plist = (int *) malloc(NP * sizeof(int));
  plist[0] = 1;
  plist[1] = 2;

  rhomax = 0.0;

  /* Initialize uu, up, id. */
  SetInitialProfile(data, uusub[0], upsub[0], idsub[0], ressub[0]);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst (ONE, constraintssub[0]);

  /* Initialize the sensitivity variables */
  SensInitZero(uu, NS);
  SensInitZero(up, NS);
  SensInitZero(constraints, NS);

  /* Identify all sensitivity variables as differential variables */
  SensSetId(id, NS);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = 0.01;
  itol = SS;    /* Specify scalar relative and absolute tolerance. */
  rtol = ZERO;
  atol = 1.0e-3;
  optIn = FALSE; 
  itask = NORMAL;

  /* Call IDAMalloc to initialize solution.
     NULL arguments are errfp and machEnv, respectively. */
  /*
  mem = IDAMalloc(NY, heatres, data, t0, uu, up, itol, &rtol, &atol,
                  id, constraints, NULL, optIn, iopt, ropt,  NULL);
		  */
  mem = SensIDAMalloc(NY, NS, Ntotal, heatres, data, t0, uu, up, 
		      itol, &rtol, &atol, id, constraints, NULL, optIn, 
		      iopt, ropt, machEnv, data->p, pbar, plist, rhomax);  
  mem = (IDAMem)mem;
  if (mem == NULL) { printf("IDAMalloc failed."); return(1); }


  /* Call SensIDABand to specify the linear solver. */
  mu = MGRID; ml = MGRID;
  retval = SensIDABand(mem, mu, ml, NULL, NULL);
  if (retval != SUCCESS) { printf("IDABand failed."); return(1); }
 
  /* Call IDACalcIC to correct the initial values. */
  retval = IDACalcIC(mem, CALC_YA_YDP_INIT, ONE, ZERO, 0, 0, 0, 0, ZERO);
  if (retval != SUCCESS) { printf("IDACalcIC failed."); return(1); }

  /* Print output heading. */

  printf("sensheatsb: Heat equation, serial example problem for SensIDA \n");
  printf("            Discretized heat equation on 2D unit square. \n");
  printf("            Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("            Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("          Total system size: %d\n\n", NY);
  printf("Number of sensitivities: Ns = %d\n", NS);
  printf("Parameter values:       p_1 = %9.2e,    p_2 = %9.2e\n", 
	 data->p[0], data->p[1]);
  printf("Scale factors:       pbar_1 = %9.2e, pbar_2 = %9.2e\n", 
	 pbar[0], pbar[1]);
  printf("Finite difference:   rhomax = %g\n", rhomax);
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: IDABAND, banded direct solver \n");
  printf("       difference quotient Jacobian, half-bandwidths = %d \n",MGRID);
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);

  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary    max(u) = max-norm of solution u \n\n");
  printf("  time     max(u)     k  nst  nni  nje   nre    h      \n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");
  umax = N_VMaxNorm(uusub[0]);

  printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %9.2e \n\n",
         t0, umax, iopt[KUSED], iopt[NST], iopt[NNI], 
         iopt[BAND_NJE], iopt[NRE], ropt[HUSED]);

  for (i = 1; i <= NS; i++){
    umax = N_VMaxNorm(uusub[i]);
    j = (plist == NULL) ? i : plist[i-1];
    printf("sensitivity s_%d:  max. norm = %11.5e\n", j, umax/pbar[j-1]);
    if (i == NS) printf("\n");
  }

  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    retval = IDASolve(mem, tout, t0, &tret, uu, up, itask);

    umax = N_VMaxNorm(uusub[0]);
    printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %9.2e \n\n",
           tret, umax, iopt[KUSED], iopt[NST], iopt[NNI], 
           iopt[BAND_NJE], iopt[NRE], ropt[HUSED]);

    for (i = 1; i <= NS; i++){
      umax = N_VMaxNorm(uusub[i]);
      j = (plist == NULL) ? i: plist[i-1];
      printf("sensitivity s_%d:  max. norm = %11.5e\n", j, umax/pbar[j-1]);
      if (i == NS) printf("\n");
    }

    if (retval < 0) { printf("IDASolve failed."); return(1); }

  } /* End of tout loop. */

  /* Print remaining counters and free memory. */
  printf("\n netf = %d,   ncfn = %d \n", iopt[NETF], iopt[NCFN]);
  /*
  IDAFree(mem);
  */
  SensIDAFree(mem);
  N_VFree(uu);
  N_VFree(up);
  N_VFree(constraints);
  N_VFree(id);
  N_VFree(res);

  if (plist != NULL) free(plist);
  free(pbar);
  free(data->p);
  free(data);

  SVecFreeSens(machEnv);

} /* End of heatsb main program. */


/*************************************************************************
 * SetInitialProfile: routine to initialize u, up, and id vectors.       */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                     N_Vector id, N_Vector res)
{
  real xfact, yfact, *udata, *updata, *iddata;
  integer mm, mm1, i, j, offset, loc;

  mm = data->mm;

  udata = N_VDATA(uu);
  updata = N_VDATA(up);
  iddata = N_VDATA(id);

  /* Initialize id to 1's. */
  N_VConst (ONE, id);

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
  N_VConst (ZERO, up);

  /* heatres sets res to negative of ODE RHS values at interior points. */
  heatres (data->ny, ZERO, uu, up, res, data);

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

int heatres(integer Neq, real tres, N_Vector uu, N_Vector up,
              N_Vector resval, void *rdata)
{
  integer mm, i, j, offset, loc;
  real *uv, *upv, *resv, coeff;
  UserData data;
  real p1, p2;

  uv = N_VDATA(uu); upv = N_VDATA(up); resv = N_VDATA(resval);

  data = (UserData)rdata;
  p1 = data->p[0];
  p2 = data->p[1];

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
           (p1*(uv[loc-1] + uv[loc+1]) + p2*(uv[loc-mm] + uv[loc+mm]) 
	    - 2.0*(p1+p2)*uv[loc]);
    }
  }

  return(SUCCESS);

} /* End of residual function heatres. */
