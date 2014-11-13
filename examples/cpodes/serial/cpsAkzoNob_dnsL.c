/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Chemical Akzo-Nobel DAE example
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_lapack.h>
#include <nvector/nvector_serial.h>

#include <sundials/sundials_math.h>

/* Problem Constants */

#define N 6

#define T0 RCONST(0.0)
#define TF RCONST(180.0)

#define RTOL RCONST(1.0e-4)
#define ATOL RCONST(1.0e-4)

#define ZERO  RCONST(0.0)
#define HALF  RCONST(0.5)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define THREE RCONST(3.0)
#define FOUR  RCONST(4.0)

typedef struct {
  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;
} *UserData;

static int f(realtype t, N_Vector y, N_Vector yd, void *f_data);
static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data);

static void PrintFinalStats(void *cpode_mem);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;
  void *cpode_mem;
  realtype reltol, abstol;
  N_Vector yy, yp, ctols;
  realtype t;
  int flag;


  /* USER DATA */

  data = (UserData) malloc(sizeof *data);
  data->k1 = RCONST(18.7);
  data->k2 = RCONST(0.58);
  data->k3 = RCONST(0.09);
  data->k4 = RCONST(0.42);
  data->K = RCONST(34.4);
  data->klA = RCONST(3.3);
  data->Ks = RCONST(115.83);
  data->pCO2 = RCONST(0.9);
  data->H = RCONST(737.0);

  /* Set the tolerances for the forward integration */
  reltol = RTOL;
  abstol = ATOL;

  /* Set projection tolerances */
  ctols = N_VNew_Serial(1);
  NV_Ith_S(ctols,0) = 1.0e-10;

  /* Initial conditions */
  yy = N_VNew_Serial(N);
  yp = N_VNew_Serial(N);
  NV_Ith_S(yy,0) = RCONST(0.444);
  NV_Ith_S(yy,1) = RCONST(0.00123);
  NV_Ith_S(yy,2) = RCONST(0.0);
  NV_Ith_S(yy,3) = RCONST(0.007);
  NV_Ith_S(yy,4) = RCONST(0.0);
  NV_Ith_S(yy,5) = data->Ks * RCONST(0.444) * RCONST(0.007);

  /* PROBLEM SOLUTION */

  printf("Problem solution\n");

  /* Create and allocate CPODES memory */
  cpode_mem = CPodeCreate(CP_BDF, CP_NEWTON);
  flag = CPodeSetUserData(cpode_mem, data);
  flag = CPodeInitExpl(cpode_mem, f, T0, yy);
  flag = CPodeSStolerances(cpode_mem, reltol, abstol);
  flag = CPLapackDense(cpode_mem, N);

  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, ctols);
  flag = CPodeSetProjTestCnstr(cpode_mem, TRUE);
  flag = CPLapackDenseProj(cpode_mem, 1, N, CPDLS_QRP);

  /* Integrate to TF */

  t = T0;
  while(t < TF) {
    flag = CPode(cpode_mem, TF, &t, yy, yp, CP_ONE_STEP);
    printf("%g  %g %g %g %g %g %g\n", t,
           NV_Ith_S(yy,0),
           NV_Ith_S(yy,1),
           NV_Ith_S(yy,2),
           NV_Ith_S(yy,3),
           NV_Ith_S(yy,4),
           NV_Ith_S(yy,5));
  }

  PrintFinalStats(cpode_mem);

  /* FREE MEMORY */

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(ctols);

  CPodeFree(&cpode_mem);

  free(data);

  return(0);
}


static int f(realtype t, N_Vector yy, N_Vector yd, void *f_data)
{
  UserData data;
  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;

  realtype y1, y2, y3, y4, y5, y6;

  realtype r1, r2, r3, r4, r5, Fin;
  realtype f1, f4;

  data = (UserData) f_data;
  k1 = data->k1;
  k2 = data->k2;
  k3 = data->k3;
  k4 = data->k4;
  K = data->K;
  klA = data->klA;
  Ks = data->Ks;
  pCO2 = data->pCO2;
  H = data->H;

  y1 = NV_Ith_S(yy,0);
  y2 = NV_Ith_S(yy,1);
  y3 = NV_Ith_S(yy,2);
  y4 = NV_Ith_S(yy,3);
  y5 = NV_Ith_S(yy,4);
  y6 = NV_Ith_S(yy,5);


  r1 = k1 * RPowerI(y1,4) * SUN_SQRT(y2);
  r2 = k2 * y3 * y4;
  r3 = k2/K * y1 * y5;
  r4 = k3 * y1 * y4 * y4;
  r5 = k4 * y6 * y6 * SUN_SQRT(y2);
  Fin = klA * ( pCO2/H - y2 );

  NV_Ith_S(yd,0) = f1 = -TWO*r1 + r2 - r3 - r4;
  NV_Ith_S(yd,1) = -HALF*r1 - r4 - HALF*r5 + Fin;
  NV_Ith_S(yd,2) = r1 - r2 + r3;
  NV_Ith_S(yd,3) = f4 = -r2 + r3 - TWO*r4;
  NV_Ith_S(yd,4) = r2 - r3 + r5;
  NV_Ith_S(yd,5) = Ks * ( y4*f1 + y1*f4 );

  return(0);
}

static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data)
{
  UserData data;
  realtype Ks;
  realtype y1, y4, y6;
  

  data = (UserData) c_data;
  Ks = data->Ks;

  y1 = NV_Ith_S(yy,0);
  y4 = NV_Ith_S(yy,3);
  y6 = NV_Ith_S(yy,5);

  NV_Ith_S(cout,0) = Ks * y1 * y4 - y6;

  return(0);
}



static void PrintFinalStats(void *cpode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  long int nproj, nce, nsetupsP, nprf;
  int flag;

  flag = CPodeGetActualInitStep(cpode_mem, &h0u);
  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  flag = CPDlsGetNumJacEvals(cpode_mem, &nje);
  flag = CPDlsGetNumFctEvals(cpode_mem, &nfeLS);

  flag = CPodeGetProjStats(cpode_mem, &nproj, &nce, &nsetupsP, &nprf);

  flag = CPodeGetNumGEvals(cpode_mem, &nge);

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld\n",
	 nst, nfe, nsetups);
  printf("nfeLS = %-6ld nje = %ld\n",
	 nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld \n",
	 nni, ncfn, netf);
  printf("nproj = %-6ld nce = %-6ld nsetupsP = %-6ld nprf = %-6ld\n",
         nproj, nce, nsetupsP, nprf);
  printf("nge = %ld\n", nge);
}
