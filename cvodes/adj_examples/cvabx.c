/************************************************************************
 * File       : cvabx.c                                                 *
 * Programmers: Radu Serban @ LLNL                                      *
 * Version of : 22 March 2002                                           *
 *----------------------------------------------------------------------*
 * Adjoint sensitivity example problem.                                 *
 * The following is a simple example problem with a banded Jacobian,    *
 * with the program for its solution by CVODES.                         *
 * The problem is the semi-discrete form of the advection-diffusion     *
 * equation in 2-D:                                                     *
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2                     *
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time              *
 * interval 0 <= t <= 1.  Homogeneous Dirichlet boundary conditions     *
 * are posed, and the initial condition is                              *
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy) .                                *
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with           *
 * central differencing, and with boundary values eliminated,           *
 * leaving an ODE system of size NEQ = MX*MY.                           *
 * This program solves the problem with the BDF method, Newton          *
 * iteration with the CVODE band linear solver, and a user-supplied     *
 * Jacobian routine.                                                    * 
 * It uses scalar relative and absolute tolerances.                     *
 * Output is printed at t = .1, .2, ..., 1.                             *
 * Run statistics (optional outputs) are printed at the end.            *
 *                                                                      *
 * Additionally, CVODES can integrate backwards in time the             *
 * the semi-discrete form of the adjoint PDE:                           *
 *   d(lambda)/dt = - d^2(lambda) / dx^2 + 0.5 d(lambda) / dx           *
 *                  - d^2(lambda) / dy^2 - 1.0                          *
 * with homogeneous Dirichlet boundary conditions and final conditions  *
 *   lambda(x,y,t=t_final) = 0.0                                        *
 * whose solution at t = 0 represents the sensitivity of                *
 *   G = int_0^t_final int_x int _y u(t,x,y) dx dy dt                   *
 * with respect to the initial conditions of the original problem.      *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"
#include "cvodea.h"
#include "cvsband.h"
#include "nvector_serial.h"
#include "band.h"

/* Problem Constants */

#define XMAX  2.0          /* domain boundaries */
#define YMAX  1.0
#define MX    40           /* mesh dimensions   */
#define MY    20
#define NEQ   MX*MY        /* number of equations       */
#define ATOL  1.e-5        /* scalar absolute tolerance */
#define RTOLB 1.e-6
#define T0    0.0          /* initial time              */
#define T1    0.1          /* first output time         */
#define DTOUT 0.1          /* output time increment     */
#define NOUT  10           /* number of output times    */
#define TOUT  1.0          /* final time */

/* Check point saved every NSTEP */
#define NSTEP 50

/* User-defined vector access macro IJth */

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. 
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the macro call vdata = N_VDATA(v),
   where v is an N_Vector. 
   The variables are ordered by the y index j, then by the x index i. */

#define IJth(vdata,i,j) (vdata[(j-1) + (i-1)*MY])

/* Type : UserData 
   contains grid constants */

typedef struct {
  real dx, dy, hdcoef, hacoef, vdcoef;
} *UserData;


/* Private Helper Functions */

static void SetIC(N_Vector u, UserData data);
static void WriteLambda(N_Vector uB);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data);

static void Jac(integer N, integer mu, integer ml, BandMat J, RhsFn f,
                void *f_data, real t, N_Vector u, N_Vector fu, N_Vector ewt,
                real h, real uround, void *jac_data, long int *nfePtr,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 

static void fB(integer NB, real tB, N_Vector u, 
               N_Vector uB, N_Vector uBdot, void *f_dataB);

static void JacB(integer NB, integer muB, integer mlB, BandMat JB, RhsFnB fB,
                 void *f_dataB, real tB, N_Vector u, 
                 N_Vector uB, N_Vector fuB, N_Vector ewtB,
                 real hB, real uroundB, void *jac_dataB, long int *nfePtrB,
                 N_Vector vtemp1B, N_Vector vtemp2B, N_Vector vtemp3B); 


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  M_Env machEnvF, machEnvB;

  UserData data;

  void *cvadj_mem;
  void *cvode_mem;

  real ropt[OPT_SIZE], dx, dy, reltol, abstol, t;
  long int iopt[OPT_SIZE];
  N_Vector u;

  real reltolB, abstolB;
  N_Vector uB;
  
  int flag, ncheck;

  /*--------------------*/
  /*----- USER DATA ----*/
  /*--------------------*/

  /* Allocate data memory */
  data = (UserData) malloc(sizeof *data);
  dx = data->dx = XMAX/(MX+1);     /* Set grid coefficients in data */
  dy = data->dy = YMAX/(MY+1);
  data->hdcoef = 1.0/(dx*dx);
  data->hacoef = 1.5/(2.0*dx);
  data->vdcoef = 1.0/(dy*dy);

  /*--------------------------*/
  /*---- ORIGINAL PROBLEM ----*/
  /*--------------------------*/

  /* Initialize serial machine environment for forward run */
  machEnvF = M_EnvInit_Serial(NEQ);

  /* Set the tolerances */
  reltol = 0.0;
  abstol = ATOL;

  /* Allocate u vector */
  u = N_VNew(NEQ, machEnvF);
  /* Initialize u vector */
  SetIC(u, data);

  /* Allocate CVODE memory for forward run */
  printf("\nAllocate CVODE memory for forward runs\n");
  cvode_mem = CVodeMalloc(NEQ, f, T0, u, BDF, NEWTON, SS, &reltol, &abstol,
                          data, NULL, FALSE, iopt, ropt, machEnvF);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return(1); }

  /* Call CVBand with  bandwidths ml = mu = MY, */
  flag = CVBand(cvode_mem, MY, MY, Jac, data);
  if (flag != SUCCESS) { printf("CVBand failed.\n"); return(1); }

  /*------------------------*/
  /*---- ADJOINT MEMORY ----*/
  /*------------------------*/

  /* Allocate global memory */
  printf("\nAllocate global memory\n");
  cvadj_mem = CVadjMalloc(cvode_mem, NSTEP);

  /*---------------------*/
  /*---- FORWARD RUN ----*/
  /*---------------------*/

  /* Perform forward run */
  printf("\nForward integration\n");
  flag = CVodeF(cvadj_mem, TOUT, u, &t, NORMAL, &ncheck);

  /* Test check point linked list */
  printf("\nList of Check Points (ncheck = %d)\n", ncheck);
  CVadjCheckPointsList(cvadj_mem);

  /*-------------------------*/
  /*---- ADJOINT PROBLEM ----*/
  /*-------------------------*/

  /* Initialize serial machine environment for backward run */
  machEnvB = M_EnvInit_Serial(NEQ);

  /* Set the tolerances */
  reltolB = RTOLB;
  abstolB = ATOL;

  /* Allocate uB */
  uB = N_VNew(NEQ, machEnvB);
  /* Initialize uB = 0 */
  N_VConst(0.0, uB);

  /* Allocate CVODE memory for backward run */
  printf("\nAllocate CVODE memory for backward run\n");
  flag = CVodeMallocB(cvadj_mem, NEQ, fB, uB, BDF, NEWTON, SS, 
                      &reltolB, &abstolB, data, NULL, 
                      FALSE, NULL, NULL, machEnvB);
  flag = CVBandB(cvadj_mem, MY, MY, JacB, data);
  if (flag != SUCCESS) { printf("CVBandB failed.\n"); return(1); }

  /*----------------------*/
  /*---- BACKWARD RUN ----*/
  /*----------------------*/
  printf("\nBackward integration\n");
  flag = CVodeB(cvadj_mem, uB);
  
  WriteLambda(uB);

  N_VFree(u);                  /* Free the u vector */
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */
  free(data);                  /* Free the user data */
  M_EnvFree_Serial(machEnvF);
  M_EnvFree_Serial(machEnvB);

  return(0);
}


/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, UserData data)
{
  int i, j;
  real x, y, dx, dy;
  real *udata;

  FILE *ff;
  real uij;

  /* Extract needed constants from data */

  dx = data->dx;
  dy = data->dy;

  /* Set pointer to data array in vector u. */

  udata = NV_DATA_S(u);

  /* Load initial profile into u vector */

  for (j=1; j <= MY; j++) {
    y = j*dy;
    for (i=1; i <= MX; i++) {
      x = i*dx;
      IJth(udata,i,j) = x*(XMAX - x)*y*(YMAX - y)*exp(5.0*x*y);
    }
  }  

  ff = fopen("cvabx.u0","w");

  for(i=1; i<=MX+2; i++) fprintf(ff, "0.0 ");
  fprintf(ff,"\n");
  for(j=1; j<= MY; j++) {
    fprintf(ff,"0.0 ");
    for(i=1; i<=MX; i++) {
      uij = IJth(udata, i, j);
      fprintf(ff,"%14.6e ",uij);
    }
    fprintf(ff,"0.0 \n");
  }
  for(i=1; i<=MX+2; i++) fprintf(ff, "0.0 ");
  fprintf(ff,"\n");
  fclose(ff);

  printf("\nInitial condition written to file cvabx.u0\n");

}

/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,u). */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data)
{
  real uij, udn, uup, ult, urt, hordc, horac, verdc, hdiff, hadv, vdiff;
  real *udata, *dudata;
  int i, j;
  UserData data;

  udata = NV_DATA_S(u);
  dudata = NV_DATA_S(udot);

  /* Extract needed constants from data */

  data = (UserData) f_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  /* Loop over all grid points. */

  for (j=1; j <= MY; j++) {

    for (i=1; i <= MX; i++) {

      /* Extract u at x_i, y_j and four neighboring points */

      uij = IJth(udata, i, j);
      udn = (j == 1)  ? 0.0 : IJth(udata, i, j-1);
      uup = (j == MY) ? 0.0 : IJth(udata, i, j+1);
      ult = (i == 1)  ? 0.0 : IJth(udata, i-1, j);
      urt = (i == MX) ? 0.0 : IJth(udata, i+1, j);

      /* Set diffusion and advection terms and load into udot */

      hdiff = hordc*(ult - 2.0*uij + urt);
      hadv = horac*(urt - ult);
      vdiff = verdc*(uup - 2.0*uij + udn);
      IJth(dudata, i, j) = hdiff + hadv + vdiff;
    }
  }
}


/* Jacobian routine. Compute J(t,u). */

static void Jac(integer N, integer mu, integer ml, BandMat J, RhsFn f,
		 void *f_data, real t, N_Vector u, N_Vector fu, N_Vector ewt,
		 real h, real uround, void *jac_data, long int *nfePtr,
		 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  integer i, j, k;
  real *kthCol, hordc, horac, verdc;
  UserData data;

  /*
    The components of f = udot that depend on u(i,j) are
    f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
      df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
      df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
      df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
      df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
      df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
  */

  data = (UserData) jac_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  for (j=1; j <= MY; j++) {
    for (i=1; i <= MX; i++) {
      k = j-1 + (i-1)*MY;
      kthCol = BAND_COL(J,k);

      /* set the kth column of J */

      BAND_COL_ELEM(kthCol,k,k) = -2.0*(verdc+hordc);
      if (i != 1)  BAND_COL_ELEM(kthCol,k-MY,k) = hordc + horac;
      if (i != MX) BAND_COL_ELEM(kthCol,k+MY,k) = hordc - horac;
      if (j != 1)  BAND_COL_ELEM(kthCol,k-1,k)  = verdc;
      if (j != MY) BAND_COL_ELEM(kthCol,k+1,k)  = verdc;
    }
  }
}

/********** Functions Called by the backward CVODE Solver ****************/

static void fB(integer NB, real tB, N_Vector u, 
               N_Vector uB, N_Vector uBdot, void *f_dataB)
{
  UserData data;
  real *uBdata, *duBdata;
  real hordc, horac, verdc;
  real uBij, uBdn, uBup, uBlt, uBrt;
  real hdiffB, hadvB, vdiffB;
  int i, j;

  uBdata = NV_DATA_S(uB);
  duBdata = NV_DATA_S(uBdot);

  /* Extract needed constants from data */

  data = (UserData) f_dataB;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  /* Loop over all grid points. */

  for (j=1; j <= MY; j++) {

    for (i=1; i <= MX; i++) {

      /* Extract u at x_i, y_j and four neighboring points */

      uBij = IJth(uBdata, i, j);
      uBdn = (j == 1)  ? 0.0 : IJth(uBdata, i, j-1);
      uBup = (j == MY) ? 0.0 : IJth(uBdata, i, j+1);
      uBlt = (i == 1)  ? 0.0 : IJth(uBdata, i-1, j);
      uBrt = (i == MX) ? 0.0 : IJth(uBdata, i+1, j);

      /* Set diffusion and advection terms and load into udot */

      hdiffB = hordc*(- uBlt + 2.0*uBij - uBrt);
      hadvB  = horac*(uBrt - uBlt);
      vdiffB = verdc*(- uBup + 2.0*uBij - uBdn);
      IJth(duBdata, i, j) = hdiffB + hadvB + vdiffB - 1.0;
    }
  }
}

static void JacB(integer NB, integer muB, integer mlB, BandMat JB, RhsFnB fB,
                 void *f_dataB, real tB, N_Vector u, 
                 N_Vector uB, N_Vector fuB, N_Vector ewtB,
                 real hB, real uroundB, void *jac_dataB, long int *nfePtrB,
                 N_Vector vtemp1B, N_Vector vtemp2B, N_Vector vtemp3B)
{
  integer i, j, k;
  real *kthCol, hordc, horac, verdc;
  UserData data;

  /* The Jacobian of the adjoint system is: JB = -J^T */

  data = (UserData) jac_dataB;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  for (j=1; j <= MY; j++) {
    for (i=1; i <= MX; i++) {
      k = j-1 + (i-1)*MY;
      kthCol = BAND_COL(JB,k);

      /* set the kth column of J */

      BAND_COL_ELEM(kthCol,k,k) = 2.0*(verdc+hordc);
      if (i != 1)  BAND_COL_ELEM(kthCol,k-MY,k) = - hordc + horac;
      if (i != MX) BAND_COL_ELEM(kthCol,k+MY,k) = - hordc - horac;
      if (j != 1)  BAND_COL_ELEM(kthCol,k-1,k)  = - verdc;
      if (j != MY) BAND_COL_ELEM(kthCol,k+1,k)  = - verdc;
    }
  }
}

static void WriteLambda(N_Vector uB)
{
  real *uBdata, uBij;
  int i, j;
  FILE *ff;

  ff = fopen("cvabx.lambda","w");

  uBdata = NV_DATA_S(uB);
  for(i=1; i<=MX+2; i++) fprintf(ff, "0.0 ");
  fprintf(ff,"\n");
  for(j=1; j<= MY; j++) {
    fprintf(ff,"0.0 ");
    for(i=1; i<=MX; i++) {
      uBij = IJth(uBdata, i, j);
      fprintf(ff,"%14.6e ",uBij);
    }
    fprintf(ff,"0.0 \n");
  }
  for(i=1; i<=MX+2; i++) fprintf(ff, "0.0 ");
  fprintf(ff,"\n");
  fclose(ff);

  printf("\nLambda written to file cvabx.lambda\n");

}
