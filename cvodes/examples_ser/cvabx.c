/************************************************************************
 * File       : cvabx.c                                                 *
 * Programmers: Radu Serban @ LLNL                                      *
 * Version of : 15 July 2003                                            *
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
#include "sundialstypes.h"
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
  realtype dx, dy, hdcoef, hacoef, vdcoef;
} *UserData;

/* Private Helper Functions */

static void SetIC(N_Vector u, UserData data);
static void WriteLambda(N_Vector uB);

/* Functions Called by the CVODE Solver */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data);

static void Jac(long int N, long int mu, long int ml, BandMat J,
                realtype t, N_Vector u, N_Vector fu, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 

static void fB(realtype tB, N_Vector u, N_Vector uB, N_Vector uBdot, void *f_dataB);

static void JacB(long int NB, long int muB, long int mlB, BandMat JB,
                 realtype tB, N_Vector u, 
                 N_Vector uB, N_Vector fuB, void *jac_dataB,
                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B); 

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  NV_Spec nvSpecF, nvSpecB;

  UserData data;

  void *cvadj_mem;
  void *cvode_mem;

  realtype dx, dy, reltol, abstol, t;
  N_Vector u;

  realtype reltolB, abstolB;
  N_Vector uB;
  
  int flag, ncheck;

  nvSpecF = nvSpecB = NULL;
  data = NULL;
  cvadj_mem = cvode_mem = NULL;
  u = uB = NULL;

  /*--------------------*/
  /*----- USER DATA ----*/
  /*--------------------*/

  /* Allocate data memory */
  data = (UserData) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2)) return(1);
  dx = data->dx = XMAX/(MX+1);     /* Set grid coefficients in data */
  dy = data->dy = YMAX/(MY+1);
  data->hdcoef = 1.0/(dx*dx);
  data->hacoef = 1.5/(2.0*dx);
  data->vdcoef = 1.0/(dy*dy);

  /*--------------------------*/
  /*---- ORIGINAL PROBLEM ----*/
  /*--------------------------*/

  /* Initialize serial vector specification for forward run */
  nvSpecF = NV_SpecInit_Serial(NEQ);
  if(check_flag((void *)nvSpecF, "NV_SpecInit", 0)) return(1);

  /* Set the tolerances */
  reltol = 0.0;
  abstol = ATOL;

  /* Allocate u vector */
  u = N_VNew(nvSpecF);
  if(check_flag((void *)u, "N_VNew", 0)) return(1);

  /* Initialize u vector */
  SetIC(u, data);

  /* Create and allocate CVODES memory for forward run */

  printf("\nCreate and allocate CVODES memory for forward runs\n");

  cvode_mem = CVodeCreate(BDF, NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeSetFdata(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetFdata", 1)) return(1);

  flag = CVodeMalloc(cvode_mem, f, T0, u, SS, &reltol, &abstol, nvSpecF);
  if(check_flag(&flag, "CVodeMalloc", 1)) return(1);

  /* Call CVBand with  bandwidths ml = mu = MY, */

  flag = CVBand(cvode_mem, NEQ, MY, MY);
  if(check_flag(&flag, "CVBand", 1)) return(1);

  flag = CVBandSetJacFn(cvode_mem, Jac);
  if(check_flag(&flag, "CVBandSetJacFn", 1)) return(1);

  flag = CVBandSetJacData(cvode_mem, data);
  if(check_flag(&flag, "CVBandSetJacData", 1)) return(1);

  /*------------------------*/
  /*---- ADJOINT MEMORY ----*/
  /*------------------------*/

  /* Allocate global memory */

  printf("\nAllocate global memory\n");

  cvadj_mem = CVadjMalloc(cvode_mem, NSTEP);
  if(check_flag((void *)cvadj_mem, "CVadjMalloc", 0)) return(1);

  /*---------------------*/
  /*---- FORWARD RUN ----*/
  /*---------------------*/

  /* Perform forward run */
  printf("\nForward integration\n");
  flag = CVodeF(cvadj_mem, TOUT, u, &t, NORMAL, &ncheck);
  if(check_flag(&flag, "CVodeF", 1)) return(1);

  /* Test check point linked list */
  printf("\nList of Check Points (ncheck = %d)\n", ncheck);
  CVadjGetCheckPointsList(cvadj_mem);

  /*-------------------------*/
  /*---- ADJOINT PROBLEM ----*/
  /*-------------------------*/

  /* Initialize serial vector specification for backward run */
  nvSpecB = NV_SpecInit_Serial(NEQ);
  if(check_flag((void *)nvSpecB, "NV_SpecInit", 0)) return(1);

  /* Set the tolerances */
  reltolB = RTOLB;
  abstolB = ATOL;

  /* Allocate uB */
  uB = N_VNew(nvSpecB);
  if(check_flag((void *)uB, "N_VNew", 0)) return(1);
  /* Initialize uB = 0 */
  N_VConst(0.0, uB);

  /* Create and allocate CVODES memory for backward run */

  printf("\nCreate and allocate CVODES memory for backward run\n");

  flag = CVodeCreateB(cvadj_mem, BDF, NEWTON);
  if(check_flag(&flag, "CVodeCreateB", 1)) return(1);

  flag = CVodeSetFdataB(cvadj_mem, data);
  if(check_flag(&flag, "CVodeSetFdataB", 1)) return(1);

  flag = CVodeMallocB(cvadj_mem, fB, TOUT, uB, SS, &reltolB, &abstolB, nvSpecB);
  if(check_flag(&flag, "CVodeMallocB", 1)) return(1);

  flag = CVBandB(cvadj_mem, NEQ, MY, MY);
  if(check_flag(&flag, "CVBandB", 1)) return(1);
  
  flag = CVBandSetJacFnB(cvadj_mem, JacB);
  if(check_flag(&flag, "CVBandSetJacFnB", 1)) return(1);

  flag = CVBandSetJacDataB(cvadj_mem, data);
  if(check_flag(&flag, "CVBandSetJacDataB", 1)) return(1);

  /*----------------------*/
  /*---- BACKWARD RUN ----*/
  /*----------------------*/

  printf("\nBackward integration\n");
  flag = CVodeB(cvadj_mem, uB);
  flag -= 1;
  if(check_flag(&flag, "CVodeB", 1)) return(1);

  WriteLambda(uB);

  N_VFree(u);                        /* Free the u vector */
  N_VFree(uB);                      /* Free the uB vector */
  CVodeFree(cvode_mem);      /* Free the CVODE problem memory */
  free(data);                     /* Free the user data */
  NV_SpecFree_Serial(nvSpecF);
  NV_SpecFree_Serial(nvSpecB);

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, UserData data)
{
  int i, j;
  realtype x, y, dx, dy;
  realtype *udata;

  FILE *ff;
  realtype uij;

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

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  realtype uij, udn, uup, ult, urt, hordc, horac, verdc, hdiff, hadv, vdiff;
  realtype *udata, *dudata;
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

static void Jac(long int N, long int mu, long int ml, BandMat J,
                realtype t, N_Vector u, N_Vector fu, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int i, j, k;
  realtype *kthCol, hordc, horac, verdc;
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

static void fB(realtype tB, N_Vector u, N_Vector uB, N_Vector uBdot, 
               void *f_dataB)
{
  UserData data;
  realtype *uBdata, *duBdata;
  realtype hordc, horac, verdc;
  realtype uBij, uBdn, uBup, uBlt, uBrt;
  realtype hdiffB, hadvB, vdiffB;
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

static void JacB(long int NB, long int muB, long int mlB, BandMat JB,
                 realtype tB, N_Vector u, 
                 N_Vector uB, N_Vector fuB, void *jac_dataB,
                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  long int i, j, k;
  realtype *kthCol, hordc, horac, verdc;
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
  realtype *uBdata, uBij;
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
