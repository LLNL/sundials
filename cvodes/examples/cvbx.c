/************************************************************************
 *                                                                      *
 * File        : cvbx.c                                                 *
 * Programmers : Scott D. Cohen and Alan C. Hindmarsh @ LLNL            *
 * Version of  : 5 March 2002                                           *
 *----------------------------------------------------------------------*
 * Modified by R. Serban to work with new serial nvector (5/3/2002)     *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem with a banded Jacobian,    *
 * with the program for its solution by CVODE.                          *
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
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* CVODE header files with a description of contents used in cvbx.c */

#include "llnltyps.h"       /* definitions of real, integer, FALSE               */
#include "cvodes.h"         /* prototypes for CVodeMalloc, CVode, and CVodeFree, */
                            /* constants OPT_SIZE, BDF, NEWTON, SS, SUCCESS,     */
                            /* NST, NFE, NSETUPS, NNI, NCFN, NETF                */
#include "cvsband.h"        /* prototype for CVBand, constant BAND_NJE           */
#include "nvector_serial.h" /* definitions of type N_Vector and macro NV_DATA_S, */
                            /* prototypes for N_VNew, N_VFree, N_VMaxNorm        */
#include "band.h"           /* definitions of type BandMat, macros               */


/* Problem Constants */

#define XMAX  2.0          /* domain boundaries */
#define YMAX  1.0
#define MX    10           /* mesh dimensions   */
#define MY     5
#define NEQ   MX*MY        /* number of equations       */
#define ATOL  1.e-5        /* scalar absolute tolerance */
#define T0    0.0          /* initial time              */
#define T1    0.1          /* first output time         */
#define DTOUT 0.1          /* output time increment     */
#define NOUT  10           /* number of output times    */

/* User-defined vector access macro IJth */

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. 
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the macro call vdata = NV_DATA_S(v),
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

static void PrintFinalStats(long int iopt[]);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data);

static void Jac(integer N, integer mu, integer ml, BandMat J, RhsFn f,
		 void *f_data, real t, N_Vector u, N_Vector fu, N_Vector ewt,
		 real h, real uround, void *jac_data, long int *nfePtr,
		 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 


/***************************** Main Program ******************************/

int main()
{
  M_Env machEnv;
  real ropt[OPT_SIZE], dx, dy, reltol, abstol, t, tout, umax;
  long int iopt[OPT_SIZE];
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, flag;

  /* Initialize serial machine environment */
  machEnv = M_EnvInit_Serial(NEQ);

  u = N_VNew(NEQ, machEnv);    /* Allocate u vector */

  reltol = 0.0;                /* Set the tolerances */
  abstol = ATOL;

  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
  dx = data->dx = XMAX/(MX+1);     /* Set grid coefficients in data */
  dy = data->dy = YMAX/(MY+1);
  data->hdcoef = 1.0/(dx*dx);
  data->hacoef = 0.5/(2.0*dx);
  data->vdcoef = 1.0/(dy*dy);

  SetIC(u, data);              /* Initialize u vector */

  /* Call CVodeMalloc to initialize CVODE: 

     NEQ     is the problem size = number of equations
     f       is the user's right hand side function in u'=f(t,u)
     T0      is the initial time
     u       is the initial dependent variable vector
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration
     SS      specifies scalar relative and absolute tolerances
     &reltol and &abstol are pointers to the scalar tolerances
     data    is the pointer to the user-defined block of coefficients
     FALSE   indicates there are no optional inputs in iopt and ropt
     iopt    and ropt arrays communicate optional integer and real input/output

     A pointer to CVODE problem memory is returned and stored in cvode_mem.  */

  cvode_mem = CVodeMalloc(NEQ, f, T0, u, BDF, NEWTON, SS, &reltol, &abstol,
                          data, NULL, FALSE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return(1); }

  /* Call CVBand to specify the CVODE band linear solver with the
     user-supplied Jacobian routine Jac, bandwidths ml = mu = MY,
     and the pointer to the user-defined block data. */

  flag = CVBand(cvode_mem, MY, MY, Jac, data);
  if (flag != SUCCESS) { printf("CVBand failed.\n"); return(1); }

  /* In loop over output points, call CVode, print results, test for error */

  printf(" \n2-D advection-diffusion equation, mesh dimensions =%3d %3d\n\n",
         MX,MY);
  umax = N_VMaxNorm(u);
  printf("At t = %4.2f    max.norm(u) =%14.6e \n", T0,umax);
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    umax = N_VMaxNorm(u);
    printf("At t = %4.2f    max.norm(u) =%14.6e   nst =%4ld \n", 
           t,umax,iopt[NST]);
    if (flag != SUCCESS) { printf("CVode failed, flag=%d.\n", flag); break; }
  }

  N_VFree(u);                  /* Free the u vector */
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */
  free(data);                  /* Free the user data */
  M_EnvFree_Serial(machEnv);   /* Free the machine environment memory */

  PrintFinalStats(iopt);       /* Print some final statistics   */
  return(0);
}


/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, UserData data)
{
  int i, j;
  real x, y, dx, dy;
  real *udata;

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
}


/* Print some final statistics located in the iopt array */

static void PrintFinalStats(long int iopt[])
{
  printf("\nFinal Statistics.. \n\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
         iopt[NST], iopt[NFE], iopt[NSETUPS], iopt[BAND_NJE]);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
         iopt[NNI], iopt[NCFN], iopt[NETF]);
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
 
