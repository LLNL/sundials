/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * simple pendulum example
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Precision specific formatting macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* Precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SIN(x)   (sin((x)))
#define COS(x)   (cos((x)))
#define SQRT(x)  (sqrt((x)))
#define ABS(x)   (fabs((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x)   (sinf((x)))
#define COS(x)   (cosf((x)))
#define SQRT(x)  (sqrtf((x)))
#define ABS(x)   (fabsf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x)   (sinl((x)))
#define COS(x)   (cosl((x)))
#define SQRT(x)  (sqrtl((x)))
#define ABS(x)   (fabsl((x)))
#endif

#define Ith(v,i)    NV_Ith_S(v,i-1)

/* Problem Constants */
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define GRAV  RCONST(13.750371636040745654980191559621114395801712)

#define TOL     RCONST(1.0e-5)
#define TOL_REF RCONST(1.0e-14)

/* Functions provided to CVODE */
static int fref(realtype t, N_Vector yy, N_Vector fy, void *f_data);

static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data);
static int proj(realtype t, N_Vector yy, N_Vector corr,
                realtype epsProj, N_Vector err, void *pdata);

/* Functions to integrate the Cartesian and reference solutions */
void GetSol(void *cpode_mem, N_Vector yy0, realtype tol, realtype tf,
            booleantype proj, N_Vector yref);

void RefSol(realtype tout, N_Vector yref);

/* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
  int         i;
  int         flag;                     /* reusable return flag    */
  realtype    tol     = RCONST(1.0e-5); /* integration tolerance   */
  realtype    tf      = RCONST(30.0);   /* final integration time  */

  void            *cpode_mem = NULL; /* CVODE memory              */
  N_Vector         yy0       = NULL; /* initial condition vector  */
  realtype        *yy0data   = NULL; /* vector data               */
  N_Vector         yref      = NULL; /* reference solution vector */


  /* Compute reference solution */
  yref = N_VNew_Serial(4);

  RefSol(tf, yref);

  /* Create serial vector to store the initial condition */
  yy0 = N_VNew_Serial(4);

  /* Set the initial condition values */
  yy0data = N_VGetArrayPointer(yy0);

  yy0data[0] = ONE;  /* x  */
  yy0data[1] = ZERO; /* y  */
  yy0data[2] = ZERO; /* xd */
  yy0data[3] = ZERO; /* yd */

  /* Create CVODE memory */
  cpode_mem = CPodeCreate(CP_BDF, CP_NEWTON);

  /* Initialize CVODE */
  flag = CPodeInitExpl(cpode_mem, f, ZERO, yy0);

  /* Set integration tolerances */
  flag = CPodeSStolerances(cpode_mem, tol, tol);

  /* Create dense SUNLinearSolver object */
  flag = CPDense(cpode_mem, 4);

  /* Set a user-supplied projection function */
  flag = CPodeProjDefine(cpode_mem, proj);

  /* Set maximum number of steps between outputs */
  flag = CPodeSetMaxNumSteps(cpode_mem, 50000);

  flag = CPodeSetStopTime(cpode_mem, tf);

  /* Compute the solution with various tolerances */
  for (i = 0; i < 5; i++) {

    /* Output tolerance and output header for this run */
    printf("\n\nTol = %8.2" ESYM"\n", tol);
    printf("Project    x         y");
    printf("         x'        y'     |     g      |    ");
    printf("nst     rhs eval    setups (J eval)  |   cf   ef\n");

    /* Compute solution with projection */
    GetSol(cpode_mem, yy0, tol, tf, TRUE, yref);

    /* Compute solution without projection */
    GetSol(cpode_mem, yy0, tol, tf, FALSE, yref);

    /* Reduce tolerance for next run */
    tol /= 10.0;
  }

  /* Free memory */
  N_VDestroy_Serial(yref);
  N_VDestroy_Serial(yy0);
  CPodeFree(&cpode_mem);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Functions to integrate the Cartesian and reference systems
 * ---------------------------------------------------------------------------*/


void GetSol(void *cpode_mem, N_Vector yy0, realtype tol, realtype tf,
            booleantype proj, N_Vector yref)
{
  N_Vector yy, yp;
  realtype t, x, y, xd, yd, g;
  int flag;

  /* Integrator stats */
  long int nst, nfe, nsetups, nje, nfeLS, ncfn, netf;

  /* Enable or disable projection */
  if (proj)
  {
    printf("  YES   ");
    CPodeSetProjFrequency(cpode_mem, 1);
  }
  else
  {
    CPodeSetProjFrequency(cpode_mem, 0);
    printf("  NO    ");
  }

  /* Create vector to store the solution */
  yy = N_VNew_Serial(4);
  yp = N_VNew_Serial(4);

  flag = CPodeReInitExpl(cpode_mem, ZERO, yy0);

  /* Set integration tolerances */
  flag = CPodeSStolerances(cpode_mem, tol, tol);

  flag = CPode(cpode_mem, tf, &t, yy, yp, CP_NORMAL);

  /* Compute the constraint violation */
  x  = Ith(yy,1);
  y  = Ith(yy,2);
  g = SUNRabs(x*x + y*y - 1.0);

  /* Compute the absolute error compared to the reference solution */
  N_VLinearSum(ONE, yy, -ONE, yref, yy);
  N_VAbs(yy, yy);

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);

  /* Output errors */
  printf("%8.2" ESYM"  %8.2" ESYM"  %8.2" ESYM"  %8.2" ESYM"  |  %8.2" ESYM"  |",
         Ith(yy,1),Ith(yy,2),Ith(yy,3),Ith(yy,4),g);

  /* Get integrator stats */
  CPodeGetNumSteps(cpode_mem, &nst);

  CPodeGetNumFctEvals(cpode_mem, &nfe);

  CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);

  CPodeGetNumErrTestFails(cpode_mem, &netf);

  CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  CPDlsGetNumJacEvals(cpode_mem, &nje);

  CPDlsGetNumFctEvals(cpode_mem, &nfeLS);

  /* Output stats */
  printf(" %6ld   %6ld+%-4ld     %4ld (%3ld)     |  %3ld  %3ld\n",
         nst, nfe, nfeLS, nsetups, nje, ncfn, netf);

  /* Free solution vector */
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);

  return;
}


/* Compute the reference system solution */
void RefSol(realtype tf, N_Vector yref)
{
  void             *cpode_mem;
  N_Vector         yy, yp;
  realtype         tol, t, th, thd;

  int       flag;

  /* Create the solution vector */
  yy = N_VNew_Serial(2);
  yp = N_VNew_Serial(2);
  Ith(yy,1) = 0.0;  /* theta */
  Ith(yy,2) = 0.0;  /* thetad */
  tol = TOL_REF;

  cpode_mem = CPodeCreate(CP_BDF, CP_NEWTON);
  flag = CPodeSetMaxNumSteps(cpode_mem, 100000);
  flag = CPodeInitExpl(cpode_mem, fref, 0.0, yy);
  flag = CPodeSStolerances(cpode_mem, tol ,tol);
  flag = CPDense(cpode_mem, 2);

  flag = CPodeSetStopTime(cpode_mem, tf);
  flag = CPode(cpode_mem, tf, &t, yy, yp, CP_NORMAL);
  th  = Ith(yy,1);
  thd = Ith(yy,2);
  Ith(yref,1) = cos(th);
  Ith(yref,2) = sin(th);
  Ith(yref,3) = -thd*sin(th);
  Ith(yref,4) =  thd*cos(th);

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  CPodeFree(&cpode_mem);

  return;
}


/* -----------------------------------------------------------------------------
 * Functions provided to CVODE
 * ---------------------------------------------------------------------------*/


/* ODE RHS function for the reference system */
static int fref(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype th, thd, g;

  g = 13.7503716373294544;

  th  = Ith(yy,1);
  thd  = Ith(yy,2);

  Ith(fy,1) = thd;
  Ith(fy,2) = -g*cos(th);

  return(0);
}


/* ODE RHS function for the Cartesian system */
static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype x, y, xd, yd, g, tmp;

  g = 13.7503716373294544;

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);

  tmp = xd*xd + yd*yd - g*y;

  Ith(fy,1) = xd;
  Ith(fy,2) = yd;
  Ith(fy,3) = -x*tmp;
  Ith(fy,4) = -y*tmp - g;

  return(0);
}


/* Projection function */
static int proj(realtype t, N_Vector yy, N_Vector corr,
                realtype epsProj, N_Vector err, void *pdata)
{
  realtype x, y, xd, yd;
  realtype x_new, y_new, xd_new, yd_new;
  realtype e1, e2, e3, e4;
  realtype e1_new, e2_new, e3_new, e4_new;
  realtype R;

  /* Extract current solution */

  x  = Ith(yy,1);
  y  = Ith(yy,2);
  xd = Ith(yy,3);
  yd = Ith(yy,4);

  /* Project onto manifold */

  R = sqrt(x*x+y*y);

  x_new = x/R;
  y_new = y/R;

  xd_new =   xd*y_new*y_new - yd*x_new*y_new;
  yd_new = - xd*x_new*y_new + yd*x_new*x_new;

  /* Return corrections */

  Ith(corr,1) = x_new  - x;
  Ith(corr,2) = y_new  - y;
  Ith(corr,3) = xd_new - xd;
  Ith(corr,4) = yd_new - yd;

  /*      +-            -+
   *      |  y*y    -x*y |
   *  P = |              |
   *      | -x*y     x*x |
   *      +-            -+
   */

  /* Return err <-  P * err */

  e1 = Ith(err,1);
  e2 = Ith(err,2);
  e3 = Ith(err,3);
  e4 = Ith(err,4);

  e1_new =  y_new*y_new * e1 - x_new*y_new * e2;
  e2_new = -x_new*y_new * e1 + x_new*x_new * e2;

  e3_new =  y_new*y_new * e3 - x_new*y_new * e4;
  e4_new = -x_new*y_new * e3 + x_new*x_new * e4;

  Ith(err,1) = e1_new;
  Ith(err,2) = e2_new;
  Ith(err,3) = e3_new;
  Ith(err,4) = e4_new;

  return(0);
}
