/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2007-04-23 23:37:26 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Modification of the CVODES example cvsdenx to illustrate treatment
 * of unphysical solution components throuh the RHS function return 
 * flag.
 *
 * Note that, to make possible negative solution components, the
 * absolute tolerances had to be loosened a bit from their values 
 * in cvsdenx.
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 *
 * NOTE: For readibility, no checks are performed on the various 
 *       function return flags.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)

#define NEQ   3
#define Y1    RCONST(1.0)
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)
#define ATOL1 RCONST(7.0e-7)
#define ATOL2 RCONST(1.0e-13)
#define ATOL3 RCONST(1.0e-5)
#define T0    RCONST(0.0)
#define T1    RCONST(0.4)
#define TMULT RCONST(10.0)
#define NOUT  12

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);
static void PrintFinalStats(void *cvode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int iout;
  booleantype check_negative;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  abstol = N_VNew_Serial(NEQ); 

  /* Initialize y */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  /* Set the tolerances */
  reltol = RTOL;
  Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* Initialize and allocate solver memory */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeInit(cvode_mem, f, T0, y);
  CVodeSVtolerances(cvode_mem, reltol, abstol);
  CVodeSetFdata(cvode_mem, &check_negative);
  /* Call CVDense to specify the CVDENSE dense linear solver */
  CVDense(cvode_mem, NEQ);

  /* Case 1: ignore negative solution components */
  printf("Ignore negative solution components\n\n");
  check_negative = FALSE;
  /* In loop, call CVode in CV_NORMAL mode */
  iout = 0;  tout = T1;
  while(1) {
    CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    PrintOutput(t, Ith(y,1), Ith(y,2), Ith(y,3));
    iout++;
    tout *= TMULT;
    if (iout == NOUT) break;
  }
  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Case 2: intercept negative solution components */
  printf("Intercept negative solution components\n\n");
  check_negative = TRUE;
  /* Reinitialize solver */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;
  CVodeReInit(cvode_mem, T0, y);
  /* In loop, call CVode in CV_NORMAL mode */
  iout = 0;  tout = T1;
  while(1) {
    CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    PrintOutput(t, Ith(y,1), Ith(y,2), Ith(y,3));
    iout++;
    tout *= TMULT;
    if (iout == NOUT) break;
  }
  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return(0);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2, y3, yd1, yd3;
  booleantype *check_negative;

  check_negative = (booleantype *)f_data;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  if ( *check_negative && (y1<0 || y2<0 || y3<0) )
    return(1);
 
  Ith(ydot,1) = yd1 = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
  Ith(ydot,3) = yd3 = RCONST(3.0e7)*y2*y2;
  Ith(ydot,2) = -yd1 - yd3;

  return(0);
}

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
  printf("At t = %0.4le   y =%14.6le  %14.6le  %14.6le\n", t, y1, y2, y3);
  return;
}

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;

  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetNumRhsEvals(cvode_mem, &nfe);
  CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  CVodeGetNumErrTestFails(cvode_mem, &netf);
  CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  CVDlsGetNumJacEvals(cvode_mem, &nje);
  CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld\n \n",
	 nni, ncfn, netf);
}

