/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2004-11-10 01:03:06 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * This simple example problem for IDA/IDAS, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.
 *
 * The problem is solved with IDA/IDAS using IDADENSE for the linear
 * solver, with a user-supplied Jacobian. Output is printed at
 * t = .4, 4, 40, ..., 4e10.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include "sundialstypes.h"
#include "sundialsmath.h"
#include "nvector_serial.h"
#include "idas.h"
#include "idadense.h"

/* Problem Constants */

#define NEQ   3
#define NOUT  12

#define ZERO RCONST(0.0);
#define ONE  RCONST(1.0);

/* Macro to define dense matrix elements, indexed from 1. */

#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Prototypes of functions called by IDA */

int resrob(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *rdata);

int jacrob(long int Neq, realtype tt, N_Vector yy, N_Vector yp,
           N_Vector resvec, realtype cj, void *jdata, DenseMat JJ,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);
static void PrintOutput(void *mem, realtype t, N_Vector y);
static void PrintFinalStats(void *mem);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, t1, tout, tret;
  int iout, retval;

  mem = NULL;
  yy = yp = avtol = NULL;
  yval = ypval = atval = NULL;

  /* Allocate N-vectors. */

  yy = N_VNew_Serial(NEQ);
  if(check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);
  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);
  avtol = N_VNew_Serial(NEQ);
  if(check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);

  /* Create and initialize  y, y', and absolute tolerance vectors. */

  yval  = NV_DATA_S(yy);
  yval[0] = ONE;
  yval[1] = ZERO;
  yval[2] = ZERO;

  ypval = NV_DATA_S(yp);
  ypval[0]  = RCONST(-0.04);
  ypval[1]  = RCONST(0.04);
  ypval[2]  = ZERO;  

  rtol = RCONST(1.0e-4);

  atval = NV_DATA_S(avtol);
  atval[0] = RCONST(1.0e-6);
  atval[1] = RCONST(1.0e-10);
  atval[2] = RCONST(1.0e-6);

  /* Integration limits */

  t0 = ZERO;
  t1 = RCONST(0.4);

  /* Call IDACreate and IDAMalloc to initialize IDA memory */

  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);
  retval = IDAMalloc(mem, resrob, t0, yy, yp, IDA_SV, &rtol, avtol);
  if(check_flag(&retval, "IDAMalloc", 1)) return(1);

  /* Call IDADense and set up the linear solver. */

  retval = IDADense(mem, NEQ);
  if(check_flag(&retval, "IDADense", 1)) return(1);
  retval = IDADenseSetJacFn(mem, jacrob);
  if(check_flag(&retval, "IDADenseSetJacFn", 1)) return(1);

  /* Loop over tout values and call IDASolve. */

  PrintHeader(rtol, avtol, yy);

  for (tout = t1, iout = 1; iout <= NOUT ; iout++, tout *= RCONST(10.0)) {
    retval=IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    if(check_flag(&retval, "IDASolve", 1)) return(1);
    PrintOutput(mem,tret,yy);
  }

  PrintFinalStats(mem);

  /* Free memory */

  IDAFree(mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(avtol);

  return(0);
  
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * Define the system residual function. 
 */

int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *rdata)
{
  realtype *yval, *ypval, *rval;

  yval = NV_DATA_S(yy); 
  ypval = NV_DATA_S(yp); 
  rval = NV_DATA_S(rr);

  rval[0]  = RCONST(-0.04)*yval[0] + RCONST(1.0e4)*yval[1]*yval[2];
  rval[1]  = -rval[0] - RCONST(3.0e7)*yval[1]*yval[1] - ypval[1];
  rval[0] -=  ypval[0];
  rval[2]  =  yval[0] + yval[1] + yval[2] - ONE;

  return(0);

}

/*
 * Define the Jacobian function. 
 */

int jacrob(long int Neq, realtype tt, N_Vector yy, N_Vector yp,
           N_Vector resvec, realtype cj, void *jdata, DenseMat JJ,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{

  realtype *yval;
  
  yval = NV_DATA_S(yy);

  IJth(JJ,1,1) = RCONST(-0.04) - cj;
  IJth(JJ,2,1) = RCONST(0.04);
  IJth(JJ,3,1) = ONE;
  IJth(JJ,1,2) = RCONST(1.0e4)*yval[2];
  IJth(JJ,2,2) = RCONST(-1.0e4)*yval[2] - RCONST(6.0e7)*yval[1] - cj;
  IJth(JJ,3,2) = ONE;
  IJth(JJ,1,3) = RCONST(1.0e4)*yval[1];
  IJth(JJ,2,3) = RCONST(-1.0e4)*yval[1];
  IJth(JJ,3,3) = ONE;

  return(0);

}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* 
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y)
{
  realtype *atval, *yval;

  atval  = NV_DATA_S(avtol);
  yval  = NV_DATA_S(y);

  printf("\nirobx: Robertson kinetics DAE serial example problem for IDA \n");
  printf("       Three equation chemical kinetics problem. \n\n");
  printf("Linear solver: IDADENSE, with user-supplied Jacobian.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg %Lg %Lg \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%Lg %Lg %Lg)\n",
         yval[0], yval[1], yval[2]);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",
         yval[0], yval[1], yval[2]);
#endif
  printf("Constraints and id not used.\n\n");
  printf("---------------------------------------------------------------------\n");
  printf("  t           y1           y2           y3");
  printf("      | nst  k      h\n");
  printf("---------------------------------------------------------------------\n");

}

/*
 * Print Output
 */

static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval  = NV_DATA_S(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.2Le %12.4Le %12.4Le %12.4Le | %3d  %1d %12.4Le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%8.2e %12.4e %12.4e %12.4e | %3d  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}

/*
 * Print final integrator statistics
 */

static void PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreD, netf, ncfn;

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDADenseGetNumJacEvals(mem, &nje);
  check_flag(&retval, "IDADenseGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDADenseGetNumResEvals(mem, &nreD);
  check_flag(&retval, "IDADenseGetNumResEvals", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreD);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);

}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
