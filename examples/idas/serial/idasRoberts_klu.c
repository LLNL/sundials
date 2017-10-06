/*
 * -----------------------------------------------------------------
 * $Revision: 4834 $
 * $Date: 2016-08-01 16:59:05 -0700 (Mon, 01 Aug 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL  (based on idasRoberts_dns.c)
 * -----------------------------------------------------------------
 * This simple example problem for IDA, due to Robertson, 
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
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01.
 *
 * The problem is solved with IDA using IDAKLU for the linear
 * solver, with a user-supplied Jacobian. Output is printed at
 * t = .4, 4, 40, ..., 4e10.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>

#include <idas/idas.h>
#include <idas/idas_klu.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Problem Constants */

#define NEQ   3
#define NOUT  12

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Prototypes of functions called by IDA */

int resrob(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data);

static int grob(realtype t, N_Vector yy, N_Vector yp,
                realtype *gout, void *user_data);

int jacrob(realtype tt,  realtype cj, 
           N_Vector yy, N_Vector yp, N_Vector resvec, 
	   SlsMat JacMat, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);
static void PrintOutput(void *mem, realtype t, N_Vector y);
static void PrintRootInfo(int root_f1, int root_f2);
static void PrintFinalStats(void *mem);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, tout1, tout, tret;
  int iout, retval, retvalr;
  int rootsfound[2];

  int nnz;

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
  yval  = N_VGetArrayPointer_Serial(yy);
  yval[0] = ONE;
  yval[1] = ZERO;
  yval[2] = ZERO;

  ypval = N_VGetArrayPointer_Serial(yp);
  ypval[0]  = RCONST(-0.04);
  ypval[1]  = RCONST(0.04);
  ypval[2]  = ZERO;  

  rtol = RCONST(1.0e-4);

  atval = N_VGetArrayPointer_Serial(avtol);
  atval[0] = RCONST(1.0e-8);
  atval[1] = RCONST(1.0e-6);
  atval[2] = RCONST(1.0e-6);

  /* Integration limits */
  t0 = ZERO;
  tout1 = RCONST(0.4);

  PrintHeader(rtol, avtol, yy);

  /* Call IDACreate and IDAMalloc to initialize IDA memory */
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);
  retval = IDAInit(mem, resrob, t0, yy, yp);
  if(check_flag(&retval, "IDAInit", 1)) return(1);
  retval = IDASVtolerances(mem, rtol, avtol);
  if(check_flag(&retval, "IDASVtolerances", 1)) return(1);

  /* Free avtol */
  N_VDestroy_Serial(avtol);

  /* Call IDARootInit to specify the root function grob with 2 components */
  retval = IDARootInit(mem, 2, grob);
  if (check_flag(&retval, "IDARootInit", 1)) return(1);

  /* Call IDAKLU and set up the linear solver. */
  nnz = NEQ * NEQ;
  retval = IDAKLU(mem, NEQ, nnz, CSC_MAT);
  if(check_flag(&retval, "IDAKLU", 1)) return(1);
  retval = IDASlsSetSparseJacFn(mem, jacrob);
  if(check_flag(&retval, "IDASlsSetSparseJacFn", 1)) return(1);

  /* In loop, call IDASolve, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached. */

  iout = 0; tout = tout1;
  while(1) {

    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

    PrintOutput(mem,tret,yy);

    if(check_flag(&retval, "IDASolve", 1)) return(1);

    if (retval == IDA_ROOT_RETURN) {
      retvalr = IDAGetRootInfo(mem, rootsfound);
      check_flag(&retvalr, "IDAGetRootInfo", 1);
      PrintRootInfo(rootsfound[0],rootsfound[1]);
    }

    if (retval == IDA_SUCCESS) {
      iout++;
      tout *= RCONST(10.0);
    }

    if (iout == NOUT) break;
  }

  PrintFinalStats(mem);

  /* Free memory */

  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);

  return(0);
  
}

/*
 *--------------------------------------------------------------------
 * Functions called by IDA
 *--------------------------------------------------------------------
 */

/*
 * Define the system residual function. 
 */

int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  realtype *yval, *ypval, *rval;

  yval = N_VGetArrayPointer_Serial(yy); 
  ypval = N_VGetArrayPointer_Serial(yp); 
  rval = N_VGetArrayPointer_Serial(rr);

  rval[0]  = RCONST(-0.04)*yval[0] + RCONST(1.0e4)*yval[1]*yval[2];
  rval[1]  = -rval[0] - RCONST(3.0e7)*yval[1]*yval[1] - ypval[1];
  rval[0] -=  ypval[0];
  rval[2]  =  yval[0] + yval[1] + yval[2] - ONE;

  return(0);
}

/*
 * Root function routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int grob(realtype t, N_Vector yy, N_Vector yp, realtype *gout,
                void *user_data)
{
  realtype *yval, y1, y3;

  yval = N_VGetArrayPointer_Serial(yy); 
  y1 = yval[0]; y3 = yval[2];
  gout[0] = y1 - RCONST(0.0001);
  gout[1] = y3 - RCONST(0.01);

  return(0);
}

/*
 * Define the Jacobian function. 
 */

int jacrob(realtype tt,  realtype cj, 
           N_Vector yy, N_Vector yp, N_Vector resvec,
	   SlsMat JacMat, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  int* colptrs;
  int* rowvals;
  realtype* data;
  
  yval = N_VGetArrayPointer_Serial(yy);
  colptrs = (*JacMat->colptrs);
  rowvals = (*JacMat->rowvals);
  data    = JacMat->data;

  SparseSetMatToZero(JacMat);

  colptrs[0] = 0;
  colptrs[1] = 3;
  colptrs[2] = 6;
  colptrs[3] = 9;

  data[0] = RCONST(-0.04) - cj;
  rowvals[0] = 0;
  data[1] = RCONST(0.04);
  rowvals[1] = 1;
  data[2] = ONE;
  rowvals[2] = 2;

  data[3] = RCONST(1.0e4)*yval[2];
  rowvals[3] = 0;
  data[4] = (RCONST(-1.0e4)*yval[2]) - (RCONST(6.0e7)*yval[1]) - cj;
  rowvals[4] = 1;
  data[5] = ONE;
  rowvals[5] = 2;

  data[6] = RCONST(1.0e4)*yval[1];
  rowvals[6] = 0;
  data[7] = RCONST(-1.0e4)*yval[1];
  rowvals[7] = 1;
  data[8] = ONE;
  rowvals[8] = 2;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */

/* 
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y)
{
  realtype *atval, *yval;

  atval  = N_VGetArrayPointer_Serial(avtol);
  yval  = N_VGetArrayPointer_Serial(y);

  printf("\nidasRoberts_klu: Robertson kinetics DAE serial example problem for IDA.\n");
  printf("               Three equation chemical kinetics problem.\n\n");
  printf("Linear solver: IDAKLU, with user-supplied Jacobian.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg %Lg %Lg \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%Lg %Lg %Lg)\n",
         yval[0], yval[1], yval[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",
         yval[0], yval[1], yval[2]);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",
         yval[0], yval[1], yval[2]);
#endif
  printf("Constraints and id not used.\n\n");
  printf("-----------------------------------------------------------------------\n");
  printf("  t             y1           y2           y3");
  printf("      | nst  k      h\n");
  printf("-----------------------------------------------------------------------\n");
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

  yval  = N_VGetArrayPointer_Serial(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.4Le %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}

static void PrintRootInfo(int root_f1, int root_f2)
{
  printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);
  return;
}

/*
 * Print final integrator statistics
 */

static void PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDASlsGetNumJacEvals(mem, &nje);
  check_flag(&retval, "IDASlsGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  /*  
      retval = IDASlsGetNumResEvals(mem, &nreLS);
      check_flag(&retval, "IDASlsGetNumResEvals", 1);
  */
  nreLS = 0;
  retval = IDAGetNumGEvals(mem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  printf("Number of root fn. evaluations     = %ld\n", nge);
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

static int check_flag(void *flagvalue, const char *funcname, int opt)
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
    errflag = (int *) flagvalue;
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
