/************************************************************************
 * File       : cvdx.c                                                  *
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh and Radu Serban @LLNL *
 * Version of : 10 July 2003                                            *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem, with the coding           *
 * needed for its solution by CVODE.  The problem is from chemical      *
 * kinetics, and consists of the following three rate equations..       *
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3                                     *
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2                        *
 *    dy3/dt = 3.e7*(y2)^2                                              *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions   *
 * y1 = 1.0, y2 = y3 = 0.  The problem is stiff.                        *
 * This program solves the problem with the BDF method, Newton          *
 * iteration with the CVODE dense linear solver, and a user-supplied    *
 * Jacobian routine.                                                    * 
 * It uses a scalar relative tolerance and a vector absolute tolerance. *
 * Output is printed in decades from t = .4 to t = 4.e10.               *
 * Run statistics (optional outputs) are printed at the end.            *
 ************************************************************************/

#include <stdio.h>


/* CVODE header files with a description of contents used in cvdx.c */

#include "sundialstypes.h"   /* definition of types realtype                  */
                             /* and the constant FALSE                        */
#include "cvodes.h"          /* prototypes for CVodeMalloc, CVode, and        */
                             /* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
                             /* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
#include "cvsdense.h"        /* prototype for CVDense, constant DENSE_NJE     */
#include "nvector_serial.h"  /* definitions of type N_Vector and macro        */
                             /* NV_Ith_S, prototypes for N_VNew, N_VFree      */
#include "dense.h"           /* definitions of type DenseMat, macro DENSE_ELEM*/


/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   3            /* number of equations  */
#define Y1    1.0          /* initial y components */
#define Y2    0.0
#define Y3    0.0
#define RTOL  1e-4         /* scalar relative tolerance            */
#define ATOL1 1e-8         /* vector absolute tolerance components */
#define ATOL2 1e-14
#define ATOL3 1e-6
#define T0    0.0          /* initial time           */
#define T1    0.4          /* first output time      */
#define TMULT 10.0         /* output time factor     */
#define NOUT  12           /* number of output times */


/* Private Helper Function */

static void PrintFinalStats(void *cvode_mem);


/* Functions Called by the CVODE Solver */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);


/***************************** Main Program ******************************/

int main()
{
  NV_Spec nvSpec;
  realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int iout, flag;

  nvSpec = NULL;
  y = abstol = NULL;
  cvode_mem = NULL;

  /* Initialize serial vector specification object */
  nvSpec = NV_SpecInit_Serial(NEQ);
  if (check_flag((void *)nvSpec, "NV_SpecInit", 0)) return(1);

  y = N_VNew(nvSpec);    /* Allocate y, abstol vectors */
  if (check_flag((void *)y, "N_VNew", 0)) return(1);
  abstol = N_VNew(nvSpec);
  if (check_flag((void *)abstol, "N_VNew", 0)) return(1);

  Ith(y,1) = Y1;               /* Initialize y */
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  reltol = RTOL;               /* Set the scalar relative tolerance */
  Ith(abstol,1) = ATOL1;       /* Set the vector absolute tolerance */
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* 
     Call CVodeCreate to create CVODES memory:

     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration

     A pointer to CVODES problem memory is returned and stored in cvode_mem.
  */
  cvode_mem = CVodeCreate(BDF, NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* 
     Call CVodeMalloc to initialize CVODES memory: 

     cvode_mem is the pointer to CVODES memory returned by CVodeCreate
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     y       is the initial dependent variable vector
     SV      specifies scalar relative and vector absolute tolerances
     &reltol is a pointer to the scalar relative tolerance
     abstol  is the absolute tolerance vector
     nvSpec  is the vector specification object 
  */
  flag = CVodeMalloc(cvode_mem, f, T0, y, SV, &reltol, abstol, nvSpec);
  if (check_flag(&flag, "CVodeMalloc", 1)) return(1);

  /* Call CVDense to specify the CVODE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine */
  flag = CVDenseSetJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDenseSetJacFn", 1)) return(1);

  /* In loop over output points, call CVode, print results, test for error */
  printf(" \n3-species kinetics problem\n\n");
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    check_flag(&flag, "CVode", 1);
    printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n",
           t, Ith(y,1), Ith(y,2), Ith(y,3));
    if (flag != SUCCESS) break;
  }

  PrintFinalStats(cvode_mem);  /* Print some final statistics   */

  /* Free y vector */
  N_VFree(y);

  /* Free abstol vector */
  N_VFree(abstol);

  /* Free CVODE problem memory */
  CVodeFree(cvode_mem);

  /* Free vector specification memory */
  NV_SpecFree_Serial(nvSpec);

  return(0);
}


/************************ Private Helper Function ************************/

/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, njeD, nfeD, nni, ncfn, netf;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDenseGetNumJacEvals(cvode_mem, &njeD);
  check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeD);
  check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeD = %-6ld njeD = %ld\n",
	 nst, nfe, nsetups, nfeD, njeD);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 nni, ncfn, netf);
}


/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,y). */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2, y3, yd1, yd3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  yd1 = Ith(ydot,1) = -0.04*y1 + 1e4*y2*y3;
  yd3 = Ith(ydot,3) = 3e7*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;
}

/* Jacobian routine. Compute J(t,y). */

static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = -0.04;  IJth(J,1,2) = 1e4*y3;          IJth(J,1,3) = 1e4*y2;
  IJth(J,2,1) =  0.04;  IJth(J,2,2) = -1e4*y3-6e7*y2;  IJth(J,2,3) = -1e4*y2;
                        IJth(J,3,2) = 6e7*y2;
}


/************************ Private Helper Function ************************/
 
/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag == SUCCESS
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag != SUCCESS */
  else if (opt == 1) {
    errflag = flagvalue;
    if (*errflag != SUCCESS) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
