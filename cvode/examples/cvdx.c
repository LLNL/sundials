/************************************************************************
 *                                                                      *
 * File: cvdx.c                                                         *
 * Programmers: Scott D. Cohen and Alan C. Hindmarsh @ LLNL             *
 * Version of 8 January 2002                                            *
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

#include "llnltyps.h" /* definitions of types real (set to double) and     */
                      /* integer (set to int), and the constant FALSE      */
#include "cvode.h"    /* prototypes for CVodeMalloc, CVode, and CVodeFree, */
                      /* constants OPT_SIZE, BDF, NEWTON, SV, SUCCESS,     */
                      /* NST, NFE, NSETUPS, NNI, NCFN, NETF                */
#include "cvdense.h"  /* prototype for CVDense, constant DENSE_NJE         */
#include "nvector.h"  /* definitions of type N_Vector and macro N_VIth,    */
                      /* prototypes for N_VNew, N_VFree                    */
#include "dense.h"    /* definitions of type DenseMat, macro DENSE_ELEM    */


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

#define Ith(v,i)    N_VIth(v,i-1)         /* Ith numbers components 1..NEQ */
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

static void PrintFinalStats(long int iopt[]);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);


/***************************** Main Program ******************************/

main()
{
  real ropt[OPT_SIZE], reltol, t, tout;
  long int iopt[OPT_SIZE];
  N_Vector y, abstol;
  void *cvode_mem;
  int iout, flag;

  y = N_VNew(NEQ, NULL);       /* Allocate y, abstol vectors */
  abstol = N_VNew(NEQ, NULL); 

  Ith(y,1) = Y1;               /* Initialize y */
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  reltol = RTOL;               /* Set the scalar relative tolerance */
  Ith(abstol,1) = ATOL1;       /* Set the vector absolute tolerance */
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* Call CVodeMalloc to initialize CVODE: 

     NEQ     is the problem size = number of equations
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     y       is the initial dependent variable vector
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration
     SV      specifies scalar relative and vector absolute tolerances
     &reltol is a pointer to the scalar relative tolerance
     abstol  is the absolute tolerance vector
     FALSE   indicates there are no optional inputs in iopt and ropt
     iopt    is an array used to communicate optional integer input and output
     ropt    is an array used to communicate optional real input and output

     A pointer to CVODE problem memory is returned and stored in cvode_mem. */

  cvode_mem = CVodeMalloc(NEQ, f, T0, y, BDF, NEWTON, SV, &reltol, abstol,
                          NULL, NULL, FALSE, iopt, ropt, NULL);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return(1); }

  /* Call CVDense to specify the CVODE dense linear solver with the
     user-supplied Jacobian routine Jac. */

  flag = CVDense(cvode_mem, Jac, NULL);
  if (flag != SUCCESS) { printf("CVDense failed.\n"); return(1); }

  /* In loop over output points, call CVode, print results, test for error */

  printf(" \n3-species kinetics problem\n\n");
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n",
           t, Ith(y,1), Ith(y,2), Ith(y,3));
    if (flag != SUCCESS) { printf("CVode failed, flag=%d.\n", flag); break; }
  }

  N_VFree(y);                  /* Free the y and abstol vectors */
  N_VFree(abstol);   
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */
  PrintFinalStats(iopt);       /* Print some final statistics   */
  return(0);
}


/************************ Private Helper Function ************************/

/* Print some final statistics located in the iopt array */

static void PrintFinalStats(long int iopt[])
{
  printf("\nFinal Statistics.. \n\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
	 iopt[NST], iopt[NFE], iopt[NSETUPS], iopt[DENSE_NJE]);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 iopt[NNI], iopt[NCFN], iopt[NETF]);
}


/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,y). */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)
{
  real y1, y2, y3, yd1, yd3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  yd1 = Ith(ydot,1) = -0.04*y1 + 1e4*y2*y3;
  yd3 = Ith(ydot,3) = 3e7*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;
}

/* Jacobian routine. Compute J(t,y). */

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  real y1, y2, y3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = -0.04;  IJth(J,1,2) = 1e4*y3;          IJth(J,1,3) = 1e4*y2;
  IJth(J,2,1) =  0.04;  IJth(J,2,2) = -1e4*y3-6e7*y2;  IJth(J,2,3) = -1e4*y2;
                        IJth(J,3,2) = 6e7*y2;
}
 
