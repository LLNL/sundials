/************************************************************************
 *                                                                      *
 * File       : chemadj.c                                               *
 * Programmers: Radu Serban @ LLNL                                      *
 * Version of : 22 March 2002                                           *
 *----------------------------------------------------------------------*
 * Adjoint sensitivity example problem.                                 *
 * The following is a simple example problem, with the coding           *
 * needed for its solution by CVODES.  The problem is from chemical     *
 * kinetics, and consists of the following three rate equations..       *
 *    dy1/dt = -p1*y1 + p2*y2*y3                                        *
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2                            *
 *    dy3/dt =  p3*(y2)^2                                               *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions   *
 * y1 = 1.0, y2 = y3 = 0.  The reaction rates are: p1=0.04, p2=1e4, and *
 * p3=3e7.  The problem is stiff.                                       *
 * This program solves the problem with the BDF method, Newton          *
 * iteration with the CVODE dense linear solver, and a user-supplied    *
 * Jacobian routine.                                                    * 
 * It uses a scalar relative tolerance and a vector absolute tolerance. *
 * Output is printed in decades from t = .4 to t = 4.e10.               *
 * Run statistics (optional outputs) are printed at the end.            *
 *                                                                      *
 * Optionally, CVODES can compute sensitivities with respect to the     *
 * problem parameters p1, p2, and p3 of the following quantity:         *
 *   G = int_t0^t1 g(t,p,y) dt                                          *
 * where                                                                *
 *   g(t,p,y) = y1 + p2 * y2 * y3                                       *
 *                                                                      *
 * The gradient dG/dp is obtained as:                                   *
 *   dG/dp = int_t0^t1 (g_p - lambda^T f_p ) dt - lambda^T(t0)*y0_p     *
 *         = - xi^T(t0) - lambda^T(t0)*y0_p                             *
 * where lambda and xi are solutions of:                                *
 *   d(lambda)/dt = - (f_y)^T * lambda + (g_y)^T                        *
 *   lambda(t1) = 0                                                     *
 * and                                                                  *
 *   d(xi)/dt = - (f_p)^T + (g_p)^T                                     *
 *   xi(t1) = 0                                                         *
 *                                                                      *
 * During the backward integration, CVODES also evaluates G as          *
 *   G = - phi(t0)                                                      *
 * where                                                                *
 *   d(phi)/dt = g(t,y,p)                                               *
 *   phi(t1) = 0                                                        *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"
#include "cvodea.h"
#include "cvsdense.h"
#include "nvector_serial.h"
#include "dense.h"


#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */
#define NEQ      3            /* number of equations                  */

#define RTOL     1e-4         /* scalar relative tolerance            */

#define ATOL1    1e-8         /* vector absolute tolerance components */
#define ATOL2    1e-14
#define ATOL3    1e-6

#define ATOLq    1e-6         /* absolute tol. for quadratures        */

#define T0       0.0          /* initial time                         */
#define TOUT     4e7          /* final time                           */

#define MAXSTP   150          /* number of steps between check points */

#define NP       3            /* number of problem parameters         */

#define ZERO     0.0

/* Type : UserData */
typedef struct {
  real p[3];
} *UserData;

/* Functions Called by the CVODE Solver */

/* f is of type RhsFn */
static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);
/* Jac is of type CVDenseJacFn */
static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);
/* fB is of type RhsFnB */
static void fB(integer NB, real t, N_Vector y, 
               N_Vector yB, N_Vector yBdot, void *f_dataB);
/* JacB is of type CVDenseJacFnB */
static void JacB(integer NB, DenseMat JB, RhsFnB fB, void *f_dataB, real t,
                 N_Vector y, N_Vector yB, N_Vector fyB, N_Vector ewtB, real hB, 
                 real uroundB, void *jac_dataB, long int *nfePtrB, N_Vector vtemp1B,
                 N_Vector vtemp2B, N_Vector vtemp3B);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  M_Env machEnvF, machEnvB;

  UserData data;

  void *cvadj_mem;
  void *cvode_mem;

  long int iopt[OPT_SIZE];
  real ropt[OPT_SIZE];
  real reltol;
  N_Vector y, abstol;

  real reltolB;
  N_Vector yB, abstolB;

  real time;
  int i, flag, ncheck;


  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = 0.04;
  data->p[1] = 1.0e4;
  data->p[2] = 3.0e7;

  /* Initialize serial machine environment for forward integration */
  machEnvF = M_EnvInit_Serial(NEQ);

  /* Initialize y */
  y = N_VNew(NEQ, machEnvF); 
  Ith(y,1) = 1.0;                
  Ith(y,2) = 0.0;
  Ith(y,3) = 0.0;

  /* Set the scalar relative tolerance reltol */
  reltol = RTOL;               
  /* Set the vector absolute tolerance abstol */
  abstol = N_VNew(NEQ, machEnvF); 
  Ith(abstol,1) = ATOL1;       
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* Optional I/O */
  for(i=0; i<OPT_SIZE; i++) {
    iopt[i] = 0;
    ropt[i] = 0.0;
  }

  /* Allocate CVODE memory for forward run */
  printf("\nAllocate CVODE memory for forward runs\n");
  cvode_mem = CVodeMalloc(NEQ, f, T0, y, BDF, NEWTON, SV, &reltol, abstol,
                          data, NULL, TRUE, iopt, ropt, machEnvF);
  flag = CVDense(cvode_mem, Jac, NULL);
  if (flag != SUCCESS) { printf("CVDense failed.\n"); return(1); }

  /* Allocate global memory */
  printf("\nAllocate global memory\n");
  cvadj_mem = CVadjMalloc(cvode_mem, MAXSTP);

  /* Perform forward run */
  printf("\nForward integration\n");
  
  flag = CVodeF(cvadj_mem, TOUT, y, &time, NORMAL, &ncheck);

  /* Test check point linked list */
  printf("\nList of Check Points (ncheck = %d)\n", ncheck);
  CVadjCheckPointsList(cvadj_mem);

  /* Initialize serial machine environment for backward run */ 
  machEnvB = M_EnvInit_Serial(NEQ+NP+1);

  /* Initialize yB */
  yB = N_VNew(NEQ+NP+1, machEnvB);
  Ith(yB,1) = 0.0;
  Ith(yB,2) = 0.0;
  Ith(yB,3) = 0.0;
  Ith(yB,4) = 0.0;
  Ith(yB,5) = 0.0;
  Ith(yB,6) = 0.0;
  Ith(yB,7) = 0.0;

  /* Set the scalar relative tolerance reltolB */
  reltolB = RTOL;               
  /* Set the vector absolute tolerance abstolB */
  abstolB = N_VNew(NEQ+NP+1, machEnvB); 
  Ith(abstolB,1) = ATOL1;       
  Ith(abstolB,2) = ATOL2;
  Ith(abstolB,3) = ATOL3;
  Ith(abstolB,4) = ATOLq;
  Ith(abstolB,5) = ATOLq;
  Ith(abstolB,6) = ATOLq;
  Ith(abstolB,7) = ATOLq;

  /* Allocate CVODE memory for backward run */
  printf("\nAllocate CVODE memory for backward run\n");
  flag = CVodeMallocB(cvadj_mem, NEQ+NP+1, fB, yB, BDF, NEWTON, SV, 
                      &reltolB, abstolB, data, NULL, 
                      FALSE, NULL, NULL, machEnvB);
  flag = CVDenseB(cvadj_mem, JacB, NULL);
  if (flag != SUCCESS) { printf("CVDenseB failed.\n"); return(1); }

  /* Backward Integration */
  flag = CVodeB(cvadj_mem, yB);

  printf("\n\n========================================================\n");
  printf("G:  %12.4e \n",-Ith(yB,7));
  printf("Gp: %12.4e %12.4e %12.4e\n", -Ith(yB,4), -Ith(yB,5), -Ith(yB,6));
  printf("========================================================\n");

  /* Free memory */
  printf("\nFree memory\n");
  CVodeFree(cvode_mem);
  CVadjFree(cvadj_mem);
  N_VFree(y);
  N_VFree(yB);
  N_VFree(abstol);
  N_VFree(abstolB);
  free(data);
  M_EnvFree_Serial(machEnvF);
  M_EnvFree_Serial(machEnvB);

  return(0);

}


/***************** Functions Called by the CVODE Solver ******************/

/*-------------------------------------------------------------------------
  f routine. Compute f(t,y). 
-------------------------------------------------------------------------*/

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)
{
  real y1, y2, y3, yd1, yd3;
  UserData data;
  real p1, p2, p3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) f_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  yd1 = Ith(ydot,1) = -p1*y1 + p2*y2*y3;
  yd3 = Ith(ydot,3) = p3*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;
}

/*-------------------------------------------------------------------------
  Jacobian routine. Compute J(t,y). 
-------------------------------------------------------------------------*/

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  real y1, y2, y3;
  UserData data;
  real p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) f_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
                      IJth(J,3,2) = 2*p3*y2;
}
 
/*-------------------------------------------------------------------------
  fB routine. Compute fB(t,y,mu). 
-------------------------------------------------------------------------*/

static void fB(integer NB, real t, N_Vector y, 
               N_Vector yB, N_Vector yBdot, void *f_dataB)
{
  UserData data;
  real y1, y2, y3;
  real p1, p2, p3;
  real l1, l2, l3;
  real x1, x2, x3;
  real phi;
  real l21, l32, y23;
  
  data = (UserData) f_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  /* The lambda vector */
  l1 = Ith(yB,1); l2 = Ith(yB,2); l3 = Ith(yB,3);

  /* The xi vector */
  x1 = Ith(yB,4); x2 = Ith(yB,5); x3 = Ith(yB,6);
  
  /* The phi variable */
  phi = Ith(yB, 7);

  /* Temporary variables */
  l21 = l2-l1;
  l32 = l3-l2;
  y23 = y2*y3;

  /* Load yBdot */
  Ith(yBdot,1) = -p1*l21 + 1.0;
  Ith(yBdot,2) = p2*y3*l21-2.0*p3*y2*l32+p2*y3;
  Ith(yBdot,3) = p2*y2*l21+p2*y2;
  Ith(yBdot,4) = -y1*l21;
  Ith(yBdot,5) = y23*l21+y23;
  Ith(yBdot,6) = -y2*y2*l32;
  Ith(yBdot,7) = y1+p2*y23;

}

/*-------------------------------------------------------------------------
  fB routine. Compute fB(t,y,mu). 
-------------------------------------------------------------------------*/

static void JacB(integer NB, DenseMat JB, RhsFnB fB, void *f_dataB, real t,
                 N_Vector y, N_Vector yB, N_Vector fyB, N_Vector ewtB, real hB, 
                 real uroundB, void *jac_dataB, long int *nfePtrB, N_Vector vtemp1B,
                 N_Vector vtemp2B, N_Vector vtemp3B)
{
  UserData data;
  real y1, y2, y3;
  real p1, p2, p3;
  
  data = (UserData) f_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  /* Load JB */
  IJth(JB,1,1) = p1;     IJth(JB,1,2) = -p1; 
  IJth(JB,2,1) = -p2*y3; IJth(JB,2,2) = p2*y3+2.0*p3*y2; IJth(JB,2,3) = -2.0*p3*y2;
  IJth(JB,3,1) = y1;     IJth(JB,3,2) = -y1; 
  IJth(JB,4,1) = -y2*y3; IJth(JB,4,2) = y2*y3; 
                         IJth(JB,5,2) = y2*y2;           IJth(JB,5,3) = -y2*y2; 
}

