/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2004-05-26 18:37:39 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem.
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of the following three rate equations.
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The reaction rates are:
 * p1=0.04, p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODE dense linear solver, and a user-supplied
 * Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * 
 * Optionally, CVODES can compute sensitivities with respect to
 * the problem parameters p1, p2, and p3 of the following quantity:
 *   G = int_t0^t1 g(t,p,y) dt
 * where
 *   g(t,p,y) = y3
 *        
 * The gradient dG/dp is obtained as:
 *   dG/dp = int_t0^t1 (g_p - lambda^T f_p ) dt - lambda^T(t0)*y0_p
 *         = - xi^T(t0) - lambda^T(t0)*y0_p
 * where lambda and xi are solutions of:
 *   d(lambda)/dt = - (f_y)^T * lambda - (g_y)^T
 *   lambda(t1) = 0
 * and
 *   d(xi)/dt = - (f_p)^T * lambda + (g_p)^T
 *   xi(t1) = 0
 * 
 * During the backward integration, CVODES also evaluates G as
 *   G = - phi(t0)
 * where
 *   d(phi)/dt = g(t,y,p)
 *   phi(t1) = 0
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "cvodes.h"
#include "cvodea.h"
#include "cvdense.h"
#include "nvector_serial.h"
#include "sundialstypes.h"

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ      3            /* number of equations                  */

#define RTOL     1e-4         /* scalar relative tolerance            */

#define ATOL1    1e-8         /* vector absolute tolerance components */
#define ATOL2    1e-14
#define ATOL3    1e-6

#define ATOLl    1e-5         /* absolute tolerance for adjoint vars. */
#define ATOLq    1e-6         /* absolute tolerance for quadratures   */

#define T0       0.0          /* initial time                         */
#define TOUT     4e7          /* final time                           */

#define TB1      4e7          /* starting point for adjoint problem   */
#define TB2      50.0         /* starting point for adjoint problem   */

#define STEPS    150          /* number of steps between check points */

#define NP       3            /* number of problem parameters         */

#define ZERO     0.0


/* Type : UserData */

typedef struct {
  realtype p[3];
} *UserData;


/* Functions Called by the CVODES Solver */

/* f is of type RhsFn */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

/* Jac is of type CVDenseJacFn */

static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* fQ is of type QuadRhsFn */

static void fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data);

/* fB is of type RhsFnB */

static void fB(realtype t, N_Vector y, 
               N_Vector yB, N_Vector yBdot, void *f_dataB);

/* JacB is of type CVDenseJacFnB */

static void JacB(long int NB, DenseMat JB, realtype t,
                 N_Vector y, N_Vector yB, N_Vector fyB, void *jac_dataB,
                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/* fQB is of type QuadRhsFnB */

static void fQB(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *fQ_dataB);


/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  NV_Spec nvSpecF, nvSpecQF, nvSpecB, nvSpecQB;

  UserData data;

  void *cvadj_mem;
  void *cvode_mem;

  realtype reltol, abstolQ;
  N_Vector y, q, abstol;

  realtype reltolB, abstolB, abstolQB;
  N_Vector yB, qB;

  realtype time;
  int flag, ncheck;

  nvSpecF = nvSpecQF = nvSpecB = nvSpecQB = NULL;
  data = NULL;
  cvadj_mem = cvode_mem = NULL;
  y = abstol = yB = qB = NULL;

  /* Print problem description */
  printf("\n\n Adjoint Sensitivity Example for Chemical Kinetics\n");
  printf(" -------------------------------------------------\n\n");
  printf("ODE: dy1/dt = -p1*y1 + p2*y2*y3\n");
  printf("     dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2\n");
  printf("     dy3/dt =  p3*(y2)^2\n\n");
  printf("Find dG/dp for\n");
  printf("     G = int_t0^tB0 g(t,p,y) dt\n");
  printf("     g(t,p,y) = y3\n\n\n");


  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2)) return(1);
  data->p[0] = 0.04;
  data->p[1] = 1.0e4;
  data->p[2] = 3.0e7;

  /* Initialize serial vector specification for forward integration */
  nvSpecF = NV_SpecInit_Serial(NEQ);
  if (check_flag((void *)nvSpecF, "NV_SpecInit", 0)) return(1);
  nvSpecQF = NV_SpecInit_Serial(1);
  if (check_flag((void *)nvSpecQF, "NV_SpecInit", 0)) return(1);

  /* Initialize y */
  y = N_VNew(nvSpecF);
  if (check_flag((void *)y, "N_VNew", 0)) return(1);
  Ith(y,1) = 1.0;
  Ith(y,2) = 0.0;
  Ith(y,3) = 0.0;

  /* Allocate q */
  q = N_VNew(nvSpecQF);
  if (check_flag((void *)q, "N_VNew", 0)) return(1);

  /* Set the scalar relative tolerance reltol */
  reltol = RTOL;    

  /* Set the vector absolute tolerance abstol */
  abstol = N_VNew(nvSpecF);
  if (check_flag((void *)abstol, "N_VNew", 0)) return(1);
  Ith(abstol,1) = ATOL1;       
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* Set the scalar absolute tolerance abstolQ */
  abstolQ = ATOLq;

  /* Create and allocate CVODES memory for forward run */
  printf("Create and allocate CVODES memory for forward runs\n");

  cvode_mem = CVodeCreate(BDF, NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeSetFdata(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetFdata", 1)) return(1);

  flag = CVodeMalloc(cvode_mem, f, T0, y, SV, &reltol, abstol, nvSpecF);
  if (check_flag(&flag, "CVodeMalloc", 1)) return(1);

  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);
  flag = CVDenseSetJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDenseSetJacFn", 1)) return(1);
  flag = CVDenseSetJacData(cvode_mem, data);
  if (check_flag(&flag, "CVDenseSetJacData", 1)) return(1);

  flag = CVodeSetQuadFdata(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetQuadFdata", 1)) return(1);
  flag = CVodeSetQuadErrCon(cvode_mem, TRUE);
  if (check_flag(&flag, "CVodeSetQuadErrCon", 1)) return(1);
  flag = CVodeSetQuadTolerances(cvode_mem, SS, &reltol, &abstolQ); 
  if (check_flag(&flag, "CVodeSetQuadTolerances", 1)) return(1);
  flag = CVodeQuadMalloc(cvode_mem, fQ, nvSpecQF);
  if (check_flag(&flag, "CVodeQuadMalloc", 1)) return(1);

  /* Allocate global memory */
  printf("Allocate global memory\n");
  cvadj_mem = CVadjMalloc(cvode_mem, STEPS);
  if (check_flag((void *)cvadj_mem, "CVadjMalloc", 0)) return(1);

  /* Perform forward run */
  printf("Forward integration ... ");
  
  flag = CVodeF(cvadj_mem, TOUT, y, &time, NORMAL, &ncheck);
  if (check_flag(&flag, "CVodeF", 1)) return(1);

  flag = CVodeGetQuad(cvode_mem, TOUT, q);
  if (check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  printf("G: %12.4e \n",Ith(q,1));

  /* Test check point linked list */
  printf("\nList of Check Points (ncheck = %d)\n", ncheck);
  CVadjGetCheckPointsList(cvadj_mem);

  /* Initialize serial nvector specification for backward run */ 
  nvSpecB  = NV_SpecInit_Serial(NEQ);   /* adjoint variables */
  if (check_flag((void *)nvSpecB, "NV_SpecInit", 0)) return(1);
  nvSpecQB = NV_SpecInit_Serial(NP);    /* quadrature variables */
  if (check_flag((void *)nvSpecQB, "NV_SpecInit", 0)) return(1);

  /* Initialize yB */
  yB = N_VNew(nvSpecB);
  if (check_flag((void *)yB, "N_VNew", 0)) return(1);
  Ith(yB,1) = 0.0;
  Ith(yB,2) = 0.0;
  Ith(yB,3) = 0.0;

  /* Initialize qB */
  qB = N_VNew(nvSpecQB);
  if (check_flag((void *)qB, "N_VNew", 0)) return(1);

  /* Set the scalar relative tolerance reltolB */
  reltolB = RTOL;               

  /* Set the scalar absolute tolerance abstolB */
  abstolB = ATOLl;

  /* Set the scalar absolute tolerance abstolQB */
  abstolQB = ATOLq;

  /* Create and allocate CVODES memory for backward run */
  printf("\nCreate and allocate CVODES memory for backward run\n");
  flag = CVodeCreateB(cvadj_mem, BDF, NEWTON);
  if (check_flag(&flag, "CVodeCreateB", 1)) return(1);
  flag = CVodeSetFdataB(cvadj_mem, data);
  if (check_flag(&flag, "CVodeSetFdataB", 1)) return(1);
  flag = CVodeMallocB(cvadj_mem, fB, TB1, yB, SS, &reltolB, &abstolB, nvSpecB);
  if (check_flag(&flag, "CVodeMallocB", 1)) return(1);

  flag = CVDenseB(cvadj_mem, NEQ);
  if (check_flag(&flag, "CVDenseB", 1)) return(1);
  flag = CVDenseSetJacFnB(cvadj_mem, JacB);
  if (check_flag(&flag, "CVDenseSetJacFnB", 1)) return(1);
  flag = CVDenseSetJacDataB(cvadj_mem, data);
  if (check_flag(&flag, "CVDenseSetJacDataB", 1)) return(1);

  flag = CVodeSetQuadFdataB(cvadj_mem, data);
  if (check_flag(&flag, "CVodeSetQuadFdataB", 1)) return(1);
  flag = CVodeSetQuadErrConB(cvadj_mem, TRUE);
  if (check_flag(&flag, "CVodeSetQuadErrConB", 1)) return(1);
  flag = CVodeSetQuadTolerancesB(cvadj_mem, SS, &reltolB, &abstolQB); 
  if (check_flag(&flag, "CVodeSetQuadTolerancesB", 1)) return(1);
  flag = CVodeQuadMallocB(cvadj_mem, fQB,nvSpecQB);
  if (check_flag(&flag, "CVodeQuadMallocB", 1)) return(1);

  /* Backward Integration */
  printf("Integrate backwards from tB0 = %12.4e\n", TB1);
  flag = CVodeB(cvadj_mem, yB);
  if (check_flag(&flag, "CVodeB", 1)) return(1);

  flag = CVodeGetQuadB(cvadj_mem, qB);
  if (check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  printf("--------------------------------------------------------\n");
  printf("tB0:        %12.4e \n",TB1);
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("--------------------------------------------------------\n\n");

  /* Reinitialize backward phase (new tB0) */
  Ith(yB,1) = 0.0;
  Ith(yB,2) = 0.0;
  Ith(yB,3) = 0.0;

  printf("Re-initialize CVODES memory for backward run\n");
  flag = CVodeReInitB(cvadj_mem, fB, TB2, yB, SS, &reltolB, &abstolB);
  if (check_flag(&flag, "CVodeReInitB", 1)) return(1);
  flag = CVodeQuadReInitB(cvadj_mem, fQB); 
  if (check_flag(&flag, "CVodeQuadReInitB", 1)) return(1);

  /* Backward Integration */
  printf("Integrate backwards from tB0 = %12.4e\n", TB2);
  flag = CVodeB(cvadj_mem, yB);
  if (check_flag(&flag, "CVodeB", 1)) return(1);

  flag = CVodeGetQuadB(cvadj_mem, qB);
  if (check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  printf("--------------------------------------------------------\n");
  printf("tB0:        %12.4e \n",TB2);
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("--------------------------------------------------------\n\n");

  /* Free memory */
  printf("Free memory\n\n");
  CVodeFree(cvode_mem);
  N_VFree(abstol);
  N_VFree(y); 
  N_VFree(q);
  N_VFree(yB);
  N_VFree(qB);
  CVadjFree(cvadj_mem);
  NV_SpecFree_Serial(nvSpecF);
  NV_SpecFree_Serial(nvSpecQF);
  NV_SpecFree_Serial(nvSpecB);
  NV_SpecFree_Serial(nvSpecQB);
  free(data);

  return(0);

}


/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,y). */
static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2, y3, yd1, yd3;
  UserData data;
  realtype p1, p2, p3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) f_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  yd1 = Ith(ydot,1) = -p1*y1 + p2*y2*y3;
  yd3 = Ith(ydot,3) = p3*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;
}

/* Jacobian routine. Compute J(t,y). */
static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) jac_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
                      IJth(J,3,2) = 2*p3*y2;
}

/* fQ routine. Compute fQ(t,y). */
static void fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data)
{

  Ith(qdot,1) = Ith(y,3);
  
}
 
/* fB routine. Compute fB(t,y,yB). */
static void fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *f_dataB)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2, l3;
  realtype l21, l32, y23;
  
  data = (UserData) f_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  /* The lambda vector */
  l1 = Ith(yB,1); l2 = Ith(yB,2); l3 = Ith(yB,3);

  /* Temporary variables */
  l21 = l2-l1;
  l32 = l3-l2;
  y23 = y2*y3;

  /* Load yBdot */
  Ith(yBdot,1) = - p1*l21;
  Ith(yBdot,2) = p2*y3*l21 - 2.0*p3*y2*l32;
  Ith(yBdot,3) = p2*y2*l21 - 1.0;
}

/* JacB routine. Compute JB(t,y,yB). */
static void JacB(long int NB, DenseMat JB, realtype t,
                 N_Vector y, N_Vector yB, N_Vector fyB, void *jac_dataB,
                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  
  data = (UserData) jac_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  /* Load JB */
  IJth(JB,1,1) = p1;     IJth(JB,1,2) = -p1; 
  IJth(JB,2,1) = -p2*y3; IJth(JB,2,2) = p2*y3+2.0*p3*y2; IJth(JB,2,3) = -2.0*p3*y2;
  IJth(JB,3,1) = -p2*y2; IJth(JB,3,2) = p2*y2;
}

/* fQB routine. Compute integrand for quadratures */
static void fQB(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *fQ_dataB)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2, l3;
  realtype l21, l32, y23;

  data = (UserData) fQ_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  /* The lambda vector */
  l1 = Ith(yB,1); l2 = Ith(yB,2); l3 = Ith(yB,3);

  /* Temporary variables */
  l21 = l2-l1;
  l32 = l3-l2;
  y23 = y2*y3;

  Ith(qBdot,1) = y1*l21;
  Ith(qBdot,2) = - y23*l21;
  Ith(qBdot,3) = y2*y2*l32;
}


/*********************** Private Helper Function ************************/

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
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

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = flagvalue;
    if (*errflag < 0) {
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
