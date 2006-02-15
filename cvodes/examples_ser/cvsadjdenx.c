/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-02-15 17:46:56 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
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
#include "nvector_serial.h"
#include "cvodes_dense.h"
#include "sundials_types.h"
#include "sundials_math.h"

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i= 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* (i,j)-th matrix component i,j = 1..NEQ */

/* Problem Constants */

#define NEQ      3             /* number of equations                  */

#define RTOL     RCONST(1e-6)  /* scalar relative tolerance            */

#define ATOL1    RCONST(1e-8)  /* vector absolute tolerance components */
#define ATOL2    RCONST(1e-14)
#define ATOL3    RCONST(1e-6)

#define ATOLl    RCONST(1e-8)  /* absolute tolerance for adjoint vars. */
#define ATOLq    RCONST(1e-6)  /* absolute tolerance for quadratures   */

#define T0       RCONST(0.0)   /* initial time                         */
#define TOUT     RCONST(4e7)   /* final time                           */

#define TB1      RCONST(4e7)   /* starting point for adjoint problem   */
#define TB2      RCONST(50.0)  /* starting point for adjoint problem   */

#define STEPS    150           /* number of steps between check points */

#define NP       3             /* number of problem parameters         */

#define ZERO     RCONST(0.0)


/* Type : UserData */

typedef struct {
  realtype p[3];
} *UserData;

/* Prototypes of user-supplied functions */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int Jac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data);
static int ewt(N_Vector y, N_Vector w, void *e_data);

static int fB(realtype t, N_Vector y, 
              N_Vector yB, N_Vector yBdot, void *f_dataB);
static int JacB(long int NB, DenseMat JB, realtype t,
                N_Vector y, N_Vector yB, N_Vector fyB, void *jac_dataB,
                N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static int fQB(realtype t, N_Vector y, N_Vector yB, 
               N_Vector qBdot, void *fQ_dataB);


/* Prototypes of private functions */

static void PrintOutput(N_Vector yB, N_Vector qB);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  void *cvadj_mem;
  void *cvode_mem;

  realtype reltolQ, abstolQ;
  N_Vector y, q;

  int steps;

  realtype reltolB, abstolB, abstolQB;
  N_Vector yB, qB;

  realtype time;
  int flag, ncheck;

  long int nst, nstB;

  CVadjCheckPointRec *ckpnt;
  int i;

  data = NULL;
  cvadj_mem = cvode_mem = NULL;
  y = yB = qB = NULL;

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
  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  /* Initialize y */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  Ith(y,1) = RCONST(1.0);
  Ith(y,2) = ZERO;
  Ith(y,3) = ZERO;

  /* Initialize q */
  q = N_VNew_Serial(1);
  if (check_flag((void *)q, "N_VNew_Serial", 0)) return(1);
  Ith(q,1) = ZERO;

  /* Set the scalar realtive and absolute tolerances reltolQ and abstolQ */
  reltolQ = RTOL;
  abstolQ = ATOLq;

  /* Create and allocate CVODES memory for forward run */
  printf("Create and allocate CVODES memory for forward runs\n");

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeMalloc(cvode_mem, f, T0, y, CV_WF, 0.0, NULL);
  if (check_flag(&flag, "CVodeMalloc", 1)) return(1);

  flag = CVodeSetEwtFn(cvode_mem, ewt, NULL);
  if (check_flag(&flag, "CVodeSetEwtFn", 1)) return(1);

  flag = CVodeSetFdata(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetFdata", 1)) return(1);

  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  flag = CVDenseSetJacFn(cvode_mem, Jac, data);
  if (check_flag(&flag, "CVDenseSetJacFn", 1)) return(1);

  flag = CVodeQuadMalloc(cvode_mem, fQ, q);
  if (check_flag(&flag, "CVodeQuadMalloc", 1)) return(1);

  flag = CVodeSetQuadFdata(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetQuadFdata", 1)) return(1);

  flag = CVodeSetQuadErrCon(cvode_mem, TRUE, CV_SS, reltolQ, &abstolQ);
  if (check_flag(&flag, "CVodeSetQuadErrCon", 1)) return(1);

  /* Allocate global memory */
  printf("Allocate global memory\n");

  steps = STEPS;
  cvadj_mem = CVadjMalloc(cvode_mem, steps, CV_HERMITE);
  /*
  cvadj_mem = CVadjMalloc(cvode_mem, steps, CV_POLYNOMIAL);
  */
  if (check_flag((void *)cvadj_mem, "CVadjMalloc", 0)) return(1);

  /* Perform forward run */
  printf("Forward integration ... ");
  
  flag = CVodeF(cvadj_mem, TOUT, y, &time, CV_NORMAL, &ncheck);
  if (check_flag(&flag, "CVodeF", 1)) return(1);
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_flag(&flag, "CVodeGetNumSteps", 1)) return(1);

  printf("done ( nst = %ld )   | ",nst);

  flag = CVodeGetQuad(cvode_mem, TOUT, q);
  if (check_flag(&flag, "CVodeGetQuad", 1)) return(1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("G: %12.4Le \n",Ith(q,1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("G: %12.4le \n",Ith(q,1));
#else
  printf("G: %12.4e \n",Ith(q,1));
#endif

  /* Test check point linked list */
  printf("\nList of Check Points (ncheck = %d)\n\n", ncheck);
  ckpnt = (CVadjCheckPointRec *) malloc ( (ncheck+1)*sizeof(CVadjCheckPointRec));
  CVadjGetCheckPointsInfo(cvadj_mem, ckpnt);
  for (i=0;i<=ncheck;i++) {
    printf("Address:       %u\n",ckpnt[i].my_addr);
    printf("Next:          %u\n",ckpnt[i].next_addr);
    printf("Time interval: %le  %le\n",ckpnt[i].t0, ckpnt[i].t1);
    printf("Step number:   %ld\n",ckpnt[i].nstep);
    printf("Order:         %d\n",ckpnt[i].order);
    printf("Step size:     %le\n",ckpnt[i].step);
    printf("\n");
  }

  /* Initialize yB */
  yB = N_VNew_Serial(NEQ);
  if (check_flag((void *)yB, "N_VNew_Serial", 0)) return(1);
  Ith(yB,1) = ZERO;
  Ith(yB,2) = ZERO;
  Ith(yB,3) = ZERO;

  /* Initialize qB */
  qB = N_VNew_Serial(NP);
  if (check_flag((void *)qB, "N_VNew", 0)) return(1);
  Ith(qB,1) = ZERO;
  Ith(qB,2) = ZERO;
  Ith(qB,3) = ZERO;

  /* Set the scalar relative tolerance reltolB */
  reltolB = RTOL;               

  /* Set the scalar absolute tolerance abstolB */
  abstolB = ATOLl;

  /* Set the scalar absolute tolerance abstolQB */
  abstolQB = ATOLq;

  /* Create and allocate CVODES memory for backward run */
  printf("\nCreate and allocate CVODES memory for backward run\n");

  flag = CVodeCreateB(cvadj_mem, CV_BDF, CV_NEWTON);
  if (check_flag(&flag, "CVodeCreateB", 1)) return(1);

  flag = CVodeMallocB(cvadj_mem, fB, TB1, yB, CV_SS, reltolB, &abstolB);
  if (check_flag(&flag, "CVodeMallocB", 1)) return(1);

  flag = CVodeSetFdataB(cvadj_mem, data);
  if (check_flag(&flag, "CVodeSetFdataB", 1)) return(1);

  flag = CVDenseB(cvadj_mem, NEQ);
  if (check_flag(&flag, "CVDenseB", 1)) return(1);

  flag = CVDenseSetJacFnB(cvadj_mem, JacB, data);
  if (check_flag(&flag, "CVDenseSetJacFnB", 1)) return(1);

  flag = CVodeQuadMallocB(cvadj_mem, fQB, qB);
  if (check_flag(&flag, "CVodeQuadMallocB", 1)) return(1);

  flag = CVodeSetQuadFdataB(cvadj_mem, data);
  if (check_flag(&flag, "CVodeSetQuadFdataB", 1)) return(1);

  flag = CVodeSetQuadErrConB(cvadj_mem, TRUE, CV_SS, reltolB, &abstolQB);
  if (check_flag(&flag, "CVodeSetQuadErrConB", 1)) return(1);

  /* Backward Integration */
  printf("Backward integration ... ");

  flag = CVodeB(cvadj_mem, T0, yB, &time, CV_NORMAL);
  if (check_flag(&flag, "CVodeB", 1)) return(1);
  CVodeGetNumSteps(CVadjGetCVodeBmem(cvadj_mem), &nstB);
  printf("done ( nst = %ld )\n", nstB);

  flag = CVodeGetQuadB(cvadj_mem, qB);
  if (check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  PrintOutput(yB, qB);

  /* Reinitialize backward phase (new tB0) */

  Ith(yB,1) = ZERO;
  Ith(yB,2) = ZERO;
  Ith(yB,3) = ZERO;

  Ith(qB,1) = ZERO;
  Ith(qB,2) = ZERO;
  Ith(qB,3) = ZERO;

  printf("Re-initialize CVODES memory for backward run\n");

  flag = CVodeReInitB(cvadj_mem, fB, TB2, yB, CV_SS, reltolB, &abstolB);
  if (check_flag(&flag, "CVodeReInitB", 1)) return(1);

  flag = CVodeQuadReInitB(cvadj_mem, fQB, qB); 
  if (check_flag(&flag, "CVodeQuadReInitB", 1)) return(1);

  printf("Backward integration ... ");

  flag = CVodeB(cvadj_mem, T0, yB, &time, CV_NORMAL);
  if (check_flag(&flag, "CVodeB", 1)) return(1);
  CVodeGetNumSteps(CVadjGetCVodeBmem(cvadj_mem), &nstB);
  printf("done ( nst = %ld )\n", nstB);

  flag = CVodeGetQuadB(cvadj_mem, qB);
  if (check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  PrintOutput(yB, qB);

  /* Free memory */
  printf("Free memory\n\n");

  CVodeFree(&cvode_mem);
  N_VDestroy_Serial(y); 
  N_VDestroy_Serial(q);
  N_VDestroy_Serial(yB);
  N_VDestroy_Serial(qB);
  CVadjFree(&cvadj_mem);
  free(data);

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y). 
*/

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
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

  return(0);
}

/* 
 * Jacobian routine. Compute J(t,y). 
*/

static int Jac(long int N, DenseMat J, realtype t,
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

  return(0);
}

/* 
 * fQ routine. Compute fQ(t,y). 
*/

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data)
{
  Ith(qdot,1) = Ith(y,3);  

  return(0);
}
 
/*
 * EwtSet function. Computes the error weights at the current solution.
 */

static int ewt(N_Vector y, N_Vector w, void *e_data)
{
  int i;
  realtype yy, ww, rtol, atol[3];

  rtol    = RTOL;
  atol[0] = ATOL1;
  atol[1] = ATOL2;
  atol[2] = ATOL3;

  for (i=1; i<=3; i++) {
    yy = Ith(y,i);
    ww = rtol * ABS(yy) + atol[i-1];  
    if (ww <= 0.0) return (-1);
    Ith(w,i) = 1.0/ww;
  }

  return(0);
}

/* 
 * fB routine. Compute fB(t,y,yB). 
*/

static int fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *f_dataB)
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
  Ith(yBdot,2) = p2*y3*l21 - RCONST(2.0)*p3*y2*l32;
  Ith(yBdot,3) = p2*y2*l21 - RCONST(1.0);

  return(0);
}

/* 
 * JacB routine. Compute JB(t,y,yB). 
*/

static int JacB(long int NB, DenseMat JB, realtype t,
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
  IJth(JB,2,1) = -p2*y3; IJth(JB,2,2) = p2*y3+2.0*p3*y2; IJth(JB,2,3) = RCONST(-2.0)*p3*y2;
  IJth(JB,3,1) = -p2*y2; IJth(JB,3,2) = p2*y2;

  return(0);
}

/*
 * fQB routine. Compute integrand for quadratures 
*/

static int fQB(realtype t, N_Vector y, N_Vector yB, 
               N_Vector qBdot, void *fQ_dataB)
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

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Print results after backward integration
 */

static void PrintOutput(N_Vector yB, N_Vector qB)
{
  printf("--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("tB0:        %12.4Le\n",TB1);
  printf("dG/dp:      %12.4Le %12.4Le %12.4Le\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
  printf("lambda(t0): %12.4Le %12.4Le %12.4Le\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("tB0:        %12.4le\n",TB1);
  printf("dG/dp:      %12.4le %12.4le %12.4le\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
  printf("lambda(t0): %12.4le %12.4le %12.4le\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
#else
  printf("tB0:        %12.4e\n",TB1);
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
#endif
  printf("--------------------------------------------------------\n\n");
}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

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
    errflag = (int *) flagvalue;
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
