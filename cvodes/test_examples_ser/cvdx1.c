/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-01-11 21:13:52 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Modification of the cvdx example to illustrate enforcing a 
 * minimum step size.
 *
 * When HMIN is reached (CV_ERR_FAILURE or CV_CONV_FAILURE return 
 * from CVode), looser tolerances are estimated such that the solver
 * can proceed with order 1 without the step size going below HMIN.
 * The absolute tolerance estimates are based on the current value 
 * of y'' (see comments for function loosen_tol).
 * While proceeding with looser tolerances, we check for sufficient
 * increase in the step size and when it reaches H_FACT*HMIN we 
 * tighten the tolerances back to their original values and reset
 * the maximum order to 5 (BDF).
 * Note that one could tighten the tolerances in stages and continue
 * with backward Euler until the tolerances reached their original
 * values.
 *
 * For readibility, no checks are performed on the various function
 * return flags.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes.h"
#include "nvector_serial.h"
#include "cvodes_dense.h"
#include "sundials_types.h"

/* User-defined vector and matrix accessor macros: Ith, IJth */
#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Problem Constants */
#define NEQ   3
#define Y1    RCONST(1.0)
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)
#define ATOL1 RCONST(1.0e-10)
#define ATOL2 RCONST(1.0e-14)
#define ATOL3 RCONST(1.0e-10)
#define T0    RCONST(0.0)
#define TOUT1 RCONST(0.4)
#define TOUT2 RCONST(400.0)

#define HMIN   RCONST(1.0e-2)
#define H_FACT RCONST(50.0)
#define A_FACT RCONST(10.0)

#define TIGHT 1
#define LOOSE 2

/* Type : UserData */
typedef struct {
  realtype p[3];        /* Problem parameters */
  N_Vector yp, yy, yyp; /* Temporary space    */
} *UserData;

/* Functions called by the solver */
static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions */
static UserData set_data();
static void free_data(UserData data);
static realtype loosen_tol(void *cvode_mem, realtype t, N_Vector y, 
                           UserData data, N_Vector atolL);
static void PrintOutput(void *cvode_mem, realtype t, realtype y1, realtype y2, realtype y3);
static void PrintSolverStats(void *cvode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  UserData data;
  realtype reltol, t, hu, alpha;
  N_Vector y, atolT, atolL;
  void *cvode_mem;
  int flag, mode;

  /* Allocate and initialize user data structure */
  data = set_data();

  /* Create serial vector of length NEQ for I.C., atolT, and atolL */
  y = N_VNew_Serial(NEQ);
  atolT = N_VNew_Serial(NEQ); 
  atolL = N_VNew_Serial(NEQ);

  /* Initialize y */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  /* Set the scalar relative tolerance and vector absolute tolerance */
  reltol = RTOL;
  Ith(atolT,1) = ATOL1;
  Ith(atolT,2) = ATOL2;
  Ith(atolT,3) = ATOL3;

  /* Call CVodeCreate and CVodeMalloc */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeMalloc(cvode_mem, f, T0, y, CV_SV, reltol, atolT);
  flag = CVodeSetFdata(cvode_mem, data);

  /* Use the CVDENSE dense linear solver with user-supplied Jacobian */
  flag = CVDense(cvode_mem, NEQ);
  flag = CVDenseSetJacFn(cvode_mem, Jac, data);

  /* Integrate to TOUT1 in NORMAL_MODE */
  flag = CVode(cvode_mem, TOUT1, y, &t, CV_NORMAL);
  PrintOutput(cvode_mem, t, Ith(y,1), Ith(y,2), Ith(y,3));

  /* Now set a minimum step size */
  printf("\nEnforce minimum step\n\n");
  flag = CVodeSetMinStep(cvode_mem, HMIN);

  flag = CVodeSetErrFile(cvode_mem, NULL);

  /* Continue to TOUT2 in ONE_STEP mode, while controlling tolerances
     to enforce the desired minimum step size */

  mode = TIGHT;

  while(t < TOUT2) {

    flag = CVode(cvode_mem, TOUT2, y, &t, CV_ONE_STEP);
    
    /* Test if the solver tries to go below HMIN */
    if (flag == CV_ERR_FAILURE || flag == CV_CONV_FAILURE) {
      PrintSolverStats(cvode_mem);
      /* Estimate relaxed tolerances */
      alpha = loosen_tol(cvode_mem, t, y, data, atolL);
      printf("Loosen tolerances:  alpha = %14.6le\n\n", alpha);
      /* Enforce backward Euler */
      CVodeSetMaxOrd(cvode_mem, 1);
      /* Reinitialize solver with loose tolerances */
      flag = CVodeReInit(cvode_mem, f, t, y, CV_SV, reltol, atolL);
      mode = LOOSE;
      /* Keep going */
      continue;
    }

    /* Check that no other error occured */
    if (flag < 0) {
      printf("\nAn error occured!!!\n\n");
      break;
    }

    PrintOutput(cvode_mem, t, Ith(y,1), Ith(y,2), Ith(y,3));

    if (mode == LOOSE) {
      flag = CVodeGetLastStep(cvode_mem, &hu);
      /* Test for sufficient increase in step size */
      if (hu > H_FACT*HMIN) {
        PrintSolverStats(cvode_mem);
        printf("Tighten tolerances\n\n");
        /* Reset maximum order */
        CVodeSetMaxOrd(cvode_mem, 5);
        /* Reinitialize solver with tight tolerances */
        flag = CVodeReInit(cvode_mem, f, t, y, CV_SV, reltol, atolT);
        mode = TIGHT;
        /* Keep going */
        continue;
      }
    }
  }

  PrintSolverStats(cvode_mem);

  /* Get solution at final time */
  t = TOUT2;
  CVodeGetDky(cvode_mem, t, 0, y);
  PrintOutput(cvode_mem, t, Ith(y,1), Ith(y,2), Ith(y,3));
  
  /* Free memory */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(atolT);
  N_VDestroy_Serial(atolL);
  free_data(data);
  CVodeFree(&cvode_mem);

  return(0);
}

/*
 * set_data
 *
 * Allocate and set user data structure
 */

static UserData set_data()
{
  UserData data;

  data = (UserData) malloc(sizeof *data);

  data->yp  = N_VNew_Serial(NEQ);
  data->yy  = N_VNew_Serial(NEQ);
  data->yyp = N_VNew_Serial(NEQ);

  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  return(data);
}

/*
 * free_data
 *
 * Deallocate user data structure
 */

static void free_data(UserData data)
{
  N_VDestroy_Serial(data->yp);
  N_VDestroy_Serial(data->yy);
  N_VDestroy_Serial(data->yyp);
  free(data);
}

/*
 * loosen_tol
 *
 * Estimate looser absolute tolerances so that the solver can proceed with
 * maximum order 1 (backward Euler) without the step size going below h_min
 *
 *                   (h_min)^2
 * atol_i = A_FACT * --------- * |y''_i|
 *                       2
 *
 * This ensures that ||y''||_WRMS <= 1 (which is essentially the error test
 * for backward Euler).
 * 
 * A_FACT is a safety factor to prevent having to loosen the tolerances
 * at two consecutive steps.
 *
 * If the current order is greater than 1, we obtain y'' from the BDF 
 * interpolant (calling CVodeGetDky). If already at order 1, we estimate
 * y'' using a difference quotient, the same way it is estimated in CVODES
 * to compute an initial step size (see CVYddNorm in cvodes.c).
 */

static realtype loosen_tol(void *cvode_mem, realtype t, N_Vector y, UserData data, 
                           N_Vector atolL)
{
  int flag, qu;
  realtype alpha;
  N_Vector yp, yy, yyp;

  flag = CVodeGetLastOrder(cvode_mem, &qu);

  if (qu > 1) {

    /* atolL <- y'' */
    flag = CVodeGetDky(cvode_mem, t, 2, atolL);

  } else {

    yp = data->yp;
    yy = data->yy;
    yyp = data->yyp;
    
    /* yp <- y'(t) */
    flag = CVodeGetDky(cvode_mem, t, 1, yp);

    /* yy <- h_min * y'(t) + y(t) */
    N_VLinearSum(HMIN, yp, 1.0, y, yy);

    /* yyp <- f( t+h_min, h*y'(t)+y(t) ) */
    f(t+HMIN, yy, yyp, data);

    /* atolL <- f(t+h, h*y'(t)+y(t)) - y'(t) */
    N_VLinearSum(1.0, yyp, -1.0, yp, atolL);

    /* atolL <- y'' */
    N_VScale(1.0/HMIN, atolL, atolL);

  }

  alpha = A_FACT * (HMIN*HMIN/2.0);

  /* atolL <- |y''| */
  N_VAbs(atolL, atolL);

  /* atolL <- alpha * |y''| */
  N_VScale(alpha, atolL, atolL);

  return(alpha);

}

/*
 * f
 *
 * Compute RHS function f(t,y). 
 */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3, yd1, yd3;

  data = (UserData) f_data;
  p1 = data->p[0]; 
  p2 = data->p[1]; 
  p3 = data->p[2];

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  yd1 = Ith(ydot,1) = -p1*y1 + p2*y2*y3;
  yd3 = Ith(ydot,3) = p3*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;
}

/*
 * Jac
 *
 * Compute Jacobian J(t,y) = df/dy. *
 */

static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;

  data = (UserData) jac_data;
  p1 = data->p[0]; 
  p2 = data->p[1]; 
  p3 = data->p[2];

  y1 = Ith(y,1);
  y2 = Ith(y,2);
  y3 = Ith(y,3);

  IJth(J,1,1) = -p1;
  IJth(J,1,2) = p2*y3;
  IJth(J,1,3) = p2*y2;

  IJth(J,2,1) =  p1;
  IJth(J,2,2) = -p2*y3-2*p3*y2;
  IJth(J,2,3) = -p2*y2;

  IJth(J,3,2) = 2*p3*y2;
}

/*
 * PrintOutput
 *
 * Print current time, order, step size, and solution
 */

static void PrintOutput(void *cvode_mem, realtype t, realtype y1, realtype y2, realtype y3)
{
  int flag, qu;
  realtype hu;

  flag = CVodeGetLastOrder(cvode_mem, &qu);
  flag = CVodeGetLastStep(cvode_mem, &hu);

  printf("At t = %0.4le    | qu = %1d  hu =%14.6le |  y =%14.6le  %14.6le  %14.6le\n", 
         t, qu, hu, y1, y2, y3);

  return;
}

/* 
 * PrinsolverStats
 *
 * Get and print some solver statistics
 */

static void PrintSolverStats(void *cvode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  int flag;

  flag = CVodeGetActualInitStep(cvode_mem, &h0u);

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  printf("\nSolver Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nni = %-6ld ncfn = %-6ld netf = %-6ld\n",
	 nst, nfe, nsetups, nni, ncfn, netf);
  printf("h0u = %14.6le\n\n",h0u);
}

