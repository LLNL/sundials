/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system from.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 4 from Section 14.1, Chapter 14: Ferraris and Tronconi
 *
 * This problem involves a blend of trigonometric and exponential terms.
 *    0.5 sin(x1 x2) - 0.25 x2/pi - 0.5 x1 = 0
 *    (1-0.25/pi) ( exp(2 x1)-e ) + e x2 / pi - 2 e x1 = 0
 * such that
 *    0.25 <= x1 <=1.0
 *    1.5 <= x2 <= 2 pi
 *
 * The treatment of the bound constraints on x1 and x2 is done using
 * the additional variables
 *    l1 = x1 - x1_min >= 0
 *    L1 = x1 - x1_max <= 0
 *    l2 = x2 - x2_min >= 0
 *    L2 = x2 - x2_max >= 0
 *
 * and using the constraint feature in KINSOL to impose
 *    l1 >= 0    l2 >= 0
 *    L1 <= 0    L2 <= 0
 *
 * The Ferraris-Tronconi test problem has two known solutions.
 * The nonlinear system is solved by KINSOL using different
 * combinations of globalization and Jacobian update strategies
 * and with different initial guesses (leading to one or the other
 * of the known solutions).
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */

/* Problem Constants */

#define NVAR   2
#define NEQ    3*NVAR

#define FTOL   RCONST(1.e-5) /* function tolerance */
#define STOL   RCONST(1.e-5) /* step tolerance     */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)

#define PI     RCONST(3.1415926)
#define E      RCONST(2.7182818)

typedef struct {
  realtype lb[NVAR];
  realtype ub[NVAR];
} *UserData;

/* Accessor macro */
#define Ith(v,i)    NV_Ith_S(v,i-1)

/* Functions Called by the KINSOL Solver */
static int func(N_Vector u, N_Vector f, void *user_data);

/* Private Helper Functions */
static void SetInitialGuess1(N_Vector u, UserData data);
static void SetInitialGuess2(N_Vector u, UserData data);
static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset);
static void PrintHeader(realtype fnormtol, realtype scsteptol);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void *kmem);
static int check_retval(void *retvalvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  SUNContext sunctx;
  UserData data;
  realtype fnormtol, scsteptol;
  N_Vector u1, u2, u, s, c;
  int glstr, mset, retval;
  void *kmem;
  SUNMatrix J;
  SUNLinearSolver LS;

  u1 = u2 = u = NULL;
  s = c = NULL;
  kmem = NULL;
  J = NULL;
  LS = NULL;
  data = NULL;

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* User data */

  data = (UserData)malloc(sizeof *data);
  data->lb[0] = PT25;       data->ub[0] = ONE;
  data->lb[1] = ONEPT5;     data->ub[1] = TWO*PI;

  /* Create serial vectors of length NEQ */
  u1 = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)u1, "N_VNew_Serial", 0)) return(1);

  u2 = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)u2, "N_VNew_Serial", 0)) return(1);

  u = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)u, "N_VNew_Serial", 0)) return(1);

  s = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)s, "N_VNew_Serial", 0)) return(1);

  c = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)c, "N_VNew_Serial", 0)) return(1);

  SetInitialGuess1(u1,data);
  SetInitialGuess2(u2,data);

  N_VConst(ONE,s); /* no scaling */

  Ith(c,1) =  ZERO;   /* no constraint on x1 */
  Ith(c,2) =  ZERO;   /* no constraint on x2 */
  Ith(c,3) =  ONE;    /* l1 = x1 - x1_min >= 0 */
  Ith(c,4) = -ONE;    /* L1 = x1 - x1_max <= 0 */
  Ith(c,5) =  ONE;    /* l2 = x2 - x2_min >= 0 */
  Ith(c,6) = -ONE;    /* L2 = x2 - x22_min <= 0 */

  fnormtol=FTOL; scsteptol=STOL;


  kmem = KINCreate(sunctx);
  if (check_retval((void *)kmem, "KINCreate", 0)) return(1);

  retval = KINSetUserData(kmem, data);
  if (check_retval(&retval, "KINSetUserData", 1)) return(1);
  retval = KINSetConstraints(kmem, c);
  if (check_retval(&retval, "KINSetConstraints", 1)) return(1);
  retval = KINSetFuncNormTol(kmem, fnormtol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);
  retval = KINSetScaledStepTol(kmem, scsteptol);
  if (check_retval(&retval, "KINSetScaledStepTol", 1)) return(1);

  retval = KINInit(kmem, func, u);
  if (check_retval(&retval, "KINInit", 1)) return(1);

  /* Create dense SUNMatrix */
  J = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if(check_retval((void *)J, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(u, J, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver to KINSOL */
  retval = KINSetLinearSolver(kmem, LS, J);
  if(check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

  /* Print out the problem size, solution parameters, initial guess. */
  PrintHeader(fnormtol, scsteptol);

  /* --------------------------- */

  printf("\n------------------------------------------\n");
  printf("\nInitial guess on lower bounds\n");
  printf("  [x1,x2] = ");
  PrintOutput(u1);

  N_VScale(ONE,u1,u);
  glstr = KIN_NONE;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale(ONE,u1,u);
  glstr = KIN_LINESEARCH;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale(ONE,u1,u);
  glstr = KIN_NONE;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale(ONE,u1,u);
  glstr = KIN_LINESEARCH;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);



  /* --------------------------- */

  printf("\n------------------------------------------\n");
  printf("\nInitial guess in middle of feasible region\n");
  printf("  [x1,x2] = ");
  PrintOutput(u2);

  N_VScale(ONE,u2,u);
  glstr = KIN_NONE;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale(ONE,u2,u);
  glstr = KIN_LINESEARCH;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale(ONE,u2,u);
  glstr = KIN_NONE;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale(ONE,u2,u);
  glstr = KIN_LINESEARCH;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);




  /* Free memory */

  N_VDestroy(u1);
  N_VDestroy(u2);
  N_VDestroy(u);
  N_VDestroy(s);
  N_VDestroy(c);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  SUNMatDestroy(J);
  free(data);
  SUNContext_Free(&sunctx);

  return(0);
}


static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset)
{
  int retval;

  printf("\n");

  if (mset==1)
    printf("Exact Newton");
  else
    printf("Modified Newton");

  if (glstr == KIN_NONE)
    printf("\n");
  else
    printf(" with line search\n");

  retval = KINSetMaxSetupCalls(kmem, mset);
  if (check_retval(&retval, "KINSetMaxSetupCalls", 1)) return(1);

  retval = KINSol(kmem, u, glstr, s, s);
  if (check_retval(&retval, "KINSol", 1)) return(1);

  printf("Solution:\n  [x1,x2] = ");
  PrintOutput(u);

  PrintFinalStats(kmem);

  return(0);

}


/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

/*
 * System function for predator-prey system
 */

static int func(N_Vector u, N_Vector f, void *user_data)
{
  realtype *udata, *fdata;
  realtype x1, l1, L1, x2, l2, L2;
  realtype *lb, *ub;
  UserData data;

  data = (UserData)user_data;
  lb = data->lb;
  ub = data->ub;

  udata = N_VGetArrayPointer(u);
  fdata = N_VGetArrayPointer(f);

  x1 = udata[0];
  x2 = udata[1];
  l1 = udata[2];
  L1 = udata[3];
  l2 = udata[4];
  L2 = udata[5];

  fdata[0] = PT5 * sin(x1*x2) - PT25 * x2 / PI - PT5 * x1;
  fdata[1] = (ONE - PT25/PI)*(exp(TWO*x1)-E) + E*x2/PI - TWO*E*x1;
  fdata[2] = l1 - x1 + lb[0];
  fdata[3] = L1 - x1 + ub[0];
  fdata[4] = l2 - x2 + lb[1];
  fdata[5] = L2 - x2 + ub[1];

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Initial guesses
 */

static void SetInitialGuess1(N_Vector u, UserData data)
{
  realtype x1, x2;
  realtype *udata;
  realtype *lb, *ub;

  udata = N_VGetArrayPointer(u);

  lb = data->lb;
  ub = data->ub;

  /* There are two known solutions for this problem */

  /* this init. guess should take us to (0.29945; 2.83693) */
  x1 = lb[0];
  x2 = lb[1];

  udata[0] = x1;
  udata[1] = x2;
  udata[2] = x1 - lb[0];
  udata[3] = x1 - ub[0];
  udata[4] = x2 - lb[1];
  udata[5] = x2 - ub[1];
}

static void SetInitialGuess2(N_Vector u, UserData data)
{
  realtype x1, x2;
  realtype *udata;
  realtype *lb, *ub;

  udata = N_VGetArrayPointer(u);

  lb = data->lb;
  ub = data->ub;

  /* There are two known solutions for this problem */

  /* this init. guess should take us to (0.5; 3.1415926) */
  x1 = PT5 * (lb[0] + ub[0]);
  x2 = PT5 * (lb[1] + ub[1]);

  udata[0] = x1;
  udata[1] = x2;
  udata[2] = x1 - lb[0];
  udata[3] = x1 - ub[0];
  udata[4] = x2 - lb[1];
  udata[5] = x2 - ub[1];
}

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype fnormtol, realtype scsteptol)
{
  printf("\nFerraris and Tronconi test problem\n");
  printf("Tolerance parameters:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  fnormtol  = %10.6Lg\n  scsteptol = %10.6Lg\n",
         fnormtol, scsteptol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("  fnormtol  = %10.6g\n  scsteptol = %10.6g\n",
         fnormtol, scsteptol);
#else
  printf("  fnormtol  = %10.6g\n  scsteptol = %10.6g\n",
         fnormtol, scsteptol);
#endif
}

/*
 * Print solution
 */

static void PrintOutput(N_Vector u)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %8.6Lg  %8.6Lg\n", Ith(u,1), Ith(u,2));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %8.6g  %8.6g\n", Ith(u,1), Ith(u,2));
#else
    printf(" %8.6g  %8.6g\n", Ith(u,1), Ith(u,2));
#endif
}

/*
 * Print final statistics contained in iopt
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nje, nfeD;
  int retval;

  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);
  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  retval = KINGetNumJacEvals(kmem, &nje);
  check_retval(&retval, "KINGetNumJacEvals", 1);
  retval = KINGetNumLinFuncEvals(kmem, &nfeD);
  check_retval(&retval, "KINGetNumLinFuncEvals", 1);

  printf("Final Statistics:\n");
  printf("  nni = %5ld    nfe  = %5ld \n", nni, nfe);
  printf("  nje = %5ld    nfeD = %5ld \n", nje, nfeD);
}

/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a retval so check if
 *             retval >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_retval(void *retvalvalue, const char *funcname, int opt)
{
  int *errretval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && retvalvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if retval < 0 */
  else if (opt == 1) {
    errretval = (int *) retvalvalue;
    if (*errretval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *errretval);
      return(1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && retvalvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}
