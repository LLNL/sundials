/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2005-05-16 18:12:40 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 * The nonlinear system is solved by KINSOL using the method
 * specified in local variable globalstrat.
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kinsol.h"
#include "kindense.h"
#include "sundialstypes.h"
#include "nvector_serial.h"
#include "sundialsmath.h"

/* Problem Constants */

#define NVAR   2
#define NEQ    3*NVAR

#define FTOL   RCONST(1.e-5)  /* function tolerance */
#define STOL   RCONST(1.e-5) /* step tolerance */

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
static void func(N_Vector u, N_Vector f, void *f_data);

/* Private Helper Functions */
static void SetInitialGuess1(N_Vector u, UserData data);
static void SetInitialGuess2(N_Vector u, UserData data);
static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset);
static void PrintHeader(int globalstrategy, realtype fnormtol, realtype scsteptol);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void *kmem);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  UserData data;
  realtype fnormtol, scsteptol;
  N_Vector u1, u2, u, s, c;
  int glstr, mset, flag;
  void *kmem;

  u1 = u2 = u = NULL;
  s = c = NULL;
  kmem = NULL;
  data = NULL;
  glstr = KIN_NONE;

  /* User data */

  data = (UserData)malloc(sizeof *data);
  data->lb[0] = PT25;       data->ub[0] = ONE;
  data->lb[1] = ONEPT5;     data->ub[1] = TWO*PI;

  /* Create serial vectors of length NEQ */
  u1 = N_VNew_Serial(NEQ);
  if (check_flag((void *)u1, "N_VNew_Serial", 0)) return(1);

  u2 = N_VNew_Serial(NEQ);
  if (check_flag((void *)u2, "N_VNew_Serial", 0)) return(1);

  u = N_VNew_Serial(NEQ);
  if (check_flag((void *)u, "N_VNew_Serial", 0)) return(1);

  s = N_VNew_Serial(NEQ);
  if (check_flag((void *)s, "N_VNew_Serial", 0)) return(1);

  c = N_VNew_Serial(NEQ);
  if (check_flag((void *)c, "N_VNew_Serial", 0)) return(1);

  SetInitialGuess1(u1,data);
  SetInitialGuess2(u2,data);

  N_VConst_Serial(ONE,s); /* no scaling */

  Ith(c,1) =  ZERO;   /* no constraint on x1 */
  Ith(c,2) =  ZERO;   /* no constraint on x2 */
  Ith(c,3) =  ONE;    /* l1 = x1 - x1_min >= 0 */
  Ith(c,4) = -ONE;    /* L1 = x1 - x1_max <= 0 */
  Ith(c,5) =  ONE;    /* l2 = x2 - x2_min >= 0 */
  Ith(c,6) = -ONE;    /* L2 = x2 - x22_min >= 0 */
  
  fnormtol=FTOL; scsteptol=STOL;


  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  flag = KINSetFdata(kmem, data);
  if (check_flag(&flag, "KINSetFdata", 1)) return(1);
  flag = KINSetConstraints(kmem, c);
  if (check_flag(&flag, "KINSetConstraints", 1)) return(1);
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
  flag = KINSetScaledStepTol(kmem, scsteptol);
  if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);

  flag = KINMalloc(kmem, func, u);
  if (check_flag(&flag, "KINMalloc", 1)) return(1);

  /* Call KINDense to specify the linear solver */

  flag = KINDense(kmem, NEQ);
  if (check_flag(&flag, "KINDense", 1)) return(1);

  /* Print out the problem size, solution parameters, initial guess. */
  PrintHeader(glstr, fnormtol, scsteptol);

  /* --------------------------- */

  N_VScale_Serial(ONE,u1,u);
  glstr = KIN_NONE;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale_Serial(ONE,u1,u);
  glstr = KIN_LINESEARCH;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale_Serial(ONE,u1,u);
  glstr = KIN_NONE;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale_Serial(ONE,u1,u);
  glstr = KIN_LINESEARCH;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);



  /* --------------------------- */

  N_VScale_Serial(ONE,u2,u);
  glstr = KIN_NONE;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale_Serial(ONE,u2,u);
  glstr = KIN_LINESEARCH;
  mset = 1;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale_Serial(ONE,u2,u);
  glstr = KIN_NONE;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);

  /* --------------------------- */

  N_VScale_Serial(ONE,u2,u);
  glstr = KIN_LINESEARCH;
  mset = 0;
  SolveIt(kmem, u, s, glstr, mset);




  /* Free memory */

  N_VDestroy_Serial(u);
  N_VDestroy_Serial(s);
  N_VDestroy_Serial(c);
  KINFree(kmem);
  free(data);

  return(0);
}


static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset)
{
  int flag;

  printf("--------------------------------\n");

  printf("glstr = %d   mset = %d\n",glstr,mset);

  printf("Initial guess:");
  PrintOutput(u);

  flag = KINSetMaxSetupCalls(kmem, mset);
  if (check_flag(&flag, "KINSetMaxSetupCalls", 1)) return(1);

  flag = KINSol(kmem, u, glstr, s, s);
  if (check_flag(&flag, "KINSol", 1)) return(1);

  printf("Solution:     ");
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

static void func(N_Vector u, N_Vector f, void *f_data)
{
  realtype *udata, *fdata;
  realtype x1, l1, L1, x2, l2, L2;
  realtype *lb, *ub;
  UserData data;
  
  data = (UserData)f_data;
  lb = data->lb;
  ub = data->ub;

  udata = NV_DATA_S(u);
  fdata = NV_DATA_S(f);

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

  udata = NV_DATA_S(u);

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

  udata = NV_DATA_S(u);

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

static void PrintHeader(int globalstrategy, realtype fnormtol, realtype scsteptol)
{
  printf("\nFerraris and Tronconi test problem --  KINSol (serial version)\n");
  printf("Total system size = %d\n", NEQ);
  printf("Linear solver is KINDENSE\n");
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf("Tolerance parameters:  fnormtol = %Lg   scsteptol = %Lg\n",
         fnormtol, scsteptol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  fnormtol = %lg   scsteptol = %lg\n",
         fnormtol, scsteptol);
#else
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n",
         fnormtol, scsteptol);
#endif

}

/* 
 * Print solution
 */

static void PrintOutput(N_Vector u)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg  %Lg\n", Ith(u,1), Ith(u,2));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %lg  %lg\n", Ith(u,1), Ith(u,2));
#else
    printf(" %g  %g\n", Ith(u,1), Ith(u,2));
#endif
}

/* 
 * Print final statistics contained in iopt 
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nje, nfeD;
  int flag;
  
  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);

  flag = KINDenseGetNumJacEvals(kmem, &nje);
  check_flag(&flag, "KINDenseGetNumJacEvals", 1);
  flag = KINDenseGetNumFuncEvals(kmem, &nfeD);
  check_flag(&flag, "KINDenseGetNumFuncEvals", 1);

  printf("Final Statistics.. \n");
  printf("nni    = %5ld    nfe   = %5ld \n", nni, nfe);
  printf("nje    = %5ld    nfeD  = %5ld \n", nje, nfeD);
}

/*
 * Check function return value...
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
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); 
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}
