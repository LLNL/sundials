/*
 * -----------------------------------------------------------------
 * $Revision: 4552 $
 * $Date: 2015-07-29 10:12:20 -0700 (Wed, 29 Jul 2015) $
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
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_klu.h>
#include <kinsol/kinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

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
  int nnz;
} *UserData;

/* Functions Called by the KINSOL Solver */
static int func(N_Vector u, N_Vector f, void *user_data);
static int jac(N_Vector y, N_Vector f, SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2);
static int jacDense(long int N,
                    N_Vector y, N_Vector f,
                    DlsMat J, void *user_data,
                    N_Vector tmp1, N_Vector tmp2);

/* Private Helper Functions */
static void SetInitialGuess1(N_Vector u, UserData data);
static void SetInitialGuess2(N_Vector u, UserData data);
static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset);
static void PrintHeader(realtype fnormtol, realtype scsteptol);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void *kmem);
static int check_flag(void *flagvalue, const char *funcname, int opt);

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

  /* User data */

  data = (UserData)malloc(sizeof *data);
  data->lb[0] = PT25;       data->ub[0] = ONE;
  data->lb[1] = ONEPT5;     data->ub[1] = TWO*PI;
  data->nnz = 12;

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

  NV_Ith_S(c,0) =  ZERO;   /* no constraint on x1 */
  NV_Ith_S(c,1) =  ZERO;   /* no constraint on x2 */
  NV_Ith_S(c,2) =  ONE;    /* l1 = x1 - x1_min >= 0  */
  NV_Ith_S(c,3) = -ONE;    /* L1 = x1 - x1_max <= 0  */
  NV_Ith_S(c,4) =  ONE;    /* l2 = x2 - x2_min >= 0  */
  NV_Ith_S(c,5) = -ONE;    /* L2 = x2 - x22_min <= 0 */
  
  fnormtol=FTOL; scsteptol=STOL;


  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  flag = KINSetUserData(kmem, data);
  if (check_flag(&flag, "KINSetUserData", 1)) return(1);
  flag = KINSetConstraints(kmem, c);
  if (check_flag(&flag, "KINSetConstraints", 1)) return(1);
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
  flag = KINSetScaledStepTol(kmem, scsteptol);
  if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);

  flag = KINInit(kmem, func, u);
  if (check_flag(&flag, "KINInit", 1)) return(1);

  /* Call KINDense to specify the linear solver */
/*   
 * flag = KINDense(kmem, NEQ); // <-- Replace with KINKLU
 * if (check_flag(&flag, "KINDense", 1)) return(1);
 * flag = KINDlsSetDenseJacFn(kmem, jacDense);
 * if (check_flag(&flag, "KINDlsSetDenseJacFn", 1)) return(1);
 */

  /* Attach KLU linear solver, which uses CSR matrix */
  flag = KINKLU(kmem, NEQ, data->nnz, CSR_MAT);
  if (check_flag(&flag, "KINKLU", 1)) return(1);

  flag = KINSlsSetSparseJacFn(kmem, jac);
  if (check_flag(&flag, "KINSlsSetSparseJacFn", 1)) return(1);

  /* Print out the problem size, solution parameters, initial guess. */
  PrintHeader(fnormtol, scsteptol);

  /* --------------------------- */

  printf("\n------------------------------------------\n");
  printf("\nInitial guess on lower bounds\n");
  printf("  [x1,x2] = ");
  PrintOutput(u1);

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

  printf("\n------------------------------------------\n");
  printf("\nInitial guess in middle of feasible region\n");
  printf("  [x1,x2] = ");
  PrintOutput(u2);

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

  N_VDestroy_Serial(u1);
  N_VDestroy_Serial(u2);
  N_VDestroy_Serial(u);
  N_VDestroy_Serial(s);
  N_VDestroy_Serial(c);
  KINFree(&kmem);
  free(data);

  return(0);
}


static int SolveIt(void *kmem, N_Vector u, N_Vector s, int glstr, int mset)
{
  int flag;

  printf("\n");

  if (mset==1)
    printf("Exact Newton");
  else
    printf("Modified Newton");

  if (glstr == KIN_NONE)
    printf("\n");
  else
    printf(" with line search\n");

  flag = KINSetMaxSetupCalls(kmem, mset);
  if (check_flag(&flag, "KINSetMaxSetupCalls", 1)) return(1);

  flag = KINSol(kmem, u, glstr, s, s);
  if (check_flag(&flag, "KINSol", 1)) return(1);

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
  UserData params;
  
  params = (UserData)user_data;
  lb = params->lb;
  ub = params->ub;

  udata = N_VGetArrayPointer_Serial(u);
  fdata = N_VGetArrayPointer_Serial(f);

  x1 = udata[0];
  x2 = udata[1];
  l1 = udata[2];
  L1 = udata[3];
  l2 = udata[4];
  L2 = udata[5];

  fdata[0] = PT5 * sin(x1*x2) - PT25 * x2 / PI - PT5 * x1;
  fdata[1] = (ONE - PT25/PI)*(SUNRexp(TWO*x1)-E) + E*x2/PI - TWO*E*x1;
  fdata[2] = l1 - x1 + lb[0];
  fdata[3] = L1 - x1 + ub[0];
  fdata[4] = l2 - x2 + lb[1];
  fdata[5] = L2 - x2 + ub[1];

  return(0);
}


/*
 * System Jacobian
 */

static int jac(N_Vector y, N_Vector f,
               SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2)
{
  realtype *yd;
  int *rowptrs;
  int *colvals;
  realtype *data;
  
  yd = N_VGetArrayPointer_Serial(y);
  rowptrs = (*J->rowptrs);
  colvals = (*J->colvals);
  data    = J->data;
  
  SparseSetMatToZero(J);

  rowptrs[0] =  0;
  rowptrs[1] =  2;
  rowptrs[2] =  4;
  rowptrs[3] =  6;
  rowptrs[4] =  8;
  rowptrs[5] = 10;
  rowptrs[6] = 12;

  /* row 0 */
  data[0] = PT5 * cos(yd[0]*yd[1]) * yd[1] - PT5;
  colvals[0] = 0;
  data[1] = PT5 * cos(yd[0]*yd[1]) * yd[0] - PT25/PI;
  colvals[1] = 1;

  /* row 1 */
  data[2] = TWO * (ONE - PT25/PI) * (SUNRexp(TWO*yd[0]) - E);
  colvals[2] = 0;
  data[3] = E/PI;
  colvals[3] = 1;
  
  /* row 2 */
  data[4] = -ONE;
  colvals[4] = 0;
  data[5] =  ONE;
  colvals[5] = 2;
  
  /* row 3 */
  data[6] = -ONE;
  colvals[6] = 0;
  data[7] =  ONE;
  colvals[7] = 3;
  
  /* row 4 */
  data[8] = -ONE;
  colvals[8] = 1;
  data[9] =  ONE;
  colvals[9] = 4;
  
  /* row 5 */
  data[10] = -ONE;
  colvals[10] = 1;
  data[11] =  ONE;
  colvals[11] = 5;
  
  return(0);

}


/*
 * System dense Jacobian
 */

static int jacDense(long int N,
                    N_Vector y, N_Vector f,
                    DlsMat J, void *user_data,
                    N_Vector tmp1, N_Vector tmp2)
{
  realtype *yd;

  yd = N_VGetArrayPointer_Serial(y);

  /* row 0 */
  DENSE_ELEM(J,0,0) = PT5 * cos(yd[0]*yd[1]) * yd[1] - PT5;
  DENSE_ELEM(J,0,1) = PT5 * cos(yd[0]*yd[1]) * yd[0] - PT25/PI;

  /* row 1 */
  DENSE_ELEM(J,1,0) = TWO * (ONE - PT25/PI) * (SUNRexp(TWO*yd[0]) - E);
  DENSE_ELEM(J,1,1) = E/PI;
  
  /* row 2 */
  DENSE_ELEM(J,2,0) = -ONE;
  DENSE_ELEM(J,2,2) =  ONE;
  
  /* row 3 */
  DENSE_ELEM(J,3,0) = -ONE;
  DENSE_ELEM(J,3,3) =  ONE;
  
  /* row 4 */
  DENSE_ELEM(J,4,1) = -ONE;
  DENSE_ELEM(J,4,4) =  ONE;
  
  /* row 5 */
  DENSE_ELEM(J,5,1) = -ONE;
  DENSE_ELEM(J,5,5) =  ONE;
  

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

  udata = N_VGetArrayPointer_Serial(u);

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

  udata = N_VGetArrayPointer_Serial(u);

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
    printf(" %8.6Lg  %8.6Lg\n", NV_Ith_S(u,0), NV_Ith_S(u,1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %8.6g  %8.6g\n", NV_Ith_S(u,0), NV_Ith_S(u,1));
#else
    printf(" %8.6g  %8.6g\n", NV_Ith_S(u,0), NV_Ith_S(u,1));
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

/*   flag = KINDlsGetNumJacEvals(kmem, &nje);
 *   check_flag(&flag, "KINDlsGetNumJacEvals", 1);
 */

  flag = KINSlsGetNumJacEvals(kmem, &nje);
  check_flag(&flag, "KINSlsGetNumJacEvals", 1);
  
  printf("Final Statistics:\n");
  printf("  nni = %5ld    nfe  = %5ld \n", nni, nfe);
  printf("  nje = %5ld    \n", nje);
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

static int check_flag(void *flagvalue, const char *funcname, int opt)
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
