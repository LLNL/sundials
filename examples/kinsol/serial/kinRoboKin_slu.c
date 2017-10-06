/*
 * -----------------------------------------------------------------
 * $Revision: 4834 $
 * $Date: 2016-08-01 16:59:05 -0700 (Mon, 01 Aug 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol S. Woodward @ LLNL.  Adapted from the file
 *    kinRoboKin_dns.c by Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * This example solves a nonlinear system from robot kinematics.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 6 from Section 14.1, Chapter 14
 * 
 * The nonlinear system is solved by KINSOL using the SUPERLU_MT linear
 * solver.
 *
 * Constraints are imposed to make all components of the solution
 * be within [-1,1].
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_superlumt.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

/* Problem Constants */

#define NVAR  8              /* variables */ 
#define NEQ   3*NVAR         /* equations + bounds */

#define FTOL   RCONST(1.e-5) /* function tolerance */
#define STOL   RCONST(1.e-5) /* step tolerance */

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)

static int func(N_Vector y, N_Vector f, void *user_data);
static int jac(N_Vector y, N_Vector f,
               SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2);
static void PrintOutput(N_Vector y);
static void PrintFinalStats(void *kmem);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  realtype fnormtol, scsteptol;
  N_Vector y, scale, constraints;
  int mset, flag, i;
  void *kmem;

  int nnz, num_threads;

  y = scale = constraints = NULL;
  kmem = NULL;

  printf("\nRobot Kinematics Example\n");
  printf("8 variables; -1 <= x_i <= 1\n");
  printf("KINSOL problem size: 8 + 2*8 = 24 \n\n");

  /* Create vectors for solution, scales, and constraints */

  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  scale = N_VNew_Serial(NEQ);
  if (check_flag((void *)scale, "N_VNew_Serial", 0)) return(1);

  constraints = N_VNew_Serial(NEQ);
  if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);

  /* Initialize and allocate memory for KINSOL */

  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  flag = KINInit(kmem, func, y); /* y passed as a template */
  if (check_flag(&flag, "KINInit", 1)) return(1);

  /* Set optional inputs */

  N_VConst_Serial(ZERO,constraints);
  for (i = NVAR+1; i <= NEQ; i++) Ith(constraints, i) = ONE;
  
  flag = KINSetConstraints(kmem, constraints);
  if (check_flag(&flag, "KINSetConstraints", 1)) return(1);

  fnormtol  = FTOL; 
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);

  scsteptol = STOL;
  flag = KINSetScaledStepTol(kmem, scsteptol);
  if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);

  /* Attach SuperLU_MT linear solver */
  /* Set the number of nonzeros in the Jacobian */
  nnz = 56;
  /* Set the number fo threads to use */
  num_threads = 2;
  flag = KINSuperLUMT(kmem, num_threads, NEQ, nnz);
  if (check_flag(&flag, "KINSuperLUMT", 1)) return(1);

  flag = KINSlsSetSparseJacFn(kmem, jac);
  if (check_flag(&flag, "KINSlsSetSparseJacFn", 1)) return(1);

  /* Indicate exact Newton */

  mset = 1;
  flag = KINSetMaxSetupCalls(kmem, mset);
  if (check_flag(&flag, "KINSetMaxSetupCalls", 1)) return(1);

  /* Initial guess */

  N_VConst_Serial(ONE, y);
  for(i = 1; i <= NVAR; i++) Ith(y,i) = SUNRsqrt(TWO)/TWO;

  printf("Initial guess:\n");
  PrintOutput(y);

  /* Call KINSol to solve problem */

  N_VConst_Serial(ONE,scale);
  flag = KINSol(kmem,           /* KINSol memory block */
                y,              /* initial guess on input; solution vector */
                KIN_LINESEARCH, /* global strategy choice */
                scale,          /* scaling vector, for the variable cc */
                scale);         /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1)) return(1);

  printf("\nComputed solution:\n");
  PrintOutput(y);

  /* Print final statistics and free memory */  

  PrintFinalStats(kmem);

  N_VDestroy_Serial(y);
  N_VDestroy_Serial(scale);
  N_VDestroy_Serial(constraints);
  KINFree(&kmem);

  return(0);
}

/* 
 * System function 
 */

static int func(N_Vector y, N_Vector f, void *user_data)
{
  realtype *yd, *fd;

  realtype x1, x2, x3, x4, x5, x6, x7, x8;
  realtype l1, l2, l3, l4, l5, l6, l7, l8;
  realtype u1, u2, u3, u4, u5, u6, u7, u8;

  realtype eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8;
  realtype lb1, lb2, lb3, lb4, lb5, lb6, lb7, lb8;
  realtype ub1, ub2, ub3, ub4, ub5, ub6, ub7, ub8;

  yd = N_VGetArrayPointer_Serial(y);
  fd = N_VGetArrayPointer_Serial(f);

  x1 = yd[0]; l1 = yd[ 8]; u1 = yd[16]; 
  x2 = yd[1]; l2 = yd[ 9]; u2 = yd[17]; 
  x3 = yd[2]; l3 = yd[10]; u3 = yd[18]; 
  x4 = yd[3]; l4 = yd[11]; u4 = yd[19]; 
  x5 = yd[4]; l5 = yd[12]; u5 = yd[20]; 
  x6 = yd[5]; l6 = yd[13]; u6 = yd[21]; 
  x7 = yd[6]; l7 = yd[14]; u7 = yd[22]; 
  x8 = yd[7]; l8 = yd[15]; u8 = yd[23]; 

  /* Nonlinear equations */

  eq1 = - 0.1238*x1 + x7 - 0.001637*x2 
    - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571;
  eq2 = 0.2638*x1 - x7 - 0.07745*x2 
    - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022;
  eq3 = 0.3578*x1 + 0.004731*x2 + x6*x8;
  eq4 = - 0.7623*x1 + 0.2238*x2 + 0.3461;
  eq5 = x1*x1 + x2*x2 - 1;
  eq6 = x3*x3 + x4*x4 - 1;
  eq7 = x5*x5 + x6*x6 - 1;
  eq8 = x7*x7 + x8*x8 - 1;

  /* Lower bounds ( l_i = 1 + x_i >= 0)*/

  lb1 = l1 - 1.0 - x1;
  lb2 = l2 - 1.0 - x2;
  lb3 = l3 - 1.0 - x3;
  lb4 = l4 - 1.0 - x4;
  lb5 = l5 - 1.0 - x5;
  lb6 = l6 - 1.0 - x6;
  lb7 = l7 - 1.0 - x7;
  lb8 = l8 - 1.0 - x8;

  /* Upper bounds ( u_i = 1 - x_i >= 0)*/

  ub1 = u1 - 1.0 + x1;
  ub2 = u2 - 1.0 + x2;
  ub3 = u3 - 1.0 + x3;
  ub4 = u4 - 1.0 + x4;
  ub5 = u5 - 1.0 + x5;
  ub6 = u6 - 1.0 + x6;
  ub7 = u7 - 1.0 + x7;
  ub8 = u8 - 1.0 + x8;

  fd[0] = eq1; fd[ 8] = lb1; fd[16] = ub1;
  fd[1] = eq2; fd[ 9] = lb2; fd[17] = ub2;
  fd[2] = eq3; fd[10] = lb3; fd[18] = ub3;
  fd[3] = eq4; fd[11] = lb4; fd[19] = ub4;
  fd[4] = eq5; fd[12] = lb5; fd[20] = ub5;
  fd[5] = eq6; fd[13] = lb6; fd[21] = ub6;
  fd[6] = eq7; fd[14] = lb7; fd[22] = ub7;
  fd[7] = eq8; fd[15] = lb8; fd[23] = ub8;

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
  realtype x1, x2, x3, x4, x5, x6, x7, x8;
  int *colptrs = *J->colptrs;
  int *rowvals = *J->rowvals;
  realtype *data = J->data;

  yd = N_VGetArrayPointer_Serial(y);

  x1 = yd[0];
  x2 = yd[1];
  x3 = yd[2];
  x4 = yd[3];
  x5 = yd[4];
  x6 = yd[5];
  x7 = yd[6];
  x8 = yd[7];

  SparseSetMatToZero(J);
  
  colptrs[0] = 0;
  colptrs[1] = 7;
  colptrs[2] = 14;
  colptrs[3] = 19;
  colptrs[4] = 24;
  colptrs[5] = 27;
  colptrs[6] = 31;
  colptrs[7] = 36;
  colptrs[8] = 40;
  colptrs[9] = 41;
  colptrs[10] = 42;
  colptrs[11] = 43;
  colptrs[12] = 44;
  colptrs[13] = 45;
  colptrs[14] = 46;
  colptrs[15] = 47;
  colptrs[16] = 48;
  colptrs[17] = 49;
  colptrs[18] = 50;
  colptrs[19] = 51;
  colptrs[20] = 52;
  colptrs[21] = 53;
  colptrs[22] = 54;
  colptrs[23] = 55;
  colptrs[24] = 56;

  /* Nonlinear equations */

  /* 
     - 0.1238*x1 + x7 - 0.001637*x2 
     - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571 
  */
  /*
  IJth(J,1,1) = - 0.1238 + 0.004731*x3;
  IJth(J,1,2) = - 0.001637 - 0.3578*x3;
  IJth(J,1,3) = 0.004731*x1 - 0.3578*x2;
  IJth(J,1,4) = - 0.9338;
  IJth(J,1,7) = 1.0;
  */

  J->data[0] = - 0.1238 + 0.004731*x3;
  rowvals[0] = 0;
  J->data[7] = - 0.001637 - 0.3578*x3;
  rowvals[7] = 0;
  J->data[14] = 0.004731*x1 - 0.3578*x2;
  rowvals[14] = 0;
  J->data[19] = - 0.9338;
  rowvals[19] = 0;
  J->data[31] = 1.0;
  rowvals[31] = 0;

  /*
    0.2638*x1 - x7 - 0.07745*x2 
    - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022
  */
  /*
  IJth(J,2,1) = 0.2638 + 0.2238*x3;
  IJth(J,2,2) = - 0.07745 + 0.7623*x3;
  IJth(J,2,3) = 0.2238*x1 + 0.7623*x2;
  IJth(J,2,4) = - 0.6734;
  IJth(J,2,7) = -1.0;
  */

  J->data[1] = 0.2638 + 0.2238*x3;
  rowvals[1] = 1;
  J->data[8] = - 0.07745 + 0.7623*x3;
  rowvals[8] = 1;
  J->data[15] = 0.2238*x1 + 0.7623*x2;
  rowvals[15] = 1;
  J->data[20] = - 0.6734;
  rowvals[20] = 1;
  J->data[32] = -1.0;
  rowvals[32] = 1;


  /*
    0.3578*x1 + 0.004731*x2 + x6*x8
  */
  /*
  IJth(J,3,1) = 0.3578;
  IJth(J,3,2) = 0.004731;
  IJth(J,3,6) = x8;
  IJth(J,3,8) = x6;
  */

  J->data[2] = 0.3578;
  rowvals[2] = 2;
  J->data[9] = 0.004731;
  rowvals[9] = 2;
  J->data[27] = x8;
  rowvals[27] = 2;
  J->data[36] = x6;
  rowvals[36] = 2;


  /*
    - 0.7623*x1 + 0.2238*x2 + 0.3461
  */
  /*
  IJth(J,4,1) = - 0.7623;
  IJth(J,4,2) = 0.2238;
  */

  J->data[3] = - 0.7623;
  rowvals[3] = 3;
  J->data[10] = 0.2238;
  rowvals[10] = 3;

  /*
    x1*x1 + x2*x2 - 1
  */
  /*
  IJth(J,5,1) = 2.0*x1;
  IJth(J,5,2) = 2.0*x2;
  */

  J->data[4] = 2.0*x1;
  rowvals[4] = 4;
  J->data[11] = 2.0*x2;
  rowvals[11] = 4;

  /*
    x3*x3 + x4*x4 - 1
  */
  /*
  IJth(J,6,3) = 2.0*x3;
  IJth(J,6,4) = 2.0*x4;
  */

  J->data[16] = 2.0*x3;
  rowvals[16] = 5;
  J->data[21] = 2.0*x4;
  rowvals[21] = 5;

  /*
    x5*x5 + x6*x6 - 1
  */
  /*
  IJth(J,7,5) = 2.0*x5;
  IJth(J,7,6) = 2.0*x6;
  */

  J->data[24] = 2.0*x5;
  rowvals[24] = 6;
  J->data[28] = 2.0*x6;
  rowvals[28] = 6;

  /*
    x7*x7 + x8*x8 - 1
  */
  /*
  IJth(J,8,7) = 2.0*x7;
  IJth(J,8,8) = 2.0*x8;
  */
  J->data[33] = 2.0*x7;
  rowvals[33] = 7;
  J->data[37] = 2.0*x8;
  rowvals[37] = 7;

  
  /*
    Lower bounds ( l_i = 1 + x_i >= 0)
    l_i - 1.0 - x_i
   */


  /*
  for(i=1;i<=8;i++) {
    IJth(J,8+i,i)   = -1.0;
    IJth(J,8+i,8+i) =  1.0;
  } 
  */

  J->data[5] = -1.0;
  rowvals[5] = 8;
  J->data[12] = -1.0;
  rowvals[12] = 9;
  J->data[17] = -1.0;
  rowvals[17] = 10;
  J->data[22] = -1.0;
  rowvals[22] = 11;
  J->data[25] = -1.0;
  rowvals[25] = 12;
  J->data[29] = -1.0;
  rowvals[29] = 13;
  J->data[34] = -1.0;
  rowvals[34] = 14;
  J->data[38] = -1.0;
  rowvals[38] = 15;

  J->data[40] = 1.0;
  rowvals[40] = 8;
  J->data[41] = 1.0;
  rowvals[41] = 9;
  J->data[42] = 1.0;
  rowvals[42] = 10;
  J->data[43] = 1.0;
  rowvals[43] = 11;
  J->data[44] = 1.0;
  rowvals[44] = 12;
  J->data[45] = 1.0;
  rowvals[45] = 13;
  J->data[46] = 1.0;
  rowvals[46] = 14;
  J->data[47] = 1.0;
  rowvals[47] = 15;


  /*
    Upper bounds ( u_i = 1 - x_i >= 0)
    u_i - 1.0 + x_i
   */
  /*
  for(i=1;i<=8;i++) {
    IJth(J,16+i,i)    = 1.0;
    IJth(J,16+i,16+i) = 1.0;
  }
  */

  J->data[6] = 1.0;
  rowvals[6] = 16;
  J->data[13] = 1.0;
  rowvals[13] = 17;
  J->data[18] = 1.0;
  rowvals[18] = 18;
  J->data[23] = 1.0;
  rowvals[23] = 19;
  J->data[26] = 1.0;
  rowvals[26] = 20;
  J->data[30] = 1.0;
  rowvals[30] = 21;
  J->data[35] = 1.0;
  rowvals[35] = 22;
  J->data[39] = 1.0;
  rowvals[39] = 23;

  J->data[48] = 1.0;
  rowvals[48] = 16;
  J->data[49] = 1.0;
  rowvals[49] = 17;
  J->data[50] = 1.0;
  rowvals[50] = 18;
  J->data[51] = 1.0;
  rowvals[51] = 19;
  J->data[52] = 1.0;
  rowvals[52] = 20;
  J->data[53] = 1.0;
  rowvals[53] = 21;
  J->data[54] = 1.0;
  rowvals[54] = 22;
  J->data[55] = 1.0;
  rowvals[55] = 23;

  return(0);

}

/* 
 * Print solution
 */

static void PrintOutput(N_Vector y)
{
  int i;

  printf("     l=x+1          x         u=1-x\n");
  printf("   ----------------------------------\n");

  for(i=1; i<=NVAR; i++) {

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %10.6Lg   %10.6Lg   %10.6Lg\n", 
           Ith(y,i+NVAR), Ith(y,i), Ith(y,i+2*NVAR));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %10.6g   %10.6g   %10.6g\n", 
           Ith(y,i+NVAR), Ith(y,i), Ith(y,i+2*NVAR));
#else
    printf(" %10.6g   %10.6g   %10.6g\n", 
           Ith(y,i+NVAR), Ith(y,i), Ith(y,i+2*NVAR));
#endif

  }

}

/* 
 * Print final statistics
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nje;
  int flag;
  
  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);

  flag = KINSlsGetNumJacEvals(kmem, &nje);
  check_flag(&flag, "KINSlsGetNumJacEvals", 1);

  printf("\nFinal Statistics.. \n");
  printf("nni    = %5ld    nfe   = %5ld \n", nni, nfe);
  printf("nje    = %5ld \n", nje);
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
