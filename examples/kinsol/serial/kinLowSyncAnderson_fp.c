/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves the linear convection-diffusion problem 
 *
 * laplacian u + 20u + 20u_x + 20u_y = -10
 *
 * on the unit square with u = 0 on the boundaries. 
 *
 * The problem was discretized using centered differences on a 128x128
 * grid. Restricted Additive Schwartz was applied on a 4x4 array of
 * subdomains with three grid lines of overlap and accelerated with
 * Anderson Acceleration, not limiting m.
 * -----------------------------------------------------------------
 * References:
 *
 * 1. Walker H and Ni P,
 *    Anderson Acceleration for Fixed-Point Iterations,
 *    SIAM Journal on Numerical Analysis,
 *    2011, vol.49(4), pp.1715-1735.
 *
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kinsol/kinsol.h"               /* access to KINSOL func., consts. */
#include "nvector/nvector_serial.h"      /* access to serial N_Vector       */
#include "sunmatrix/sunmatrix_sparse.h"  /* access to Sparse SunMatrix      */

/* precision specific formatting macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define ABS(x)  (fabs((x)))
#define SQRT(x) (sqrt((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define ABS(x)  (fabsf((x)))
#define SQRT(x) (sqrtf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define ABS(x)  (fabsl((x)))
#define SQRT(x) (sqrtl((x)))
#endif

/* problem constants */
#define NUMREGIONS   4                       /* number of subdomains for RAS             */
#define OVERLAPWIDTH 3                       /* number of overlapping grid lines for RAS */
#define N            128                     /* dimension of one side of grid            */
//#define N            10                     /* dimension of one side of grid            */

#define ZERO         RCONST(0.0)             /* real 0.0   */
#define C            RCONST(20.0)            /* real 20.0  */
#define D            RCONST(20.0)            /* real 20.0  */
#define ONE          RCONST(1.0)             /* real 1.0   */
#define NEGONE       RCONST(-1.0)            /* real -1.0  */

/* UserData */ 
typedef struct {
  SUNMatrix A;
  SUNMatrix L;
  SUNMatrix R;
  N_Vector v;
  N_Vector accum;
} *UserData;

/* Linear Convection-Diffusion fixed point function */
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

/* Function to return the 2D Convection Diffusion Matrix */
static int ConvectionDiffusionMatrix2D(SUNMatrix Q);

/* Function to return the matrix for overlapping regions in RAS */
static int OverlappingRestrictionMatrix2D(SUNMatrix Q);

/* Function to return the matrix for nonoverlapping regions in RAS*/
static int NonOverlappingRestrictionMatrix2D(SUNMatrix Q);

/* Check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Check the system solution */
static int check_ans(N_Vector u, realtype tol);

/* -----------------------------------------------------------------------------
 * Main program
 * ---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int       retval  = 0;
  N_Vector  u       = NULL;
  N_Vector  scale   = NULL;
  realtype  tol     = 100 * SQRT(UNIT_ROUNDOFF);
  long int  mxiter  = 75;
  long int  maa     = 0;           /* no acceleration          */
  realtype  damping = RCONST(1.0); /* no damping               */
  int       orth    = 0;           /* standard MGS in Anderson */
  long int  nni, nfe;
  realtype* data;
  void*     kmem;
  UserData  udata;                 /* contain matrices for function evals */ 

  /* Check if a acceleration/damping/orthogonaliztion values were provided */
  if (argc > 1) maa     = (long int) atoi(argv[1]);
  if (argc > 2) damping = (realtype) atof(argv[2]);
  if (argc > 3) orth = (int) atof(argv[3]);

  /* -------------------------
   * Print problem description
   * ------------------------- */

  printf("Solve the linear system:\n");
  printf("    laplacian u + 20u + 20u_x + 20u_y = -10\n");
  printf("                                    u = 0\n");
  printf("     on [0,1] x [0,1] with boundaries = 0\n");
  printf("Solution method: Anderson accelerated RAS posed as fixed point iteration.\n");
  printf("    tolerance    = %"GSYM"\n", tol);
  printf("    max iters    = %ld\n", mxiter);
  printf("    accel vec    = %ld\n", maa);
  printf("    damping      = %"GSYM"\n", damping);
  printf("    orth routine = %ld\n", orth);
  
  /* ---------------------------------------------------------------------------
   * Allocate user data -- containing matrices required for function evaluations 
   * --------------------------------------------------------------------------- */
  udata = (UserData)malloc(sizeof *udata);

  /* create 2d convection diffusion matrix here */
  /* allocate space for v and accum here */
  udata->v     = N_VNew_Serial(N*N);
  if (check_retval((void *)udata->v, "N_VNew_Serial", 0)) return(1);

  udata->accum = N_VNew_Serial(N*N);
  if (check_retval((void *)udata->accum, "N_VNew_Serial", 0)) return(1);
  
  /* allocate space for A, and all L and R matrices, used in problems domain decomposition */
  printf("Before convection diffusion matrix created\n");
  ConvectionDiffusionMatrix2D(udata->A);
  printf("After convection diffusion matrix created\n");

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  u = N_VNew_Serial(N*N);
  if (check_retval((void *)u, "N_VNew_Serial", 0)) return(1);

  scale = N_VClone(u);
  if (check_retval((void *)scale, "N_VClone", 0)) return(1);

  /* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * ----------------------------------------- */

  kmem = KINCreate();
  if (check_retval((void *)kmem, "KINCreate", 0)) return(1);

  /* Set number of prior residuals used in Anderson acceleration */
  retval = KINSetMAA(kmem, maa);
  
  /* Set orthogonalization routine used in Anderson acceleration */
  retval = KINSetOrthAA(kmem, orth);

  /* Initialize FixedPoint Iteration as the function we're using */
  retval = KINInit(kmem, FPFunction, u);
  if (check_retval(&retval, "KINInit", 1)) return(1);

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Specify stopping tolerance based on residual */
  retval = KINSetFuncNormTol(kmem, tol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);

  /* Set maximum number of iterations */
  retval = KINSetNumMaxIters(kmem, mxiter);
  if (check_retval(&retval, "KINSetNumMaxItersFuncNormTol", 1)) return(1);

  /* Set Anderson acceleration damping parameter */
  retval = KINSetDampingAA(kmem, damping);
  if (check_retval(&retval, "KINSetDampingAA", 1)) return(1);

  /* -------------
   * Initial guess
   * ------------- */

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return(1);
  N_VConst(ZERO, u);

  /* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- */

  /* No scaling used */
  N_VConst(ONE, scale);

  /* Call main solver */
  retval = KINSol(kmem,         /* KINSol memory block */
                  u,            /* initial guess on input; solution vector */
                  KIN_FP,       /* global strategy choice */
                  scale,        /* scaling vector, for the variable cc */
                  scale);       /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) return(1);

  /* ------------------------------------
   * Get solver statistics
   * ------------------------------------ */

  /* get solver stats */
  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);

  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("Number of solver iterations: %6ld\n", nni);
  printf("Number of function evaluations: %6ld\n", nfe);

  /* ------------------------------------
   * Print solution and check error
   * ------------------------------------ */

  /* check solution */
  retval = check_ans(u, tol);

  /* -----------
   * Free memory
   * ----------- */

  N_VDestroy(u);
  N_VDestroy(scale);
  KINFree(&kmem);
  free(udata);

  return(retval);
}
    

/* -----------------------------------------------------------------------------
 * Linear convection-diffusion problem 
 *
 * laplacian u + 20u + 20u_x + 20u_y = -10
 * on the unit square with u = 0 on the boundaries 
 *
 * Linear fixed point function 
 *
 * ** PUT FUNCTION DESCRIPTION HERE ** 
 *
 * ---------------------------------------------------------------------------*/
int FPFunction(N_Vector u, N_Vector g, void* user_data)
{
  SUNMatrix subA, A, R, L;
  N_Vector  v, accum, temp; 
  UserData udata;
  int i, j;
  
  /* Need to setup size of subA and temp */
 
  /* Grab matrices and vectors from user_data */ 
  /*udata = (UserData)user_data;
  A = udata->A;
  R = udata->R;
  L = udata->L;
  v = udata->v;
  accum = udata->accum;*/

  /* v = A * u */ 
  /*SUNMatMatvec(A, u, v);*/

  /* v = f - v */

  /* accum = 0 */
  /*N_VConst(ZERO, accum);*/

  for (i = 0; i < NUMREGIONS; i++)
  {
    for (j = 0; j < NUMREGIONS; j++)
    {
      /* setup R = OverlappingRestrictionMatrix2D */
      /* setup L = NonOverlappingRestrictionMatrix2D */
      /* subA = R * A * R^T */
      /* temp = R * v */
      /* solve subA * temp = temp (overwrite old temp) */
      /* temp = L^T * temp */
      /* accum = accum + temp */
      continue; /* just for now */ 
    }  
  }
  /* v = u + accum */
  /* copy v into u and return */

  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to return the 2D Convection Diffusion Matrix to be used in
 * the fixed point function call
 * ---------------------------------------------------------------------------*/
int ConvectionDiffusionMatrix2D(SUNMatrix Q)
{
  SUNMatrix I, Laplacian, CenteredDiff;
  realtype  h, h2, h2_inv;
  realtype  *colj, *matdata;
  sunindextype *rowptrs, *colindices;
  int       i, j, nnz, nnz_ctr, fail;

  /* Scaled Identity matrix */
  nnz = N*N;
  I = SUNSparseMatrix(N*N, N*N, nnz, CSR_MAT);
  rowptrs = SUNSparseMatrix_IndexPointers(I);
  rowptrs[0] = 0; /* first row starts at index 0 */
  colindices = SUNSparseMatrix_IndexValues(I);
  matdata = SUNSparseMatrix_Data(I);
  for (j = 0; j < N*N; j++)
  {    
    matdata[i] = C;
    colindices[i] = i;
    rowptrs[i+1] = rowptrs[i] + 1;
  }

  /* Laplacian Matrix */
  h = ONE / (N + 1);
  h2_inv = ONE / (h*h);
  
  nnz = 81408;
  //nnz = 460;
  Laplacian = SUNSparseMatrix(N*N, N*N, nnz, CSR_MAT);
  rowptrs = SUNSparseMatrix_IndexPointers(Laplacian);
  rowptrs[0] = 0; /* first row starts at index 0 */
  colindices = SUNSparseMatrix_IndexValues(Laplacian);
  matdata = SUNSparseMatrix_Data(Laplacian);
  for (i = 0; i < N*N; i++)
  {
    nnz_ctr = rowptrs[i];
    if (i == 0)
    {
      /* FIRST ROW OF MATRIX */
      /* vals = h2_inv * [-4, 1, 1] */
      /* cols = [i, i+1, i+N] */
      matdata[nnz_ctr] = -4; matdata[nnz_ctr+1] = 1; matdata[nnz_ctr+2] = 1;
      colindices[nnz_ctr] = i; colindices[nnz_ctr+1] = i+1; colindices[nnz_ctr+2] = i+N;
      rowptrs[i+1] = rowptrs[i] + 3;
    }
    else if (i == N*N-1)
    {
      /* LAST ROW OF MATRIX */
      /* vals = h2_inv * [1, 1, -4] */
      /* cols = [i-N, i-1, i] */
      matdata[nnz_ctr] = 1; matdata[nnz_ctr+1] = 1; matdata[nnz_ctr+2] = -4;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i;
      rowptrs[i+1] = rowptrs[i] + 3;
    }
    else if ((i+1) % N == 0)
    {
      /* vals = h2_inv * [1, 1, -4, 1] */
      /* cols = [i-N, i-1, i, i+N] */
      matdata[nnz_ctr] = 1; matdata[nnz_ctr+1] = 1; matdata[nnz_ctr+2] = -4; matdata[nnz_ctr+3] = 1;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i; colindices[nnz_ctr+3] = i+N;
      rowptrs[i+1] = rowptrs[i] + 4;
    }
    else if (i % N == 0)
    {
      /* vals = h2_inv * [1, -4, 1, 1] */
      /* cols = [i-N, i, i+1, i+N] */
      matdata[nnz_ctr] = 1; matdata[nnz_ctr+1] = -4; matdata[nnz_ctr+2] = 1; matdata[nnz_ctr+3] = 1;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i; colindices[nnz_ctr+2] = i+1; colindices[nnz_ctr+3] = i+N;
      rowptrs[i+1] = rowptrs[i] + 4;
    }
    else if (i < N)
    {
      /* vals = h2_inv * [1, -4, 1, 1] */
      /* cols = [i-1, i, i+1, i+N] */
      matdata[nnz_ctr] = 1; matdata[nnz_ctr+1] = -4; matdata[nnz_ctr+2] = 1; matdata[nnz_ctr+3] = 1;
      colindices[nnz_ctr] = i-1; colindices[nnz_ctr+1] = i; colindices[nnz_ctr+2] = i+1; colindices[nnz_ctr+3] = i+N;
      rowptrs[i+1] = rowptrs[i] + 4;
    }
    else if ((i+N) > (N*N-1)) 
    {
      /* vals = h2_inv * [1, 1, -4, 1] */
      /* cols = [i-N, i-1, i, i+1] */
      matdata[nnz_ctr] = 1; matdata[nnz_ctr+1] = 1; matdata[nnz_ctr+2] = -4; matdata[nnz_ctr+3] = 1;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i; colindices[nnz_ctr+3] = i+1;
      rowptrs[i+1] = rowptrs[i] + 4;
    }
    else 
    {
      /* vals = h2_inv * [1, 1, -4, 1, 1] */
      /* cols = [i-N, i-1, i, i+1, i+N] */
      matdata[nnz_ctr] = 1; matdata[nnz_ctr+1] = 1; matdata[nnz_ctr+2] = -4; matdata[nnz_ctr+3] = 1;
        matdata[nnz_ctr+4] = 1;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i; colindices[nnz_ctr+3] = i+1;
        colindices[nnz_ctr+4] = i + N;
      rowptrs[i+1] = rowptrs[i] + 5;
    }
  }
  
  /* Centered Difference Matrix */
  h2 = ONE / (2*h);

  nnz = 65024;
  //nnz = 360;
  CenteredDiff = SUNSparseMatrix(N*N, N*N, nnz, CSR_MAT);
  rowptrs = SUNSparseMatrix_IndexPointers(CenteredDiff);
  rowptrs[0] = 0; /* first row starts at index 0 */
  colindices = SUNSparseMatrix_IndexValues(CenteredDiff);
  matdata = SUNSparseMatrix_Data(CenteredDiff);
  for (i = 0; i < N*N; i++)
  {
    nnz_ctr = rowptrs[i];
    if (i == 0)
    {
        /* FIRST ROW OF MATRIX */
        /* vals = h2 * [-1, -1] */
        /* cols = [i+1, i+N] */
        matdata[nnz_ctr] = NEGONE*h2; matdata[nnz_ctr+1] = NEGONE*h2;
        colindices[nnz_ctr] = i+1; colindices[nnz_ctr+1] = i+N;
        rowptrs[i+1] = rowptrs[i] + 2;
    }
    else if (i == N*N-1)
    {
      /* LAST ROW OF MATRIX */
      /* vals = h2 * [1, 1] */
      /* cols = [i-N, i-1] */
      matdata[nnz_ctr] = ONE*h2; matdata[nnz_ctr+1] = ONE*h2;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1;
      rowptrs[i+1] = rowptrs[i] + 2;
    }
    else if ((i+1) % N == 0)
    {
      /* vals = h2 * [1, 1, -1] */
      /* cols = [i-N, i-1, i+N] */
      matdata[nnz_ctr] = ONE*h2; matdata[nnz_ctr+1] = ONE*h2; matdata[nnz_ctr+2] = NEGONE*h2;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i+N;
      rowptrs[i+1] = rowptrs[i] + 3;
    }
    else if (i % N == 0)
    {
      /* vals = h2 * [1, -1, -1] */
      /* cols = [i-N, i+1, i+N] */
      matdata[nnz_ctr] = ONE*h2; matdata[nnz_ctr+1] = NEGONE*h2; matdata[nnz_ctr+2] = NEGONE*h2;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i+1; colindices[nnz_ctr+2] = i+N;
      rowptrs[i+1] = rowptrs[i] + 3;
    }
    else if (i < N)
    {
      /* vals = h2 * [1, -1, -1] */
      /* cols = [i-1, i+1, i+N] */
      matdata[nnz_ctr] = ONE*h2; matdata[nnz_ctr+1] = NEGONE*h2; matdata[nnz_ctr+2] = NEGONE*h2;
      colindices[nnz_ctr] = i-1; colindices[nnz_ctr+1] = i+1; colindices[nnz_ctr+2] = i+N;
      rowptrs[i+1] = rowptrs[i] + 3;
    }
    else if ((i+N) > (N*N-1)) 
    {
      /* vals = h2 * [1, 1, -1] */
      /* cols = [i-N, i-1, i+1] */
      matdata[nnz_ctr] = ONE*h2; matdata[nnz_ctr+1] = ONE*h2; matdata[nnz_ctr+2] = NEGONE*h2;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i+1;
      rowptrs[i+1] = rowptrs[i] + 3;
    }
    else 
    {
      /* vals = h2 * [1, 1, -1, -1] */
      /* cols = [i-N, i-1, i+1, i+N] */
      matdata[nnz_ctr] = ONE*h2; matdata[nnz_ctr+1] = ONE*h2; matdata[nnz_ctr+2] = NEGONE*h2; matdata[nnz_ctr+3] = NEGONE*h2;
      colindices[nnz_ctr] = i-N; colindices[nnz_ctr+1] = i-1; colindices[nnz_ctr+2] = i+1; colindices[nnz_ctr+3] = i+N;
      rowptrs[i+1] = rowptrs[i] + 4;
    }
  }
  
  /* Q = h2_inv * Laplacian + C * I + h2 * CenteredDiff */
  Q = SUNMatClone(Laplacian);
  SUNMatCopy(Laplacian, Q);             /* Q = Laplacian */
  /* CORRECT UP TO HERE */
  fail = SUNMatScaleAdd(h2_inv, Q, I);         /* Q = h2_inv * Laplacian + C * I */
  if (fail) {
    printf(">>> FAILED -- SUNMatScaleAdd returned %d \n", fail);
    SUNMatDestroy(Q); return(1);
  }
  //fail = SUNMatScaleAdd(ONE, Q, CenteredDiff); /* Q = h2_inv * Laplacian + C * I  + h2 * CenteredDiff*/
  /*if (fail) {
    printf(">>> FAILED -- SUNMatScaleAdd returned %d \n", fail);
    SUNMatDestroy(Q); return(1);
  }*/

  /* Clean up matrices */
  SUNMatDestroy(Laplacian);
  SUNMatDestroy(CenteredDiff);
  SUNMatDestroy(I);

  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to return the matrix for overlapping regions in RAS
 * ---------------------------------------------------------------------------*/
int OverlappingRestrictionMatrix2D(SUNMatrix Q)
{
  /*C = SUNSparseMatrix(5, 6, 9, CSR_MAT);*/
  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to return the matrix for nonoverlapping regions in RAS
 * ---------------------------------------------------------------------------*/
int NonOverlappingRestrictionMatrix2D(SUNMatrix Q)
{
  /*C = SUNSparseMatrix(5, 6, 9, CSR_MAT);*/
  return(0);
}


/* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 *
 * NEEDS TO BE UPDATED
 * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector u, realtype tol)
{
  /* Get vector data array */

  /* print the solution */

  /* solution error */

  /* print the solution error */

  printf("PASS\n");
  return(0);
}

/* -----------------------------------------------------------------------------
 * Check function return value
 *   opt == 0 check if returned NULL pointer
 *   opt == 1 check if returned a non-zero value
 * ---------------------------------------------------------------------------*/
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if the function returned a NULL pointer -- no memory allocated */
  if (opt == 0) {
    if (returnvalue == NULL) {
      fprintf(stderr, "\nERROR: %s() failed -- returned NULL\n\n", funcname);
      return(1);
    } else {
      return(0);
    }
  }

  /* Check if the function returned an non-zero value -- internal failure */
  if (opt == 1) {
    errflag = (int *) returnvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nERROR: %s() failed -- returned %d\n\n", funcname, *errflag);
      return(1);
    } else {
      return(0);
    }
  }

  /* if we make it here then opt was not 0 or 1 */
  fprintf(stderr, "\nERROR: check_retval failed -- Invalid opt value\n\n");
  return(1);
}
