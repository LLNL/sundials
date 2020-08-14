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

#include "kinsol/kinsol.h"               /* access to KINSOL func., consts.   */
#include "nvector/nvector_serial.h"      /* access to serial N_Vector         */
#include "sunmatrix/sunmatrix_band.h"    /* access to banded SUNMatrix        */
#include "sunlinsol/sunlinsol_band.h"    /* access to banded SUNMatrix solver */
#include "sundials/sundials_math.h"      /* access to SUNRSqrt function       */

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

#define TRUE         RCONST(26.651793729865151050262284) /* Norm of true solution */ 

/* UserData */ 
typedef struct {
  SUNMatrix subA;
  N_Vector v;
  N_Vector accum;
  N_Vector temp;
} *UserData;

/* Linear Convection-Diffusion fixed point function */
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

/* Function to calculate the residual for the current iterations guess u */
static int ConvectionDiffusionResidual2D(N_Vector v, N_Vector u);

/* Function to apply the matrix for overlapping regions in RAS */
static int OverlappingRestrictionMatrix2D(N_Vector v, N_Vector temp, int i, int j, int *subX, int *subY);

/* Function to apply the matrix for nonoverlapping regions in RAS*/
static int NonOverlappingRestrictionMatrix2D(N_Vector temp, N_Vector out, int i, int j);

/* Function to setup the subdomain system for overlapping regions in RAS*/
static int SubDomain2D(SUNMatrix subA, int bandwidth);

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
  realtype  tol     = SQRT(UNIT_ROUNDOFF);
  long int  mxiter  = 75;
  long int  maa     = 0;           /* no acceleration          */
  realtype  damping = RCONST(1.0); /* no damping               */
  int       orth    = 0;           /* standard MGS in Anderson */
  long int  nni, nfe;
  realtype* data;
  void*     kmem;
  //UserData  udata;                 /* contain matrices for function evals */ 

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
  /*udata = NULL;
  udata = (UserData)malloc(sizeof *udata);*/

  /* create 2d convection diffusion matrix here */
  /* allocate space for v and accum here */
  /*udata->v     = N_VNew_Serial(N*N);
  if (check_retval((void *)udata->v, "N_VNew_Serial", 0)) return(1);

  udata->accum = N_VNew_Serial(N*N);
  if (check_retval((void *)udata->accum, "N_VNew_Serial", 0)) return(1);*/

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
  
  /* Set userdata to be used in function evaluation */
  /*retval = KINSetUserData(kmem, udata);
  if (check_retval(&retval, "KINSetUserData", 1)) return(1);*/

  /* Set number of prior residuals used in Anderson acceleration */
  retval = KINSetMAA(kmem, maa);
  if (check_retval(&retval, "KINSetMAA", 1)) return(1);
  
  /* Set orthogonalization routine used in Anderson acceleration */
  retval = KINSetOrthAA(kmem, orth);
  if (check_retval(&retval, "KINSetOrthAA", 1)) return(1);

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
  /*free(udata);*/

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
  SUNMatrix subA;
  SUNLinearSolver LS;
  N_Vector  v, accum, temp, temp2; 
  UserData udata;
  int i, j, k, retval, subX, subY;
  realtype* vdata;
 
  /* Grab matrices and vectors from user_data */ 
  /*udata = (UserData)user_data;*/
  /*v     = udata->v;
  printf("after grabbing v\n");
  accum = udata->accum;*/
  v     = N_VNew_Serial(N*N);
  if (check_retval((void *)v, "N_VNew_Serial", 0)) return(1);
  accum = N_VNew_Serial(N*N);
  if (check_retval((void *)accum, "N_VNew_Serial", 0)) return(1);

  /* v = f - A * u */
  retval = ConvectionDiffusionResidual2D(v, u);
  if (check_retval(&retval, "ConvectionDiffusionResidual2D", 1)) return(1);
  if (check_retval((void *)v, "v NULL", 0)) return(1);
  /*vdata = N_VGetArrayPointer(v);
  for (i = 10; i < 15; i++) printf("v[%d] %e\n", i, vdata[i]);*/
  printf("v norm %e\n", SUNRsqrt(N_VDotProd(v, v)));

  /* accum = 0 */
  N_VConst(ZERO, accum);

  for (i = 0; i < NUMREGIONS; i++)
  {
    for (j = 0; j < NUMREGIONS; j++)
    {
      /* setup R = OverlappingRestrictionMatrix2D */
      /* temp = R*v */
      temp  = N_VNew_Serial(N);
      if (check_retval((void *)temp, "N_VNew_Serial", 0)) return(1);
      retval = OverlappingRestrictionMatrix2D(v, temp, i, j, &subX, &subY);
      if (check_retval(&retval, "OverlappingRestrictionMatrix2D", 1)) return(1);
      if (check_retval((void *)temp, "temp NULL", 0)) return(1);

      /* setup subA = R * A * R^T */
      k = N_VGetLength_Serial(temp);
      subA = SUNBandMatrix(k, subX, subX);
      retval = SubDomain2D(subA, subX);
      if (check_retval(&retval, "SubDomain2D", 1)) return(1);

      /* solve subA * temp = temp (overwrite old temp) */
      LS = SUNLinSol_Band(temp, subA);
      retval = SUNLinSolInitialize(LS);
      check_retval(&retval, "SUNLinSolInitialize", 1);
      retval = SUNLinSolSetup(LS, subA);
      check_retval(&retval, "SUNLinSolSetup", 1);
      retval = SUNLinSolSolve(LS, subA, temp, temp, UNIT_ROUNDOFF);
      check_retval(&retval, "SUNLinSolSolve", 1);
      printf("temp norm %e\n", SUNRsqrt(N_VDotProd(temp, temp)));

      /* setup L^T = NonOverlappingRestrictionMatrix2D */
      /* temp2 = L^T * temp */
      temp2  = N_VNew_Serial(N*N);
      if (check_retval((void *)temp2, "N_VNew_Serial", 0)) return(1);
      retval = NonOverlappingRestrictionMatrix2D(temp, temp2, i, j);
      if (check_retval(&retval, "NonOverlappingRestrictionMatrix2D", 1)) return(1);
      if (check_retval((void *)temp2, "temp2 NULL", 0)) return(1);

      /* accum = accum + temp */
      N_VLinearSum(ONE, accum, ONE, temp2, accum);

      N_VDestroy(temp); /* Destroying and reallocating for each subdomain
                           -- inefficient, but it's serial */
      //N_VDestroy(temp2);
      SUNMatDestroy(subA); /* For now while we're not allocating them beforehand 
                              -- destroy and recreate after each iteration */
      SUNLinSolFree(LS);   /* Need to create new Linear Solver object next iteration
                              -- with new subA */ 
    }  
  }
  /* u = u + accum (overwrite old u) */
  /*printf("accum norm %e\n", SUNRsqrt(N_VDotProd(accum, accum)));*/
  N_VLinearSum(ONE, accum, ONE, u, u);

  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to calculate the residual for the current iterate u 
 * ---------------------------------------------------------------------------*/
int ConvectionDiffusionResidual2D(N_Vector v, N_Vector u)
{
  realtype *vdata, *uvals;
  realtype h, h2_inv, hhalf_inv, a, b, c, d, e;
  int i;
  
  /* Define some constants */
  h = 1.0 / (N+1);
  h2_inv    = 1.0 / (h*h); /* 1/h^2 */
  hhalf_inv = 1.0 / (2*h);  /* 1/(2*h) */
  a = h2_inv * -4.0 + C;
  b = h2_inv + D * hhalf_inv;
  c = h2_inv + D * hhalf_inv;
  d = h2_inv + D * hhalf_inv;
  e = h2_inv + D * hhalf_inv;

  /* v = f - A * u */
  vdata = N_VGetArrayPointer(v);
  uvals = N_VGetArrayPointer(u);

  for (i=0; i<N*N; i++) {
    vdata[i] = -10;

    if (i == 0) {
      vdata[i] -= (a*uvals[i] + c*uvals[i+1] + e*uvals[i+N]);
    }
    else if (i == 1) {
      vdata[i] -= (a*uvals[i] + b*uvals[i-1] + c*uvals[i+1] + e*uvals[i+N]);
    }
    else if (i == (N*N-1)) {
      vdata[i] -= (a*uvals[i] + b*uvals[i-1] + c*uvals[i+1] + d*uvals[i-N]);
    }
    else if (i == (N*N-2)) {
      vdata[i] -= (a*uvals[i] + b*uvals[i-1] + d*uvals[i-N]);
    }
    else {
      vdata[i] -= a * uvals[i];
      if ((i%N) != 0)     vdata[i] -= b * uvals[i-1];
      if (((i+1)%N) != 0) vdata[i] -= c * uvals[i+1];
      if ((i-N) > 0)      vdata[i] -= d * uvals[i-N];
      if ((i+N) < (N*N))  vdata[i] -= e * uvals[i+N];
    }
  }

  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to return the matrix for overlapping regions in RAS
 * ---------------------------------------------------------------------------*/
static int OverlappingRestrictionMatrix2D(N_Vector v, N_Vector temp, int i, int j, int *subX, int *subY)
{
  int ind_i, ind_j;
  int subDomainLengthX, subDomainLengthY;
  int x, y, xstart, xend, ystart, yend;
  realtype *vdata, *temp_data;

  subDomainLengthX = N / NUMREGIONS;
  subDomainLengthY = subDomainLengthX;
  xstart = i * subDomainLengthX;
  xend = xstart + subDomainLengthX - 1;
  ystart = j * subDomainLengthY;
  yend = ystart + subDomainLengthY - 1;

  if (i > 0) {
    xstart = xstart - OVERLAPWIDTH;
    subDomainLengthX = subDomainLengthX + OVERLAPWIDTH; 
  }
  if (i < (NUMREGIONS-1)) {
    xend = xend + OVERLAPWIDTH;
    subDomainLengthX = subDomainLengthX + OVERLAPWIDTH; 
  }
  
  if (j > 0) {
    ystart = ystart - OVERLAPWIDTH;
    subDomainLengthY = subDomainLengthY + OVERLAPWIDTH; 
  }
  if (j < (NUMREGIONS-1)) {
    yend = yend + OVERLAPWIDTH;
    subDomainLengthY = subDomainLengthY + OVERLAPWIDTH; 
  }
 
  /* temp = R * v */ 
  //temp  = N_VNew_Serial(subDomainLengthX * subDomainLengthY);
  /* Create temp data array of correct length and attach to temp vector */
  free(NV_DATA_S(temp));
  temp_data = NULL;
  temp_data = (realtype *) malloc(subDomainLengthX * subDomainLengthY * sizeof(realtype));
  NV_DATA_S(temp) = temp_data; 
  NV_LENGTH_S(temp) = subDomainLengthX * subDomainLengthY;
  N_VConst(ZERO, temp);
 
  vdata = N_VGetArrayPointer(v); 
  /* Double nested for loop for application of R*v */
  for (x = xstart; x < xend+1; x++) {
    for (y = ystart; y < yend+1; y++) {
      ind_i = x-xstart + (y-ystart)*subDomainLengthX;
      ind_j = y * N + x;
      temp_data[ind_i] = vdata[ind_j];
    }
  }

  *subX = subDomainLengthX;
  *subY = subDomainLengthY;

  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to return the matrix for nonoverlapping regions in RAS
 * ---------------------------------------------------------------------------*/
int NonOverlappingRestrictionMatrix2D(N_Vector temp, N_Vector out, int i, int j)
{
  int ind_i, ind_j, xoffset, yoffset;
  int subDomainLengthX, subDomainLengthY;
  int x, y, xstart, xend, ystart, yend;
  realtype *out_data, *temp_data; 

  subDomainLengthX = N / NUMREGIONS;
  subDomainLengthY = subDomainLengthX;
  xstart = i * subDomainLengthX;
  xend = xstart + subDomainLengthX - 1;
  ystart = j * subDomainLengthY;
  yend = ystart + subDomainLengthY - 1;

  if (i > 0) {
    subDomainLengthX = subDomainLengthX + OVERLAPWIDTH; 
  }
  if (i < (NUMREGIONS-1)) {
    subDomainLengthX = subDomainLengthX + OVERLAPWIDTH; 
  }
  
  if (j > 0) {
    subDomainLengthY = subDomainLengthY + OVERLAPWIDTH; 
  }
  if (j < (NUMREGIONS-1)) {
    subDomainLengthY = subDomainLengthY + OVERLAPWIDTH; 
  }

  /* out = L^T * temp */
  //out       = N_VNew_Serial(N*N);
  //if (check_retval((void *)out, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, out);
  out_data  = N_VGetArrayPointer(out);
  temp_data = N_VGetArrayPointer(out);

  /* Double nested for loop for application of L^T*temp */
  for (x = xstart; x < xend; x++) {
    for (y = ystart; y < yend; y++) {
      xoffset = 0;
      if (i > 0) xoffset = OVERLAPWIDTH;
      yoffset = 0;
      if (j > 0) yoffset = OVERLAPWIDTH; 

      ind_i = x-xstart + xoffset + (y-ystart+yoffset)*subDomainLengthX;
      ind_j = y * N + x - 1; 
      out_data[ind_j] = temp_data[ind_i]; 
    }
  }

  return(0);
}

/* -----------------------------------------------------------------------------
 * Function to form the subdomain matrix for overlapping regions in RAS
 * ---------------------------------------------------------------------------*/
int SubDomain2D(SUNMatrix subA, int bandwidth)
{
  realtype h, h2_inv, hhalf_inv, a, b;
  realtype *colj;
  int i, j, k, kstart, kend, cols; 

  /* Define some constants */
  h = 1.0 / (N+1);
  h2_inv    = 1.0 / (h*h); /* 1/h^2 */
  hhalf_inv = 1.0 / (2*h);  /* 1/(2*h) */
  a = h2_inv * -4.0 + C;
  b = h2_inv + D * hhalf_inv;
 
  cols = SUNBandMatrix_Columns(subA); 
  for (j=0; j<cols; j++) {
    colj = SUNBandMatrix_Column(subA, j); /* Grab column */

    /* Zero out the whole column */
    kstart = (j<bandwidth) ? -j : -bandwidth;
    kend = (j>cols-1-bandwidth) ? cols-1-j : bandwidth;
    for (k=kstart; k < kend; k++) {
      colj[k] = ZERO;
    }
   
    /* Add the nonzero values */ 
    colj[0] = a;                              /* Store diagonal */
    if ((j % bandwidth) != 0) colj[-1] = b;   /* Store above diagonal */
    if (((j+1)%bandwidth) != 0) colj[1] = b;  /* Store below diagonal */
    if (j < bandwidth) colj[bandwidth] = b;   /* Store below diagonal bandwidth value */  
    if (j > bandwidth) colj[-bandwidth] = b;  /* Store above diagonal bandwidth value */  

  }

  return(0);
}


/* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 *
 * NEEDS TO BE UPDATED
 * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector u, realtype tol)
{
  realtype* data = NULL;
  realtype  norm, norm_error;

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* print the norm of solution */
  norm = SUNRsqrt(N_VDotProd(u, u));
  printf("Computed solution norm:\n");
  printf("    ||u||_2 = %"GSYM"\n", norm);

  /* solution error */
  norm_error = ABS(norm - TRUE);

  /* print the solution error */
  printf("Solution error:\n");
  printf("    ex = %"GSYM"\n", norm_error);

  if (norm_error > tol) {
    printf("FAIL\n");
    return(1);
  }

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
