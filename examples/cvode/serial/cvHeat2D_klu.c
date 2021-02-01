/* -----------------------------------------------------------------
 * Programmer(s): Chris Nguyen
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem for CVODE: 2D heat equation, serial, sparse.
 * Based on idaHeat2D_klu.c and cvRoberts_klu.c
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the KLU sparse solver.
 *
 * The PDE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MGRID x MGRID
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MGRID^2. Here MGRID = 10.
 *
 * The system is solved with CVODE using the direct sparse linear system
 * solver, half-bandwidths equal to M, and default
 * difference-quotient Jacobian.
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 * -----------------------------------------------------------------*/

//////////////////////////////////////////////////////////////////////////
// Note: was still trying to get case MGRID=3 to work amd convert to 
// CVODE format. Any comments using // were left in for future work notes
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>                /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>    /* access to KLU sparse direct solver   */
#include <sundials/sundials_types.h>    /* defs. of realtype, sunindextype      */

/* User-defined vector and matrix accessor macro: Ith */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.
*/

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */

/* Problem Constants */
/*
#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.0)
*/
#define MGRID 3
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.0)
#define NEQ   MGRID*MGRID      /* number of equations    */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.1)      /* first output time      */
#define TMULT RCONST(2.0)      /* output time factor     */
#define NOUT  12               /* number of output times */
#define TOTAL 4*MGRID+8*(MGRID-2)+(MGRID-4)*(MGRID+4*(MGRID-2)) /* total num of nonzero elements */

/* Type: UserData */

typedef struct {
  sunindextype mm;
  realtype dx;
  realtype coeff;
} *UserData;

/* Prototypes of functions called by CVODE */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static int jacHeat(realtype tt, N_Vector yy, N_Vector fy, SUNMatrix JacMat, 
                   void *user_data, N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Exact same setup as jacHeat. Function needed for special case MGRID=3  */
static int jacHeat3(realtype tt, N_Vector yy, N_Vector fy, SUNMatrix JacMat, 
                    void *user_data, N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */

static void PrintHeader(realtype rtol, realtype atol);
static void PrintOutput(void *mem, realtype t, N_Vector uu);
static void PrintFinalStats(void *cvode_mem);
static int SetInitialProfile(UserData data, N_Vector uu, N_Vector res);

static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *cvode_mem;
  UserData data;
  N_Vector uu, constraints, res; 
  SUNMatrix A;
  SUNLinearSolver LS;
  int retval, iout;
  sunindextype mu, ml, netf, ncfn;
  realtype rtol, atol, t0, t1, tout, tret;

  int nnz; /* number of non-zeroes  */
  
  cvode_mem = NULL;
  data = NULL;
  uu = constraints = res = NULL;
  A = NULL;
  LS = NULL;

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew_Serial(NEQ);
  if(check_retval((void *)uu, "N_VNew_Serial", 0)) return(1);
  constraints = N_VNew_Serial(NEQ);
  if(check_retval((void *)constraints, "N_VNew_Serial", 0)) return(1);
  res = N_VNew_Serial(NEQ);
  if(check_retval((void *)res, "N_VNew_Serial", 0)) return(1);
 
  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  if(check_retval((void *)data, "malloc", 2)) return(1);
  data->mm = MGRID;
  data->dx = ONE/(MGRID - ONE);
  data->coeff = ONE/( (data->dx) * (data->dx) );

  /* Initialize uu, up. */
  SetInitialProfile(data, uu, res);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = RCONST(0.01);
  rtol = ZERO;
  atol = RCONST(1.0e-8);

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /*
  retval = IDASetConstraints(cvode_mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1)) return(1);
  */
  N_VDestroy(constraints);

  retval = CVodeInit(cvode_mem, f, t0, uu);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  retval = CVodeSStolerances(cvode_mem, rtol, atol);
  if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create sparse SUNMatrix for use in linear solves */
  nnz = NEQ * NEQ;
  A = SUNSparseMatrix(NEQ, NEQ, nnz, CSC_MAT);
  if(check_retval((void *)A, "SUNSparseMatrix", 0)) return(1);

  /* Create KLU solver object for use by CVode */
  LS = SUNLinSol_KLU(y, A);
  if(check_retval((void *)LS, "SUNLinSol_KLU", 0)) return(1);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine Jac based pm size of Jac */
  if(MGRID >= 4){
    retval = CVodeSetJacFn(cvode_mem, jacHeat);
  }
  /* special case MGRID=3  */
  else if(MGRID==3){
    retval = CVodeSetJacFn(cvode_mem, jacHeat3);
  }
  /* MGRID<=2 is pure boundary points, nothing to solve  */
  else{
    printf("MGRID size is too small to run.\n");
    return(1);
  }
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* Print output heading. */
  PrintHeader(rtol, atol);
  
  PrintOutput(cvode_mem, t0, uu);

  /* Loop over output times, call CVode, and print results. */
  
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {
    
    retval = CVode(cvode_mem, tout, uu, &tret, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1)) return(1);

    PrintOutput(cvode_mem, tret, uu);
  
  }
  
  /* Print remaining counters and free memory. */
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);

  PrintFinalStats(cvode_mem);

  /* Free memory  */

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(uu);
  N_VDestroy(res);
  free(data);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODE
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute function f(t,u). 
 */
/////////////////////////////////////////
// How to define f for heat equation like
// in Roberts problem?
////////////////////////////////////////
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype u1, ud1;

  u1 = Ith(u,1);

  ud1 = Ith(udot,1) = RCONST(4)*u1;

  return(0);
}

/*
int heatres(realtype tres, N_Vector uu, N_Vector resval, 
            void *user_data)
{
  sunindextype mm, i, j, offset, loc;
  realtype *uv, *resv, coeff;
  UserData data;
  
  uv = N_VGetArrayPointer(uu); resv = N_VGetArrayPointer(resval);

  data = (UserData)user_data;
  mm = data->mm;
  coeff = data->coeff;
  
  //Initialize resval to uu, to take care of boundary equations.
  N_VScale(ZERO, uu, resval);
  
  //Loop over interior points; set res = up - (central difference). 
  for (j = 1; j < mm-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      resv[loc] = upv[loc] - coeff * 
	  (uv[loc-1] + uv[loc+1] + uv[loc-mm] + uv[loc+mm] - RCONST(4.0)*uv[loc]);
    }
  }
  
  return(0);
  
  }
*/

/* Jacobian matrix setup for MGRID=3  */
static int jacHeat3(realtype tt, N_Vector yy, N_Vector fy, SUNMatrix JacMat, 
                    void *user_data, N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(JacMat);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(JacMat);
  realtype *data = SUNSparseMatrix_Data(JacMat);
  realtype dx =  ONE/(MGRID - ONE);
  realtype beta = RCONST(4.0)/(dx*dx);

  yval = N_VGetArrayPointer(yy);

  SUNMatZero(JacMat); // initialize Jacobian matrix
  
  // set up number of elements in each column 
  
  colptrs[0] = 0;
  colptrs[1] = 1;
 
  /*
  JacMat -> colptrs[0]  = 0;
  JacMat -> colptrs[1]  = 1;
  JacMat -> colptrs[2]  = 3;
  JacMat -> colptrs[3]  = 4;
  JacMat -> colptrs[4]  = 6;
  JacMat -> colptrs[5]  = 7;
  JacMat -> colptrs[6]  = 9;
  JacMat -> colptrs[7]  = 10;
  JacMat -> colptrs[8]  = 12;
  JacMat -> colptrs[9]  = 13;
  */
  //set up data and row values stored 

  data[0] = TWO*TWO;
  rowvals[0] = 0;

  /*
  JacMat -> data[0] = ONE;
  JacMat -> rowvals[0] = 0;  
  JacMat -> data[1] = ONE;
  JacMat -> rowvals[1] = 1;
  JacMat -> data[2] = -ONE/(dx*dx);
  JacMat -> rowvals[2] = 4;  
  JacMat -> data[3] = ONE;
  JacMat -> rowvals[3] = 2;
  JacMat -> data[4] = ONE;
  JacMat -> rowvals[4] = 3;
  JacMat -> data[5] = -ONE/(dx*dx);
  JacMat -> rowvals[5] = 4;  
  JacMat -> data[6] = beta;
  JacMat -> rowvals[6] = 4;
  JacMat -> data[7] = -ONE/(dx*dx);
  JacMat -> rowvals[7] = 4;
  JacMat -> data[8] = ONE;
  JacMat -> rowvals[8] = 5;
  JacMat -> data[9] = ONE;
  JacMat -> rowvals[9] = 6; 
  JacMat -> data[10] = -ONE/(dx*dx);
  JacMat -> rowvals[10] = 4;
  JacMat -> data[11] = ONE;
  JacMat -> rowvals[11] = 7;
  JacMat -> data[12] = ONE;
  JacMat -> rowvals[12] = 8;
  */

  return(0);
}

//////////////////////////////////////////////////////////////////
// Remove boundary points in matrix, redefine without
//////////////////////////////////////////////////////////////////

/* Jacobian matrix setup for MGRID>=4  */
static int jacHeat(realtype tt, N_Vector yy, N_Vector fy, SUNMatrix JacMat, 
                   void *user_data, N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(J);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(J);
  realtype *data = SUNSparseMatrix_Data(J);
  realtype dx =  ONE/(MGRID - ONE);
  realtype beta = RCONST(4.0)/(dx*dx);
  int i,j, repeat=0;

  yval = N_VGetArrayPointer(yy);

  SUNMatZero(JacMat); /* initialize Jacobian matrix  */

  /* 
   *-----------------------------------------------
   * set up number of elements in each column 
   *-----------------------------------------------
   */
  
  /**** first column block ****/
  colptrs[0] = 0;
  colptrs[1] = 1;
  /* count by twos in the middle  */
  for(i=2;i<MGRID;i++) colptrs[i] = colptrs[i-1]+2;
  colptrs[MGRID] = 2*MGRID-2;

  /**** second column block ****/
  colptrs[MGRID+1] = 2*MGRID;
  colptrs[MGRID+2] = 2*MGRID+3;
  /* count by fours in the middle */
  for(i=0;i<MGRID-4;i++) colptrs[MGRID+3+i] = colptrs[MGRID+3+i-1]+4;
  colptrs[2*MGRID-1] = 2*MGRID+4*(MGRID-2)-2;
  colptrs[2*MGRID]   = 2*MGRID+4*(MGRID-2);
  
  /**** repeated (MGRID-4 times) middle column blocks ****/
  for(i=0;i<MGRID-4;i++){
    colptrs[2*MGRID+1+repeat]   = (colptrs[2*MGRID+1+repeat-1])+2;
    colptrs[2*MGRID+1+repeat+1] = (colptrs[2*MGRID+1+repeat])+4;
    
    /* count by fives in the middle */
    for(j=0;j<MGRID-4;j++) colptrs[2*MGRID+1+repeat+2+j] = colptrs[2*MGRID+1+repeat+1+j]+5;
   
    colptrs[2*MGRID+1+repeat+(MGRID-4)+2] = colptrs[2*MGRID+1+repeat+(MGRID-4)+1]+4;
    colptrs[2*MGRID+1+repeat+(MGRID-4)+3] = colptrs[2*MGRID+1+repeat+(MGRID-4)+2]+2;  

    repeat+=MGRID; /* shift that accounts for accumulated number of columns */
  }
  
  /**** last-1 column block ****/
  colptrs[MGRID*MGRID-2*MGRID+1] = TOTAL-2*MGRID-4*(MGRID-2)+2;
  colptrs[MGRID*MGRID-2*MGRID+2] = TOTAL-2*MGRID-4*(MGRID-2)+5;
  /* count by fours in the middle */
  for(i=0;i<MGRID-4;i++) colptrs[MGRID*MGRID-2*MGRID+3+i] = colptrs[MGRID*MGRID-2*MGRID+3+i-1]+4;
  colptrs[MGRID*MGRID-MGRID-1] = TOTAL-2*MGRID;
  colptrs[MGRID*MGRID-MGRID]   = TOTAL-2*MGRID+2;

  /**** last column block ****/
  colptrs[MGRID*MGRID-MGRID+1] = TOTAL-MGRID-(MGRID-2)+1;
  /* count by twos in the middle */
  for(i=0;i<MGRID-2;i++) colptrs[MGRID*MGRID-MGRID+2+i] = colptrs[MGRID*MGRID-MGRID+2+i-1]+2;
  colptrs[MGRID*MGRID-1] = TOTAL-1;
  colptrs[MGRID*MGRID]   = TOTAL;
  

  /*
   *-----------------------------------------------
   * set up data stored
   *-----------------------------------------------
   */

  /**** first column block ****/
  data[0] = ONE;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=1;i<MGRID+(MGRID-2)  ;i+=2) data[i] = ONE;
  for(i=2;i<MGRID+(MGRID-2)-1;i+=2) data[i] = -ONE/(dx*dx);

  /**** second column block ****/
  data[MGRID+MGRID-2] = ONE;
  data[MGRID+MGRID-1] = -ONE/(dx*dx);
  data[MGRID+MGRID]   = beta;
  data[MGRID+MGRID+1] = -ONE/(dx*dx);
  data[MGRID+MGRID+2] = -ONE/(dx*dx);
  /* middle data elements */
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+3+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+4+4*i] = beta;
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+5+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+6+4*i] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-5] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-4] = beta;
  data[2*MGRID+4*(MGRID-2)-3] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-2] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-1] = ONE;
    
  /**** repeated (MGRID-4 times) middle column blocks ****/
  repeat=0;
  for(i=0;i<MGRID-4;i++){
    data[2*MGRID+4*(MGRID-2)+repeat]   = ONE;
    data[2*MGRID+4*(MGRID-2)+repeat+1] = -ONE/(dx*dx);
    
    data[2*MGRID+4*(MGRID-2)+repeat+2] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+3] = beta;
    data[2*MGRID+4*(MGRID-2)+repeat+4] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+5] = -ONE/(dx*dx);

    /* 5 in 5*j chosen since there are 5 elements in each column */
    /* this column loops MGRID-4 times within the outer loop */
    for(j=0;j<MGRID-4;j++){
      data[2*MGRID+4*(MGRID-2)+repeat+6+5*j]  = -ONE/(dx*dx);
      data[2*MGRID+4*(MGRID-2)+repeat+7+5*j]  = -ONE/(dx*dx);
      data[2*MGRID+4*(MGRID-2)+repeat+8+5*j]  = beta;
      data[2*MGRID+4*(MGRID-2)+repeat+9+5*j]  = -ONE/(dx*dx);
      data[2*MGRID+4*(MGRID-2)+repeat+10+5*j] = -ONE/(dx*dx);
    }
    
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+6] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+7] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+8] = beta;
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+9] = -ONE/(dx*dx);
    
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+10] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+11] = ONE;
    
    repeat+=MGRID+4*(MGRID-2); /* shift that accounts for accumulated columns and elements */
  }
  
  /**** last-1 column block ****/
  data[TOTAL-6*(MGRID-2)-4] = ONE;
  data[TOTAL-6*(MGRID-2)-3] = -ONE/(dx*dx);
  data[TOTAL-6*(MGRID-2)-2] = -ONE/(dx*dx);
  data[TOTAL-6*(MGRID-2)-1] = beta;
  data[TOTAL-6*(MGRID-2)  ] = -ONE/(dx*dx);
  /* middle data elements */
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+1+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+2+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+3+4*i] = beta;
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+4+4*i] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-7] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-6] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-5] = beta;
  data[TOTAL-2*(MGRID-2)-4] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-3] = ONE;

  /**** last column block ****/
  data[TOTAL-2*(MGRID-2)-2] = ONE;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=TOTAL-2*(MGRID-2)-1;i<TOTAL-2;i+=2) data[i] = -ONE/(dx*dx);
  for(i=TOTAL-2*(MGRID-2)  ;i<TOTAL-1;i+=2) data[i] = ONE;
  data[TOTAL-1] = ONE;
  
  /*
   *-----------------------------------------------
   * row values 
   *-----------------------------------------------
   */

  /**** first block ****/
  rowvals[0] = 0;
  /* alternating pattern in data, separate loop for each pattern */
  for(i=1;i<MGRID+(MGRID-2)  ;i+=2) rowvals[i] = (i+1)/2;
  for(i=2;i<MGRID+(MGRID-2)-1;i+=2) rowvals[i] = i/2+MGRID; // i+1 unnecessary here
  
  /**** second column block ****/
  rowvals[MGRID+MGRID-2] = MGRID;
  rowvals[MGRID+MGRID-1] = MGRID+1;
  rowvals[MGRID+MGRID]   = MGRID+1;
  rowvals[MGRID+MGRID+1] = MGRID+2;
  rowvals[MGRID+MGRID+2] = 2*MGRID+1;
  /* middle row values */
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+3+4*i] = MGRID+1+i;
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+4+4*i] = MGRID+2+i;
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+5+4*i] = MGRID+3+i;
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+6+4*i] = 2*MGRID+2+i;
  rowvals[2*MGRID+4*(MGRID-2)-5] = MGRID+(MGRID-2)-1;
  rowvals[2*MGRID+4*(MGRID-2)-4] = MGRID+(MGRID-2); // starting from here, add two diag patterns
  rowvals[2*MGRID+4*(MGRID-2)-3] = 2*MGRID+(MGRID-2);
  rowvals[2*MGRID+4*(MGRID-2)-2] = MGRID+(MGRID-2);
  rowvals[2*MGRID+4*(MGRID-2)-1] = MGRID+(MGRID-2)+1;
  
  /**** repeated (MGRID-4 times) middle column blocks ****/
  repeat=0;
  for(i=0;i<MGRID-4;i++){
    rowvals[2*MGRID+4*(MGRID-2)+repeat]   = MGRID+(MGRID-2)+2+MGRID*i;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+1] = MGRID+(MGRID-2)+2+MGRID*i+1;
    
    rowvals[2*MGRID+4*(MGRID-2)+repeat+2] = MGRID+(MGRID-2)+2+MGRID*i+1-MGRID;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+3] = MGRID+(MGRID-2)+2+MGRID*i+1;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+4] = MGRID+(MGRID-2)+2+MGRID*i+2; //*this
    rowvals[2*MGRID+4*(MGRID-2)+repeat+5] = MGRID+(MGRID-2)+2+MGRID*i+1+MGRID;
    
    /* 5 in 5*j chosen since there are 5 elements in each column */
    /* column repeats MGRID-4 times within the outer loop */
    for(j=0;j<MGRID-4;j++){
      rowvals[2*MGRID+4*(MGRID-2)+repeat+6+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+1-MGRID+1+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+7+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+1+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+8+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+2+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+9+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+2+1+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+10+5*j] = MGRID+(MGRID-2)+2+MGRID*i+1+MGRID+1+j;
    }
    
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+6] = MGRID+(MGRID-2)+2+MGRID*i-2;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+7] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID-1;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+8] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID;//*this+MGRID
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+9] = MGRID+(MGRID-2)+2+MGRID*i-2+2*MGRID;

    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+10] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+11] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID+1;
    
    repeat+=MGRID+4*(MGRID-2); /* shift that accounts for accumulated columns and elements */
  }  

  
  /**** last-1 column block ****/
  rowvals[TOTAL-6*(MGRID-2)-4] = MGRID*MGRID-1-2*(MGRID-1)-1;
  rowvals[TOTAL-6*(MGRID-2)-3] = MGRID*MGRID-1-2*(MGRID-1); // starting with this as base
  rowvals[TOTAL-6*(MGRID-2)-2] = MGRID*MGRID-1-2*(MGRID-1)-MGRID;
  rowvals[TOTAL-6*(MGRID-2)-1] = MGRID*MGRID-1-2*(MGRID-1);
  rowvals[TOTAL-6*(MGRID-2)  ] = MGRID*MGRID-1-2*(MGRID-1)+1;
  /* middle row values */
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+1+4*i] = MGRID*MGRID-1-2*(MGRID-1)-MGRID+1+i;
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+2+4*i] = MGRID*MGRID-1-2*(MGRID-1)+i;
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+3+4*i] = MGRID*MGRID-1-2*(MGRID-1)+1+i;//copied above
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+4+4*i] = MGRID*MGRID-1-2*(MGRID-1)+2+i;
  rowvals[TOTAL-2*(MGRID-2)-7] = MGRID*MGRID-2*MGRID-2;
  rowvals[TOTAL-2*(MGRID-2)-6] = MGRID*MGRID-MGRID-3;
  rowvals[TOTAL-2*(MGRID-2)-5] = MGRID*MGRID-MGRID-2;
  rowvals[TOTAL-2*(MGRID-2)-4] = MGRID*MGRID-MGRID-2;
  rowvals[TOTAL-2*(MGRID-2)-3] = MGRID*MGRID-MGRID-1;

  /* last column block */
  rowvals[TOTAL-2*(MGRID-2)-2] = MGRID*MGRID-MGRID;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=0;i<(MGRID-2);i++) rowvals[TOTAL-2*(MGRID-2)-1+2*i] = MGRID*MGRID-2*MGRID+1+i;
  for(i=0;i<(MGRID-2);i++) rowvals[TOTAL-2*(MGRID-2)  +2*i] = MGRID*MGRID-MGRID+1+i;  
  rowvals[TOTAL-1] = MGRID*MGRID-1;

  //SUNSparseMatrix_Print(JacMat);
 
  return(0);
}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize u vector.       
 */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector res)
{
  realtype xfact, yfact, *udata;
  sunindextype mm, mm1, i, j, offset, loc;
  
  mm = data->mm;
  mm1 = mm - 1;
  
  udata = N_VGetArrayPointer(uu);

  /* Initialize uu on all grid points. */ 
  for (j = 0; j < mm; j++) {
    yfact = data->dx * j;
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      xfact = data->dx * i;
      loc = offset + i;
      udata[loc] = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }

  /* heatres sets res to negative of PDE RHS values at interior points. */
  heatres(ZERO, uu, res, data);
  
  /* Finally, set values of u, up, and id at boundary points. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) {
        udata[loc] = BVAL;  
      }
    }
  }
  
  return(0);

}

/* 
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, realtype atol)
{
  printf("\ncvHeat2D_klu: Heat equation, serial example problem for CVODE\n");
  printf("          Discretized heat equation on 2D unit square.\n");
  printf("          Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("          Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Linear solver: CVKLU, sparse direct solver \n");
  printf("       difference quotient Jacobian \n");
  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time       umax     k  nst  nni  nje   nre     h       \n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n");
}

/*
 * Print Output
 */

static void PrintOutput(void *mem, realtype t, N_Vector uu)
{
  int retval;
  realtype umax, hused;
  long int nst, nni, nje, nre, nreLS;
  int kused;

  umax = N_VMaxNorm(uu);
  
  retval = CVodeGetLastOrder(mem, &kused);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetNumSteps(mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumNonlinSolvIters(mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  //retval = CVodeGetNumResEvals(mem, &nre);
  //check_retval(&retval, "CVodeGetNumResEvals", 1);
  retval = CVodeGetLastStep(mem, &hused);
  check_retval(&retval, "CVodeGetLastStep", 1);
  retval = CVodeGetNumJacEvals(mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);

/*
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %9.2Le \n",
         t, umax, kused, nst, nni, nje, nre, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, hused);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, hused);
#endif
*/

#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %9.2Le \n",
         t, umax, kused, nst, nni, nje, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, hused);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, hused);
#endif

}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nni, ncfn, netf, nge;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/* 
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", 
              funcname, *retval);
      return(1); 
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
