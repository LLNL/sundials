/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Chris Nguyen
 * -----------------------------------------------------------------
 * Example problem for CVODE: 2D heat equation, serial, sparse. 
 * Based on idaHeat2D_klu.c and cvRoberts_klu.c
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the sparse solver CVKLU.
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
 * difference-quotient Jacobian. For purposes of illustration,
 * The constraints u >= 0 are posed for all components. Output is 
 * taken at t = 0, .01, .02, .04, ..., 10.24.
 * -----------------------------------------------------------------
 */

//////////////////////////////////////////////////////////////////////////
// Note: was still trying to get case MGRID=3 to work amd convert to 
// CVODE format. Any comments using // were left in for future work notes
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>
#include <cvode/cvode_klu.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

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
#define NEQ   MGRID*MGRID      /* number of equations  */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.1)      /* first output time      */
#define TMULT RCONST(2.0)     /* output time factor     */
#define NOUT  12               /* number of output times */
#define TOTAL 4*MGRID+8*(MGRID-2)+(MGRID-4)*(MGRID+4*(MGRID-2)) /* total num of nonzero elements */

/* Type: UserData */

typedef struct {
  long int mm;
  realtype dx;
  realtype coeff;
} *UserData;

/* Prototypes of functions called by CVODE */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static int jacHeat(realtype tt,
	    N_Vector yy, N_Vector fy, N_Vector resvec,
	    SlsMat JacMat, void *user_data,
	    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Exact same setup as jacHeat. Function needed for special case MGRID=3  */
static int jacHeat3(realtype tt, 
	    N_Vector yy, N_Vector fy, N_Vector resvec,
	    SlsMat JacMat, void *user_data,
	    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */

static void PrintHeader(realtype rtol, realtype atol);
static void PrintOutput(void *mem, realtype t, N_Vector uu);
static void PrintFinalStats(void *cvode_mem);
static int SetInitialProfile(UserData data, N_Vector uu, N_Vector res);

static int check_flag(void *flagvalue, char *funcname, int opt);

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
  int flag, iout;
  long int mu, ml, netf, ncfn;
  realtype rtol, atol, t0, t1, tout, tret;

  int nnz; /* number of non-zeroes  */
  
  cvode_mem = NULL;
  data = NULL;
  uu = constraints = res = NULL;

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew_Serial(NEQ);
  if(check_flag((void *)uu, "N_VNew_Serial", 0)) return(1);
  constraints = N_VNew_Serial(NEQ);
  if(check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
  res = N_VNew_Serial(NEQ);
  if(check_flag((void *)res, "N_VNew_Serial", 0)) return(1);
 
  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2)) return(1);
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

  /* Call CVodeCreate to initialize solution */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /*
  flag = IDASetConstraints(cvode_mem, constraints);
  if(check_flag(&flag, "IDASetConstraints", 1)) return(1);
  */
  N_VDestroy_Serial(constraints);

  flag = CVodeInit(cvode_mem, f, t0, uu);
  if(check_flag(&flag, "CVodeInit", 1)) return(1);

  flag = CVodeSStolerances(cvode_mem, rtol, atol);
  if(check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  /* Call CVKLU and set up the linear solver  */
  nnz = NEQ*NEQ;
  flag = CVKLU(cvode_mem, NEQ, nnz);
  if(check_flag(&flag, "CVKLU", 1)) return(1);
  /* check size of Jacobian matrix  */
  if(MGRID >= 4){
    flag = CVSlsSetSparseJacFn(cvode_mem, jacHeat);
  }
  /* special case MGRID=3  */
  else if(MGRID==3){
    flag = CVSlsSetSparseJacFn(cvode_mem, jacHeat3);
  }
  /* MGRID<=2 is pure boundary points, nothing to solve  */
  else{
    printf("MGRID size is too small to run.\n");
    return(1);
  }
  if(check_flag(&flag, "CVSlsSetSparseJacFn", 1)) return(1);

  /* Print output heading. */
  PrintHeader(rtol, atol);
  
  PrintOutput(cvode_mem, t0, uu);

  /* Loop over output times, call CVode, and print results. */
  
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {
    
    flag = CVode(cvode_mem, tout, uu, &tret, CV_NORMAL);
    if(check_flag(&flag, "CVode", 1)) return(1);

    PrintOutput(cvode_mem, tret, uu);
  
  }
  
  /* Print remaining counters and free memory. */
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);

  PrintFinalStats(cvode_mem);

  /* Free memory  */

  CVodeFree(&cvode_mem);
  N_VDestroy_Serial(uu);
  N_VDestroy_Serial(res);
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
  long int mm, i, j, offset, loc;
  realtype *uv, *resv, coeff;
  UserData data;
  
  uv = N_VGetArrayPointer_Serial(uu); resv = N_VGetArrayPointer_Serial(resval);

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
static int jacHeat3(realtype tt, 
           N_Vector yy, N_Vector fy,
	   SlsMat JacMat, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  yval = N_VGetArrayPointer_Serial(yy);
  realtype dx =  ONE/(MGRID - ONE);
  realtype beta = RCONST(4.0)/(dx*dx);

  SparseSetMatToZero(JacMat); // initialize Jacobian matrix
  
  //set up number of elements in each column 
  
  JacMat -> colptrs[0] = 0;
  JacMat -> colptrs[1] = 1;
 
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

  JacMat -> data[0] = TWO*TWO;
  JacMat -> rowvals[0] = 0;
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
static int jacHeat(realtype tt, 
           N_Vector yy, N_Vector fy, N_Vector resvec,
	   SlsMat JacMat, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  yval = N_VGetArrayPointer_Serial(yy);
  realtype dx =  ONE/(MGRID - ONE);
  realtype beta = RCONST(4.0)/(dx*dx);
  int i,j, repeat=0;

  SparseSetMatToZero(JacMat); /* initialize Jacobian matrix  */

  /* 
   *-----------------------------------------------
   * set up number of elements in each column 
   *-----------------------------------------------
   */
  
  /**** first column block ****/
  JacMat -> colptrs[0] = 0;
  JacMat -> colptrs[1] = 1;
  /* count by twos in the middle  */
  for(i=2;i<MGRID;i++) JacMat -> colptrs[i] = (JacMat -> colptrs[i-1])+2;
  JacMat -> colptrs[MGRID] = 2*MGRID-2;

  /**** second column block ****/
  JacMat -> colptrs[MGRID+1] = 2*MGRID;
  JacMat -> colptrs[MGRID+2] = 2*MGRID+3;
  /* count by fours in the middle */
  for(i=0;i<MGRID-4;i++) JacMat -> colptrs[MGRID+3+i] = (JacMat -> colptrs[MGRID+3+i-1])+4;
  JacMat -> colptrs[2*MGRID-1] = 2*MGRID+4*(MGRID-2)-2;
  JacMat -> colptrs[2*MGRID] = 2*MGRID+4*(MGRID-2);
  
  /**** repeated (MGRID-4 times) middle column blocks ****/
  for(i=0;i<MGRID-4;i++){
    JacMat -> colptrs[2*MGRID+1+repeat]   = (JacMat -> colptrs[2*MGRID+1+repeat-1])+2;
    JacMat -> colptrs[2*MGRID+1+repeat+1] = (JacMat -> colptrs[2*MGRID+1+repeat])+4;
    
    /* count by fives in the middle */
    for(j=0;j<MGRID-4;j++) JacMat -> colptrs[2*MGRID+1+repeat+2+j] = 
			  (JacMat -> colptrs[2*MGRID+1+repeat+1+j])+5;
   
    JacMat -> colptrs[2*MGRID+1+repeat+(MGRID-4)+2] = (JacMat -> colptrs[2*MGRID+1+repeat+(MGRID-4)+1])+4;
    JacMat -> colptrs[2*MGRID+1+repeat+(MGRID-4)+3] = (JacMat -> colptrs[2*MGRID+1+repeat+(MGRID-4)+2])+2;  

    repeat+=MGRID; /* shift that accounts for accumulated number of columns */
  }
  
  /**** last-1 column block ****/
  JacMat -> colptrs[MGRID*MGRID-2*MGRID+1] = TOTAL-2*MGRID-4*(MGRID-2)+2;
  JacMat -> colptrs[MGRID*MGRID-2*MGRID+2] = TOTAL-2*MGRID-4*(MGRID-2)+5;
  /* count by fours in the middle */
  for(i=0;i<MGRID-4;i++) JacMat -> colptrs[MGRID*MGRID-2*MGRID+3+i] = 
			(JacMat -> colptrs[MGRID*MGRID-2*MGRID+3+i-1])+4;
  JacMat -> colptrs[MGRID*MGRID-MGRID-1] = TOTAL-2*MGRID;
  JacMat -> colptrs[MGRID*MGRID-MGRID]   = TOTAL-2*MGRID+2;

  /**** last column block ****/
  JacMat -> colptrs[MGRID*MGRID-MGRID+1] = TOTAL-MGRID-(MGRID-2)+1;
  /* count by twos in the middle */
  for(i=0;i<MGRID-2;i++) JacMat -> colptrs[MGRID*MGRID-MGRID+2+i] = 
			(JacMat -> colptrs[MGRID*MGRID-MGRID+2+i-1])+2;
  JacMat -> colptrs[MGRID*MGRID-1] = TOTAL-1;
  JacMat -> colptrs[MGRID*MGRID]   = TOTAL;
  

  /*
   *-----------------------------------------------
   * set up data stored
   *-----------------------------------------------
   */

  /**** first column block ****/
  JacMat -> data[0] = ONE;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=1;i<MGRID+(MGRID-2)  ;i+=2) JacMat -> data[i] = ONE;
  for(i=2;i<MGRID+(MGRID-2)-1;i+=2) JacMat -> data[i] = -ONE/(dx*dx);

  /**** second column block ****/
  JacMat -> data[MGRID+MGRID-2] = ONE;
  JacMat -> data[MGRID+MGRID-1] = -ONE/(dx*dx);
  JacMat -> data[MGRID+MGRID]   = beta;
  JacMat -> data[MGRID+MGRID+1] = -ONE/(dx*dx);
  JacMat -> data[MGRID+MGRID+2] = -ONE/(dx*dx);
  /* middle data elements */
  for(i=0;i<(MGRID-4);i++) JacMat -> data[MGRID+MGRID+3+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) JacMat -> data[MGRID+MGRID+4+4*i] = beta;
  for(i=0;i<(MGRID-4);i++) JacMat -> data[MGRID+MGRID+5+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) JacMat -> data[MGRID+MGRID+6+4*i] = -ONE/(dx*dx);
  JacMat -> data[2*MGRID+4*(MGRID-2)-5] = -ONE/(dx*dx);
  JacMat -> data[2*MGRID+4*(MGRID-2)-4] = beta;
  JacMat -> data[2*MGRID+4*(MGRID-2)-3] = -ONE/(dx*dx);
  JacMat -> data[2*MGRID+4*(MGRID-2)-2] = -ONE/(dx*dx);
  JacMat -> data[2*MGRID+4*(MGRID-2)-1] = ONE;
    
  /**** repeated (MGRID-4 times) middle column blocks ****/
  repeat=0;
  for(i=0;i<MGRID-4;i++){
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat]   = ONE;
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+1] = -ONE/(dx*dx);
    
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+2] = -ONE/(dx*dx);
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+3] = beta;
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+4] = -ONE/(dx*dx);
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+5] = -ONE/(dx*dx);

    /* 5 in 5*j chosen since there are 5 elements in each column */
    /* this column loops MGRID-4 times within the outer loop */
    for(j=0;j<MGRID-4;j++){
      JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+6+5*j]  = -ONE/(dx*dx);
      JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+7+5*j]  = -ONE/(dx*dx);
      JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+8+5*j]  = beta;
      JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+9+5*j]  = -ONE/(dx*dx);
      JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+10+5*j] = -ONE/(dx*dx);
    }
    
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+6] = -ONE/(dx*dx);
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+7] = -ONE/(dx*dx);
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+8] = beta;
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+9] = -ONE/(dx*dx);
    
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+10] = -ONE/(dx*dx);
    JacMat -> data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+11] = ONE;
    
    repeat+=MGRID+4*(MGRID-2); /* shift that accounts for accumulated columns and elements */
  }
  
  /**** last-1 column block ****/
  JacMat -> data[TOTAL-6*(MGRID-2)-4] = ONE;
  JacMat -> data[TOTAL-6*(MGRID-2)-3] = -ONE/(dx*dx);
  JacMat -> data[TOTAL-6*(MGRID-2)-2] = -ONE/(dx*dx);
  JacMat -> data[TOTAL-6*(MGRID-2)-1] = beta;
  JacMat -> data[TOTAL-6*(MGRID-2)  ] = -ONE/(dx*dx);
  /* middle data elements */
  for(i=0;i<(MGRID-4);i++) JacMat -> data[TOTAL-6*(MGRID-2)+1+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) JacMat -> data[TOTAL-6*(MGRID-2)+2+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) JacMat -> data[TOTAL-6*(MGRID-2)+3+4*i] = beta;
  for(i=0;i<(MGRID-4);i++) JacMat -> data[TOTAL-6*(MGRID-2)+4+4*i] = -ONE/(dx*dx);
  JacMat -> data[TOTAL-2*(MGRID-2)-7] = -ONE/(dx*dx);
  JacMat -> data[TOTAL-2*(MGRID-2)-6] = -ONE/(dx*dx);
  JacMat -> data[TOTAL-2*(MGRID-2)-5] = beta;
  JacMat -> data[TOTAL-2*(MGRID-2)-4] = -ONE/(dx*dx);
  JacMat -> data[TOTAL-2*(MGRID-2)-3] = ONE;

  /**** last column block ****/
  JacMat -> data[TOTAL-2*(MGRID-2)-2] = ONE;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=TOTAL-2*(MGRID-2)-1;i<TOTAL-2;i+=2) JacMat -> data[i] = -ONE/(dx*dx);
  for(i=TOTAL-2*(MGRID-2)  ;i<TOTAL-1;i+=2) JacMat -> data[i] = ONE;
  JacMat -> data[TOTAL-1] = ONE;
  
  /*
   *-----------------------------------------------
   * row values 
   *-----------------------------------------------
   */

  /**** first block ****/
  JacMat -> rowvals[0] = 0;
  /* alternating pattern in data, separate loop for each pattern */
  for(i=1;i<MGRID+(MGRID-2)  ;i+=2) JacMat -> rowvals[i] = (i+1)/2;
  for(i=2;i<MGRID+(MGRID-2)-1;i+=2) JacMat -> rowvals[i] = i/2+MGRID; // i+1 unnecessary here
  
  /**** second column block ****/
  JacMat -> rowvals[MGRID+MGRID-2] = MGRID;
  JacMat -> rowvals[MGRID+MGRID-1] = MGRID+1;
  JacMat -> rowvals[MGRID+MGRID]   = MGRID+1;
  JacMat -> rowvals[MGRID+MGRID+1] = MGRID+2;
  JacMat -> rowvals[MGRID+MGRID+2] = 2*MGRID+1;
  /* middle row values */
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[MGRID+MGRID+3+4*i] = MGRID+1+i;
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[MGRID+MGRID+4+4*i] = MGRID+2+i;
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[MGRID+MGRID+5+4*i] = MGRID+3+i;
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[MGRID+MGRID+6+4*i] = 2*MGRID+2+i;
  JacMat -> rowvals[2*MGRID+4*(MGRID-2)-5] = MGRID+(MGRID-2)-1;
  JacMat -> rowvals[2*MGRID+4*(MGRID-2)-4] = MGRID+(MGRID-2); // starting from here, add two diag patterns
  JacMat -> rowvals[2*MGRID+4*(MGRID-2)-3] = 2*MGRID+(MGRID-2);
  JacMat -> rowvals[2*MGRID+4*(MGRID-2)-2] = MGRID+(MGRID-2);
  JacMat -> rowvals[2*MGRID+4*(MGRID-2)-1] = MGRID+(MGRID-2)+1;
  
  /**** repeated (MGRID-4 times) middle column blocks ****/
  repeat=0;
  for(i=0;i<MGRID-4;i++){
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat]   = MGRID+(MGRID-2)+2+MGRID*i;
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+1] = MGRID+(MGRID-2)+2+MGRID*i+1;
    
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+2] = MGRID+(MGRID-2)+2+MGRID*i+1-MGRID;
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+3] = MGRID+(MGRID-2)+2+MGRID*i+1;
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+4] = MGRID+(MGRID-2)+2+MGRID*i+2; //*this
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+5] = MGRID+(MGRID-2)+2+MGRID*i+1+MGRID;
    
    /* 5 in 5*j chosen since there are 5 elements in each column */
    /* column repeats MGRID-4 times within the outer loop */
    for(j=0;j<MGRID-4;j++){
      JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+6+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+1-MGRID+1+j;
      JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+7+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+1+j;
      JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+8+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+2+j;
      JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+9+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+2+1+j;
      JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+10+5*j] = MGRID+(MGRID-2)+2+MGRID*i+1+MGRID+1+j;
    }
    
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+6] = MGRID+(MGRID-2)+2+MGRID*i-2;
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+7] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID-1;
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+8] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID;//*this+MGRID
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+9] = MGRID+(MGRID-2)+2+MGRID*i-2+2*MGRID;

    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+10] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID;
    JacMat -> rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+11] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID+1;
    
    repeat+=MGRID+4*(MGRID-2); /* shift that accounts for accumulated columns and elements */
  }  

  
  /**** last-1 column block ****/
  JacMat -> rowvals[TOTAL-6*(MGRID-2)-4] = MGRID*MGRID-1-2*(MGRID-1)-1;
  JacMat -> rowvals[TOTAL-6*(MGRID-2)-3] = MGRID*MGRID-1-2*(MGRID-1); // starting with this as base
  JacMat -> rowvals[TOTAL-6*(MGRID-2)-2] = MGRID*MGRID-1-2*(MGRID-1)-MGRID;
  JacMat -> rowvals[TOTAL-6*(MGRID-2)-1] = MGRID*MGRID-1-2*(MGRID-1);
  JacMat -> rowvals[TOTAL-6*(MGRID-2)  ] = MGRID*MGRID-1-2*(MGRID-1)+1;
  /* middle row values */
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[TOTAL-6*(MGRID-2)+1+4*i] = MGRID*MGRID-1-2*(MGRID-1)-MGRID+1+i;
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[TOTAL-6*(MGRID-2)+2+4*i] = MGRID*MGRID-1-2*(MGRID-1)+i;
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[TOTAL-6*(MGRID-2)+3+4*i] = MGRID*MGRID-1-2*(MGRID-1)+1+i;//copied above
  for(i=0;i<(MGRID-4);i++) JacMat -> rowvals[TOTAL-6*(MGRID-2)+4+4*i] = MGRID*MGRID-1-2*(MGRID-1)+2+i;
  JacMat -> rowvals[TOTAL-2*(MGRID-2)-7] = MGRID*MGRID-2*MGRID-2;
  JacMat -> rowvals[TOTAL-2*(MGRID-2)-6] = MGRID*MGRID-MGRID-3;
  JacMat -> rowvals[TOTAL-2*(MGRID-2)-5] = MGRID*MGRID-MGRID-2;
  JacMat -> rowvals[TOTAL-2*(MGRID-2)-4] = MGRID*MGRID-MGRID-2;
  JacMat -> rowvals[TOTAL-2*(MGRID-2)-3] = MGRID*MGRID-MGRID-1;

  /* last column block */
  JacMat -> rowvals[TOTAL-2*(MGRID-2)-2] = MGRID*MGRID-MGRID;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=0;i<(MGRID-2);i++) JacMat -> rowvals[TOTAL-2*(MGRID-2)-1+2*i] = MGRID*MGRID-2*MGRID+1+i;
  for(i=0;i<(MGRID-2);i++) JacMat -> rowvals[TOTAL-2*(MGRID-2)  +2*i] = MGRID*MGRID-MGRID+1+i;  
  JacMat -> rowvals[TOTAL-1] = MGRID*MGRID-1;

  //SparsePrintMat(JacMat);
 
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
  long int mm, mm1, i, j, offset, loc;
  
  mm = data->mm;
  mm1 = mm - 1;
  
  udata = N_VGetArrayPointer_Serial(uu);

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
  //  printf("Constraints set to force all solution components >= 0. \n");
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
  int flag;
  realtype umax, hused;
  long int nst, nni, nje, nre, nreLS;
  int kused;

  umax = N_VMaxNorm(uu);
  
  flag = CVodeGetLastOrder(mem, &kused);
  check_flag(&flag, "IDAGetLastOrder", 1);
  flag = CVodeGetNumSteps(mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = CVodeGetNumNonlinSolvIters(mem, &nni);
  check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);
  //flag = CVodeGetNumResEvals(mem, &nre);
  //check_flag(&flag, "IDAGetNumResEvals", 1);
  flag = CVodeGetLastStep(mem, &hused);
  check_flag(&flag, "IDAGetLastStep", 1);
  flag = CVodeSlsGetNumJacEvals(mem, &nje);
  check_flag(&flag, "IDASlsGetNumJacEvals", 1);

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
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVSlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVSlsGetNumJacEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

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
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
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
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
