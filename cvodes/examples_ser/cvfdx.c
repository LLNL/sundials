/************************************************************************
 *                                                                      *
 * File       : cvfdx.c                                                 *
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and Radu Serban      * 
 *              @ LLNL                                                  *
 * Version of : 11 February 2004                                        *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem, with the coding           *
 * needed for its solution by CVODES.  The problem is from chemical     *
 * kinetics, and consists of the following three rate equations..       *
 *    dy1/dt = -p1*y1 + p2*y2*y3                                        *
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2                            *
 *    dy3/dt =  p3*(y2)^2                                               *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions   *
 * y1 = 1.0, y2 = y3 = 0.  The reaction rates are: p1=0.04, p2=1e4, and *
 * p3=3e7.  The problem is stiff.                                       *
 * This program solves the problem with the BDF method, Newton          *
 * iteration with the CVODES dense linear solver, and a user-supplied   *
 * Jacobian routine.                                                    * 
 * It uses a scalar relative tolerance and a vector absolute tolerance. *
 * Output is printed in decades from t = .4 to t = 4.e10.               *
 * Run statistics (optional outputs) are printed at the end.            *
 *                                                                      *
 * Optionally, CVODES can compute sensitivities with respect to the     *
 * problem parameters p1, p2, and p3.                                   *
 * The sensitivity right hand side is given analytically through the    *
 * user routine fS (of type SensRhs1Fn).                                *
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and       *
 * STAGGERED1) can be used and sensitivities may be included in the     *
 * error test or not (error control set on FULL or PARTIAL,             *
 * respectively).                                                       *
 *                                                                      *
 * Execution:                                                           *
 *                                                                      *
 * If no sensitivities are desired:                                     *
 *    % cvsdx -nosensi                                                  *
 * If sensitivities are to be computed:                                 *
 *    % cvsdx -sensi sensi_meth err_con                                 *
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of    *
 * {full, partial}.                                                     * 
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sundialstypes.h"   /* definition of types realtype (set to double)      */
                             /* and the constant FALSE                            */
#include "cvodes.h"          /* prototypes for CVodeMalloc, CVode, and CVodeFree, */
                             /* constants OPT_SIZE, BDF, NEWTON, SV, SUCCESS,     */
                             /* NST, NFE, NSETUPS, NNI, NCFN, NETF                */
#include "cvsdense.h"        /* prototype for CVDense, constant DENSE_NJE         */
#include "nvector_serial.h"  /* definitions of type N_Vector and macro NV_Ith_S,  */
                             /* prototypes for N_VNew, N_VFree                    */
#include "dense.h"           /* definitions of type DenseMat, macro DENSE_ELEM    */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   3            /* number of equations  */
#define Y1    1.0          /* initial y components */
#define Y2    0.0
#define Y3    0.0
#define RTOL  1e-4         /* scalar relative tolerance            */
#define ATOL1 1e-8         /* vector absolute tolerance components */
#define ATOL2 1e-14
#define ATOL3 1e-6
#define T0    0.0          /* initial time           */
#define T1    0.4          /* first output time      */
#define TMULT 10.0         /* output time factor     */
#define NOUT  12           /* number of output times */

#define NP    3
#define NS    3

#define ZERO  0.0

/* Type : UserData */

typedef struct {
  realtype p[3];
} *UserData;

/* Private Helper Function */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth, int *err_con);
static void WrongArgs(char *name);
static void PrintFinalStats(void *cvode_mem, booleantype sensi, int sensi_meth, int err_con);
static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS);


/* Functions Called by the CVODES Solver */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
               int iS, N_Vector yS, N_Vector ySdot, 
               void *fS_data, N_Vector tmp1, N_Vector tmp2);


/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  NV_Spec nvSpec;
  void *cvode_mem;
  UserData data;
  realtype reltol, t, tout;
  N_Vector y, abstol;
  int iout, flag;

  realtype pbar[NP];
  int is, *plist; 
  N_Vector *yS;
  booleantype sensi;
  int sensi_meth, err_con;

  nvSpec = NULL;
  cvode_mem = NULL;
  data = NULL;
  y = abstol = NULL;
  yS = NULL;
  plist = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* Initialize serial vector specification */
  nvSpec = NV_SpecInit_Serial(NEQ);
  if (check_flag((void *)nvSpec, "NV_SpecInit", 0)) return(1);

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2)) return(1);
  data->p[0] = 0.04;
  data->p[1] = 1.0e4;
  data->p[2] = 3.0e7;

  /* INITIAL STATES */
  y = N_VNew(nvSpec);
  if (check_flag((void *)y, "N_VNew", 0)) return(1);
  abstol = N_VNew(nvSpec);
  if (check_flag((void *)abstol, "N_VNew", 0)) return(1);

  /* Initialize y */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  /* TOLERANCES */
  /* Set the scalar relative tolerance */
  reltol = RTOL;               

  /* Set the vector absolute tolerance */
  Ith(abstol,1) = ATOL1;       
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* CVODE_CREATE */
  cvode_mem = CVodeCreate(BDF, NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeSetFdata(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetFdata", 1)) return(1);

  /* CVODE_MALLOC */
  flag = CVodeMalloc(cvode_mem, f, T0, y, SV, &reltol, abstol, nvSpec);
  if (check_flag(&flag, "CVodeMalloc", 1)) return(1);

  /* CVDENSE */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  flag = CVDenseSetJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDenseSetJacFn", 1)) return(1);

  flag = CVDenseSetJacData(cvode_mem, data);
  if (check_flag(&flag, "CVDenseSetJacData", 1)) return(1);

  printf("\n3-species chemical kinetics problem\n");

  /* SENSITIVITY */
  if (sensi) {

    pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];
    plist = (int *) malloc(NS * sizeof(int));
    if (check_flag((void *)plist, "malloc", 2)) return(1);
    for (is=0; is<NS; is++) plist[is] = is+1;

    yS = N_VNew_S(NS, nvSpec);
    if (check_flag((void *)yS, "N_VNew", 0)) return(1);
    for (is=0;is<NS;is++)
      N_VConst(0.0, yS[is]);

    flag = CVodeSetSensRhs1Fn(cvode_mem, fS);
    if (check_flag(&flag, "CVodeSetSensRhs1Fn", 1)) return(1);
    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    if (check_flag(&flag, "CVodeSetSensFdata", 1)) return(1);
    flag = CVodeSetSensFdata(cvode_mem, data);
    if (check_flag(&flag, "CVodeSetSensFdata", 1)) return(1);
    flag = CVodeSetSensPbar(cvode_mem, pbar);
    if (check_flag(&flag, "CVodeSetSensPbar", 1)) return(1);

    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, plist, yS);
    if(check_flag(&flag, "CVodeSensMalloc", 1)) return(1);

    printf("Sensitivity: YES ");
    if(sensi_meth == SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      if(sensi_meth == STAGGERED) printf("( STAGGERED +");
      else                        printf("( STAGGERED1 +");   
    if(err_con == FULL) printf(" FULL ERROR CONTROL )");
    else                printf(" PARTIAL ERROR CONTROL )");

  } else {

    printf("Sensitivity: NO ");

  }
  
  /* In loop over output points, call CVode, print results, test for error */
  
  printf("\n\n");
  printf("===================================================");
  printf("==================================\n");
  printf("     T     Q       H      NST                    y1");
  printf("           y2           y3    \n");
  printf("===================================================");
  printf("==================================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    if (check_flag(&flag, "CVode", 1)) break;
    PrintOutput(cvode_mem, t, y);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, yS);
      if (check_flag(&flag, "CVodeGetSens", 1)) break;
      PrintOutputS(yS);
    } 
    printf("-------------------------------------------------");
    printf("------------------------------------\n");
  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi,sensi_meth,err_con);

  /* Free memory */

  /* Free y vector */
  N_VFree(y);

  /* Free abstol vector */
  N_VFree(abstol);   

  /* Free yS vector */
  if (sensi) N_VFree_S(NS, yS);

  /* Free user data */
  free(data);

  /* Free CVODES problem memory */
  CVodeFree(cvode_mem);

  /* Free vector specification */
  NV_SpecFree_Serial(nvSpec);

  /* Free plist */
  if (sensi) free(plist);

  return(0);
}


/************************ Private Helper Function ************************/

/* ======================================================================= */
/* Exit if arguments are incorrect */

static void ProcessArgs(int argc, char *argv[], 
                        booleantype *sensi, int *sensi_meth, int *err_con)
{
  *sensi = FALSE;
  *sensi_meth = -1;
  *err_con = -1;

  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = TRUE;
  else
    WrongArgs(argv[0]);
  
  if (*sensi) {

    if (argc != 4)
      WrongArgs(argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      *sensi_meth = STAGGERED1;
    else 
      WrongArgs(argv[0]);

    if (strcmp(argv[3],"full") == 0)
      *err_con = FULL;
    else if (strcmp(argv[3],"partial") == 0)
      *err_con = PARTIAL;
    else
      WrongArgs(argv[0]);
  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = full or partial\n");
    
    exit(0);
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and solution  */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, flag;
  realtype hu, *udata;
  
  udata = NV_DATA_S(u);

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetLastOrder(cvode_mem, &qu);
  check_flag(&flag, "CVodeGetLastOrder", 1);
  flag = CVodeGetLastStep(cvode_mem, &hu);
  check_flag(&flag, "CVodeGetLastStep", 1);

  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
  printf("                                Solution       ");
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
}

/* ======================================================================= */
/* Print sensitivities */

static void PrintOutputS(N_Vector *uS)
{
  realtype *sdata;

  sdata = NV_DATA_S(uS[0]);
  printf("                                Sensitivity 1  ");
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
  
  sdata = NV_DATA_S(uS[1]);
  printf("                                Sensitivity 2  ");
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);

  sdata = NV_DATA_S(uS[2]);
  printf("                                Sensitivity 3  ");
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
}

/* ======================================================================= */
/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem, booleantype sensi, int sensi_meth, int err_con)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int njeD, nfeD;
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

  if (sensi) {
    flag = CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
    check_flag(&flag, "CVodeGetNumSensRhsEvals", 1);
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1);
    flag = CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
    check_flag(&flag, "CVodeGetNumSensLinSolvSetups", 1);
    flag = CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
    check_flag(&flag, "CVodeGetNumSensErrTestFails", 1);
    flag = CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
    check_flag(&flag, "CVodeGetNumSensNonlinSolvIters", 1);
    flag = CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);
    check_flag(&flag, "CVodeGetNumSensNonlinSolvConvFails", 1);
  }

  flag = CVDenseGetNumJacEvals(cvode_mem, &njeD);
  check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeD);
  check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

  printf("\n");
  printf("njeD    = %5ld    nfeD     = %5ld\n", njeD, nfeD);

}


/***************** Functions Called by the CVODES Solver ******************/

/* ======================================================================= */
/* f routine. Compute f(t,y). */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
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
}

/* ======================================================================= */
/* Jacobian routine. Compute J(t,y). */

static void Jac(long int N, DenseMat J, realtype t,
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
}
 
/* ======================================================================= */
/* fS routine. Compute sensitivity r.h.s. */

static void fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
               int iS, N_Vector yS, N_Vector ySdot, 
               void *fS_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2, sd3;

  data = (UserData) fS_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
  s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);

  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
  sd3 = 2*p3*y2*s2;
  sd2 = -sd1-sd3;

  switch (iS) {
  case 0:
    sd1 += -y1;
    sd2 +=  y1;
    break;
  case 1:
    sd1 +=  y2*y3;
    sd2 += -y2*y3;
    break;
  case 2:
    sd2 += -y2*y2;
    sd3 +=  y2*y2;
    break;
  }
  
  Ith(ySdot,1) = sd1;
  Ith(ySdot,2) = sd2;
  Ith(ySdot,3) = sd3;
}


/************************ Private Helper Function ************************/

/* ======================================================================= */
/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag == SUCCESS
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag != SUCCESS */
  else if (opt == 1) {
    errflag = flagvalue;
    if (*errflag != SUCCESS) {
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
