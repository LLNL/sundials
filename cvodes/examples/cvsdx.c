/************************************************************************
 *                                                                      *
 * File       : cvdx.c                                                  *
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and Radu Serban      * 
 *              @ LLNL                                                  *
 * Version of : 20 March 2002                                           *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem, with the coding           *
 * needed for its solution by CVODE.  The problem is from chemical      *
 * kinetics, and consists of the following three rate equations..       *
 *    dy1/dt = -p1*y1 + p2*y2*y3                                        *
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2                            *
 *    dy3/dt =  p3*(y2)^2                                               *
 * on the interval from t = 0.0 to t = 4.e10, with initial conditions   *
 * y1 = 1.0, y2 = y3 = 0.  The reaction rates are: p1=0.04, p2=1e4, and *
 * p3=3e7.  The problem is stiff.                                       *
 * This program solves the problem with the BDF method, Newton          *
 * iteration with the CVODE dense linear solver, and a user-supplied    *
 * Jacobian routine.                                                    * 
 * It uses a scalar relative tolerance and a vector absolute tolerance. *
 * Output is printed in decades from t = .4 to t = 4.e10.               *
 * Run statistics (optional outputs) are printed at the end.            *
 *                                                                      *
 * Optionally, CVODES can compute sensitivities with respect to the     *
 * problem parameters p1, p2, and p3.                                   *
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
#include "llnltyps.h"        /* definitions of types real (set to double) and     */
                             /* integer (set to int), and the constant FALSE      */
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
  real p[3];
} *UserData;

/* Private Helper Function */

static void WrongArgs(char *argv[]);
static void PrintFinalStats(boole sensi, int sensi_meth, int err_con, long int iopt[]);
static void PrintOutput(long int iopt[], real ropt[], real t, N_Vector u);
static void PrintOutputS(N_Vector *uS);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);
static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  M_Env machEnv;
  UserData data;
  real ropt[OPT_SIZE], reltol, t, tout;
  long int iopt[OPT_SIZE];
  N_Vector y, abstol;
  void *cvode_mem;
  int iout, flag;

  real pbar[NP], rhomax;
  integer is, *plist; 
  N_Vector *yS;
  boole sensi;
  int sensi_meth, err_con, ifS;

  /* Process arguments */
  if (argc < 2)
    WrongArgs(argv);

  if (strcmp(argv[1],"-nosensi") == 0)
    sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    sensi = TRUE;
  else
    WrongArgs(argv);

  if (sensi) {

    if (argc != 4)
      WrongArgs(argv);

    if (strcmp(argv[2],"sim") == 0)
      sensi_meth = SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      sensi_meth = STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      sensi_meth = STAGGERED1;
    else 
      WrongArgs(argv);

    if (strcmp(argv[3],"full") == 0)
      err_con = FULL;
    else if (strcmp(argv[3],"partial") == 0)
      err_con = PARTIAL;
    else
      WrongArgs(argv);

  }

  /* Initialize serial machine environment */
  machEnv = M_EnvInit_Serial(NEQ);

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = 0.04;
  data->p[1] = 1.0e4;
  data->p[2] = 3.0e7;

  /* INITIAL STATES */
  y = N_VNew(NEQ, machEnv);
  abstol = N_VNew(NEQ, machEnv);

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

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(NEQ, f, T0, y, BDF, NEWTON, SV, &reltol, abstol,
                          data, NULL, FALSE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { 
    printf("CVodeMalloc failed.\n"); 
    return(1); 
  }

  /* CVDENSE */
  flag = CVDense(cvode_mem, Jac, NULL);
  if (flag != SUCCESS) { printf("CVDense failed.\n"); return(1); }

  /* SENSITIVITY */
  if(sensi) {
    pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];
    plist = (integer *) malloc(NS * sizeof(integer));
    for(is=0;is<NS;is++) plist[is] = is+1;

    yS = N_VNew_S(NS, NEQ, machEnv);
    for(is=0;is<NS;is++)
      N_VConst(0.0, yS[is]);

    ifS = ALLSENS;
    if(sensi_meth==STAGGERED1) ifS = ONESENS;

    rhomax = ZERO;
    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, pbar, plist,
                           ifS, NULL, err_con, rhomax, yS, NULL, NULL);
    if (flag != SUCCESS) {
      printf("CVodeSensMalloc failed, flag=%d\n",flag);
      return(1);
    }
  }
  
  /* In loop over output points, call CVode, print results, test for error */

  printf("\n3-species chemical kinetics problem\n\n");
  printf("=====================================================================================\n");
  printf("     T     Q       H      NST                    y1           y2           y3    \n");
  printf("=====================================================================================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    if (flag != SUCCESS) {
      printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }
    PrintOutput(iopt, ropt, t, y);
    if (sensi) {
      flag = CVodeSensExtract(cvode_mem, t, yS);
      if (flag != SUCCESS) { 
        printf("CVodeSensExtract failed, flag=%d.\n", flag); 
        break; 
      }
      PrintOutputS(yS);
    } 
    printf("-------------------------------------------------------------------------------------\n");
  }

  /* Print final statistics */
  PrintFinalStats(sensi,sensi_meth,err_con,iopt);

  /* Free memory */
  N_VFree(y);                  /* Free the y and abstol vectors       */
  N_VFree(abstol);   
  if(sensi) N_VFree_S(NS, yS); /* Free the yS vectors                 */
  free(data);                  /* Free user data                      */
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory       */
  M_EnvFree_Serial(machEnv);   /* Free the machine environment memory */

  return(0);
}


/************************ Private Helper Function ************************/

/* ======================================================================= */
/* Exit if arguments are incorrect */

static void WrongArgs(char *argv[])
{
  printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",argv[0]);
  printf("         sensi_meth = sim, stg, or stg1\n");
  printf("         err_con    = full or partial\n");
  
  exit(0);
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and solution  */

static void PrintOutput(long int iopt[], real ropt[], real t, N_Vector u)
{

  real *udata;
  
  udata = NV_DATA_S(u);

  printf("%8.3e %2ld  %8.3e %5ld\n", t,iopt[QU],ropt[HU],iopt[NST]);
  printf("                                Solution       ");
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
  
}
/* ======================================================================= */
/* Print sensitivities */

static void PrintOutputS(N_Vector *uS)
{

  real *sdata;

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

static void PrintFinalStats(boole sensi, int sensi_meth, int err_con, long int iopt[])
{

  printf("\n\n========================================================");
  printf("\nFinal Statistics");
  printf("\nSensitivity: ");

  if(sensi) {
    printf("YES ");
    if(sensi_meth == SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      if(sensi_meth == STAGGERED) printf("( STAGGERED +");
      else                        printf("( STAGGERED1 +");   
    if(err_con == FULL) printf(" FULL ERROR CONTROL )");
    else                printf(" PARTIAL ERROR CONTROL )");
  } else {
    printf("NO");
  }

  printf("\n\n");
  /*
  printf("lenrw   = %5ld    leniw = %5ld\n", iopt[LENRW], iopt[LENIW]);
  printf("llrw    = %5ld    lliw  = %5ld\n", iopt[SPGMR_LRW], iopt[SPGMR_LIW]);
  */
  printf("nst     = %5ld                \n\n", iopt[NST]);
  printf("nfe     = %5ld    nfSe  = %5ld  \n", iopt[NFE],  iopt[NFSE]);
  printf("nni     = %5ld    nniS  = %5ld  \n", iopt[NNI],  iopt[NNIS]);
  printf("ncfn    = %5ld    ncfnS = %5ld  \n", iopt[NCFN], iopt[NCFNS]);
  printf("netf    = %5ld    netfS = %5ld\n\n", iopt[NETF], iopt[NETFS]);
  printf("nsetups = %5ld                  \n", iopt[NSETUPS]);
  printf("nje     = %5ld                  \n", iopt[DENSE_NJE]);

  printf("========================================================\n");

}


/***************** Functions Called by the CVODE Solver ******************/

/* ======================================================================= */
/* f routine. Compute f(t,y). */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)
{
  real y1, y2, y3, yd1, yd3;
  UserData data;
  real p1, p2, p3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) f_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  yd1 = Ith(ydot,1) = -p1*y1 + p2*y2*y3;
  yd3 = Ith(ydot,3) = p3*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;
}

/* ======================================================================= */
/* Jacobian routine. Compute J(t,y). */

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  real y1, y2, y3;
  UserData data;
  real p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) f_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
                      IJth(J,3,2) = 2*p3*y2;
}
 
