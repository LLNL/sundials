#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"
#include "cvodes.h"
#include "iterativ.h"
#include "cvsspgmr.h"
#include "cvsbandpre.h"
#include "cvsband.h"
#include "band.h"
#include "nvector.h"
#include "llnlmath.h"

/* Problem Constants */

#define T0           0.0          /* initial time */
#define NOUT         10           /* number of output times */
#define TOUT         0.5          /* time between two outputs */

#define XMIN         0.0          /* grid boundaries in x */
#define XMAX         1.0           
#define YMIN         0.0          /* grid boundaries in y */
#define YMAX         1.0

#define M            40           /* M = number of interior points */
#define NEQ          M*M          /* NEQ = number of equations */

#define NP           2

/* CVodeMalloc Constants */

#define RTOL    1.0e-4            /* scalar relative tolerance */
#define ATOL    1.0e-4            /* scalar absolute tolerance */

#define BAND  0
#define SPGMR 1

/* Sensitivity Constants */

#define NS    2

#define ZERO  RCONST(0.0)

/* Type : UserData */ 
typedef struct {
  real *p;
} *UserData;

/* User-defined matrix accessor macro: IJth */
#define IJth(vdata,i,j) ( vdata[ (i-1)*M + (j-1) ] )

/* Private Helper Functions */
static void SetInitialProfiles(N_Vector u);
static void PrintOutput(long int iopt[], real ropt[], real t, 
                        N_Vector y, N_Vector *uS);
static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[], int lin_solv);

/* Functions Called by the CVODE Solver */
static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);

/***************************** Main Program ******************************/

main(int argc, char *argv[])
{
  real abstol, reltol, t, tout, ropt[OPT_SIZE];
  long int iopt[OPT_SIZE];
  N_Vector u;
  UserData data;
  CVBandPreData bpdata;
  void *cvode_mem;
  int ml, mu, iout, flag, i;
  machEnvType machEnv;

  real *pbar, rhomax;
  integer is, *plist;
  N_Vector *uS;

  int lin_solv, sensi, sensi_meth, err_con;

  /* MACHINE ENVIRONMENT */
  machEnv = NULL;

  /* PROBLEM PARAMETERS */
  data = (UserData)malloc(sizeof *data);
  data->p = (real *)malloc(NS*sizeof(real));
  data->p[0] = 1.0;
  data->p[1] = 1.0;

  /* INITIAL STATES */
  u = N_VNew(NEQ, machEnv);
  SetInitialProfiles(u);
  
  /* TOLERANCES */
  abstol=ATOL; 
  reltol=RTOL;

  /* OPTIONAL INPUT/OUTPUT */
  for (i = 0; i < OPT_SIZE; i++) {
    iopt[i] = 0;
    ropt[i] = 0.0;
  }

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(NEQ, f, T0, u, BDF, NEWTON, SS, &reltol,
                          &abstol, data, NULL, TRUE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed."); return(1); }

  /* LINEAR SOLVER */
  printf("\nLinear solver type (%d:BAND , %d:SPGMR): ",BAND,SPGMR);
  scanf("%d",&lin_solv);
  if(lin_solv == BAND) {
    /* CVBAND */
    flag = CVBand(cvode_mem, M, M, NULL, NULL);
    if (flag != SUCCESS) { printf("CVBand failed.\n"); return(1); }
  } else if(lin_solv == SPGMR) {
    /* INITIALIZE PRECONDITIONER */
    ml = mu = M;
    bpdata = CVBandPreAlloc (NEQ, f, data, mu, ml);
    /* CVSPGMR */
    flag = CVSpgmr(cvode_mem, LEFT, MODIFIED_GS, 0, 0.0, CVBandPrecond, CVBandPSolve, 
                   bpdata, NULL, NULL);
    if (flag != SUCCESS) { printf("CVSpgmr failed.\n"); return(1); }
  } else {
    printf("\n\n WRONG linear solver type \n\n");
    exit(0);
  }

  /* SENSITIVTY */
  printf("\nPerform sensitivity analysis? (0:NO , 1:YES): ");scanf("%d",&sensi);
  if(sensi) {
    pbar = (real *) malloc(NP*sizeof(real));
    for(is=0; is<NP; is++) pbar[is] = data->p[is];
    plist = (integer *) malloc(NS * sizeof(integer));
    for(is=0; is<NS; is++) plist[is] = is+1;

    uS = N_VNew_S(NS,NEQ,machEnv);
    for(is=0;is<NS;is++)
      N_VConst(ZERO,uS[is]);

    rhomax = ZERO;
   
    printf("\nSensitivity method (%d:SIMULTANEOUS , %d:STAGGERED): ",SIMULTANEOUS,STAGGERED);
    scanf("%d",&sensi_meth);
    printf("\nError control (%d:FULL , %d:PARTIAL): ",FULL,PARTIAL);
    scanf("%d",&err_con);

    flag = CVodeSensMalloc(cvode_mem,NS,sensi_meth,data->p,pbar,plist,
                           ALLSENS,NULL,err_con,rhomax,uS,NULL,NULL);
    if (flag != SUCCESS) {printf("CVodeSensMalloc failed, flag=%d\n",flag);return(1);}

  }

  /* In loop over output points, call CVode, print results, test for error */
  printf("\n  NST      T        H      Q   NCF  NCFS   NEF  NEFS\n");
  printf(  "----------------------------------------------------\n");
  for (iout=0,tout=TOUT; iout < NOUT; iout++,tout+=TOUT) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    if (flag != SUCCESS) { printf("CVode failed, flag=%d.\n", flag); break; }

    printf("%5ld %8.3e %8.3e %2d %5ld %5ld %5ld %5ld\n",
           iopt[NST],t,ropt[HU],iopt[QU],iopt[NCFN],iopt[NCFNS],iopt[NETF],iopt[NETFS]);
  }

  /* Free memory and print final statistics */  

  N_VFree(u);
  if(sensi) N_VFree_S(NS, uS);
  CVodeFree(cvode_mem);

  PrintFinalStats(sensi,sensi_meth,err_con,iopt,lin_solv);
  return(0);
}


/*********************** Private Helper Functions ************************/

/* ======================================================================= */
/* Set initial conditions in y */

static void SetInitialProfiles(N_Vector u)
{
  int i, j;
  real del, x, y, *udata;

  /* Set pointer to data array in vector y. */

  udata = N_VDATA(u);

  /* Load initial profiles of u into y vector */

  del = 1.0/(M+1);
  for (i=1; i<=M; i++) {
    x = i*del;
    for (j=1; j<=M; j++) {
      y = j*del;
      IJth(udata,i,j) = 16.0 * x * (1.0-x) * y * (1.0-y);
    }
  }

}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(long int iopt[], real ropt[], real t, 
                        N_Vector u, N_Vector *uS)
{
  real *udata, *uSdata;

  udata = N_VDATA(u);

  printf("t = %.2e q=%d h=%.3e nst=%d\n\n", t,iopt[QU],ropt[HU],iopt[NST]);
  printf("=================================================\n");

}

/* ======================================================================= */
/* Print final statistics contained in iopt */

static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[], int lin_solv)
{
  
  printf("\n\n========================================================");
  printf("\nFinal Statistics");
  printf("\nSensitivity: ");

  if(sensi) {
    printf("YES ");
    if(sensi_meth == SIMULTANEOUS) printf("( SIMULTANEOUS +");
    else                           printf("( STAGGERED +");
    if(err_con == FULL) printf(" FULL ERROR CONTROL )");
    else                printf(" PARTIAL ERROR CONTROL )");
  } else {
    printf("NO");
  }

  printf("\n\n");
  printf("nst     = %5ld                \n\n", iopt[NST]);
  printf("nfe     = %5ld    nfSe  = %5ld  \n", iopt[NFE],  iopt[NFSE]);
  printf("nni     = %5ld    nniS  = %5ld  \n", iopt[NNI],  iopt[NNIS]);
  printf("ncfn    = %5ld    ncfnS = %5ld  \n", iopt[NCFN], iopt[NCFNS]);
  printf("netf    = %5ld    netfS = %5ld\n\n", iopt[NETF], iopt[NETFS]);
  printf("nsetups = %5ld                  \n", iopt[NSETUPS]);
  if(lin_solv == SPGMR) {
    printf("nli     = %5ld    ncfl  = %5ld  \n", iopt[SPGMR_NLI], iopt[SPGMR_NCFL]);
    printf("npe     = %5ld    nps   = %5ld  \n", iopt[SPGMR_NPE], iopt[SPGMR_NPS]);
  } else {
    printf("nje     = %5ld                  \n", iopt[BAND_NJE]);
  }

  printf("========================================================\n");

}


/***************** Functions Called by the CVODE Solver ******************/

/* ======================================================================= */
/* f routine. Compute f(t,y). */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data)
{
  real *udata, *udotdata;
  real px, py, del, rdel2, uxx, uyy;
  int i, j;
  UserData data;

  data     = (UserData) f_data;
  udata    = N_VDATA(u);
  udotdata = N_VDATA(udot);

  px = data->p[0];
  py = data->p[1];

  del  = 1.0/(M+1);
  rdel2 = 1.0/(del*del);

  /* Internal points */
  for(i=2; i<M; i++) {
    for(j=2; j<M; j++) {
      uxx = IJth(udata,i-1,j) - 2.0*IJth(udata,i,j) + IJth(udata,i+1,j);
      uyy = IJth(udata,i,j-1) - 2.0*IJth(udata,i,j) + IJth(udata,i,j+1);
      IJth(udotdata,i,j) = rdel2*(px*uxx+py*uyy);
    }
  }

  /* Top boundary (i=1) */
  for(j=2; j<M; j++) {
    uxx =                   - 2.0*IJth(udata,1,j) + IJth(udata,2,j);
    uyy = IJth(udata,1,j-1) - 2.0*IJth(udata,1,j) + IJth(udata,1,j+1);
    IJth(udotdata,1,j) = rdel2*(px*uxx+py*uyy);
  }

  /* Bottom boundary (i=M) */
  for(j=2; j<M; j++) {
    uxx = IJth(udata,M-1,j) - 2.0*IJth(udata,M,j);
    uyy = IJth(udata,M,j-1) - 2.0*IJth(udata,M,j) + IJth(udata,M,j+1);
    IJth(udotdata,M,j) = rdel2*(px*uxx+py*uyy);
  }

  /* Left boundary (j=1) */
  for(i=2; i<M; i++) {
    uxx = IJth(udata,i-1,1) - 2.0*IJth(udata,i,1) + IJth(udata,i+1,1);
    uyy =                   - 2.0*IJth(udata,i,1) + IJth(udata,i,2);
    IJth(udotdata,i,1) = rdel2*(px*uxx+py*uyy);
  }

  /* Right boundary (j=M) */
  for(i=2; i<M; i++) {
    uxx = IJth(udata,i-1,M) - 2.0*IJth(udata,i,M) + IJth(udata,i+1,M);
    uyy = IJth(udata,i,M-1) - 2.0*IJth(udata,i,M);
    IJth(udotdata,i,M) = rdel2*(px*uxx+py*uyy);
  }

  /* Top-left corner (i=1, j=1) */
  uxx =                   - 2.0*IJth(udata,1,1) + IJth(udata,2,1);
  uyy =                   - 2.0*IJth(udata,1,1) + IJth(udata,1,2);
  IJth(udotdata,1,1) = rdel2*(px*uxx+py*uyy);

  /* Top-right corner (i=1, j=M) */
  uxx =                   - 2.0*IJth(udata,1,M) + IJth(udata,2,M);
  uyy = IJth(udata,1,M-1) - 2.0*IJth(udata,1,M);
  IJth(udotdata,1,M) = rdel2*(px*uxx+py*uyy);

  /* Bottom-left corner (i=M, j=1) */
  uxx = IJth(udata,M-1,1) - 2.0*IJth(udata,M,1);
  uyy =                   - 2.0*IJth(udata,M,1) + IJth(udata,M,2);
  IJth(udotdata,M,1) = rdel2*(px*uxx+py*uyy);

  /* Bottom-right corner (i=M, j=M) */
  uxx = IJth(udata,M-1,M) - 2.0*IJth(udata,M,M);
  uyy = IJth(udata,M,M-1) - 2.0*IJth(udata,M,M);
  IJth(udotdata,M,M) = rdel2*(px*uxx+py*uyy);

}

 
