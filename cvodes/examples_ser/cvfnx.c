/************************************************************************
 *                                                                      *
 * File       : cvfnx.c                                                 *
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne, and *
 *              Radu Serban @ LLNL                                      *
 * Version of : 11 July 2003                                            *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem, with the program for its  *
 * solution by CVODES.  The problem is the semi-discrete form of the    *
 * advection-diffusion equation in 1-D:                                 *
 *   du/dt = q1 * d^2 u / dx^2 + q2 * du/dx                             *
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.      *
 * Homogeneous Dirichlet boundary conditions are posed, and the         *
 * initial condition is                                                 *
 *   u(x,y,t=0) = x(2-x)exp(2x) .                                       *
 * The PDE is discretized on a uniform grid of size MX+2 with           *
 * central differencing, and with boundary values eliminated,           *
 * leaving an ODE system of size NEQ = MX.                              *
 * This program solves the problem with the option for nonstiff systems:*
 * ADAMS method and functional iteration.                               *
 * It uses scalar relative and absolute tolerances.                     *
 * Output is printed at t = .5, 1.0, ..., 5.                            *
 * Run statistics (optional outputs) are printed at the end.            *
 *                                                                      *
 * Optionally, CVODES can compute sensitivities with respect to the     *
 * problem parameters q1 and q2.                                        *
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and       *
 * STAGGERED1) can be used and sensitivities may be included in the     *
 * error test or not (error control set on FULL or PARTIAL,             *
 * respectively).                                                       *
 *                                                                      *
 * Execution:                                                           *
 *                                                                      *
 * If no sensitivities are desired:                                     *
 *    % cvsnx -nosensi                                                  *
 * If sensitivities are to be computed:                                 *
 *    % cvsnx -sensi sensi_meth err_con                                 *
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of    *
 * {full, partial}.                                                     * 
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sundialstypes.h"
#include "cvodes.h"
#include "nvector_serial.h"

/* Problem Constants */
#define XMAX  2.0          /* domain boundary           */
#define MX    10           /* mesh dimension            */
#define NEQ   MX           /* number of equations       */
#define ATOL  1.e-5        /* scalar absolute tolerance */
#define T0    0.0          /* initial time              */
#define T1    0.5          /* first output time         */
#define DTOUT 0.5          /* output time increment     */
#define NOUT  10           /* number of output times    */

#define NP    2
#define NS    2

#define ZERO  RCONST(0.0)

/* Type : UserData 
   contains problem parameters, grid constants, work array. */

typedef struct {
  realtype *p;
  realtype dx, hdcoef, hacoef;
} *UserData;


/* Private Helper Functions */

static void WrongArgs(char *argv[]);
static void SetIC(N_Vector u, realtype dx);
static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi, int sensi_meth, int err_con);

/* Functions Called by the CVODES Solver */
static void f(realtype t, N_Vector u, N_Vector udot, void *f_data);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  NV_Spec nvSpec;
  void *cvode_mem;
  UserData data;
  realtype dx, reltol, abstol, t, tout;
  N_Vector u;
  int iout, flag;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS=NULL;
  booleantype sensi=FALSE;
  int sensi_meth=-1, err_con=-1;

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

  /* Initialize serial vector specification */
  nvSpec = NV_SpecInit_Serial(NEQ);

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data); /* Allocate data memory */
  data->p = (realtype *) malloc(NP * sizeof(realtype));
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->p[0] = 1.0;
  data->p[1] = 0.5;

  /* INITIAL STATES */
  u = N_VNew(nvSpec);    /* Allocate u vector */
  SetIC(u, dx);           /* Initialize u vector */

  /* TOLERANCES */
  reltol = 0.0;                /* Set the tolerances */
  abstol = ATOL;

  /* CVODE_CREATE */
  cvode_mem = CVodeCreate(ADAMS, FUNCTIONAL);
  if (cvode_mem == NULL) { printf("CVodeCreate failed.\n"); return(1); }

  flag = CVodeSetFdata(cvode_mem, data);
  if (flag != SUCCESS) { printf("CVodeSetFdata failed.\n"); return(1); }

  /* CVODE_MALLOC */
  flag = CVodeMalloc(cvode_mem, f, T0, u, SS, &reltol, &abstol, nvSpec);
  if (flag != SUCCESS) { printf("CVodeMalloc failed.\n"); return(1); }
  
  /* SENSITIVTY */
  if(sensi) {
    pbar  = (realtype *) malloc(NP * sizeof(realtype));
    pbar[0] = 1.0;
    pbar[1] = 0.5;
    plist = (int *) malloc(NS * sizeof(int));
    for(is=0; is<NS; is++)
      plist[is] = is+1; /* sensitivity w.r.t. i-th parameter */

    uS = N_VNew_S(NS, nvSpec);
    for(is=0;is<NS;is++)
      N_VConst(0.0, uS[is]);

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    flag = CVodeSetSensRho(cvode_mem, ZERO);
    flag = CVodeSetSensPbar(cvode_mem, pbar);

    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, plist, uS);
    if (flag != SUCCESS) { printf("CVodeSensMalloc failed, flag=%d\n",flag); return(1); }
  }

  /* In loop over output points, call CVode, print results, test for error */

  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n\n", MX);
  printf("============================================================\n");
  printf("     T     Q       H      NST                    Max norm   \n");
  printf("============================================================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    if (flag != SUCCESS) { 
      printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }
    PrintOutput(cvode_mem, t, u);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, uS);
      if (flag != OKAY) { printf("CVodeGetSens failed, flag=%d.\n", flag); break; }
      PrintOutputS(uS);
    } 
    printf("------------------------------------------------------------\n");
  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi,sensi_meth,err_con);

  /* Free memory */
  N_VFree(u);                  /* Free the u vector              */
  if(sensi) N_VFree_S(NS, uS); /* Free the uS vectors            */
  free(data->p);               /* Free the p vector              */
  free(data);                  /* Free block of UserData         */
  CVodeFree(cvode_mem);        /* Free the CVODES problem memory */
  NV_SpecFree_Serial(nvSpec);  /* Free the vector specification  */

  return(0);
}


/************************ Private Helper Functions ***********************/

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
/* Set initial conditions in u vector */

static void SetIC(N_Vector u, realtype dx)
{
  int i;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  udata = NV_DATA_S(u);

  /* Load initial profile into u vector */
  for (i=0; i<NEQ; i++) {
    x = (i+1)*dx;
    udata[i] = x*(XMAX - x)*exp(2.0*x);
  }  
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and max norm of solution  */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
  int nst, qu;
  realtype hu;
  
  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetLastOrder(cvode_mem, &qu);
  CVodeGetLastStep(cvode_mem, &hu);

  printf("%8.3e %2d  %8.3e %5d\n", t, qu, hu ,nst);
  printf("                                Solution       ");
  printf("%12.4e \n", N_VMaxNorm(u));
  
}
/* ======================================================================= */
/* Print max norm of sensitivities */

static void PrintOutputS(N_Vector *uS)
{

  printf("                                Sensitivity 1  ");
  printf("%12.4e \n", N_VMaxNorm(uS[0]));
  printf("                                Sensitivity 2  ");
  printf("%12.4e \n", N_VMaxNorm(uS[1]));

}

/* ======================================================================= */
/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem, booleantype sensi, int sensi_meth, int err_con)
{
  int nst;
  int nfe, nsetups, nni, ncfn, netf;
  int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;

  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetNumRhsEvals(cvode_mem, &nfe);
  CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  CVodeGetNumErrTestFails(cvode_mem, &netf);
  CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  if (sensi) {
    CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
    CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
    CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
    CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
    CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);
  }

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
  printf("nst     = %5d\n\n", nst);
  printf("nfe     = %5d\n",   nfe);
  printf("netf    = %5d    nsetups  = %5d\n", netf, nsetups);
  printf("nni     = %5d    ncfn     = %5d\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5d    nfeS     = %5d\n", nfSe, nfeS);
    printf("netfs   = %5d    nsetupsS = %5d\n", netfS, nsetupsS);
    printf("nniS    = %5d    ncfnS    = %5d\n", nniS, ncfnS);
  }

  printf("========================================================\n");


}

/***************** Function Called by the CVODES Solver ********************/

/* ======================================================================= */
/* f routine. Compute f(t,u). */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  realtype ui, ult, urt, hordc, horac, hdiff, hadv;
  realtype dx;
  realtype *udata, *dudata;
  int i;
  UserData data;

  udata = NV_DATA_S(u);
  dudata = NV_DATA_S(udot);

  /* Extract needed problem constants from data */
  data = (UserData) f_data;
  dx    = data->dx;
  hordc = data->p[0]/(dx*dx);
  horac = data->p[1]/(2.0*dx);

  /* Loop over all grid points. */
  for (i=0; i<NEQ; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = udata[i];
    if(i!=0) 
      ult = udata[i-1];
    else
      ult = 0.0;
    if(i!=NEQ-1)
      urt = udata[i+1];
    else
      urt = 0.0;

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - 2.0*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i] = hdiff + hadv;
  }
}
