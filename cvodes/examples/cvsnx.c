/************************************************************************
 *                                                                      *
 * File: pvnx.c                                                         *
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne, and *
 *              Radu Serban @ LLNL                                      *
 * Version of 14 November 2001                                          *
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
#include <math.h>
#include "llnltyps.h"
#include "cvodes.h"
#include "nvector.h"

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
  real *p;
  real dx, hdcoef, hacoef;
} *UserData;


/* Private Helper Functions */

static void WrongArgs(char *argv[]);
static void SetIC(N_Vector u, real dx, integer N);
static void PrintOutput(long int iopt[], real ropt[], real t, N_Vector u);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(boole sensi, int sensi_meth, int err_con, long int iopt[]);

/* Functions Called by the CVODES Solver */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data);


/***************************** Main Program ******************************/

main(int argc, char *argv[])
{
  real ropt[OPT_SIZE], dx, reltol, abstol, t, tout;
  long int iopt[OPT_SIZE];
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, flag;
  machEnvType machEnv;

  real *pbar, rhomax;
  integer is, *plist;
  N_Vector *uS;
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

  machEnv = NULL;

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data); /* Allocate data memory */
  data->p = (real *) malloc(NP * sizeof(real));
  dx = data->dx = XMAX/((real)(MX+1));
  data->p[0] = 1.0;
  data->p[1] = 0.5;

  /* INITIAL STATES */
  u = N_VNew(NEQ, machEnv);    /* Allocate u vector */
  SetIC(u, dx, NEQ);           /* Initialize u vector */

  /* TOLERANCES */
  reltol = 0.0;                /* Set the tolerances */
  abstol = ATOL;

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(NEQ, f, T0, u, ADAMS, FUNCTIONAL, SS, &reltol,
                          &abstol, data, NULL, FALSE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { 
    printf("CVodeMalloc failed.\n");
    return(1);
  }

  /* SENSITIVTY */
  if(sensi) {
    pbar  = (real *) malloc(NP * sizeof(real));
    pbar[0] = 1.0;
    pbar[1] = 0.5;
    plist = (integer *) malloc(NS * sizeof(integer));
    for(is=0; is<NS; is++)
      plist[is] = is+1; /* sensitivity w.r.t. i-th parameter */

    uS = N_VNew_S(NS,NEQ,machEnv);
    for(is=0;is<NS;is++)
      N_VConst(0.0,uS[is]);

    rhomax = ZERO;
    ifS = ALLSENS;
    if(sensi_meth==STAGGERED1) ifS = ONESENS;
    flag = CVodeSensMalloc(cvode_mem,NS,sensi_meth,data->p,pbar,plist,
                           ifS,NULL,err_con,rhomax,uS,NULL,NULL);
    if (flag != SUCCESS) {
      printf("CVodeSensMalloc failed, flag=%d\n",flag);
      return(1);
    }
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
    PrintOutput(iopt, ropt, t, u);
    if (sensi) {
      flag = CVodeSensExtract(cvode_mem, t, uS);
      if (flag != SUCCESS) { 
        printf("CVodeSensExtract failed, flag=%d.\n", flag); 
        break; 
      }
      PrintOutputS(uS);
    } 
    printf("------------------------------------------------------------\n");
  }

  /* Print final statistics */
  PrintFinalStats(sensi,sensi_meth,err_con,iopt);

  /* Free memory */
  N_VFree(u);                  /* Free the u vector             */
  if(sensi) N_VFree_S(NS, uS); /* Free the uS vectors           */
  free(data->p);               /* Free the p vector             */
  free(data);                  /* Free block of UserData        */
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */

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

static void SetIC(N_Vector u, real dx, integer N)
{
  int i;
  real x;
  real *udata;

  /* Set pointer to data array and get local length of u. */
  udata = N_VDATA(u);

  /* Load initial profile into u vector */
  for (i=0; i<N; i++) {
    x = (i+1)*dx;
    udata[i] = x*(XMAX - x)*exp(2.0*x);
  }  
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and max norm of solution  */

static void PrintOutput(long int iopt[], real ropt[], real t, N_Vector u)
{
  
  printf("%8.3e %2d  %8.3e %5ld\n", t,iopt[QU],ropt[HU],iopt[NST]);
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

  printf("========================================================\n");

}

/***************** Function Called by the CVODES Solver ********************/

/* ======================================================================= */
/* f routine. Compute f(t,u). */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data)
{
  real ui, ult, urt, hordc, horac, hdiff, hadv;
  real dx;
  real *udata, *dudata;
  int i;
  UserData data;

  udata = N_VDATA(u);
  dudata = N_VDATA(udot);

  /* Extract needed problem constants from data */
  data = (UserData) f_data;
  dx    = data->dx;
  hordc = data->p[0]/(dx*dx);
  horac = data->p[1]/(2.0*dx);

  /* Loop over all grid points. */
  for (i=0; i<N; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = udata[i];
    if(i!=0) 
      ult = udata[i-1];
    else
      ult = 0.0;
    if(i!=N-1)
      urt = udata[i+1];
    else
      urt = 0.0;

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - 2.0*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i] = hdiff + hadv;
  }
}
