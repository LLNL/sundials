#include <stdio.h>
#include <math.h>
#include "llnltyps.h"
#include "cvodes.h"
#include "nvector.h"
#include "cvsdense.h"

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

static void SetIC(N_Vector u, real dx, integer N);

static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[], int iter);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data);


/***************************** Main Program ******************************/

main()
{
  real ropt[OPT_SIZE], dx, reltol, abstol, t, tout;
  long int iopt[OPT_SIZE];
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, flag;
  machEnvType machEnv;
  int lmm, iter;

  real *pbar, rhomax;
  integer is, *plist;
  N_Vector *uS;
  int sensi, sensi_meth, err_con, ifS;

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

  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
  
  printf("\nLinear multistep method (%d:ADAMS , %d:BDF): ",ADAMS,BDF);
  scanf("%d",&lmm);
  printf("\nNonlinear iteration (%d:FUNCTIONAL , %d:NEWTON): ",FUNCTIONAL,NEWTON);
  scanf("%d",&iter);

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(NEQ, f, T0, u, lmm, iter, SS, &reltol,
                          &abstol, data, NULL, FALSE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return(1); }

  if(iter==NEWTON) {
    flag = CVDense(cvode_mem,NULL,NULL);
    if (flag != SUCCESS) { printf("CVDense failed.\n"); return(1); }
  }

  /* SENSITIVTY DATA */
  printf("\nPerform sensitivity analysis? (0:NO , 1:YES): ");scanf("%d",&sensi);
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

    printf("\nSensitivity method (%d:SIMULTANEOUS , %d:STAGGERED, %d:STAGGERED1): ",
           SIMULTANEOUS,STAGGERED,STAGGERED1);
    scanf("%d",&sensi_meth);
    printf("\nError control (%d:FULL , %d:PARTIAL): ",FULL,PARTIAL);
    scanf("%d",&err_con);

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

  printf("At t = %4.2f    max.norm(u) =%14.6e \n", T0,N_VMaxNorm(u));
  
  if(sensi) {
    for(is=0;is<NS;is++)
      printf("sensitivity s_%d:  max.norm =%14.6e \n", is, N_VMaxNorm(uS[is]));
    printf("\n");
  }

  /* In loop over output points, call CVode, print results, test for error */

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    printf("At t = %4.2f    max.norm(u) =%14.6e   nst =%4d \n", 
           t,N_VMaxNorm(u),iopt[NST]);
    if (flag != SUCCESS) { printf("CVode failed, flag=%d.\n", flag); break; }

    if(sensi) {
      flag = CVodeSensExtract(cvode_mem, t, uS);
      for(is=0;is<NS;is++)
        printf("sensitivity s_%d:  max.norm =%14.6e \n", is, N_VMaxNorm(uS[is]));
      printf("\n");
    }
  }

  PrintFinalStats(sensi,sensi_meth,err_con,iopt,iter);

  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */
  N_VFree(u);                  /* Free the u vector */

  free(data->p);
  free(data);                  /* Free block of UserData */

  return(0);
}


/************************ Private Helper Functions ***********************/

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


/* Print some final statistics located in the iopt array */

static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[], int iter)
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
  if(iter == NEWTON) {
    printf("nje     = %5ld                  \n", iopt[DENSE_NJE]);
  }

  printf("========================================================\n");

}


/***************** Function Called by the CVODE Solver ******************/

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
