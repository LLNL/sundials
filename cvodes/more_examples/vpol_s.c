/************************************************************************
 * Van der Pol oscillator..                                             *
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.         *
 * This second-order ODE is converted to a first-order system by        *
 * defining y0 = x and y1 = xdot.                                       *
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "llnltyps.h" /* contains the definition for real, integer    */
#include "cvodes.h"   /* main CVODE header file                       */
#include "nvector.h"  /* contains the definition of type N_Vector     */
#include "llnlmath.h" /* contains the macros ABS, SQR                 */

/* Problem Constants */
#define NEQ        2
#define NOUT       4
#define T0         0.0
#define T1         1.39283880203
#define DTOUT      2.214773875
#define TOL_FACTOR 1.0e4

#define ATOL  1.0e-6
#define RTOL  0.0
#define ITOL  SS
#define ERRFP stdout
#define OPTIN FALSE

/* Private Helper Functions */
static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[]);

/* Functions Called by the CVODE Solver */
static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real tn,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);

/*=====================================================================*/

main()
{
  real ropt[OPT_SIZE], reltol=RTOL, abstol=ATOL, t, tout;
  long int iopt[OPT_SIZE];
  int flag, iout, ii;
  N_Vector y;
  void *cvode_mem;
  boole firstrun;

  real p[1], pbar[1], rhomax;
  integer NS, is, plist[1];
  N_Vector *yS;
  int sensi, sensi_meth, err_con;

  p[0] = 3.0; /* ETA */

  y = N_VNew(NEQ, NULL);

  N_VIth(y,0) = 2.0;
  N_VIth(y,1) = 0.0;
  
  cvode_mem = CVodeMalloc(NEQ, f, T0, y, ADAMS, FUNCTIONAL, SS,
                          &reltol, &abstol, p, ERRFP, OPTIN, iopt, ropt, NULL);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed."); return(1); }
  
  printf("\nPerform sensitivity analysis? (0:NO , 1:YES): ");scanf("%d",&sensi);
  if(sensi) {
    pbar[0] = p[0];
    plist[0] = 1; /* sensitivity w.r.t. 1st parameter */

    NS = 1;
    yS = N_VNew_S(NS,NEQ,NULL);
    N_VConst(0.0,yS[0]);

    printf("\nSensitivity method (%d:SIMULTANEOUS , %d:STAGGERED): ",SIMULTANEOUS,STAGGERED);
    scanf("%d",&sensi_meth);
    printf("\nError control (%d:FULL , %d:PARTIAL): ",FULL,PARTIAL);
    scanf("%d",&err_con);

    rhomax = 0.0;
    flag = CVodeSensMalloc(cvode_mem,NS,sensi_meth,p,pbar,plist,
                           ALLSENS,NULL,err_con,rhomax,yS,NULL,NULL);
    if (flag != SUCCESS) {
      printf("CVodeSensMalloc failed, flag=%d\n",flag);
      return(1);
    }
  }

  ii = 0;
  for(iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    do {
      flag = CVode(cvode_mem, tout, y, &t, ONE_STEP);
      if(ii%20==0) {
        printf("\n  NST      T        H      Q   NCF  NCFS   NEF  NEFS\n");
        printf(  "----------------------------------------------------\n");
      }
      printf("%5ld %8.3e %8.3e %2d %5ld %5ld %5ld %5ld\n",
             iopt[NST],t,ropt[HU],iopt[QU],iopt[NCFN],
             iopt[NCFNS],iopt[NETF],iopt[NETFS]);
      if (flag != SUCCESS) {printf("\n\n CVode error = %d \n\n", flag);return(1);}
      ii++;
    } while(t<tout);
  }

  PrintFinalStats(sensi,sensi_meth,err_con,iopt);
  
  CVodeFree(cvode_mem);
  N_VFree(y);

}

/*=====================================================================*/

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)
{
  real y0, y1;
  real *p, ETA;

  p = (real *) f_data;
  ETA = p[0];

  y0 = N_VIth(y,0);
  y1 = N_VIth(y,1);

  N_VIth(ydot,0) = y1;
  N_VIth(ydot,1) = (1.0 - SQR(y0))* ETA * y1 - y0;
} 

/*=====================================================================*/

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real tn,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{ 
  real y0, y1;
  real *p, ETA;

  p = (real *) f_data;
  ETA = p[0];

  y0 = N_VIth(y,0);
  y1 = N_VIth(y,1);

  DENSE_ELEM(J,0,1) = 1.0;
  DENSE_ELEM(J,1,0) = -2.0 * ETA * y0 * y1 - 1.0;
  DENSE_ELEM(J,1,1) = ETA * (1.0 - SQR(y0));
}

/*=====================================================================*/

static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[])
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



