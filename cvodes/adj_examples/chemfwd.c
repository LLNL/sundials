#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"
#include "cvodes.h"
#include "cvsdense.h"
#include "nvector.h"
#include "dense.h"

#define Ith(v,i)    N_VIth(v,i-1)         /* Ith numbers components 1..NEQ */
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
#define ATOLq 1e-6         /* absolute tol. for quadrature variable */
#define T0    0.0          /* initial time */
#define T1    4e7          /* final time   */
#define NP    3

#define ZERO  0.0

/* Type : UserData */
typedef struct {
  real p[3];
} *UserData;

/* Private Helper Function */

static void PrintFinalStats(long int iopt[]);
static void PrintOutput(real t, N_Vector u, N_Vector *uS);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);
static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  UserData data;
  real pbar[NP], rhomax, ropt[OPT_SIZE], reltol, t;
  long int iopt[OPT_SIZE];
  N_Vector y, abstol, *yS;
  void *cvode_mem;
  int i, is, flag, *plist;

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = 0.04;
  data->p[1] = 1.0e4;
  data->p[2] = 3.0e7;

  /* INITIAL STATES */
  y = N_VNew(NEQ+1, NULL); 
  abstol = N_VNew(NEQ+1, NULL); 

  /* Initialize y */
  Ith(y,1) = Y1;                
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;
  Ith(y,4) = ZERO;

  /* TOLERANCES */
  /* Set the scalar relative tolerance */
  reltol = RTOL;               
  /* Set the vector absolute tolerance */
  Ith(abstol,1) = ATOL1;       
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;
  Ith(abstol,4) = ATOLq;

  /* OPTIONAL I/O */
  for(i=0; i<OPT_SIZE; i++) {
    iopt[i] = 0;
    ropt[i] = 0.0;
  }
  iopt[MXSTEP] = 1000;

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(NEQ+1, f, T0, y, BDF, NEWTON, SV, &reltol, abstol,
                          data, NULL, TRUE, iopt, ropt, NULL);
  if (cvode_mem == NULL) { 
    printf("CVodeMalloc failed.\n"); 
    return(1); 
  }

  /* CVDENSE */
  flag = CVDense(cvode_mem, Jac, NULL);
  if (flag != SUCCESS) { printf("CVDense failed.\n"); return(1); }

  /* SENSITIVITY */
  pbar[0] = data->p[0];
  pbar[1] = data->p[1];
  pbar[2] = data->p[2];
  plist = (integer *) malloc(NP * sizeof(integer));
  for(is=0; is<NP; is++) plist[is] = is+1;
  
  yS = N_VNew_S(NP, NEQ+1, NULL);
  for(is=0; is<NP; is++)
    N_VConst(0.0,yS[is]);

  rhomax = ZERO;
  flag = CVodeSensMalloc(cvode_mem,NP,STAGGERED,data->p,pbar,plist,
                         ALLSENS,NULL,FULL,rhomax,yS,NULL,NULL);
  if (flag != SUCCESS) {
    printf("CVodeSensMalloc failed, flag=%d\n",flag);
    return(1);
  }

  flag = CVode(cvode_mem, T1, y, &t, NORMAL);

  if (flag == SUCCESS) {
    flag = CVodeSensExtract(cvode_mem, T1, yS);
    PrintOutput(t, y, yS);
  } else {
    printf("CVode failed, flag=%d.\n", flag); 
  }

  /* Print final statistics */
  PrintFinalStats(iopt);       /* Print some final statistics   */

  /* Free memory */
  N_VFree(y);                  /* Free the y and abstol vectors */
  N_VFree(abstol);   
  free(data);
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */

  return(0);
}


/************************ Private Helper Function ************************/
/* ======================================================================= */
/* Print t, solution, and gradient of derived function */

static void PrintOutput(real t, N_Vector u, N_Vector *uS)
{

  real *udata, Gp1, Gp2, Gp3;
  
  printf("\n\n========================================================\n");

  printf("t:  %8.3e\n", t);

  udata = N_VDATA(u);

  printf("y:  %12.4e %12.4e %12.4e\n", udata[0], udata[1], udata[2]);
  printf("G:  %12.4e\n", udata[3]);
  
  udata = N_VDATA(uS[0]);
  Gp1 = udata[3];
  udata = N_VDATA(uS[1]);
  Gp2 = udata[3];
  udata = N_VDATA(uS[2]);
  Gp3 = udata[3];

  printf("Gp: %12.4e %12.4e %12.4e\n", Gp1, Gp2, Gp3);

}
/* ======================================================================= */
/* Print some final statistics located in the iopt array */

static void PrintFinalStats(long int iopt[])
{

  printf("\n========================================================");
  printf("\nFinal Statistics");

  printf("\n\n");
  /*
  printf("lenrw   = %5ld    leniw = %5ld\n", iopt[LENRW], iopt[LENIW]);
  printf("llrw    = %5ld    lliw  = %5ld\n", iopt[SPGMR_LRW], iopt[SPGMR_LIW]);
  */
  printf("nst     = %5ld\n\n", iopt[NST]);
  printf("nfe     = %5ld\n", iopt[NFE]);
  printf("nni     = %5ld\n", iopt[NNI]);
  printf("ncfn    = %5ld\n", iopt[NCFN]);
  printf("netf    = %5ld\n\n", iopt[NETF]);
  printf("nsetups = %5ld\n", iopt[NSETUPS]);
  printf("nje     = %5ld\n", iopt[DENSE_NJE]);

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

        Ith(ydot,4) = y1 + p2*y2*y3;
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
  IJth(J,4,1) = 1.0;  IJth(J,4,2) = p2*y3;          IJth(J,4,3) = p2*y2;
}
 
