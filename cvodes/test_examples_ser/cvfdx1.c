/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2005-07-01 00:00:30 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Modification of the cvfdx to illustrate switching on and off
 * sensitivity computations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sundialstypes.h"
#include "cvodes.h"
#include "cvdense.h"
#include "nvector_serial.h"
#include "dense.h"

/* Accessor macros */
#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Problem Constants */
#define NEQ   3
#define Y1    RCONST(1.0)
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1e-6)
#define ATOL1 RCONST(1e-8)
#define ATOL2 RCONST(1e-14)
#define ATOL3 RCONST(1e-6)
#define T0    RCONST(0.0)
#define T1    RCONST(0.4)
#define TMULT RCONST(10.0)
#define NOUT  12

#define NP    3
#define NS    3

#define ZERO  RCONST(0.0)

/* Type : UserData */
typedef struct {
  realtype p[3];
} *UserData;

/* Prototypes of functions by CVODES */
static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static void Jac(long int N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static void fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
               int iS, N_Vector yS, N_Vector ySdot, 
               void *fS_data, N_Vector tmp1, N_Vector tmp2);

/* Prototypes of private functions */
int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, int sensi);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *cvode_mem;
  UserData data;
  realtype reltol;
  N_Vector y0, y, abstol;
  int flag;

  realtype pbar[NP];
  int is, *plist; 
  N_Vector *yS0, *yS;
  booleantype sensi;

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  /* Initial conditions */
  y0 = N_VNew_Serial(NEQ);
  Ith(y0,1) = Y1;
  Ith(y0,2) = Y2;
  Ith(y0,3) = Y3;

  /* Solution vector */
  y = N_VNew_Serial(NEQ);
  
  /* Tolerances: scalar relative tolerance, vector absolute tolerance */
  reltol = RTOL;               
  abstol = N_VNew_Serial(NEQ);
  Ith(abstol,1) = ATOL1;       
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* Create, set, and allocate CVODES object*/
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeSetFdata(cvode_mem, data);
  flag = CVodeMalloc(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);

  /* Attach linear solver */
  flag = CVDense(cvode_mem, NEQ);
  flag = CVDenseSetJacFn(cvode_mem, Jac, data);

  /* Sensitivity-related settings */

  pbar[0] = data->p[0];
  pbar[1] = data->p[1];
  pbar[2] = data->p[2];
  plist = (int *) malloc(NS * sizeof(int));
  for (is=0; is<NS; is++) plist[is] = is+1;

  yS0 = N_VCloneVectorArray_Serial(NS, y);
  for (is=0;is<NS;is++) N_VConst(ZERO, yS0[is]);

  yS = N_VCloneVectorArray_Serial(NS, y);
  
  flag = CVodeSensMalloc(cvode_mem, NS, CV_SIMULTANEOUS, yS0);

  flag = CVodeSetSensRhs1Fn(cvode_mem, fS);
  flag = CVodeSetSensErrCon(cvode_mem, TRUE);
  flag = CVodeSetSensFdata(cvode_mem, data);
  flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);

  /* Deactivate sensitivity and run CVODE */

  sensi = FALSE;
  flag = CVodeSensToggle(cvode_mem, sensi);

  printf("Sensitivity: NO\n");
  flag = runCVode(cvode_mem, y, yS, sensi);

  /* Reactivate sensitivity and run CVODE with sensitivity */

  sensi = TRUE;
  flag = CVodeSensToggle(cvode_mem, sensi);
  
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);

  printf("Sensitivity: YES (SIMULTANEOUS + FULL ERROR CONTROL)\n");
  flag = runCVode(cvode_mem, y, yS, sensi);

  /* Deactivate sensitivity and run CVODE twice with different params. */

  sensi = FALSE;
  flag = CVodeSensToggle(cvode_mem, sensi);

  data->p[0] = RCONST(0.05);
  data->p[1] = RCONST(2.0e4);
  data->p[2] = RCONST(2.9e7);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);

  printf("Sensitivity: NO\n");
  flag = runCVode(cvode_mem, y, yS, sensi);

  data->p[0] = RCONST(0.06);
  data->p[1] = RCONST(3.0e4);
  data->p[2] = RCONST(2.8e7);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);

  printf("Sensitivity: NO\n");
  flag = runCVode(cvode_mem, y, yS, sensi);

  /* Reactivate sensitivity and run CVODE again using a different method 
     NOTE: CVodeSensToggle is not needed, as CVodeSensReInit does its job */

  sensi = TRUE;

  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);
  flag = CVodeSensReInit(cvode_mem, CV_STAGGERED, yS0);
  flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
  flag = CVodeSetSensErrCon(cvode_mem, FALSE);

  printf("Sensitivity: YES (STAGGERED + PARTIAL ERROR CONTROL)\n");
  flag = runCVode(cvode_mem, y, yS, sensi);

  /* Free sensitivity-related memory and do one more CVODE run 
     NOTE: CVodeSensToggle is not needed, as CVodeSensFree does its job */
  
  sensi = FALSE;  
  CVodeSensFree(cvode_mem);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);

  printf("Sensitivity: NO\n");
  flag = runCVode(cvode_mem, y, yS, sensi);
  
  /* Free memory */

  N_VDestroy_Serial(y0);                 /* Free y0 vector */
  N_VDestroy_Serial(y);                  /* Free y vector */
  N_VDestroy_Serial(abstol);             /* Free abstol vector */
  N_VDestroyVectorArray_Serial(yS0, NS); /* Free yS0 vector */
  N_VDestroyVectorArray_Serial(yS, NS);  /* Free yS vector */
  free(plist);                           /* Free plist */

  free(data);                            /* Free user data */
  CVodeFree(cvode_mem);                  /* Free CVODES memory */

  return(0);

}

int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, int sensi)
{

  realtype t, tout;
  int iout, flag;

  /* In loop over output points, call CVode, print results, test for error */
  
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, yS);
    }

  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi);

  printf("\n\n");

  return(flag);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y). 
 */

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


/* 
 * Jacobian routine. Compute J(t,y). 
 */

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
 
/* 
 * fS routine. Compute sensitivity r.h.s. 
 */

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

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* 
 * Print some final statistics from the CVODES memory.
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int njeD, nfeD;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  if (sensi) {
    flag = CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    flag = CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
    flag = CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
    flag = CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
    flag = CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);
  }

  flag = CVDenseGetNumJacEvals(cvode_mem, &njeD);
  flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeD);

  printf("nst     = %5ld\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  printf("njeD    = %5ld    nfeD     = %5ld\n", njeD, nfeD);

  if(sensi) {
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }


}

