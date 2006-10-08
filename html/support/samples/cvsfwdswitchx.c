/*
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Modification of the cvsfwddendx to illustrate switching on and off
 * sensitivity computations.
 *
 * NOTE: For readibility, no checks are performed on the various 
 *       function return flags.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>

/* Accessor macros */
#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Problem Constants */
#define MXSTEPS 2000
#define T0      RCONST(0.0)
#define T1      RCONST(4.0e10)

#define ZERO  RCONST(0.0)

/* Type : UserData */
typedef struct {
  booleantype sensi;
  booleantype errconS;
  booleantype fsDQ;
  int meth;
  realtype p[3];
} *UserData;

/* Prototypes of functions by CVODES */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int Jac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *fS_data, N_Vector tmp1, N_Vector tmp2);

/* Prototypes of private functions */
static int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, UserData data);
static void PrintHeader(UserData data);
static void PrintFinalStats(void *cvode_mem, UserData data);

/* Readibility replacements */
#define sensi   data->sensi
#define errconS data->errconS
#define fsDQ    data->fsDQ
#define meth    data->meth

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;
  void *cvode_mem;

  long int Neq;
  realtype reltol;
  N_Vector y0, y, abstol;

  int Ns;
  realtype *pbar;
  int is, *plist; 
  N_Vector *yS0, *yS;

  int flag;

  /* Allocate and initialize parameters in user data structure */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  /* Problem size */
  Neq = 3;

  /* Allocate vectors */
  y0 = N_VNew_Serial(Neq);      /* initial conditions */
  y = N_VNew_Serial(Neq);       /* solution vector */
  abstol = N_VNew_Serial(Neq);  /* absolute tolerances */

  /* Set initial conditions */
  Ith(y0,1) = RCONST(1.0);
  Ith(y0,2) = RCONST(0.0);
  Ith(y0,3) = RCONST(0.0);

  /* Set integration tolerances */
  reltol = RCONST(1e-6);
  Ith(abstol,1) = RCONST(1e-8);
  Ith(abstol,2) = RCONST(1e-14);
  Ith(abstol,3) = RCONST(1e-6);

  /* Create, set, and allocate CVODES object*/
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeSetFdata(cvode_mem, data);
  flag = CVodeSetMaxNumSteps(cvode_mem, MXSTEPS);
  flag = CVodeMalloc(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);

  /* Attach linear solver */
  flag = CVDense(cvode_mem, Neq);
  flag = CVDenseSetJacFn(cvode_mem, Jac, data);

  /* Sensitivity-related settings */

  sensi = TRUE;           /* sensitivity ON */
  meth = CV_SIMULTANEOUS; /* simultaneous corrector method */
  errconS = TRUE;         /* full error control */
  fsDQ = FALSE;           /* user-provided sensitvity RHS */

  Ns = 3;

  pbar = (realtype *) malloc(Ns * sizeof(realtype));
  pbar[0] = data->p[0];
  pbar[1] = data->p[1];
  pbar[2] = data->p[2];

  plist = (int *) malloc(Ns * sizeof(int));
  for (is=0; is<Ns; is++) plist[is] = is;

  yS0 = N_VCloneVectorArray_Serial(Ns, y);
  for (is=0;is<Ns;is++) N_VConst(ZERO, yS0[is]);

  yS = N_VCloneVectorArray_Serial(Ns, y);
  
  flag = CVodeSensMalloc(cvode_mem, Ns, meth, yS0);

  flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);

  /*
    Sensitivities are enabled
    Set full error control
    Set user-provided sensitivity RHS
    Run CVODES
  */

  flag = CVodeSetSensErrCon(cvode_mem, errconS);
  flag = CVodeSetSensRhs1Fn(cvode_mem, fS, data);

  flag = runCVode(cvode_mem, y, yS, data);

  /*
    Change parameters
    Toggle sensitivities OFF
    Reinitialize and run CVODES
  */

  data->p[0] = RCONST(0.05);
  data->p[1] = RCONST(2.0e4);
  data->p[2] = RCONST(2.9e7);

  sensi = FALSE;

  flag = CVodeSensToggleOff(cvode_mem);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);
  flag = runCVode(cvode_mem, y, yS, data);

  /*
    Change parameters
    Switch to internal DQ sensitivity RHS function
    Toggle sensitivities ON (reinitialize sensitivities)
    Reinitialize and run CVODES
  */
  
  data->p[0] = RCONST(0.06);
  data->p[1] = RCONST(3.0e4);
  data->p[2] = RCONST(2.8e7);

  sensi = TRUE;
  fsDQ = TRUE;

  flag = CVodeSetSensRhs1Fn(cvode_mem, NULL, NULL);
  flag = CVodeSensReInit(cvode_mem, meth, yS0);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);
  flag = runCVode(cvode_mem, y, yS, data);

  /*
    Switch to partial error control
    Switch back to user-provided sensitivity RHS
    Toggle sensitivities ON (reinitialize sensitivities)
    Change method to staggered
    Reinitialize and run CVODES
  */

  sensi = TRUE;
  errconS = FALSE;
  fsDQ = FALSE;
  meth = CV_STAGGERED;

  flag = CVodeSetSensErrCon(cvode_mem, errconS);
  flag = CVodeSetSensRhs1Fn(cvode_mem, fS, data);
  flag = CVodeSensReInit(cvode_mem, meth, yS0);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);
  flag = runCVode(cvode_mem, y, yS, data);

  /*
    Free sensitivity-related memory
    (CVodeSensToggle is not needed, as CVodeSensFree toggles sensitivities OFF)
    Reinitialize and run CVODES
  */
  
  sensi = FALSE;

  CVodeSensFree(cvode_mem);
  flag = CVodeReInit(cvode_mem, f, T0, y0, CV_SV, reltol, abstol);
  flag = runCVode(cvode_mem, y, yS, data);
  
  /* Free memory */

  N_VDestroy_Serial(y0);                 /* Free y0 vector */
  N_VDestroy_Serial(y);                  /* Free y vector */
  N_VDestroy_Serial(abstol);             /* Free abstol vector */
  N_VDestroyVectorArray_Serial(yS0, Ns); /* Free yS0 vector */
  N_VDestroyVectorArray_Serial(yS, Ns);  /* Free yS vector */
  free(plist);                           /* Free plist */

  free(data);                            /* Free user data */
  CVodeFree(&cvode_mem);                 /* Free CVODES memory */

  return(0);

}

static int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, UserData data)
{
  realtype t;
  int flag;

  /* Print header for current run */
  PrintHeader(data);

  /* Call CVode in CV_NORMAL mode */  
  flag = CVode(cvode_mem, T1, y, &t, CV_NORMAL);

  /* Print final statistics */
  PrintFinalStats(cvode_mem, data);

  printf("\n");

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

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
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

  return(0);
}


/* 
 * Jacobian routine. Compute J(t,y). 
 */

static int Jac(long int N, DenseMat J, realtype t,
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

  return(0);
}
 
/* 
 * fS routine. Compute sensitivity r.h.s. 
 */

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
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

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

static void PrintHeader(UserData data)
{
  /* Print sensitivity control flags */
  printf("Sensitivity: ");
  if (sensi) {
    printf("YES (");
    switch (meth) {
    case CV_SIMULTANEOUS:
      printf("SIMULTANEOUS + ");
      break;
    case CV_STAGGERED:
      printf("STAGGERED + ");
      break;
    case CV_STAGGERED1:
      printf("STAGGERED-1 + ");
      break;
    }
    if (errconS) printf("FULL ERROR CONTROL + ");
    else         printf("PARTIAL ERROR CONTROL + ");
    if (fsDQ)    printf("DQ sensitivity RHS)\n");
    else         printf("user-provided sensitivity RHS)\n");
  } else {
    printf("NO\n");
  }

  /* Print current problem parameters */
  printf("Parameters: [%8.4e  %8.4e  %8.4e]\n",data->p[0], data->p[1], data->p[2]);
}

/* 
 * Print some final statistics from the CVODES memory.
 */

static void PrintFinalStats(void *cvode_mem, UserData data)
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

  printf("Run statistics:\n");

  printf("   nst     = %5ld\n", nst);
  printf("   nfe     = %5ld\n",   nfe);
  printf("   netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("   nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  printf("   njeD    = %5ld    nfeD     = %5ld\n", njeD, nfeD);

  if(sensi) {
    printf("   -----------------------------------\n");
    printf("   nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("   netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("   nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }


}

