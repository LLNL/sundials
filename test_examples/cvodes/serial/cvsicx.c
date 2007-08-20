/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2007-08-20 20:56:24 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example of using CVODES to compute sensitivities with respect to
 * initial conditions.
 *
 * The system consists of 3 ODEs modeling transport of lead through
 * the body:
 *   y1' = -(k01+k21+k31) y1 + k12 y2 + k13 y3 + b1
 *   y2' = k21 y1 - (k02+k12) y2
 *   y3' = k31 y1 - k13 y3
 *
 *   with nominal IC:
 *     y1(0) = y2(0) = y3(0) = 0
 *
 * We denote:
 *   a11 = -(k01+k21+k31)   a12 = k12          a13 = k13
 *   a21 = k21              a22 = - (k02+k12)  
 *   a31 = k31                                 a33 = -k13
 *  
 * NOTE: For readibility, no checks are performed on the various 
 *       function return flags.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>

/* Problem Constants */
#define  k01  RCONST(0.0211)
#define  k02  RCONST(0.0162)
#define  k21  RCONST(0.0111)
#define  k12  RCONST(0.0124)
#define  k31  RCONST(0.0039)
#define  k13  RCONST(0.000035)
#define  b1   RCONST(49.3)

#define y10   RCONST(0.0)
#define y20   RCONST(0.0)
#define y30   RCONST(0.0)

#define T0 RCONST(0.0)
#define TF RCONST(400.0)

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Accessor macros */
#define Ith(v,i) NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Type : UserData */
typedef struct {
  realtype a11, a12, a13;
  realtype a21, a22;
  realtype a31, a33;
  realtype p[3];
} *UserData;

/* Prototypes of functions by CVODES */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int Jac(int N, realtype t,
               N_Vector y, N_Vector fy, 
               DlsMat J, void *jac_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *fS_data, N_Vector tmp1, N_Vector tmp2);

/* Prototypes for private functions */
static int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, UserData data, FILE *f);
static void PrintFinalStats(void *cvode_mem);

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
  realtype reltol, abstol;
  N_Vector y0, y;

  int Ns, is;
  N_Vector *yS0, *yS;

  int flag;

  FILE *f1;

  /* 
   * Allocate and initialize user data structure 
   */

  data = (UserData) malloc(sizeof *data);
  data->a11 = -(k01+k21+k31);
  data->a12 = k12;
  data->a13 = k13;
  data->a21 = k21;
  data->a22 = - (k02+k12);
  data->a31 = k31;
  data->a33 = -k13;

  /*
   * Integration settings
   */

  /* Problem size */
  Neq = 3;

  /* Allocate and set vector of initial conditions */
  y0 = N_VNew_Serial(Neq);
  Ith(y0,1) = y10;
  Ith(y0,2) = y20;
  Ith(y0,3) = y30;

  /* Allocate solution vector */
  y = N_VNew_Serial(Neq);

  /* Set integration tolerances */
  reltol = RCONST(1e-7);
  abstol = RCONST(1e-6);

  /* Create, set, and allocate CVODES object*/
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeSetUserData(cvode_mem, data);
  flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
  flag = CVodeInit(cvode_mem, f, T0, y0);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);

  /* Attach linear solver */
  flag = CVDense(cvode_mem, Neq);
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);

  /*
   * Sensitivity-related settings
   */

  Ns = Neq;

  yS0 = N_VCloneVectorArray_Serial(Ns, y);
  for (is=0;is<Ns;is++) {
    N_VConst(ZERO, yS0[is]);
    Ith(yS0[is],is+1) = ONE;
  }

  yS = N_VCloneVectorArray_Serial(Ns, y);

  flag = CVodeSensInit1(cvode_mem, Ns, CV_SIMULTANEOUS, fS, yS0);
  flag = CVodeSensEEtolerances(cvode_mem);
  flag = CVodeSetSensErrCon(cvode_mem, TRUE);

  /* Note that, by not calling CVodeSetSensParams, we will use
     the default pbar values which means that CVODES will use 
     the same tolerances for sensitivities as for the state variables */

  /* Run CVODES */
  printf("Use analitycal sensitivity RHS (output written to cvsic.dat1)\n\n");
  f1 = fopen("cvsic.dat1","w");
  flag = runCVode(cvode_mem, y, yS, data, f1);
  fclose(f1);

  /* Free memory */

  N_VDestroy_Serial(y);
  N_VDestroy_Serial(y0);
  N_VDestroyVectorArray_Serial(yS, Ns);
  N_VDestroyVectorArray_Serial(yS0, Ns);

  /*free(plist);*/

  free(data);

  CVodeFree(&cvode_mem);

  return(0);

}


static int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, UserData data, FILE *f)
{
  realtype t;
  int flag;

  /* Call CVode in CV_ONE_STEP mode */  
  t = T0;
  while(t<TF) {
    flag = CVode(cvode_mem, TF, y, &t, CV_ONE_STEP);
    flag = CVodeGetSens(cvode_mem, &t, yS);
    fprintf(f, "%lf   ", t);
    fprintf(f, "%lf %lf %lf   ",Ith(y,1), Ith(y,2), Ith(y,3));
    fprintf(f, "%lf %lf %lf   ",Ith(yS[0],1), Ith(yS[0],2), Ith(yS[0],3));
    fprintf(f, "%lf %lf %lf   ",Ith(yS[1],1), Ith(yS[1],2), Ith(yS[1],3));
    fprintf(f, "%lf %lf %lf   ",Ith(yS[2],1), Ith(yS[2],2), Ith(yS[2],3));
    fprintf(f, "\n");
  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem);

  printf("\n");

  return(flag);

}


static void PrintFinalStats(void *cvode_mem)
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

  flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
  flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
  flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
  flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
  flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
  flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);

  flag = CVDlsGetNumJacEvals(cvode_mem, &njeD);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeD);

  printf("Run statistics:\n");

  printf("   nst     = %5ld\n", nst);
  printf("   nfe     = %5ld\n",   nfe);
  printf("   netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("   nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  printf("   njeD    = %5ld    nfeD     = %5ld\n", njeD, nfeD);

  printf("   -----------------------------------\n");
  printf("   nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
  printf("   netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
  printf("   nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);


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
  UserData data;
  realtype a11, a12, a13;
  realtype a21, a22;
  realtype a31, a33;
  realtype y1, y2, y3;

  /* extract constants from user data */

  data = (UserData) f_data;

  a11 = data->a11;
  a12 = data->a12;
  a13 = data->a13;

  a21 = data->a21;
  a22 = data->a22;
  
  a31 = data->a31;
  a33 = data->a33;

  /* extract solution components */

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  /* compute rhs */

  Ith(ydot,1) = a11 * y1 + a12 * y2  + a13 * y3  + b1;
  Ith(ydot,2) = a21 * y1 + a22 * y2;
  Ith(ydot,3) = a31 * y1             + a33 * y3;     

  return(0);
}

/* 
 * Jacobian routine. Compute J(t,y). 
 */

static int Jac(int N, realtype t,
               N_Vector y, N_Vector fy, 
               DlsMat J, void *jac_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype a11, a12, a13;
  realtype a21, a22;
  realtype a31, a33;

  IJth(J,1,1) = a11;  IJth(J,1,2) = a12;  IJth(J,1,3) = a13;
  IJth(J,2,1) = a21;  IJth(J,2,2) = a22;
  IJth(J,3,1) = a31;                      IJth(J,3,3) = a33; 

  return(0);
}

/* 
 * Compute RHS of the sensitivity equations with respect to the iS parameter.
 * When the parameters are IC (i.e. the ODE RHS does not depend on any of the 
 * parameters), the sensitivity RHS has the same form for all iS.
 */

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *fS_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype a11, a12, a13;
  realtype a21, a22;
  realtype a31, a33;

  realtype s1, s2, s3;

  /* extract constants from user data */

  data = (UserData) fS_data;

  a11 = data->a11;
  a12 = data->a12;
  a13 = data->a13;

  a21 = data->a21;
  a22 = data->a22;
  
  a31 = data->a31;
  a33 = data->a33;

  /* extract components of current sensitivity vector */

  s1 = Ith(yS,1); 
  s2 = Ith(yS,2); 
  s3 = Ith(yS,3);

  Ith(ySdot,1) = a11 * s1 + a12 * s2  + a13 * s3;
  Ith(ySdot,2) = a21 * s1 + a22 * s2;
  Ith(ySdot,3) = a31 * s1             + a33 * s3;  

  return(0);
}
