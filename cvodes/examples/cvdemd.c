/**************************************************************************
 *                                                                        *
 * File        : cvdemd.c                                                 *
 * Programmers : Scott D. Cohen and Alan C. Hindmarsh @ LLNL              *
 * Version of  : 26 June 2002                                             *
 *------------------------------------------------------------------------*
 * Modified by R. Serban to work with new serial nvector (5/3/2002)       *
 *------------------------------------------------------------------------*
 *                                                                        *
 * Demonstration program for CVODE - direct linear solvers. Two           *
 * separate problems are solved using both the ADAMS and BDF linear       *
 * multistep methods in combination with FUNCTIONAL and NEWTON            *
 * iterations :                                                           *
 *                                                                        *
 * Problem 1.. Van der Pol oscillator..                                   *
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.           *
 * This second-order ODE is converted to a first-order system by          *
 * defining y0 = x and y1 = xdot.                                         *
 * The NEWTON iteration cases use the following types of Jacobian         *
 * approximation: (1) dense, user-supplied, (2) dense, difference         *
 * quotient approximation, (3) diagonal approximation.                    *
 *                                                                        * 
 * Problem 2.. ydot = A * y, where A is a banded lower triangular         *
 * matrix derived from 2-D advection PDE.                                 *
 * The NEWTON iteration cases use the following types of Jacobian         *
 * approximation: (1) band, user-supplied, (2) band, difference           *
 * quotient approximation, (3) diagonal approximation.                    *
 *                                                                        *
 * For each problem, in the series of eight runs, CVodeMalloc is called   *
 * only once, for the first run, whereas CVReInit is called for each of   *
 * the remaining seven runs.                                              *
 *                                                                        *
 * Notes.. This program demonstrates the usage of the sequential CVODE    *
 * macros NV_Ith_S, NV_DATA_S, DENSE_ELEM, BAND_COL, and BAND_COL_ELEM.   *
 * The NV_Ith_S macro is used to reference the components of an N_Vector. *
 * It works for any size N=NEQ, but due to efficiency concerns it should  *
 * only by used when the problem size is small. The Problem 1 right       *
 * hand side and Jacobian functions f1 and Jac1 both use NV_Ith_S. The    *
 * NV_DATA_S macro gives the user access to the memory used for the       *
 * component storage of an N_Vector. In the sequential case, the user     *
 * may assume that this is one contiguous array of reals. The NV_DATA_S   *
 * macro gives a more efficient (than the NV_Ith_S macro) means to access *
 * the components of an N_Vector and should be used when the problem      *
 * size is large. The Problem 2 right hand side function f2 uses the      *
 * NV_DATA_S macro. The DENSE_ELEM macro used in Jac1 gives access to an  *
 * element of a dense matrix of type DenseMat. It should be used only     *
 * when the problem size is small (the size of a DenseMat is NEQ x NEQ)   *
 * due to efficiency concerns. For larger problem sizes, the macro        *
 * DENSE_COL can be used in order to work directly with a column of a     *
 * DenseMat. The BAND_COL and BAND_COL_ELEM allow efficient columnwise    *
 * access to the the elements of a band matrix of type BandMat. These     *
 * macros are used in the Jac2 function.                                  *
 **************************************************************************/

#include <stdio.h>
#include <math.h>
#include "sundialstypes.h"  /* definitions for realtype, integertype        */
#include "cvodes.h"         /* main CVODE header file                       */
#include "cvsdense.h"       /* use CVDENSE linear solver each internal step */
#include "cvsband.h"        /* use CVBAND linear solver each internal step  */
#include "cvsdiag.h"        /* use CVDIAG linear solver each internal step  */
#include "nvector_serial.h" /* contains the definition of type N_Vector     */
#include "sundialsmath.h"   /* contains the macros ABS, SQR                 */

/* Shared Problem Constants */

#define ATOL  1.0e-6
#define RTOL  0.0
#define ITOL  SS
#define ERRFP stdout
#define OPTIN FALSE

/* Problem #1 Constants */

#define P1_NEQ        2
#define P1_ETA        3.0
#define P1_NOUT       4
#define P1_T0         0.0
#define P1_T1         1.39283880203
#define P1_DTOUT      2.214773875
#define P1_TOL_FACTOR 1.0e4

/* Problem #2 Constants */

#define P2_MESHX        5
#define P2_MESHY        5
#define P2_NEQ          (P2_MESHX * P2_MESHY)
#define P2_ALPH1        1.0
#define P2_ALPH2        1.0
#define P2_NOUT         5
#define P2_ML           5
#define P2_MU           0
#define P2_T0           0.0
#define P2_T1           0.01
#define P2_TOUT_MULT   10.0
#define P2_TOL_FACTOR   1.0e3

/* Linear Solver Options */

enum {FUNC, DENSE_USER, DENSE_DQ, DIAG, BAND_USER, BAND_DQ};

/* Private Helper Functions */

static int  Problem1(void);
static void PrintIntro1(void);
static int  Problem2(void);
static void PrintIntro2(void);
static realtype MaxError(N_Vector y, realtype t);
static int PrepareNextRun(void *cvode_mem, int lmm, int miter, integertype mu,
                          integertype ml);
static void PrintFinalStats(long int iopt[], int miter, realtype ero);

/* Functions Called by the CVODE Solver */

static void f1(integertype N, realtype t, N_Vector y, N_Vector ydot, void *f_data);

static void Jac1(integertype N, DenseMat J, RhsFn f, void *f_data, realtype tn,
                 N_Vector y, N_Vector fy, N_Vector ewt, realtype h, realtype uround,
                 void *jac_data, long int *nfePtr, N_Vector vtemp1,
                 N_Vector vtemp2, N_Vector vtemp3);

static void f2(integertype N, realtype t, N_Vector y, N_Vector ydot, void *f_data);

static void Jac2(integertype N, integertype mu, integertype ml, BandMat J, RhsFn f,
                 void *f_data, realtype tn, N_Vector y, N_Vector fy, N_Vector ewt,
                 realtype h, realtype uround, void *jac_data, long int *nfePtr,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 

/* Implementation */

int main()
{
  int nerr;

  nerr = Problem1();
  nerr += Problem2();
  printf("\n\n-------------------------------------------------------------");
  printf("\n-------------------------------------------------------------");
  printf("\n\n Number of errors encountered = %d \n", nerr);
  return(0);
}

static int Problem1(void)
{
  M_Env machEnv;
  realtype ropt[OPT_SIZE], reltol=RTOL, abstol=ATOL, t, tout, ero, er;
  long int iopt[OPT_SIZE];
  int lmm, miter, flag, iter, iout, nerr=0;
  N_Vector y;
  void *cvode_mem = NULL;
  booleantype firstrun;

  machEnv = M_EnvInit_Serial(P1_NEQ);

  y = N_VNew(P1_NEQ, machEnv);
  PrintIntro1();

  for (lmm=ADAMS; lmm <= BDF; lmm++) {
    for (miter=FUNC; miter <= DIAG; miter++) {
      ero = 0.0;
      NV_Ith_S(y,0) = 2.0;
      NV_Ith_S(y,1) = 0.0;
      iter = (miter == FUNC) ? FUNCTIONAL : NEWTON;
      
      firstrun = (lmm==ADAMS) && (miter==FUNC);
      if (firstrun) {
        cvode_mem = CVodeMalloc(P1_NEQ, f1, P1_T0, y, lmm, iter, ITOL,
                                &reltol, &abstol, NULL, ERRFP, OPTIN, iopt, ropt, machEnv);
        if (cvode_mem == NULL) { printf("CVodeMalloc failed."); return(1); }
      } else {
        flag = CVReInit(cvode_mem, f1, P1_T0, y, lmm, iter,
                        ITOL, &reltol, &abstol, NULL, ERRFP, OPTIN, iopt, ropt, machEnv);
        if (flag != SUCCESS) { printf("CVReInit failed."); return(1); }
      }
      
      flag = PrepareNextRun(cvode_mem, lmm, miter, 0, 0);     
      if (flag != SUCCESS) { printf("PrepareNextRun failed."); return(1); }
      
      printf("\n     t           x              xdot         qu     hu \n");
      
      for(iout=1, tout=P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT) {
        flag = CVode(cvode_mem, tout, y, &t, NORMAL);
        printf("%10.5f    %12.5e   %12.5e   %2ld    %6.4e\n",
               t, NV_Ith_S(y,0), NV_Ith_S(y,1), iopt[QU], ropt[HU]);
        if (flag != SUCCESS) {
          nerr++;
          printf("\n\n CVode returned error flag = %d \n\n", flag);
          break;
        }
        if (iout%2 == 0) {
          er = ABS(NV_Ith_S(y,0)) / abstol;
          if (er > ero) ero = er;
          if (er > P1_TOL_FACTOR) {
            nerr++;
            printf("\n\n Error exceeds %g * tolerance \n\n", P1_TOL_FACTOR);
          }
        }
      }
      
      PrintFinalStats(iopt, miter, ero);
    }
  }

  CVodeFree(cvode_mem);
  N_VFree(y);
  M_EnvFree_Serial(machEnv);

  return(nerr);
}

static void PrintIntro1(void)
{
  printf("Demonstration program for CVODE package - direct linear solvers\n");
  printf("\n\n");
  
  printf("Problem 1.. Van der Pol oscillator..\n");
  printf(" xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0\n");
  printf(" neq = %d,  itol = %s,  reltol = %.2g,  abstol = %.2g",
	 P1_NEQ, (ITOL == SS) ? "SS" : "SV", RTOL, ATOL);
}


static void f1(integertype N, realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y0, y1;
  
  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  NV_Ith_S(ydot,0) = y1;
  NV_Ith_S(ydot,1) = (1.0 - SQR(y0))* P1_ETA * y1 - y0;
} 

static void Jac1(integertype N, DenseMat J, RhsFn f, void *f_data, realtype tn,
                 N_Vector y, N_Vector fy, N_Vector ewt, realtype h, realtype uround,
                 void *jac_data, long int *nfePtr, N_Vector vtemp1,
                 N_Vector vtemp2, N_Vector vtemp3)
{ 
  realtype y0, y1;

  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  DENSE_ELEM(J,0,1) = 1.0;
  DENSE_ELEM(J,1,0) = -2.0 * P1_ETA * y0 * y1 - 1.0;
  DENSE_ELEM(J,1,1) = P1_ETA * (1.0 - SQR(y0));
}

static int Problem2(void)
{
  M_Env machEnv;
  realtype ropt[OPT_SIZE], reltol=RTOL, abstol=ATOL, t, tout, er, erm, ero;
  long int iopt[OPT_SIZE];
  int lmm, miter, iter, flag, iout, nerr=0;
  N_Vector y;
  void *cvode_mem = NULL;
  booleantype firstrun;
 
  machEnv = M_EnvInit_Serial(P2_NEQ);

  y = N_VNew(P2_NEQ, machEnv);
  PrintIntro2();
  
  for (lmm=ADAMS; lmm <= BDF; lmm++) {
    for (miter=FUNC; miter <= BAND_DQ; miter++) {     
      if ((miter==DENSE_USER) || (miter==DENSE_DQ)) continue;
      ero = 0.0;
      N_VConst(0.0, y);
      NV_Ith_S(y,0) = 1.0;
      iter = (miter == FUNC) ? FUNCTIONAL : NEWTON;

      firstrun = (lmm==ADAMS) && (miter==FUNC);

      if (firstrun) {
        cvode_mem = CVodeMalloc(P2_NEQ, f2, P2_T0, y, lmm, iter, ITOL,
                                &reltol, &abstol, NULL, ERRFP, OPTIN, 
                                iopt, ropt, machEnv);
        if (cvode_mem == NULL) { printf("CVodeMalloc failed."); continue; }
      } else {
        flag = CVReInit(cvode_mem, f2, P2_T0, y, lmm, iter, ITOL,
                        &reltol, &abstol, NULL, ERRFP, OPTIN, 
                        iopt, ropt, machEnv);
        if (flag != SUCCESS) { printf("CVReInit failed."); return(1); }
      }
      
      flag = PrepareNextRun(cvode_mem, lmm, miter, P2_MU, P2_ML);
      if (flag != SUCCESS) { printf("PrepareNextRun failed."); return(1); }
      
      printf("\n      t        max.err      qu     hu \n");
      
      for(iout=1, tout=P2_T1; iout <= P2_NOUT; iout++, tout*=P2_TOUT_MULT) {
        flag = CVode(cvode_mem, tout, y, &t, NORMAL);
        erm = MaxError(y, t);
        printf("%10.3f  %12.4e   %2ld   %12.4e\n", t, erm, iopt[QU], ropt[HU]);
        if (flag != SUCCESS) {
          nerr++;
          printf("\n\n CVode returned error flag = %d \n\n", flag);
          break;
        }
        er = erm / abstol;
        if (er > ero) ero = er;
        if (er > P2_TOL_FACTOR) {
          nerr++;
          printf("\n\n Error exceeds %g * tolerance \n\n", P2_TOL_FACTOR);
        }
      }
      
      PrintFinalStats(iopt, miter, ero);
    }
  }
  
  CVodeFree(cvode_mem);
  N_VFree(y);
  M_EnvFree_Serial(machEnv);

  return(nerr);
}

static void PrintIntro2(void)
{
  printf("\n\n-------------------------------------------------------------");
  printf("\n-------------------------------------------------------------");
  printf("\n\nProblem 2.. ydot = A * y, where A is a banded lower\n");
  printf("triangular matrix derived from 2-D advection PDE\n\n");
  printf(" neq = %d, ml = %d, mu = %d\n", P2_NEQ, P2_ML, P2_MU);
  printf(" itol = %s, reltol = %.2g, abstol = %.2g",
         (ITOL == SS) ? "SS" : "SV", RTOL, ATOL);
}


static void f2(integertype N, realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  integertype i, j, k;
  realtype d, *ydata, *dydata;
  
  ydata = NV_DATA_S(y);
  dydata = NV_DATA_S(ydot);

  /*
     Excluding boundaries, 

     ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
         i,j    i,j       i,j             i-1,j             i,j-1
  */

  for (j=0; j < P2_MESHY; j++) {
    for (i=0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      d = -2.0*ydata[k];
      if (i != 0) d += P2_ALPH1 * ydata[k-1];
      if (j != 0) d += P2_ALPH2 * ydata[k-P2_MESHX];
      dydata[k] = d;
    }
  }
}

static void Jac2(integertype N, integertype mu, integertype ml, BandMat J, RhsFn f,
                 void *f_data, realtype tn, N_Vector y, N_Vector fy, N_Vector ewt,
                 realtype h, realtype uround, void *jac_data, long int *nfePtr,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  integertype i, j, k;
  realtype *kthCol;

  /*
     The components of f(t,y) which depend on y    are
                                               i,j
     f    , f      , and f      : 
      i,j    i+1,j        i,j+1

     f    = -2 y    + alpha1 * y      + alpha2 * y
      i,j       i,j             i-1,j             i,j-1

     f      = -2 y      + alpha1 * y    + alpha2 * y
      i+1,j       i+1,j             i,j             i+1,j-1

     f      = -2 y      + alpha1 * y        + alpha2 * y
      i,j+1       i,j+1             i-1,j+1             i,j
  */

  for (j=0; j < P2_MESHY; j++) {
    for (i=0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      kthCol = BAND_COL(J,k);
      BAND_COL_ELEM(kthCol,k,k) = -2.0;
      if (i != P2_MESHX-1) BAND_COL_ELEM(kthCol,k+1,k) = P2_ALPH1;
      if (j != P2_MESHY-1) BAND_COL_ELEM(kthCol,k+P2_MESHX,k) = P2_ALPH2;
    }
  }
}

static realtype MaxError(N_Vector y, realtype t)
{
  integertype i, j, k;
  realtype *ydata, er, ex=0.0, yt, maxError=0.0, ifact_inv, jfact_inv=1.0;
  
  if (t == 0.0) return(0.0);

  ydata = NV_DATA_S(y);
  if (t <= 30.0) ex = exp(-2.0*t); 
  
  for (j = 0; j < P2_MESHY; j++) {
    ifact_inv = 1.0;
    for (i = 0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      yt = RPowerI(t,i+j) * ex * ifact_inv * jfact_inv;
      er = ABS(ydata[k] - yt);
      if (er > maxError) maxError = er;
      ifact_inv /= (i+1);
    }
    jfact_inv /= (j+1);
  }
  return(maxError);
}

static int PrepareNextRun(void *cvode_mem, int lmm, int miter, integertype mu,
                          integertype ml)
{
  int flag = SUCCESS;
  
  printf("\n\n-------------------------------------------------------------");
  
  printf("\n\nLinear Multistep Method : ");
  if (lmm == ADAMS) {
    printf("ADAMS\n");
  } else {
    printf("BDF\n");
  }
  
  printf("Iteration               : ");
  if (miter == FUNC) {
    printf("FUNCTIONAL\n");
  } else {
    printf("NEWTON\n");
    printf("Linear Solver           : ");
    switch(miter) {
    case DENSE_USER : 
      printf("Dense, User-Supplied Jacobian\n");
      flag = CVDense(cvode_mem, Jac1, NULL);
      break;
    case DENSE_DQ   : 
      printf("Dense, Difference Quotient Jacobian\n");
      flag = CVReInitDense(cvode_mem, NULL, NULL);
      break;
    case DIAG       : 
      printf("Diagonal Jacobian\n");
      flag = CVDiag(cvode_mem);
      break;
    case BAND_USER  : 
      printf("Band, User-Supplied Jacobian\n");
      flag = CVBand(cvode_mem, mu, ml, Jac2, NULL);
      break;
    case BAND_DQ  :   
      printf("Band, Difference Quotient Jacobian\n");
      flag = CVReInitBand(cvode_mem, mu, ml, NULL, NULL);
      break;    
    }
  }							
  return(flag);
  
}

static void PrintFinalStats(long int iopt[], int miter, realtype ero)
{
  int nje = 0, llrw = 0, lliw = 0;
  
  printf("\n Final statistics for this run..\n\n");
  printf(" CVode real workspace length              = %4ld \n", iopt[LENRW]);
  printf(" CVode integer workspace length           = %4ld \n", iopt[LENIW]);
  printf(" Number of steps                          = %4ld \n", iopt[NST]);
  printf(" Number of f-s                            = %4ld \n", iopt[NFE]);
  printf(" Number of setups                         = %4ld \n", iopt[NSETUPS]);
  printf(" Number of nonlinear iterations           = %4ld \n", iopt[NNI]);
  printf(" Number of nonlinear convergence failures = %4ld \n", iopt[NCFN]);
  printf(" Number of error test failures            = %4ld \n\n", iopt[NETF]);
  
  if (miter != FUNC) {
    switch(miter) {
    case DENSE_USER :
    case DENSE_DQ   : 
      nje  = DENSE_NJE;
      llrw = DENSE_LRW;
      lliw = DENSE_LIW;
      break;
    case BAND_USER  :
    case BAND_DQ    : 
      nje  = BAND_NJE;
      llrw = BAND_LRW;
      lliw = BAND_LIW;
      break;  
    case DIAG       : 
      nje = NSETUPS;
      llrw = DIAG_LRW;
      lliw = DIAG_LIW;
      break;
    }
    printf(" Linear solver real workspace length      = %4ld \n", iopt[llrw]);
    printf(" Linear solver integer workspace length   = %4ld \n", iopt[lliw]);
    printf(" Number of Jacobian evaluations           = %4ld \n\n", iopt[nje]);
  }
  
  printf(" Error overrun = %.3f \n", ero);
  
}







