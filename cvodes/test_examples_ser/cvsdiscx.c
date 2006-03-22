/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-03-22 00:28:14 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Simple 1D example to illustrate integrating over discontinuities:
 * 
 * A) Discontinuity in solution
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       y' = -y   ; y(1) = 1    ; t = [1,2]
 *
 * B) Discontinuity in RHS (y')
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       z' = -5*z ; z(1) = y(1) ; t = [1,2]
 *    This case is solved twice, first by explicitly treating the 
 *    discontinuity point and secondly by letting the integrator 
 *    deal with the discontinuity.
 *
 * NOTE: For readibility, no checks are performed on the various 
 *       function return flags.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include "cvodes.h"
#include "nvector_serial.h"
#include "cvodes_dense.h"

#define Ith(v,i) NV_Ith_S(v,i-1)

#define RHS1 1
#define RHS2 2

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int main()
{
  void *cvode_mem;
  N_Vector y;
  int neq, flag;
  realtype reltol, abstol, t0, t1, t2, t;
  long int nst1, nst2, nst;

  neq = 1;

  reltol = RCONST(1.0e-3);  
  abstol = RCONST(1.0e-4);

  t0 = RCONST(0.0);
  t1 = RCONST(1.0);
  t2 = RCONST(2.0);

  y = N_VNew_Serial(neq);
  
  /*
   * ---------------------------------------------------------------
   * Discontinuity in solution 
   * Note that it is not required to set TSTOP. 
   * ---------------------------------------------------------------
   */

  printf("\nDiscontinuity in solution\n\n");

  /* Initialize solver */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  Ith(y,1) = 1.0;
  CVodeMalloc(cvode_mem, f, t0, y, CV_SS, reltol, &abstol);
  CVodeSetFdata(cvode_mem, &flag);
  CVDense(cvode_mem, neq);

  /* Integrate to the point of discontinuity */
  CVodeSetStopTime(cvode_mem, t1);
  flag = RHS1;
  t = t0;
  printf("%12.8e  %12.8e\n",t,Ith(y,1));
  while (t<t1) {
    CVode(cvode_mem, t1, y, &t, CV_ONE_STEP_TSTOP);
    printf("%12.8e  %12.8e\n",t,Ith(y,1));
  }
  CVodeGetNumSteps(cvode_mem, &nst1);

  /* Include discontinuity and reinitialize solver */
  Ith(y,1) = 1.0;
  CVodeReInit(cvode_mem, f, t1, y, CV_SS, reltol, &abstol); 

  /* Integrate from discontinuity */
  CVodeSetStopTime(cvode_mem, t2);
  flag = RHS1;
  t = t1;
  printf("%12.8e  %12.8e\n",t,Ith(y,1));
  while (t<t2) {
    CVode(cvode_mem, t2, y, &t, CV_ONE_STEP_TSTOP);
    printf("%12.8e  %12.8e\n",t,Ith(y,1));
  }
  CVodeGetNumSteps(cvode_mem, &nst2);
  
  /* Print statistics and free solver memory */
  nst = nst1 + nst2;
  printf("\nNumber of steps: %d + %d = %d\n",nst1, nst2, nst);

  CVodeFree(&cvode_mem);

  /*
   * ---------------------------------------------------------------
   * Discontinuity in RHS: Case 1 - explicit treatment 
   * Note that it is not required to set TSTOP, but without it
   * we would have to find y(t1) to reinitialize the solver. 
   * ---------------------------------------------------------------
   */

  printf("\nDiscontinuity in RHS: Case 1 - explicit treatment\n\n");

  /* Initialize solver */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  Ith(y,1) = 1.0;
  CVodeMalloc(cvode_mem, f, t0, y, CV_SS, reltol, &abstol);
  CVodeSetFdata(cvode_mem, &flag);
  CVDense(cvode_mem, neq);

  /* Integrate to the point of discontinuity */
  CVodeSetStopTime(cvode_mem, t1);
  flag = RHS1;
  t = t0;
  printf("%12.8e  %12.8e\n",t,Ith(y,1));
  while (t<t1) {
    CVode(cvode_mem, t1, y, &t, CV_ONE_STEP_TSTOP);
    printf("%12.8e  %12.8e\n",t,Ith(y,1));
  }
  CVodeGetNumSteps(cvode_mem, &nst1);

  /* If TSTOP was not set, we'd need to find y(t1): */ 
  /* CVodeGetDky(cvode_mem, t1, 0, y); */

  /* Reinitialize solver */
  CVodeReInit(cvode_mem, f, t1, y, CV_SS, reltol, &abstol); 

  /* Integrate from discontinuity with new RHS */
  CVodeSetStopTime(cvode_mem, t2);
  flag = RHS2;
  t = t1;
  printf("%12.8e  %12.8e\n",t,Ith(y,1));
  while (t<t2) {
    CVode(cvode_mem, t2, y, &t, CV_ONE_STEP_TSTOP);
    printf("%12.8e  %12.8e\n",t,Ith(y,1));
  }
  CVodeGetNumSteps(cvode_mem, &nst2);

  /* Print statistics and free solver memory */
  nst = nst1 + nst2;
  printf("\nNumber of steps: %d + %d = %d\n",nst1, nst2, nst);

  CVodeFree(&cvode_mem);
  
  /*
   * ---------------------------------------------------------------
   * Discontinuity in RHS: Case 2 - let CVODES deal with it 
   * Note that here we MUST set TSTOP to ensure that the 
   * change in the RHS happens at the appropriate time 
   * ---------------------------------------------------------------
   */

  printf("\nDiscontinuity in RHS: Case 2 - let CVODES deal with it\n\n");

  /* Initialize solver */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  Ith(y,1) = 1.0;
  CVodeMalloc(cvode_mem, f, t0, y, CV_SS, reltol, &abstol);
  CVodeSetFdata(cvode_mem, &flag);
  CVDense(cvode_mem, neq);

  /* Integrate to the point of discontinuity */
  CVodeSetStopTime(cvode_mem, t1);
  flag = RHS1;
  t = t0;
  printf("%12.8e  %12.8e\n",t,Ith(y,1));
  while (t<t1) {
    CVode(cvode_mem, t1, y, &t, CV_ONE_STEP_TSTOP);
    printf("%12.8e  %12.8e\n",t,Ith(y,1));
  }
  CVodeGetNumSteps(cvode_mem, &nst1);

  /* Continue integration from discontinuity with new RHS */
  CVodeSetStopTime(cvode_mem, t2);
  flag = RHS2;
  t = t1;
  printf("%12.8e  %12.8e\n",t,Ith(y,1));
  while (t<t2) {
    CVode(cvode_mem, t2, y, &t, CV_ONE_STEP_TSTOP);
    printf("%12.8e  %12.8e\n",t,Ith(y,1));
  }
  CVodeGetNumSteps(cvode_mem, &nst);

  /* Print statistics and free solver memory */
  nst2 = nst - nst1;
  printf("\nNumber of steps: %d + %d = %d\n",nst1, nst2, nst);

  CVodeFree(&cvode_mem);


  /* Free y vector and return */
  N_VDestroy(y);
  return(0);
}

/*
 * RHS function
 * The form of the RHS function is controlled by the flag passed as f_data:
 *   flag = RHS1 -> y' = -y
 *   flag = RHS2 -> y' = -5*y
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  int *flag;

  flag = (int *) f_data;

  switch(*flag) {
  case RHS1:
    Ith(ydot,1) = -Ith(y,1);
    break;
  case RHS2:
    Ith(ydot,1) = -5.0*Ith(y,1);
    break;
  }

  return(0);
}

