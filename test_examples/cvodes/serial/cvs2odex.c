/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Modification of cvsdenx demonstrating how to integrate two 
 * different but closely related ODE systems at once.
 *
 * First system:
 *
 *  y1' = -.04*y1 + 1.e4*y2*y3
 *  y2' = 3.e7*(y2)^2
 *  y3' = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *
 * y1(0)=1 y2(0)=y3(0)=0
 *
 *
 * Second system:
 *  y1' = -.04*y1 + 1.e4*y2*y3
 *  y2' = 3.e7*(y2)^2
 *
 * where y3 = -y1-y2+1
 *
 * y1(0)=1 y2(0)=0
 *
 * NOTE: For readibility, no checks are performed on the various 
 *       function return flags.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>

#define Ith(v,i) NV_Ith_S(v,i-1)

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int main()
{
  int neq1, neq2;
  int flag1, flag2;
  N_Vector y1, y2, abstol1, abstol2;
  void *cvode_mem1, *cvode_mem2;
  realtype reltol, t0, t1, t;

  reltol = 1.0e-4;  
  t0 = 0.0;
  t1 = 0.4;

  flag1 = 1;
  flag2 = 2;

  neq1 = 3;
  y1 = N_VNew_Serial(neq1);
  abstol1 = N_VNew_Serial(neq1); 
  Ith(y1,1) = 1.0;  Ith(y1,2) = 0.0;  Ith(y1,3) = 0.0;
  Ith(abstol1,1) = 1.0e-8;  Ith(abstol1,2) = 1.0e-6;  Ith(abstol1,3) = 1.0e-14;

  neq2 = 2;
  y2 = N_VNew_Serial(neq2);
  abstol2 = N_VNew_Serial(neq2); 
  Ith(y2,1) = 1.0;  Ith(y2,2) = 0.0;
  Ith(abstol2,1) = 1.0e-8;  Ith(abstol2,2) = 1.0e-6;
  
  cvode_mem1 = CVodeCreate(CV_BDF, CV_NEWTON);
  cvode_mem2 = CVodeCreate(CV_BDF, CV_NEWTON);
  
  CVodeSetUserData(cvode_mem1, &flag1);
  CVodeSetUserData(cvode_mem2, &flag2);

  CVodeInit(cvode_mem1, f, t0, y1);
  CVodeSVtolerances(cvode_mem1, reltol, abstol1);
  CVodeInit(cvode_mem2, f, t0, y2);
  CVodeSVtolerances(cvode_mem2, reltol, abstol2);
  
  CVDense(cvode_mem1, neq1);
  CVDense(cvode_mem2, neq2);

  CVode(cvode_mem1, t1, y1, &t, CV_NORMAL);
  printf("t = %g  y1=%g  y2=%g  y3=%g\n",t, Ith(y1,1), Ith(y1,2), Ith(y1,3));

  CVode(cvode_mem2, t1, y2, &t, CV_NORMAL);
  printf("t = %g  y1=%g  y2=%g\n",t, Ith(y2,1), Ith(y2,2));
  
  N_VDestroy_Serial(y1);
  N_VDestroy_Serial(abstol1);
  CVodeFree(&cvode_mem1);

  N_VDestroy_Serial(y2);
  N_VDestroy_Serial(abstol2);
  CVodeFree(&cvode_mem2);

  return(0);
}


static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  int *flag;

  realtype y1, y2, y3, yd1, yd2;

  flag = (int *) f_data;

    
  y1 = Ith(y,1); y2 = Ith(y,2); 

  if (*flag == 1)
    y3 = Ith(y,3);
  else
    y3 = -y1-y2+1.0;
  
  yd1 = Ith(ydot,1) = -0.04*y1 + 1.0e4*y2*y3;
  yd2 = Ith(ydot,2) = 3.0e7*y3*y3;

  if (*flag == 1)
    Ith(ydot,3) = -yd1-yd2;


  return(0);
}

