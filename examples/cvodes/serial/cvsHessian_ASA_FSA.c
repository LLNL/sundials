/*
 * -----------------------------------------------------------------
 * $Revision: 4956 $
 * $Date: 2016-09-23 11:15:59 -0700 (Fri, 23 Sep 2016) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 *
 * Hessian through adjoint sensitivity example problem.
 *
 *        [ - p1 * y1^2 - y3 ]           [ 1 ]
 *   y' = [    - y2          ]    y(0) = [ 1 ]
 *        [ -p2^2 * y2 * y3  ]           [ 1 ]
 *
 *   p1 = 1.0
 *   p2 = 2.0
 *
 *           2
 *          /
 *   G(p) = |  0.5 * ( y1^2 + y2^2 + y3^2 ) dt
 *          /
 *          0
 *
 * Compute the gradient (ASA) and Hessian (FSA over ASA) of G(p).
 *
 * See D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

typedef struct {
  realtype p1, p2;
} *UserData;

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data);
static int fS(int Ns, realtype t,
              N_Vector y, N_Vector ydot,
              N_Vector *yS, N_Vector *ySdot,
              void *user_data,
              N_Vector tmp1, N_Vector tmp2);
static int fQS(int Ns, realtype t,
               N_Vector y, N_Vector *yS, 
               N_Vector yQdot, N_Vector *yQSdot,
               void *user_data,
               N_Vector tmp, N_Vector tmpQ);

static int fB1(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB);
static int fQB1(realtype t, N_Vector y, N_Vector *yS, 
                N_Vector yB, N_Vector qBdot, void *user_dataB);


static int fB2(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB);
static int fQB2(realtype t, N_Vector y, N_Vector *yS,
                N_Vector yB, N_Vector qBdot, void *user_dataB);

void PrintFwdStats(void *cvode_mem);
void PrintBckStats(void *cvode_mem, int idx);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  void *cvode_mem;

  long int Neq, Np2;
  int Np;

  realtype t0, tf;

  realtype reltol;
  realtype abstol, abstolQ, abstolB, abstolQB;

  N_Vector y, yQ;
  N_Vector *yS, *yQS;
  N_Vector yB1, yB2;
  N_Vector yQB1, yQB2;

  int steps, ncheck;
  int indexB1, indexB2;

  int flag;
  realtype time;

  realtype dp;
  realtype G, Gp, Gm;
  realtype grdG_fwd[2], grdG_bck[2], grdG_cntr[2];
  realtype H11, H22;

  /* User data structure */

  data = (UserData) malloc(sizeof *data);
  data->p1 = RCONST(1.0);
  data->p2 = RCONST(2.0);

  /* Problem size, integration interval, and tolerances */

  Neq = 3;
  Np  = 2;
  Np2 = 2*Np;

  t0 = 0.0;
  tf = 2.0;

  reltol = 1.0e-8;

  abstol = 1.0e-8;
  abstolQ = 1.0e-8;

  abstolB = 1.0e-8;
  abstolQB = 1.0e-8;

  /* Initializations for forward problem */

  y = N_VNew_Serial(Neq);
  N_VConst(ONE, y);

  yQ = N_VNew_Serial(1);
  N_VConst(ZERO, yQ);

  yS = N_VCloneVectorArray_Serial(Np, y);
  N_VConst(ZERO, yS[0]);
  N_VConst(ZERO, yS[1]);

  yQS = N_VCloneVectorArray_Serial(Np, yQ);
  N_VConst(ZERO, yQS[0]);
  N_VConst(ZERO, yQS[1]);

  /* Create and initialize forward problem */

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  flag = CVodeInit(cvode_mem, f, t0, y);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  flag = CVodeSetUserData(cvode_mem, data);

  flag = CVDense(cvode_mem, Neq);

  flag = CVodeQuadInit(cvode_mem, fQ, yQ);
  flag = CVodeQuadSStolerances(cvode_mem, reltol, abstolQ);
  flag = CVodeSetQuadErrCon(cvode_mem, TRUE);

  flag = CVodeSensInit(cvode_mem, Np, CV_SIMULTANEOUS, fS, yS);
  flag = CVodeSensEEtolerances(cvode_mem);
  flag = CVodeSetSensErrCon(cvode_mem, TRUE);

  flag = CVodeQuadSensInit(cvode_mem, fQS, yQS);

  flag = CVodeQuadSensEEtolerances(cvode_mem);
  flag = CVodeSetQuadSensErrCon(cvode_mem, TRUE);

  /* Initialize ASA */

  steps = 100;
  flag = CVodeAdjInit(cvode_mem, steps, CV_POLYNOMIAL);

  /* Forward integration */

  printf("-------------------\n");
  printf("Forward integration\n");
  printf("-------------------\n\n");

  flag = CVodeF(cvode_mem, tf, y, &time, CV_NORMAL, &ncheck);

  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  G = Ith(yQ,1);

  flag = CVodeGetSens(cvode_mem, &time, yS);

  flag = CVodeGetQuadSens(cvode_mem, &time, yQS);

  printf("ncheck = %d\n", ncheck);
  printf("\n");
  printf("     y:    %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:    %12.4e\n", Ith(yQ,1));
  printf("\n");
  printf("     yS1:  %12.4e %12.4e %12.4e\n", Ith(yS[0],1), Ith(yS[0],2), Ith(yS[0],3));
  printf("     yS2:  %12.4e %12.4e %12.4e\n", Ith(yS[1],1), Ith(yS[1],2), Ith(yS[1],3));
  printf("\n");
  printf("   dG/dp:  %12.4e %12.4e\n", Ith(yQS[0],1), Ith(yQS[1],1));
  printf("\n");

  printf("Final Statistics for forward pb.\n");
  printf("--------------------------------\n");
  PrintFwdStats(cvode_mem);


  /* Initializations for backward problems */

  yB1 = N_VNew_Serial(2*Neq);
  N_VConst(ZERO, yB1);

  yQB1 = N_VNew_Serial(Np2);
  N_VConst(ZERO, yQB1);

  yB2 = N_VNew_Serial(2*Neq);
  N_VConst(ZERO, yB2);

  yQB2 = N_VNew_Serial(Np2);
  N_VConst(ZERO, yQB2);

  /* Create and initialize backward problems (one for each column of the Hessian) */

  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB1);
  flag = CVodeInitBS(cvode_mem, indexB1, fB1, tf, yB1);
  flag = CVodeSStolerancesB(cvode_mem, indexB1, reltol, abstolB);
  flag = CVodeSetUserDataB(cvode_mem, indexB1, data);
  flag = CVodeQuadInitBS(cvode_mem, indexB1, fQB1, yQB1);
  flag = CVodeQuadSStolerancesB(cvode_mem, indexB1, reltol, abstolQB);
  flag = CVodeSetQuadErrConB(cvode_mem, indexB1, TRUE);
  flag = CVDenseB(cvode_mem, indexB1, 2*Neq);

  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB2);
  flag = CVodeInitBS(cvode_mem, indexB2, fB2, tf, yB2);
  flag = CVodeSStolerancesB(cvode_mem, indexB2, reltol, abstolB);
  flag = CVodeSetUserDataB(cvode_mem, indexB2, data);
  flag = CVodeQuadInitBS(cvode_mem, indexB2, fQB2, yQB2);
  flag = CVodeQuadSStolerancesB(cvode_mem, indexB2, reltol, abstolQB);
  flag = CVodeSetQuadErrConB(cvode_mem, indexB2, TRUE);
  flag = CVDenseB(cvode_mem, indexB2, 2*Neq);

  /* Backward integration */

  printf("---------------------------------------------\n");
  printf("Backward integration ... (2 adjoint problems)\n");
  printf("---------------------------------------------\n\n");

  flag = CVodeB(cvode_mem, t0, CV_NORMAL);

  flag = CVodeGetB(cvode_mem, indexB1, &time, yB1);
  flag = CVodeGetQuadB(cvode_mem, indexB1, &time, yQB1);

  flag = CVodeGetB(cvode_mem, indexB2, &time, yB2);
  flag = CVodeGetQuadB(cvode_mem, indexB2, &time, yQB2);

  printf("   dG/dp:  %12.4e %12.4e   (from backward pb. 1)\n", -Ith(yQB1,1), -Ith(yQB1,2));
  printf("           %12.4e %12.4e   (from backward pb. 2)\n", -Ith(yQB2,1), -Ith(yQB2,2));
  printf("\n");
  printf("   H = d2G/dp2:\n");
  printf("        (1)            (2)\n");
  printf("  %12.4e   %12.4e\n", -Ith(yQB1,3) , -Ith(yQB2,3));
  printf("  %12.4e   %12.4e\n", -Ith(yQB1,4) , -Ith(yQB2,4));
  printf("\n");

  printf("Final Statistics for backward pb. 1\n");
  printf("-----------------------------------\n");
  PrintBckStats(cvode_mem, indexB1);

  printf("Final Statistics for backward pb. 2\n");
  printf("-----------------------------------\n");
  PrintBckStats(cvode_mem, indexB2);

  /* Free CVODES memory */

  CVodeFree(&cvode_mem);

  /* Finite difference tests */

  dp = RCONST(1.0e-2);

  printf("-----------------------\n");
  printf("Finite Difference tests\n");
  printf("-----------------------\n\n");

  printf("del_p = %g\n\n",dp);

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);
  flag = CVodeInit(cvode_mem, f, t0, y);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  flag = CVodeSetUserData(cvode_mem, data);
  flag = CVDense(cvode_mem, Neq);
  flag = CVodeQuadInit(cvode_mem, fQ, yQ);
  flag = CVodeQuadSStolerances(cvode_mem, reltol, abstolQ);
  flag = CVodeSetQuadErrCon(cvode_mem, TRUE);

  data->p1 += dp;

  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  Gp = Ith(yQ,1);
  printf("p1+  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));

  data->p1 -= 2.0*dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);
  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);
  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  Gm = Ith(yQ,1);
  printf("p1-  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));
 
  data->p1 += dp;

  grdG_fwd[0] = (Gp-G)/dp;
  grdG_bck[0] = (G-Gm)/dp;
  grdG_cntr[0] = (Gp-Gm)/(2.0*dp);
  H11 = (Gp - 2.0*G + Gm) / (dp*dp);

  data->p2 += dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);
  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);
  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  Gp = Ith(yQ,1);
  printf("p2+  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));
 
  data->p2 -= 2.0*dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);
  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);
  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  Gm = Ith(yQ,1);
  printf("p2-  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));

  data->p2 += dp;

  grdG_fwd[1] = (Gp-G)/dp;
  grdG_bck[1] = (G-Gm)/dp;
  grdG_cntr[1] = (Gp-Gm)/(2.0*dp);
  H22 = (Gp - 2.0*G + Gm) / (dp*dp);

  printf("\n");

  printf("   dG/dp:  %12.4e %12.4e   (fwd FD)\n", grdG_fwd[0], grdG_fwd[1]);
  printf("           %12.4e %12.4e   (bck FD)\n", grdG_bck[0], grdG_bck[1]);
  printf("           %12.4e %12.4e   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4e\n", H11);
  printf("  H(2,2):  %12.4e\n", H22);


  /* Free memory */

  CVodeFree(&cvode_mem);

  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yQ);

  N_VDestroyVectorArray_Serial(yS, Np);
  N_VDestroyVectorArray_Serial(yQS, Np);

  N_VDestroy_Serial(yB1);
  N_VDestroy_Serial(yQB1);
  N_VDestroy_Serial(yB2);
  N_VDestroy_Serial(yQB2);

  free(data);

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */


static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2;

  data = (UserData) user_data;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  Ith(ydot,1) = -p1*y1*y1 - y3;
  Ith(ydot,2) = -y2;
  Ith(ydot,3) = -p2*p2*y2*y3;

  return(0);
}

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data)
{
  realtype y1, y2, y3;

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  Ith(qdot,1) = 0.5 * ( y1*y1 + y2*y2 + y3*y3 );

  return(0);
}

static int fS(int Ns, realtype t,
              N_Vector y, N_Vector ydot,
              N_Vector *yS, N_Vector *ySdot,
              void *user_data,
              N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype fys1, fys2, fys3;
  realtype p1, p2;

  data = (UserData) user_data;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  /* 1st sensitivity RHS */

  s1 = Ith(yS[0],1);
  s2 = Ith(yS[0],2);
  s3 = Ith(yS[0],3);

  fys1 = - 2.0*p1*y1 * s1 - s3;
  fys2 = - s2;
  fys3 = - p2*p2*y3 * s2 - p2*p2*y2 * s3;

  Ith(ySdot[0],1) = fys1 - y1*y1;
  Ith(ySdot[0],2) = fys2;
  Ith(ySdot[0],3) = fys3;

  /* 2nd sensitivity RHS */

  s1 = Ith(yS[1],1);
  s2 = Ith(yS[1],2);
  s3 = Ith(yS[1],3);

  fys1 = - 2.0*p1*y1 * s1 - s3;
  fys2 = - s2;
  fys3 = - p2*p2*y3 * s2 - p2*p2*y2 * s3;

  Ith(ySdot[1],1) = fys1;
  Ith(ySdot[1],2) = fys2;
  Ith(ySdot[1],3) = fys3 - 2.0*p2*y2*y3;

  return(0);
}

static int fQS(int Ns, realtype t,
               N_Vector y, N_Vector *yS, 
               N_Vector yQdot, N_Vector *yQSdot,
               void *user_data,
               N_Vector tmp, N_Vector tmpQ)
{
  realtype y1, y2, y3;
  realtype s1, s2, s3;

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);


  /* 1st sensitivity RHS */

  s1 = Ith(yS[0],1);
  s2 = Ith(yS[0],2);
  s3 = Ith(yS[0],3);

  Ith(yQSdot[0],1) = y1*s1 + y2*s2 + y3*s3;


  /* 1st sensitivity RHS */

  s1 = Ith(yS[1],1);
  s2 = Ith(yS[1],2);
  s3 = Ith(yS[1],3);

  Ith(yQSdot[1],1) = y1*s1 + y2*s2 + y3*s3;

  return(0);
}

static int fB1(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 1 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */
  
  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  s1 = Ith(yS[0],1); 
  s2 = Ith(yS[0],2); 
  s3 = Ith(yS[0],3);

  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  
  Ith(yBdot,1) = 2.0*p1*y1 * l1     - y1;
  Ith(yBdot,2) = l2 + p2*p2*y3 * l3 - y2;
  Ith(yBdot,3) = l1 + p2*p2*y2 * l3 - y3;

  Ith(yBdot,4) = 2.0*p1*y1 * m1     + l1 * 2.0*(y1 + p1*s1) - s1;
  Ith(yBdot,5) = m2 + p2*p2*y3 * m3 + l3 * p2*p2*s3         - s2;
  Ith(yBdot,6) = m1 + p2*p2*y2 * m3 + l3 * p2*p2*s2         - s3;

  return(0);
}

static int fQB1(realtype t, N_Vector y, N_Vector *yS,
                N_Vector yB, N_Vector qBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 1 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */
  
  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  s1 = Ith(yS[0],1); 
  s2 = Ith(yS[0],2); 
  s3 = Ith(yS[0],3);
  
  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  Ith(qBdot,1) = -y1*y1 * l1;
  Ith(qBdot,2) = -2.0*p2*y2*y3 * l3;

  Ith(qBdot,3) = -y1*y1 * m1        - l1 * 2.0*y1*s1;
  Ith(qBdot,4) = -2.0*p2*y2*y3 * m3 - l3 * 2.0*(p2*y3*s2 + p2*y2*s3);

  return(0);
}




static int fB2(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 2 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */

  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  s1 = Ith(yS[1],1); 
  s2 = Ith(yS[1],2); 
  s3 = Ith(yS[1],3);
  
  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  Ith(yBdot,1) = 2.0*p1*y1 * l1     - y1;
  Ith(yBdot,2) = l2 + p2*p2*y3 * l3 - y2;
  Ith(yBdot,3) = l1 + p2*p2*y2 * l3 - y3;

  Ith(yBdot,4) = 2.0*p1*y1 * m1     + l1 * 2.0*p1*s1              - s1;
  Ith(yBdot,5) = m2 + p2*p2*y3 * m3 + l3 * (2.0*p2*y3 + p2*p2*s3) - s2;
  Ith(yBdot,6) = m1 + p2*p2*y2 * m3 + l3 * (2.0*p2*y2 + p2*p2*s2) - s3;


  return(0);
}


static int fQB2(realtype t, N_Vector y, N_Vector *yS,
                N_Vector yB, N_Vector qBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 2 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */
  
  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  s1 = Ith(yS[1],1); 
  s2 = Ith(yS[1],2); 
  s3 = Ith(yS[1],3);

  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  Ith(qBdot,1) = -y1*y1 * l1;
  Ith(qBdot,2) = -2.0*p2*y2*y3 * l3;

  Ith(qBdot,3) = -y1*y1 * m1        - l1 * 2.0*y1*s1;
  Ith(qBdot,4) = -2.0*p2*y2*y3 * m3 - l3 * 2.0*(p2*y3*s2 + p2*y2*s3 + y2*y3);

  return(0);
}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

void PrintFwdStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nfQSe, netfQS;

  int qlast, qcur;
  realtype h0u, hlast, hcur, tcur;

  int flag;


  flag = CVodeGetIntegratorStats(cvode_mem, &nst, &nfe, &nsetups, &netf, 
                                 &qlast, &qcur,
                                 &h0u, &hlast, &hcur,
                                 &tcur);

  flag = CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn);

  flag = CVodeGetQuadStats(cvode_mem, &nfQe, &netfQ);

  flag = CVodeGetSensStats(cvode_mem, &nfSe, &nfeS, &netfS, &nsetupsS);

  flag = CVodeGetSensNonlinSolvStats(cvode_mem, &nniS, &ncfnS);

  flag = CVodeGetQuadSensStats(cvode_mem, &nfQSe, &netfQS);


  printf(" Number steps: %5ld\n\n", nst);
  printf(" Function evaluations:\n");
  printf("  f:        %5ld\n  fQ:       %5ld\n  fS:       %5ld\n  fQS:      %5ld\n",
         nfe, nfQe, nfSe, nfQSe);
  printf(" Error test failures:\n");
  printf("  netf:     %5ld\n  netfQ:    %5ld\n  netfS:    %5ld\n  netfQS:   %5ld\n",
         netf, netfQ, netfS, netfQS);
  printf(" Linear solver setups:\n");
  printf("  nsetups:  %5ld\n  nsetupsS: %5ld\n", nsetups, nsetupsS);
  printf(" Nonlinear iterations:\n");
  printf("  nni:      %5ld\n  nniS:     %5ld\n", nni, nniS);
  printf(" Convergence failures:\n");
  printf("  ncfn:     %5ld\n  ncfnS:    %5ld\n", ncfn, ncfnS);

  printf("\n");

}


void PrintBckStats(void *cvode_mem, int idx)
{
  void *cvode_mem_bck;

  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;

  int qlast, qcur;
  realtype h0u, hlast, hcur, tcur;

  int flag;

  cvode_mem_bck = CVodeGetAdjCVodeBmem(cvode_mem, idx);

  flag = CVodeGetIntegratorStats(cvode_mem_bck, &nst, &nfe, &nsetups, &netf, 
                                 &qlast, &qcur,
                                 &h0u, &hlast, &hcur,
                                 &tcur);

  flag = CVodeGetNonlinSolvStats(cvode_mem_bck, &nni, &ncfn);

  flag = CVodeGetQuadStats(cvode_mem_bck, &nfQe, &netfQ);

  printf(" Number steps: %5ld\n\n", nst);
  printf(" Function evaluations:\n");
  printf("  f:        %5ld\n  fQ:       %5ld\n", nfe, nfQe);
  printf(" Error test failures:\n");
  printf("  netf:     %5ld\n  netfQ:    %5ld\n", netf, netfQ);
  printf(" Linear solver setups:\n");
  printf("  nsetups:  %5ld\n", nsetups);
  printf(" Nonlinear iterations:\n");
  printf("  nni:      %5ld\n", nni);
  printf(" Convergence failures:\n");
  printf("  ncfn:     %5ld\n", ncfn);

  printf("\n");


}
