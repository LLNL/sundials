/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2007-08-01 01:26:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal 
 * force applied to the connecting rod and a spring-damper connecting 
 * the crank and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 *
 * IDAS also computes sensitivities with respect to the problem 
 * parameters k (spring constant) and c (damper constant) of the 
 * kinetic energy:
 *   E = int_t0^tend g(t,y,p) dt, 
 * where
 *   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 
 *              0.5*m1*a^2/4*v1^2 + 0.5*m2*v2^2
 *              
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>

#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>


#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i= 1..NEQ */

/* Problem Constants */

#define NEQ   10
#define NP     2

#define TBEGIN  RCONST(0.0)
#define TEND    RCONST(5.000)

#define NOUT  40

#define RTOLF   RCONST(1.0e-06)
#define ATOLF   RCONST(1.0e-07)

#define RTOLQ   RCONST(1.0e-06)
#define ATOLQ   RCONST(1.0e-08)

#define RTOLB   RCONST(1.0e-06)
#define ATOLB   RCONST(1.0e-08)

#define RTOLFD  RCONST(1.0e-06)
#define ATOLFD  RCONST(1.0e-08)


#define ZERO     RCONST(0.00)
#define QUARTER  RCONST(0.25)
#define HALF     RCONST(0.50)
#define ONE      RCONST(1.00)
#define TWO      RCONST(2.00)
#define FOUR     RCONST(4.00)

typedef struct {
  realtype a;
  realtype J1, J2, m1, m2;
  realtype l0;
  realtype params[2];
  realtype F;
} *UserData;

static int ressc(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data);
static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data);

static int rhsQS(int Ns, realtype t, N_Vector yy, N_Vector yp, 
                 N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsQS,
                 void *user_data,  N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS);

static int resscB(realtype tt, 
                  N_Vector yy, N_Vector yp,
                  N_Vector yyB, N_Vector ypB, N_Vector rrB,
                  void *user_dataB);

static int rhsQB(realtype tt, 
                 N_Vector yy, N_Vector yp, 
                 N_Vector yyB, N_Vector ypB, 
                 N_Vector rrQB, void *user_dataB);


static void setIC(N_Vector yy, N_Vector yp, UserData data);
static void setICB(N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, UserData data);
static void force(N_Vector yy, realtype *Q, UserData data);

static void PrintOutput(void *mem, realtype t, N_Vector y);
static void PrintFinalStats(void *mem);
static void PrintFinalStatsB(void *mem, int bckPb);
static int check_flag(void *flagvalue, char *funcname, int opt);
/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  UserData data;

  void *mem;
  N_Vector yy, yp, yyB, ypB, id, q, *yyS, *ypS, *qS, qB;
  realtype rtol, atol, *yval, *ypval;
  realtype t0, tf, tout, dt, tret;
  realtype pbar[2];
  realtype dp, G, Gm[2], Gp[2];
  int flag, iout, is;
  int indexB, ncheck;
  realtype atolS[NP];

  id = N_VNew_Serial(NEQ);
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  q = N_VNew_Serial(1);

  yyS= N_VCloneVectorArray(NP,yy);
  ypS= N_VCloneVectorArray(NP,yp);
  qS = N_VCloneVectorArray_Serial(NP, q);

  data = (UserData) malloc(sizeof *data);

  data->a = 0.5;   /* half-length of crank */
  data->J1 = 1.0;  /* crank moment of inertia */
  data->m2 = 1.0;  /* mass of connecting rod */
  data->m1 = 1.0;
  data->J2 = 2.0;  /* moment of inertia of connecting rod */
  data->params[0] = 1.0;   /* spring constant */
  data->params[1] = 1.0;   /* damper constant */
  data->l0 = 1.0;  /* spring free length */
  data->F = 1.0;   /* external constant force */

  N_VConst(ONE, id);
  NV_Ith_S(id, 9) = ZERO;
  NV_Ith_S(id, 8) = ZERO;
  NV_Ith_S(id, 7) = ZERO;
  NV_Ith_S(id, 6) = ZERO;
  
  printf("\n");

  /* Consistent IC*/
  setIC(yy, yp, data);

  for (is=0;is<NP;is++) {
    N_VConst(ZERO, yyS[is]);
    N_VConst(ZERO, ypS[is]);
  }

  /* IDA initialization */
  mem = IDACreate();
  flag = IDAInit(mem, ressc, TBEGIN, yy, yp);
  flag = IDASStolerances(mem, RTOLF, ATOLF);
  flag = IDASetUserData(mem, data);
  flag = IDASetId(mem, id);
  flag = IDASetSuppressAlg(mem, TRUE);
  flag = IDASetMaxNumSteps(mem, 20000);

  /* Call IDADense and set up the linear solver. */
  flag = IDADense(mem, NEQ);

  //flag = IDASensInit(mem, NP, IDA_STAGGERED, NULL, yyS, ypS);
  flag = IDASensInit(mem, NP, IDA_SIMULTANEOUS, NULL, yyS, ypS);
  pbar[0] = data->params[0];pbar[1] = data->params[1];
  flag = IDASetSensParams(mem, data->params, pbar, NULL);
  flag = IDASensEEtolerances(mem);
  IDASetSensErrCon(mem, TRUE);
  
  N_VConst(ZERO, q);
  flag = IDAQuadInit(mem, rhsQ, q);
  flag = IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  flag = IDASetQuadErrCon(mem, TRUE);
  
  N_VConst(ZERO, qS[0]);
  flag = IDAQuadSensInit(mem, rhsQS, qS);
  atolS[0] = atolS[1] = ATOLQ;
  flag = IDAQuadSensSStolerances(mem, RTOLQ, atolS);
  flag = IDASetQuadSensErrCon(mem, TRUE);  
  
  //PrintOutput(mem,TBEGIN,yy);

  flag = IDAAdjInit(mem, 100, IDA_HERMITE);
  //flag = IDAAdjInit(mem, 1000, IDA_POLYNOMIAL);

  /* Perform forward run */
  printf("Forward integration ... ");

  flag = IDASolveF(mem, TEND, &tret, yy, yp, IDA_NORMAL, &ncheck);
  //flag = IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  if (check_flag(&flag, "IDASolve", 1)) return(1);

  //PrintOutput(mem,tret,yy);

  PrintFinalStats(mem);

  IDAGetQuad(mem, &tret, q);
  printf("--------------------------------------------\n");
  printf("  G = %24.16f\n", Ith(q,1));
  printf("--------------------------------------------\n\n");
  
  IDAGetQuadSens(mem, &tret, qS);
  printf("-------------F O R W A R D------------------\n");
  printf("   dG/dp:  %12.4le %12.4le\n", Ith(qS[0],1), Ith(qS[1],1));
  printf("--------------------------------------------\n\n");


  /***********************************************************
   *****************   A D J O I N T   ***********************
   *********************************************************/

  /* Integrate backward problem. */
  yyB = N_VNew_Serial(NEQ);
  ypB = N_VNew_Serial(NEQ);

  setICB(yy, yp, yyB, ypB, data);

  qB = N_VNew_Serial(NP);
  N_VConst(ZERO, qB);

  N_Vector resB;
  resB = N_VNew_Serial(NEQ);
  
  resscB(TEND, yy, yp, yyB, ypB, resB, data);
  printf("norm of IC for backward problem: %12.4e\n", N_VMaxNorm_Serial(resB));
  //N_VPrint_Serial(resB);
  N_VDestroy_Serial(resB);

  flag = IDACreateB(mem, &indexB);
  if (check_flag(&flag, "IDACreateB", 1)) return(1);

  flag = IDAInitB(mem, indexB, resscB, TEND, yyB, ypB);
  if (check_flag(&flag, "IDAInitB", 1)) return(1);

  flag = IDASStolerancesB(mem, indexB, RTOLB, ATOLB);
  if (check_flag(&flag, "IDASStolerancesB", 1)) return(1);

  flag = IDASetUserDataB(mem, indexB, data);
  if (check_flag(&flag, "IDASetUserDataB", 1)) return(1);

  flag = IDASetMaxNumStepsB(mem, indexB, 10000);


  IDASetIdB(mem, indexB, id);
  IDASetSuppressAlgB(mem, indexB, TRUE);

  flag = IDAQuadInitB(mem, indexB, rhsQB, qB);
  IDAQuadSStolerancesB(mem, indexB, RTOLB, ATOLB);
  IDASetQuadErrConB(mem, indexB, TRUE);

  flag = IDADenseB(mem, indexB, NEQ);
  if (check_flag(&flag, "IDADenseB", 1)) return(1);

  printf("Backward integration ...\n");
  flag = IDASolveB(mem, TBEGIN, IDA_NORMAL);
  if (check_flag(&flag, "IDASolveB", 1)) return(1);

  IDAGetDky(IDAGetAdjIDABmem(mem, indexB), TBEGIN, 0, yyB);
  IDAGetDky(IDAGetAdjIDABmem(mem, indexB), TBEGIN, 1, ypB);
  
  setIC(yy, yp, data);
  
  resscB(TBEGIN, yy, yp, yyB, ypB, resB, data);
  printf("Norm of residual for backward pb sln:%12.4e\n\n", N_VMaxNorm_Serial(resB));
  //N_VPrint_Serial(resB);

  IDAGetQuadB(mem, indexB, &tret, qB);
  printf("-------------A D J O I N T-----------------\n");
  printf("   dG/dp:  %12.4le %12.4le\n", Ith(qB,1), Ith(qB,2));
  printf("-------------------------------------------\n\n");


  PrintFinalStatsB(mem,indexB);

  IDAFree(&mem);



  /* Finite differences for dG/dp */
  dp = 0.00001;
  data->params[0] = ONE;
  data->params[1] = ONE;

  mem = IDACreate();

  setIC(yy, yp, data);
  flag = IDAInit(mem, ressc, TBEGIN, yy, yp);
  flag = IDASStolerances(mem, RTOLFD, ATOLFD);
  flag = IDASetUserData(mem, data);
  flag = IDASetId(mem, id);
  flag = IDASetSuppressAlg(mem, TRUE);
  /* Call IDADense and set up the linear solver. */
  flag = IDADense(mem, NEQ);

  N_VConst(ZERO, q);
  IDAQuadInit(mem, rhsQ, q);
  IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  IDASetQuadErrCon(mem, TRUE);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);

  IDAGetQuad(mem,&tret,q);
  G = Ith(q,1);
  printf("  G  =%12.6e\n", Ith(q,1));

  /******************************
  * BACKWARD for k
  ******************************/
  data->params[0] -= dp;
  setIC(yy, yp, data);

  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gm[0] = Ith(q,1);
  printf("Gm[0]=%12.6e\n", Ith(q,1));

  /****************************
  * FORWARD for k *
  ****************************/
  data->params[0] += (TWO*dp);
  setIC(yy, yp, data);
  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gp[0] = Ith(q,1);
  printf("Gp[0]=%12.6e\n", Ith(q,1));


  /* Backward for c */
  data->params[0] = ONE;
  data->params[1] -= dp;
  setIC(yy, yp, data);
  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gm[1] = Ith(q,1);

  /* Forward for c */
  data->params[1] += (TWO*dp);
  setIC(yy, yp, data);
  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gp[1] = Ith(q,1);

  IDAFree(&mem);


  printf("---------------BACKWARD------------------\n");
  printf("   dG/dp:  %12.4le %12.4le\n", (G-Gm[0])/dp, (G-Gm[1])/dp);
  printf("-----------------------------------------\n\n");

  printf("---------------FORWARD--------------------\n");
  printf("   dG/dp:  %12.4le %12.4le\n", (Gp[0]-G)/dp, (Gp[1]-G)/dp);
  printf("-----------------------------------------\n\n");

  printf("--------------CENTERED-------------------\n");
  printf("   dG/dp:  %12.4le %12.4le\n", (Gp[0]-Gm[0])/(TWO*dp) ,(Gp[1]-Gm[1])/(TWO*dp));
  printf("-----------------------------------------\n\n");


  /* Free memory */
  free(data);

  N_VDestroy(id);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(q);
  /*  N_VDestroyVectorArray(yyS, NP);
  N_VDestroyVectorArray(ypS, NP);
  N_VDestroy_Serial(yyB);
  N_VDestroy_Serial(ypB);
  N_VDestroy_Serial(qB);
  */
  return(0);
  
}

static void setIC(N_Vector yy, N_Vector yp, UserData data)
{
  realtype pi;
  realtype a, J1, m2, J2;
  realtype q, p, x;
  realtype Q[3];

  N_VConst(ZERO, yy);
  N_VConst(ZERO, yp);

  pi = FOUR*atan(ONE);

  a = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;
  
  q = pi/TWO;
  p = asin(-a);
  x = cos(p);

  NV_Ith_S(yy,0) = q;
  NV_Ith_S(yy,1) = x;
  NV_Ith_S(yy,2) = p;
  
  force(yy, Q, data);

  NV_Ith_S(yp,3) = Q[0]/J1;
  NV_Ith_S(yp,4) = Q[1]/m2;
  NV_Ith_S(yp,5) = Q[2]/J2;

}

static void setICB(N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, UserData data)
{
  realtype a, J1, J2, m1, m2;

  N_VConst(ZERO, yyB);
  N_VConst(ZERO, ypB);

  a = data->a;
  J1 = data->J1;
  J2 = data->J2;
  m1 = data->m1;
  m2 = data->m2;

  Ith(ypB,4) = -(J1-QUARTER*HALF*a*a*m1)*Ith(yy,4)/J1;
  Ith(ypB,5) = -Ith(yy,5); 
  Ith(ypB,6) = -Ith(yy,6);
}

static void force(N_Vector yy, realtype *Q, UserData data)
{
  realtype a, k, c, l0, F;
  realtype q, x, p;
  realtype qd, xd, pd;  
  realtype s1, c1, s2, c2, s21, c21;
  realtype l2, l, ld;
  realtype f, fl;

  a = data->a;
  k = data->params[0];
  c = data->params[1];
  l0 = data->l0;
  F = data->F;

  q = NV_Ith_S(yy,0);
  x = NV_Ith_S(yy,1);
  p = NV_Ith_S(yy,2);

  qd = NV_Ith_S(yy,3);
  xd = NV_Ith_S(yy,4);
  pd = NV_Ith_S(yy,5);

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);
  s21 = s2*c1 - c2*s1;
  c21 = c2*c1 + s2*s1;

  l2 = x*x - x*(c2+a*c1) + (ONE + a*a)/FOUR + a*c21/TWO;
  l = RSqrt(l2);
  ld = TWO*x*xd - xd*(c2+a*c1) + x*(s2*pd+a*s1*qd) - a*s21*(pd-qd)/TWO;
  ld /= TWO*l;

  f = k*(l-l0) + c*ld;
  fl = f/l;

  Q[0] = - fl * a * (s21/TWO + x*s1) / TWO;
  Q[1] = fl * (c2/TWO - x + a*c1/TWO) + F;
  Q[2] = - fl * (x*s2 - a*s21/TWO) / TWO - F*s2;

}

static int ressc(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  UserData data;
  realtype Q[3];
  realtype a, J1, m2, J2;
  realtype *yval, *ypval, *rval;
  realtype q, x, p;
  realtype qd, xd, pd;  
  realtype lam1, lam2, mu1, mu2;
  realtype s1, c1, s2, c2;

  data = (UserData) user_data;

  a  = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  yval = NV_DATA_S(yy); 
  ypval = NV_DATA_S(yp); 
  rval = NV_DATA_S(rr);

  q = yval[0];
  x = yval[1];
  p = yval[2];

  qd = yval[3];
  xd = yval[4];
  pd = yval[5];

  lam1 = yval[6];
  lam2 = yval[7];

  mu1 = yval[8];
  mu2 = yval[9];

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);

  force(yy, Q, data);

  rval[0] = ypval[0] - qd + a*s1*mu1 - a*c1*mu2;
  rval[1] = ypval[1] - xd + mu1;
  rval[2] = ypval[2] - pd + s2*mu1 - c2*mu2; 

  rval[3] = J1*ypval[3] - Q[0] + a*s1*lam1 - a*c1*lam2;
  rval[4] = m2*ypval[4] - Q[1] + lam1;
  rval[5] = J2*ypval[5] - Q[2] + s2*lam1 - c2*lam2; 

  rval[6] = x - c2 - a*c1;
  rval[7] = -s2 - a*s1;

  rval[8] = a*s1*qd + xd + s2*pd;
  rval[9] = -a*c1*qd - c2*pd;

  return(0);
}

static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data)
{
  realtype v1, v2, v3;
  realtype m1, J1, J1hat, m2, J2, a;
  UserData data;
  
  data = (UserData) user_data;
  J1 = data->J1;
  m1 = data->m1;
  m2 = data->m2;
  J2 = data->J2;
  a  = data->a;

  J1hat = J1 - HALF*QUARTER*a*a*m1;

  v1 = Ith(yy,4); v2 = Ith(yy,5); v3 = Ith(yy,6);

  Ith(qdot,1) = HALF*(J1hat*v1*v1 + m2*v2*v2 + J2*v3*v3);

  return(0);
}

static int rhsQS(int Ns, realtype t, N_Vector yy, N_Vector yp, 
                 N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsQS,
                 void *user_data,  N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS)
{
  realtype v1, v2, v3;
  realtype m1, J1, J1hat, m2, J2, a;
  UserData data;
  realtype s1, s2, s3;
  
  data = (UserData) user_data;
  J1 = data->J1;
  m1 = data->m1;
  m2 = data->m2;
  J2 = data->J2;
  a  = data->a;

  J1hat = J1 - HALF*QUARTER*a*a*m1;

  v1 = Ith(yy,4); v2 = Ith(yy,5); v3 = Ith(yy,6);
  
  /* Sensitivities of v. */
  s1 = Ith(yyS[0],4);
  s2 = Ith(yyS[0],5);
  s3 = Ith(yyS[0],6);

  Ith(rhsQS[0], 1) = J1hat*v1*s1 + m2*v2*s2 + J2*v3*s3;

  s1 = Ith(yyS[1],4);
  s2 = Ith(yyS[1],5);
  s3 = Ith(yyS[1],6);

  Ith(rhsQS[1], 1) = J1hat*v1*s1 + m2*v2*s2 + J2*v3*s3;

  return(0);
}
static int resscB(realtype tt, 
                  N_Vector yy, N_Vector yp,
                  N_Vector yyB, N_Vector ypB, N_Vector rrB,
                  void *user_dataB)
{
  UserData data;
  realtype y1, y2, y3, y4, y5, y6, y7, y8, y9, y10;
  realtype l1, l2, l3, l4, l5, l6, l7, l8, l9, l10;
  realtype lp1, lp2, lp3, lp4, lp5, lp6;
  realtype f, l, lp, llp, dldy1, dlpdy1, dfldy1, dldy2, dlpdy2, dfldy2;
  realtype dldy3, dlpdy3, dfldy3, dfldy4, dfldy5, dfldy6;
  realtype m1, J1, m2, J2, a, k, c, l0, F; //problem's parameters
  realtype cy1, sy1, cy3, sy3, cy13, sy13; //cosines and sines
  realtype J1hat, lsq;

  data = (UserData) user_dataB;
  a  = data->a;
  J1 = data->J1;
  m1 = data->m1;
  m2 = data->m2;
  J2 = data->J2;
  F = data->F;
  l0 = data->l0;
  c = data->params[1];
  k = data->params[0];

  /* The y  vector */
  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3); y4 = Ith(yy,4);  y5 = Ith(yy, 5); 
  y6 = Ith(yy,6); y7 = Ith(yy,7); y8 = Ith(yy,8); y9 = Ith(yy,9); y10 = Ith(yy,10);

  /* The lambda vector */
  l1 = Ith(yyB,1); l2 = Ith(yyB,2); l3 = Ith(yyB,3); l4 = Ith(yyB,4);  l5 = Ith(yyB, 5); 
  l6 = Ith(yyB,6); l7 = Ith(yyB,7); l8 = Ith(yyB,8); l9 = Ith(yyB,9); l10 = Ith(yyB,10);

  /* The lambda dot vector */
  lp1 = Ith(ypB,1); lp2 = Ith(ypB,2); lp3 = Ith(ypB,3);
  lp4 = Ith(ypB,4); lp5 = Ith(ypB,5); lp6 = Ith(ypB,6);

  cy1=cos(y1); sy1=sin(y1); cy3=cos(y3); sy3=sin(y3);
  cy13=cos(y1-y3);sy13 = sin(y1-y3);

  /* l, l' and f*/
  lsq = y2*y2 - y2*(cy3+a*cy1) + (ONE + a*a)/FOUR + a*cy13/TWO;
  l = RSqrt(lsq);
  llp = TWO*y2*y5 - y5*(cy3+a*cy1) + y2*(sy3*y6+a*sy1*y4) - a*sy13*(y4-y6)/TWO;
  lp = llp/(TWO*l);
  f = k*(l-l0) + c*lp;

  /* d(f/l)/dy1 */
  dldy1  = HALF*a*(y2*sy1-HALF*a*sy13) / l; 
  dlpdy1 = HALF*( (a*y5*sy1 + a*y2*y4*cy1 + HALF*a*cy13*(y6-y4))*l - llp*dldy1) / lsq; 
  dfldy1 = ( (k*dldy1+c*dlpdy1)*l - f*dldy1 ) / lsq;

  /* d(f/l)/dy2 */
  dldy2  = HALF*(TWO*y2-cy3-a*cy1)/l;
  dlpdy2 = HALF*( (TWO*y5+sy3*y6+a*sy1*y4)*l - llp*dldy2 ) / lsq;
  dfldy2 = ( (k*dldy2+c*dlpdy2)*l - f*dldy2 ) / lsq;

  /* d(f/l)/dy3 */
  dldy3  = HALF*(y2*sy3+HALF*a*sy13)/l;
  dlpdy3 = HALF*( (y5*sy3 + y2*cy3*y6 - HALF*a*cy13*(y6-y4))*l - llp*dldy3 ) / lsq;
  dfldy3 = ( (k*dldy3+c*dlpdy3)*l - f*dldy3 ) /lsq;

  /* d(f/l)/dy4..6 */
  dfldy4 = HALF*c*a*(y2*sy1-HALF*sy13)/lsq;
  dfldy5 = HALF*c*(TWO*y2-cy3-a*cy1)/lsq;
  dfldy6 = HALF*c*(y2*sy3+HALF*a*sy13)/lsq;

  J1hat = J1-QUARTER*HALF*a*a*m1;

  Ith(rrB, 1) = lp1 -
    a*(y9*cy1+y10*sy1)*l1 + 
    ( HALF*f/l*a*(HALF*cy13-y2*cy1) + HALF*a*(HALF*sy13-y2*sy1)*dfldy1 - a*(y7*cy1+y8*sy1)  )*l4 - 
    ( HALF*f/l*a*sy1-(HALF*cy3-y2+HALF*a*cy1)*dfldy1 )*l5 - 
    HALF*(HALF*f/l*a*cy13 + (y2*sy3+HALF*a*sy13)*dfldy1)*l6 - 
    //( HALF*f/l*(y2*cy3-HALF*a*cy13) + (y2*sy3+HALF*a*sy13)*dfldy
    a*(sy1*l7 - cy1*l8) - a*y4*(cy1*l9+sy1*l10) ;

  Ith(rrB, 2) = lp2     - 
    (HALF*f/l*a*sy1 - HALF*a*(HALF*sy13-y2*sy1)*dfldy2)*l4 - 
    (f/l - (HALF*cy3-y2+HALF*a*cy1)*dfldy2)*l5 - 
    HALF*(f/l*sy3+(y2*sy3+HALF*a*sy13)*dfldy2)*l6 - 
    l7;

  Ith(rrB, 3) = lp3     - 
    (y9*cy3+y10*sy3)*l3 - 
    HALF*a*(HALF*f/l*cy13 - (HALF*sy13-y2*sy1)*dfldy3)*l4 - 
    (HALF*f/l*sy3 - (HALF*cy3-y2+HALF*a*cy1)*dfldy3 )*l5 - 
    ( HALF*f/l*(y2*cy3-HALF*a*cy13) + HALF*(y2*sy3+HALF*a*sy13)*dfldy3   + (F+y7)*cy3 + y8*sy3 )*l6 - 
    sy3*l7 + cy3*l8 - y6*cy3*l9 - y6*sy3*l10;

  Ith(rrB, 4) = J1*lp4 + l1 +
    HALF*a*(HALF*sy13-y2*sy1)*dfldy4 * l4 +
    (HALF*cy3-y2+HALF*a*cy1)*dfldy4 * l5 -
    HALF*(y2*sy3+HALF*a*sy13)*dfldy4 * l6 -
    a*(sy1*l9-cy1*l10)+
    J1hat*y4;

  Ith(rrB, 5) = m2*lp5  + l2 +
    HALF*a*(HALF*sy13-y2*sy1)*dfldy5 * l4 +
    (HALF*cy3-y2+HALF*a*cy1)*dfldy5 * l5 -
    HALF*(y2*sy3+HALF*a*sy13)*dfldy5 * l6 -
    l9+
    m2*y5;

  Ith(rrB, 6) = J2*lp6  + l3 +
    HALF*a*(HALF*sy13-y2*sy1)*dfldy6 * l4 +
    (HALF*cy3-y2+HALF*a*cy1)*dfldy6 * l5 -
    HALF*(y2*sy3+HALF*a*sy13)*dfldy6 * l6 -
    sy3*l9 +cy3*l10+
    J2*y6;

  Ith(rrB, 7) =         - a*sy1*l4 - l5 - sy3*l6;
  Ith(rrB, 8) =         + a*cy1*l4      + cy3*l6 ;
  Ith(rrB, 9) =         - a*sy1*l1 - l2 - sy3*l3;
  Ith(rrB,10) =         + a*cy1*l1      + cy3*l3;
  
  return(0);
}

static int rhsQB(realtype tt, 
                 N_Vector yy, N_Vector yp, 
                 N_Vector yyB, N_Vector ypB, 
                 N_Vector rrQB, void *user_dataB)
{
  realtype y1, y2, y3, y4, y5, y6, y7, y8, y9, y10;
  realtype l1, l2, l3, l4, l5, l6, l7, l8, l9, l10;
  realtype cy1, cy3, sy1, sy3, cy13, sy13, t1, t2, t3;
  realtype l, lp, lp_l, ll0_l;
  realtype m1, J1, m2, J2, a, k, c, l0, F; //problem's parameters
  UserData data;

  data = (UserData) user_dataB;
  a  = data->a;
  J1 = data->J1;
  m1 = data->m1;
  m2 = data->m2;
  J2 = data->J2;
  F = data->F;
  l0 = data->l0;
  c = data->params[1];
  k = data->params[0];

  /* The y  vector */
  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3); y4 = Ith(yy,4);  y5 = Ith(yy, 5); 
  y6 = Ith(yy,6); y7 = Ith(yy,7); y8 = Ith(yy,8); y9 = Ith(yy,9); y10 = Ith(yy,10);

  /* The lambda vector */
  l1 = Ith(yyB,1); l2 = Ith(yyB,2); l3 = Ith(yyB,3); l4 = Ith(yyB,4);  l5 = Ith(yyB, 5); 
  l6 = Ith(yyB,6); l7 = Ith(yyB,7); l8 = Ith(yyB,8); l9 = Ith(yyB,9); l10 = Ith(yyB,10);

  cy1=cos(y1); sy1=sin(y1); cy3=cos(y3); sy3=sin(y3);
  cy13=cos(y1-y3);sy13 = sin(y1-y3);

  l = y2*y2 - y2*(cy3+a*cy1) + QUARTER*(ONE+a*a) + HALF*a*cy13;
  l = RSqrt(l);

  lp = (TWO*y2*y5 - y5*(cy3+a*cy1) + y2*(sy3*y6+a*sy1*y4) - HALF*a*sy13*(y4-y6))/(TWO*l);

  ll0_l = (l-l0)/l;
  lp_l  = lp/l;
  
  t1 = HALF*a*(HALF*sy13-y2*sy1);
  t2 = HALF*cy3-y2+HALF*a*cy1;
  t3 = -HALF*(y2*sy3+HALF*a*sy13);

  Ith(rrQB,1) = t1*l4 + t2*l5 + t3*l6;


  Ith(rrQB,2) = - lp_l  * Ith(rrQB,1);
  Ith(rrQB,1) = - ll0_l * Ith(rrQB,1);

  return(0);
}

static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int flag, kused;
  long int nst;
  realtype hused;

  yval  = NV_DATA_S(y);

  flag = IDAGetLastOrder(mem, &kused);
  flag = IDAGetNumSteps(mem, &nst);
  flag = IDAGetLastStep(mem, &hused);

  printf("%10.4le %12.4le %12.4le %12.4le %3ld  %1d %12.4le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
}


static void PrintFinalStats(void *mem)
{
  int flag;
  long int nst, nni, nje, nre, nreLS, netf, ncfn;

  flag = IDAGetNumSteps(mem, &nst);
  flag = IDAGetNumResEvals(mem, &nre);
  flag = IDADlsGetNumJacEvals(mem, &nje);
  flag = IDAGetNumNonlinSolvIters(mem, &nni);
  flag = IDAGetNumErrTestFails(mem, &netf);
  flag = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  flag = IDADlsGetNumResEvals(mem, &nreLS);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
}

static void PrintFinalStatsB(void *mem, int bckPb)
{
  int flag;
  long int nst, nni, nje, nre, nreLS, netf, ncfn;

  flag = IDAGetNumSteps(IDAGetAdjIDABmem(mem,bckPb), &nst);
  flag = IDAGetNumResEvals(IDAGetAdjIDABmem(mem,bckPb), &nre);
  flag = IDADlsGetNumJacEvals(IDAGetAdjIDABmem(mem,bckPb), &nje);
  flag = IDAGetNumNonlinSolvIters(IDAGetAdjIDABmem(mem,bckPb), &nni);
  flag = IDAGetNumErrTestFails(IDAGetAdjIDABmem(mem,bckPb), &netf);
  flag = IDAGetNumNonlinSolvConvFails(IDAGetAdjIDABmem(mem,bckPb), &ncfn);
  flag = IDADlsGetNumResEvals(IDAGetAdjIDABmem(mem,bckPb), &nreLS);

  printf("\nFinal Run Statistics (BACKWARD PB): \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
}



static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}


