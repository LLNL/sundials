/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008-04-18 19:42:43 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Slider-crank example for CPODES
 *
 * The multibody system consists of two bodies (crank and
 * connecting rod) with a translational-spring-damper (TSDA)
 * and a constant force acting on the connecting rod.
 *
 * The system has a single degree of freedom. It is modeled with 3
 * generalized coordinates (crank angle, horizontal position of the
 * translational joint, and angle of the connecting rod) and
 * therefore has 2 constraints.
 *
 * Example 6.1.8, pp. 271 in
 * Ed. Haug - Intermediate Dynamics, Prentiss Hall, 1992
 *
 * For its solution with CPODES, the resulting index-3 DAE is reformulated
 * as an ODE with invariants, considering the underlying ODE (obtained by
 * eliminating the Lagrange multipliers from the differential EOM and
 * the acceleration-level constraints) together with the position and 
 * velocity constraints.
 *
 *  |                                  |
 *  |       /\                         |       /\
 *  |      /  \                        |      /  \
 *  | a/2 /    \ 1/2                   |     /    \
 *  |    /      \                      |    /      \
 *  |   /--TSDA--\                     |   /        \
 *  |  /          \                    |  /          \
 *  | / a/2        \ 1/2               | /            \---
 *  |/              \                  |/  \ q1   q3 / \  \
 *  +-----------------------------     +----------------------------
 *                    \                |             \  |\
 *                     \               |               -|-\
 *                      \ 1            |                |  \
 *                       \             |<----- q2 ----->|   \
 *                        \                                  \
 *                         \                                  \
 *                          \ --> F                            \
 *
 * The local reference frame on the crank is positioned at the
 * revolute joint on the ground. The crank has length a, mass m1, and
 * intertia (with respect to the local frame) J1.
 * The local reference frame on the conncting rod is positioned at the
 * translational joint. The connecting rod has length 2, mass m2, and
 * inertia J2.
 * The TSDA has spring constant k, damping constant c, and free length l0.
 * A constant horizontal force F acts on the connecting rod.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <limits.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

#include <sundials/sundials_math.h>

/* Problem Constants */

#define NS 3
#define NC 2

#define CTOL RCONST(1.0e-10)

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FOUR RCONST(4.0)
#define FIVE RCONST(5.0)
#define TEN  RCONST(10.0)

#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)


typedef struct {

  realtype a;
  realtype J1, J2, m2;
  realtype k, c, l0;
  realtype F;

  realtype **M;
  realtype *b;
  int *piv;

} *UserData;

static int f(realtype t, N_Vector y, N_Vector yd, void *f_data);
static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data);
static int cjac(int Nc, int Ny, 
                realtype t, N_Vector y, N_Vector cy,
                DlsMat Jac, void *jac_data,
                N_Vector tmp1, N_Vector tmp2);

void setIC(N_Vector yy, UserData data);
void get_acc(realtype *pos, realtype *vel, realtype *acc, UserData data);

static void PrintOutput(void *cpode_mem, FILE *fout, realtype t, N_Vector y);
static void PrintFinalStats(void *cpode_mem);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  FILE *fout;

  void *cpode_mem;
  realtype rtol, atol;
  N_Vector yy, yp, ctols;
  realtype t, t0, tf;
  int flag;

  /* User data */

  data = (UserData) malloc(sizeof *data);

  data->M = newDenseMat(NS+NC, NS+NC);
  data->b = newRealArray(NS+NC);
  data->piv = newIntArray(NS+NC);

  data->a = 0.5;
  data->J1 = 1.0;
  data->m2 = 1.0;
  data->J2 = 2.0;
  data->k = 1.0;
  data->c = 1.0;
  data->l0 = 1.0;
  data->F = 1.0;

  /* Tolerances */
  rtol = RCONST(1.0e-12);
  atol = RCONST(1.0e-8);

  /* Integration limits */
  t0 = ZERO;
  tf = RCONST(10.0);

  /* Set projection tolerances */
  ctols = N_VNew_Serial(2*NC);
  N_VConst(CTOL, ctols);

  /* Initial conditions */
  yy = N_VNew_Serial(2*NS);
  yp = N_VNew_Serial(2*NS);

  setIC(yy, data);

  /* Create and allocate CPODES memory */
  cpode_mem = CPodeCreate(CP_BDF, CP_NEWTON);
  flag = CPodeSetUserData(cpode_mem, data);
  flag = CPodeInitExpl(cpode_mem, f, t0, yy);
  flag = CPodeSStolerances(cpode_mem, rtol, atol);
  flag = CPodeSetMaxNumSteps(cpode_mem, 2000);

  /* Use the CPODES dense direct linear solver with DQ Jacobian */
  flag = CPDense(cpode_mem, 2*NS);

  /* Use the CPODES projection function using the CPODES 
   * dense direct linear algebra module with user-provided
   * constraint Jacobian function */
  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, ctols);
  flag = CPDenseProj(cpode_mem, 2*NC, 2*NS, CPDLS_QR);
  flag = CPDlsProjSetDenseJacFn(cpode_mem, cjac);

  /* In loop, call CPode, print results, and test for error. */

  fout = fopen("cpsSlCrank_dns.sol", "w");

  t = t0;
  while( t < tf ) {

    flag = CPode(cpode_mem, tf, &t, yy, yp, CP_ONE_STEP);

    if (flag < 0) break;

    PrintOutput(cpode_mem, fout, t, yy);

  }

  fclose(fout);

  printf("\nSolution file written to cpsSlCrank_dns.sol\n\n");
  PrintFinalStats(cpode_mem);

  /* FREE MEMORY */

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(ctols);

  CPodeFree(&cpode_mem);

  free(data);

  return(0);
}


void setIC(N_Vector yy, UserData data)
{
  realtype pi;
  realtype a, J1, m2, J2;
  realtype q, p, x;

  N_VConst(ZERO, yy);

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
  
}


void get_acc(realtype *pos, realtype *vel, realtype *acc, UserData data)
{
  realtype **M, *b;
  int *piv;

  realtype a, a2, k, c, l0, F;
  realtype J1, m2, J2;
  realtype q, x, p;
  realtype qd, xd, pd;  
  realtype s1, c1, s2, c2, s21, c21;
  realtype l2, l, ld;
  realtype f, fl;

  realtype M11, M12, M22, det;
  realtype b1, b2;

  realtype Q1, Q2, Q3;
  realtype g1, g2;
  
  realtype lam1, lam2;

  M = data->M;
  b = data->b;
  piv = data->piv;

  a = data->a;
  a2 = a*a;

  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  k = data->k;
  c = data->c;
  l0 = data->l0;
  F = data->F;

  q = pos[0];
  x = pos[1];
  p = pos[2];

  qd = vel[0];
  xd = vel[1];
  pd = vel[2];

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

  Q1 = - fl * a * (s21/TWO + x*s1) / TWO;
  Q2 = fl * (c2/TWO - x + a*c1/TWO) + F;
  Q3 = - fl * (x*s2 - a*s21/TWO) / TWO - F*s2;

  g1 = -(a*c1*qd*qd + c2*pd*pd);
  g2 = -(a*s1*qd*qd + s2*pd*pd);

  M11 = -a2*s1*s1/J1-ONE/m2-s2*s2/J2;
  M12 = a2*s1*c1/J1+c2*s2/J2;
  /*M21 = M12;*/
  M22 = -a2*c1*c1/J1-c2*c2/J2;

  det = M11*M22-M12*M12;

  b1 = g1-Q1*a*s1/J1-Q2/m2-Q3*s2/J2;
  b2 = g2+Q1*a*c1/J1+Q3*c2/J2; 

  /* Solve M*lam = b */
  lam1 = (M22*b1-M12*b2)/det;
  lam2 = (M11*b2-M12*b1)/det;

  acc[0] = (Q1-a*s1*lam1+a*c1*lam2)/J1;
  acc[1] = (Q2-lam1)/m2;
  acc[2] = (Q3-s2*lam1+c2*lam2)/J2;

  return;
}

static int f(realtype t, N_Vector yy, N_Vector yd, void *f_data)
{
  UserData data;
  realtype *yyval, *ydval;
  realtype *pos, *vel, acc[3];
  
  data = (UserData) f_data;

  /* Extract pos and vel */

  yyval =  NV_DATA_S(yy);
  pos = yyval;
  vel = &(yyval[3]);

  /* Compute acc */

  get_acc(pos, vel, acc, data);

  /* Load yd */

  ydval =  NV_DATA_S(yd);

  ydval[0] = vel[0];
  ydval[1] = vel[1];
  ydval[2] = vel[2];

  ydval[3] = acc[0];
  ydval[4] = acc[1];
  ydval[5] = acc[2];

  return(0);
}



static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data)
{
  UserData data;
  realtype a;
  realtype *yval, *cval;
  realtype q, x, p;
  realtype qd, xd, pd;  
  realtype s1, c1, s2, c2;

  data = (UserData) c_data;

  a  = data->a;

  yval = NV_DATA_S(yy); 
  cval = NV_DATA_S(cout);

  q = yval[0];
  x = yval[1];
  p = yval[2];

  qd = yval[3];
  xd = yval[4];
  pd = yval[5];

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);

  cval[0] = x - c2 - a*c1;
  cval[1] = -s2 - a*s1;

  cval[2] = a*s1*qd + xd + s2*pd;
  cval[3] = -a*c1*qd - c2*pd;

  return(0);
}


static int cjac(int Nc, int Ny, 
                realtype t, N_Vector y, N_Vector cy,
                DlsMat Jac, void *jac_data,
                N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype *yval;
  realtype q, x, p;
  realtype qd, xd, pd;    
  realtype a, as1, ac1, s2, c2;

  data = (UserData) jac_data;

  a  = data->a;

  yval = NV_DATA_S(y); 

  q = yval[0];
  x = yval[1];
  p = yval[2];

  qd = yval[3];
  xd = yval[4];
  pd = yval[5];

  as1 = a*sin(q);
  ac1 = a*cos(q);
  s2 = sin(p);
  c2 = cos(p);


  IJth(Jac,1,1) = as1;
  IJth(Jac,2,1) = ONE;
  IJth(Jac,3,1) = s2;
  IJth(Jac,4,1) = ZERO;
  IJth(Jac,5,1) = ZERO;
  IJth(Jac,6,1) = ZERO;

  IJth(Jac,1,2) = -ac1;
  IJth(Jac,2,2) = ZERO;
  IJth(Jac,3,2) = -c2;
  IJth(Jac,4,2) = ZERO;
  IJth(Jac,5,2) = ZERO;
  IJth(Jac,6,2) = ZERO;

  IJth(Jac,1,3) = ac1*qd*qd;
  IJth(Jac,2,3) = ZERO;
  IJth(Jac,3,3) = c2*pd*pd;
  IJth(Jac,4,3) = as1;
  IJth(Jac,5,3) = ONE;
  IJth(Jac,6,3) = s2;

  IJth(Jac,1,4) = as1*qd*qd;
  IJth(Jac,2,4) = ZERO;
  IJth(Jac,3,4) = s2*pd*pd;
  IJth(Jac,4,4) = -ac1;
  IJth(Jac,5,4) = ZERO;
  IJth(Jac,6,4) = -c2;

  return(0);
}


static void PrintOutput(void *mem, FILE *fout, realtype t, N_Vector y)
{
  realtype *yval;

  yval  = NV_DATA_S(y);

  fprintf(fout, "%10.4e %12.4e %12.4e %12.4e\n",
          t, yval[0], yval[1], yval[2]);
}


static void PrintFinalStats(void *cpode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;
  long int nproj, nce, nsetupsP, nprf;
  int flag;

  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  flag = CPDlsGetNumJacEvals(cpode_mem, &nje);
  flag = CPDlsGetNumFctEvals(cpode_mem, &nfeLS);

  flag = CPodeGetProjStats(cpode_mem, &nproj, &nce, &nsetupsP, &nprf);

  printf("nst = %ld\nnfe = %ld\nnsetups = %ld\n\n", nst, nfe, nsetups);
  printf("nfeLs = %ld  nje = %ld\n", nfeLS, nje);
  printf("nni = %ld  ncfn =  %ld  netf = %ld\n", nni, ncfn, netf);
  printf("nproj = %ld  nce = %ld  nsetupsP = %ld nprf = %ld\n",
         nproj, nce, nsetupsP, nprf);
}
