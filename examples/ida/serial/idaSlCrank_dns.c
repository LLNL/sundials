/*
 * -----------------------------------------------------------------
 * $Revision: 4809 $
 * $Date: 2016-07-18 11:16:25 -0700 (Mon, 18 Jul 2016) $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal force
 * applied to the connecting rod and a spring-damper connecting the crank
 * and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 *
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Problem Constants */

#define NEQ   10

#define TEND  RCONST(10.0)
#define NOUT  41


#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FOUR RCONST(4.0)

typedef struct {
  realtype a;
  realtype J1, J2, m2;
  realtype k, c, l0;
  realtype F;
} *UserData;

int ressc(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data);

void setIC(N_Vector yy, N_Vector yp, UserData data);
void force(N_Vector yy, realtype *Q, UserData data);

static void PrintHeader(realtype rtol, realtype atol, N_Vector y);
static void PrintOutput(void *mem, realtype t, N_Vector y);
static void PrintFinalStats(void *mem);

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  UserData data;

  void *mem;
  N_Vector yy, yp, id;
  realtype rtol, atol;
  realtype t0, tf, tout, dt, tret;
  int flag, iout;

  /* User data */

  data = (UserData) malloc(sizeof *data);

  data->a = 0.5;   /* half-length of crank */
  data->J1 = 1.0;  /* crank moment of inertia */
  data->m2 = 1.0;  /* mass of connecting rod */
  data->J2 = 2.0;  /* moment of inertia of connecting rod */
  data->k = 1.0;   /* spring constant */
  data->c = 1.0;   /* damper constant */
  data->l0 = 1.0;  /* spring free length */
  data->F = 1.0;   /* external constant force */

  /* Create N_Vectors */
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  id = N_VNew_Serial(NEQ);

  /* Consistent IC */
  setIC(yy, yp, data);

  /* ID array */
  N_VConst(ONE, id);
  NV_Ith_S(id,6) = ZERO;
  NV_Ith_S(id,7) = ZERO;
  NV_Ith_S(id,8) = ZERO;
  NV_Ith_S(id,9) = ZERO;

  /* Tolerances */
  rtol = RCONST(1.0e-6);
  atol = RCONST(1.0e-6);

  /* Integration limits */
  t0 = ZERO;
  tf = TEND;
  dt = (tf-t0)/(NOUT-1);

  /* IDA initialization */
  mem = IDACreate();
  flag = IDAInit(mem, ressc, t0, yy, yp);
  flag = IDASStolerances(mem, rtol, atol);
  flag = IDASetUserData(mem, data);
  flag = IDASetId(mem, id);
  flag = IDASetSuppressAlg(mem, TRUE);

  /* Call IDADense and set up the linear solver. */
  flag = IDADense(mem, NEQ);

  PrintHeader(rtol, atol, yy);

  /* In loop, call IDASolve, print results, and test for error. */

  PrintOutput(mem,t0,yy);

  tout = dt;
  for (iout=1; iout<NOUT; iout++) {
    tout = iout*dt;
    flag = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    if (flag < 0) break;

    PrintOutput(mem,tret,yy);

  }

  PrintFinalStats(mem);

  /* Free memory */

  free(data);
  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(id);

  return(0);
  
}

void setIC(N_Vector yy, N_Vector yp, UserData data)
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

void force(N_Vector yy, realtype *Q, UserData data)
{
  realtype a, k, c, l0, F;
  realtype q, x, p;
  realtype qd, xd, pd;  
  realtype s1, c1, s2, c2, s21, c21;
  realtype l2, l, ld;
  realtype f, fl;

  a = data->a;
  k = data->k;
  c = data->c;
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
  l = SUNRsqrt(l2);
  ld = TWO*x*xd - xd*(c2+a*c1) + x*(s2*pd+a*s1*qd) - a*s21*(pd-qd)/TWO;
  ld /= TWO*l;

  f = k*(l-l0) + c*ld;
  fl = f/l;

  Q[0] = - fl * a * (s21/TWO + x*s1) / TWO;
  Q[1] = fl * (c2/TWO - x + a*c1/TWO) + F;
  Q[2] = - fl * (x*s2 - a*s21/TWO) / TWO - F*s2;

}

int ressc(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
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

  yval = N_VGetArrayPointer_Serial(yy); 
  ypval = N_VGetArrayPointer_Serial(yp); 
  rval = N_VGetArrayPointer_Serial(rr);

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

static void PrintHeader(realtype rtol, realtype atol, N_Vector y)
{
  printf("\nidaSlCrank_dns: Slider-Crank DAE serial example problem for IDAS\n");
  printf("Linear solver: IDADENSE, Jacobian is computed by IDAS.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n",
         rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n",
         rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n",
         rtol, atol);
#endif
  printf("-----------------------------------------------------------------------\n");
  printf("  t            y1          y2           y3");
  printf("      | nst  k      h\n");
  printf("-----------------------------------------------------------------------\n");
}

static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int flag, kused;
  long int nst;
  realtype hused;

  yval  = N_VGetArrayPointer_Serial(y);

  flag = IDAGetLastOrder(mem, &kused);
  flag = IDAGetNumSteps(mem, &nst);
  flag = IDAGetLastStep(mem, &hused);

  printf("%10.4e %12.4e %12.4e %12.4e %3ld  %1d %12.4e\n", 
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

