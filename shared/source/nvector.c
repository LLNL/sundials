/****************************************************************
 *                                                              *
 * File          : nvector.c                                    *
 * Programmers   : Radu Serban, LLNL                            *
 * Version of    : 26 February 2002                             *
 *                                                              *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the implementation file for a generic NVECTOR        *
 * package. It contains the implementation of the N_Vector      *
 * kernels listed in nvector.h.                                 *
 *                                                              *
 ****************************************************************/

#include "nvector.h"    /* generic M_Env and N_Vector */

N_Vector N_VNew(integer n, M_Env machEnv)
{
  N_Vector v_new;
  v_new = machEnv->ops->nvnew(n, machEnv);
  return(v_new);
}

N_Vector_S N_VNew_S(integer ns, integer n, M_Env machEnv)
{
  N_Vector_S vs_new;
  vs_new = machEnv->ops->nvnewS(ns, n, machEnv);
  return(vs_new);
}

void N_VFree(N_Vector v)
{
  v->menv->ops->nvfree(v);
}

void N_VFree_S(integer ns, N_Vector_S vs) 
{
  (*vs)->menv->ops->nvfreeS(ns, vs);
}

N_Vector N_VMake(integer n, real *v_data, M_Env machEnv)
{
  N_Vector v_new;
  v_new = machEnv->ops->nvmake(n, v_data, machEnv);
  return(v_new);
}

void N_VDispose(N_Vector v)
{
  v->menv->ops->nvdispose(v);
}

real *N_VGetData(N_Vector v)
{
  real *data;
  data = v->menv->ops->nvgetdata(v);
  return(data);
}

void N_VSetData(real *v_data, N_Vector v)
{
  v->menv->ops->nvsetdata(v_data, v);
}

void N_VLinearSum(real a, N_Vector x, real b, N_Vector y, N_Vector z)
{
  z->menv->ops->nvlinearsum(a, x, b, y, z);
}

void N_VConst(real c, N_Vector z)
{
  z->menv->ops->nvconst(c, z);
}

void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  z->menv->ops->nvprod(x, y, z);
}

void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  z->menv->ops->nvdiv(x, y, z);
}

void N_VScale(real c, N_Vector x, N_Vector z) 
{
  z->menv->ops->nvscale(c, x, z);
}

void N_VAbs(N_Vector x, N_Vector z)
{
  z->menv->ops->nvabs(x, z);
}

void N_VInv(N_Vector x, N_Vector z)
{
  z->menv->ops->nvinv(x, z);
}

void N_VAddConst(N_Vector x, real b, N_Vector z)
{
  z->menv->ops->nvaddconst(x, b, z);
}

real N_VDotProd(N_Vector x, N_Vector y)
{
  real prod;
  prod = y->menv->ops->nvdotprod(x, y);
  return(prod);
}

real N_VMaxNorm(N_Vector x)
{
  real norm;
  norm = x->menv->ops->nvmaxnorm(x);
  return(norm);
}

real N_VWrmsNorm(N_Vector x, N_Vector w)
{
  real norm;
  norm = x->menv->ops->nvwrmsnorm(x, w);
  return(norm);
}

real N_VMin(N_Vector x)
{
  real minval;
  minval = x->menv->ops->nvmin(x);
  return(minval);
}

real N_VWL2Norm(N_Vector x, N_Vector w)
{
  real norm;
  norm = x->menv->ops->nvwl2norm(x, w);
  return(norm);
}

real N_VL1Norm(N_Vector x)
{
  real norm;
  norm = x->menv->ops->nvl1norm(x);
  return(norm);
}

void N_VOneMask(N_Vector x)
{
  x->menv->ops->nvonemask(x);
}

void N_VCompare(real c, N_Vector x, N_Vector z)
{
  z->menv->ops->nvcompare(c, x, z);
}

boole N_VInvTest(N_Vector x, N_Vector z)
{
  boole flag;
  flag = z->menv->ops->nvinvtest(x, z);
  return(flag);
}

boole N_VConstrProdPos(N_Vector c, N_Vector x)
{
  boole flag;
  flag = x->menv->ops->nvconstrprodpos(c, x);
  return(flag);
}

boole N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  boole flag;
  flag = x->menv->ops->nvconstrmask(c, x, m);
  return(flag);
}

real N_VMinQuotient(N_Vector num, N_Vector denom)
{
  real quotient;
  quotient = num->menv->ops->nvminquotient(num, denom);
  return(quotient);
}

void N_VPrint(N_Vector x)
{
  x->menv->ops->nvprint(x);
}
