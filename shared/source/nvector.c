/*******************************************************************
 *                                                                 *
 * File          : nvector.c                                       *
 * Programmers   : Radu Serban, LLNL                               *
 * Version of    : 26 March 2003                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for a generic NVECTOR           *
 * package. It contains the implementation of the N_Vector         *
 * kernels listed in nvector.h.                                    *
 *                                                                 *
 *******************************************************************/

#include <stdlib.h>
#include "nvector.h"    /* generic M_Env and N_Vector */

N_Vector N_VNew(M_Env machEnv)
{
  N_Vector v_new;
  v_new = machEnv->ops->nvnew(machEnv);
  return(v_new);
}

void N_VSpace(M_Env machEnv, long int *lrw, long int *liw)
{
  machEnv->ops->nvspace(machEnv, lrw, liw);
}

void N_VFree(N_Vector v)
{
  v->menv->ops->nvfree(v);
}

N_Vector_S N_VNew_S(integertype ns, M_Env machEnv)
{
  N_Vector_S vs_new;
  integertype is, js;

  if (ns <= 0) return(NULL);
  if (machEnv == NULL) return(NULL);
  vs_new = (N_Vector_S) malloc(ns * sizeof(N_Vector *));

  for (is=0; is<ns; is++) {
    vs_new[is] = N_VNew(machEnv);
    if (vs_new[is] == NULL) {
      for (js=0; js<is; js++) N_VFree(vs_new[js]);
      free(vs_new);
      return(NULL);
    }
  }

  return(vs_new);
}

void N_VFree_S(integertype ns, N_Vector_S vs) 
{
  integertype is;
  
  for (is=0; is<ns; is++) N_VFree(vs[is]);
  free(vs);
}

N_Vector N_VMake(realtype *v_data, M_Env machEnv)
{
  N_Vector v_new;
  v_new = machEnv->ops->nvmake(v_data, machEnv);
  return(v_new);
}

void N_VDispose(N_Vector v)
{
  v->menv->ops->nvdispose(v);
}

realtype *N_VGetData(N_Vector v)
{
  realtype *data;
  data = v->menv->ops->nvgetdata(v);
  return(data);
}

void N_VSetData(realtype *v_data, N_Vector v)
{
  v->menv->ops->nvsetdata(v_data, v);
}

void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  z->menv->ops->nvlinearsum(a, x, b, y, z);
}

void N_VConst(realtype c, N_Vector z)
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

void N_VScale(realtype c, N_Vector x, N_Vector z) 
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

void N_VAddConst(N_Vector x, realtype b, N_Vector z)
{
  z->menv->ops->nvaddconst(x, b, z);
}

realtype N_VDotProd(N_Vector x, N_Vector y)
{
  realtype prod;
  prod = y->menv->ops->nvdotprod(x, y);
  return(prod);
}

realtype N_VMaxNorm(N_Vector x)
{
  realtype norm;
  norm = x->menv->ops->nvmaxnorm(x);
  return(norm);
}

realtype N_VWrmsNorm(N_Vector x, N_Vector w)
{
  realtype norm;
  norm = x->menv->ops->nvwrmsnorm(x, w);
  return(norm);
}

realtype N_VMin(N_Vector x)
{
  realtype minval;
  minval = x->menv->ops->nvmin(x);
  return(minval);
}

realtype N_VWL2Norm(N_Vector x, N_Vector w)
{
  realtype norm;
  norm = x->menv->ops->nvwl2norm(x, w);
  return(norm);
}

realtype N_VL1Norm(N_Vector x)
{
  realtype norm;
  norm = x->menv->ops->nvl1norm(x);
  return(norm);
}

void N_VOneMask(N_Vector x)
{
  x->menv->ops->nvonemask(x);
}

void N_VCompare(realtype c, N_Vector x, N_Vector z)
{
  z->menv->ops->nvcompare(c, x, z);
}

booleantype N_VInvTest(N_Vector x, N_Vector z)
{
  booleantype flag;
  flag = z->menv->ops->nvinvtest(x, z);
  return(flag);
}

booleantype N_VConstrProdPos(N_Vector c, N_Vector x)
{
  booleantype flag;
  flag = x->menv->ops->nvconstrprodpos(c, x);
  return(flag);
}

booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  booleantype flag;
  flag = x->menv->ops->nvconstrmask(c, x, m);
  return(flag);
}

realtype N_VMinQuotient(N_Vector num, N_Vector denom)
{
  realtype quotient;
  quotient = num->menv->ops->nvminquotient(num, denom);
  return(quotient);
}

void N_VPrint(N_Vector x)
{
  x->menv->ops->nvprint(x);
}
