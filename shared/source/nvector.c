/*******************************************************************
 *                                                                 *
 * File          : nvector.c                                       *
 * Programmers   : Radu Serban, LLNL                               *
 * Version of    : 07 February 2004                                *
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
#include "nvector.h"    /* generic NV_Spec and N_Vector */

N_Vector N_VNew(NV_Spec nvSpec)
{
  N_Vector v_new;
  v_new = nvSpec->ops->nvnew(nvSpec);
  return(v_new);
}

void N_VSpace(NV_Spec nvSpec, long int *lrw, long int *liw)
{
  nvSpec->ops->nvspace(nvSpec, lrw, liw);
}

void N_VFree(N_Vector v)
{
  v->nvspec->ops->nvfree(v);
}

N_Vector_S N_VNew_S(int ns, NV_Spec nvSpec)
{
  N_Vector_S vs_new;
  int is, js;

  if (ns <= 0) return(NULL);
  if (nvSpec == NULL) return(NULL);
  vs_new = (N_Vector_S) malloc(ns * sizeof(N_Vector *));
  if(vs_new == NULL) return(NULL);

  for (is=0; is<ns; is++) {
    vs_new[is] = N_VNew(nvSpec);
    if (vs_new[is] == NULL) {
      for (js=0; js<is; js++) N_VFree(vs_new[js]);
      free(vs_new);
      return(NULL);
    }
  }

  return(vs_new);
}

void N_VFree_S(int ns, N_Vector_S vs) 
{
  int is;
  
  for (is=0; is<ns; is++) N_VFree(vs[is]);
  free(vs);
}

N_Vector N_VMake(void *v_data, NV_Spec nvSpec)
{
  N_Vector v_new;
  v_new = nvSpec->ops->nvmake(v_data, nvSpec);
  return(v_new);
}

void N_VDispose(N_Vector v)
{
  v->nvspec->ops->nvdispose(v);
}

void *N_VGetData(N_Vector v)
{
  void *data;
  data = v->nvspec->ops->nvgetdata(v);
  return(data);
}

void N_VSetData(void *v_data, N_Vector v)
{
  v->nvspec->ops->nvsetdata(v_data, v);
}

void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  z->nvspec->ops->nvlinearsum(a, x, b, y, z);
}

void N_VConst(realtype c, N_Vector z)
{
  z->nvspec->ops->nvconst(c, z);
}

void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  z->nvspec->ops->nvprod(x, y, z);
}

void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  z->nvspec->ops->nvdiv(x, y, z);
}

void N_VScale(realtype c, N_Vector x, N_Vector z) 
{
  z->nvspec->ops->nvscale(c, x, z);
}

void N_VAbs(N_Vector x, N_Vector z)
{
  z->nvspec->ops->nvabs(x, z);
}

void N_VInv(N_Vector x, N_Vector z)
{
  z->nvspec->ops->nvinv(x, z);
}

void N_VAddConst(N_Vector x, realtype b, N_Vector z)
{
  z->nvspec->ops->nvaddconst(x, b, z);
}

realtype N_VDotProd(N_Vector x, N_Vector y)
{
  realtype prod;
  prod = y->nvspec->ops->nvdotprod(x, y);
  return(prod);
}

realtype N_VMaxNorm(N_Vector x)
{
  realtype norm;
  norm = x->nvspec->ops->nvmaxnorm(x);
  return(norm);
}

realtype N_VWrmsNorm(N_Vector x, N_Vector w)
{
  realtype norm;
  norm = x->nvspec->ops->nvwrmsnorm(x, w);
  return(norm);
}

realtype N_VWrmsNormMask(N_Vector x, N_Vector w, N_Vector id)
{
  realtype norm;
  norm = x->nvspec->ops->nvwrmsnormmask(x, w, id);
  return(norm);
}

realtype N_VMin(N_Vector x)
{
  realtype minval;
  minval = x->nvspec->ops->nvmin(x);
  return(minval);
}

realtype N_VWL2Norm(N_Vector x, N_Vector w)
{
  realtype norm;
  norm = x->nvspec->ops->nvwl2norm(x, w);
  return(norm);
}

realtype N_VL1Norm(N_Vector x)
{
  realtype norm;
  norm = x->nvspec->ops->nvl1norm(x);
  return(norm);
}

void N_VCompare(realtype c, N_Vector x, N_Vector z)
{
  z->nvspec->ops->nvcompare(c, x, z);
}

booleantype N_VInvTest(N_Vector x, N_Vector z)
{
  booleantype flag;
  flag = z->nvspec->ops->nvinvtest(x, z);
  return(flag);
}

booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  booleantype flag;
  flag = x->nvspec->ops->nvconstrmask(c, x, m);
  return(flag);
}

realtype N_VMinQuotient(N_Vector num, N_Vector denom)
{
  realtype quotient;
  quotient = num->nvspec->ops->nvminquotient(num, denom);
  return(quotient);
}

void N_VPrint(N_Vector x)
{
  x->nvspec->ops->nvprint(x);
}
