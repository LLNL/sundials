/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2004-07-22 22:20:43 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban, LLNL                               
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/shared/LICENSE
 * -----------------------------------------------------------------
 * This is the implementation file for a generic NVECTOR package. 
 * It contains the implementation of the N_Vector kernels listed in 
 * nvector.h.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "nvector.h"

/*
 * Functions in the 'ops' structure
 */

N_Vector N_VClone(N_Vector w)
{
  N_Vector v;
  v = w->ops->nvclone(w);
  return(v);
}

void N_VDestroy(N_Vector v)
{
  v->ops->nvdestroy(v);
}

N_Vector N_VCloneEmpty(N_Vector w)
{
  N_Vector v;
  v = w->ops->nvcloneempty(w);
  return(v);
}

void N_VDestroyEmpty(N_Vector v)
{
  v->ops->nvdestroyempty(v);
}

void N_VSpace(N_Vector v, long int *lrw, long int *liw)
{
  v->ops->nvspace(v, lrw, liw);
}

realtype *N_VGetArrayPointer(N_Vector v)
{
  realtype *data;
  data = v->ops->nvgetarraypointer(v);
  return(data);
}

void N_VSetArrayPointer(realtype *v_data, N_Vector v)
{
  v->ops->nvsetarraypointer(v_data, v);
}

void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  z->ops->nvlinearsum(a, x, b, y, z);
}

void N_VConst(realtype c, N_Vector z)
{
  z->ops->nvconst(c, z);
}

void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  z->ops->nvprod(x, y, z);
}

void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  z->ops->nvdiv(x, y, z);
}

void N_VScale(realtype c, N_Vector x, N_Vector z) 
{
  z->ops->nvscale(c, x, z);
}

void N_VAbs(N_Vector x, N_Vector z)
{
  z->ops->nvabs(x, z);
}

void N_VInv(N_Vector x, N_Vector z)
{
  z->ops->nvinv(x, z);
}

void N_VAddConst(N_Vector x, realtype b, N_Vector z)
{
  z->ops->nvaddconst(x, b, z);
}

realtype N_VDotProd(N_Vector x, N_Vector y)
{
  realtype prod;
  prod = y->ops->nvdotprod(x, y);
  return(prod);
}

realtype N_VMaxNorm(N_Vector x)
{
  realtype norm;
  norm = x->ops->nvmaxnorm(x);
  return(norm);
}

realtype N_VWrmsNorm(N_Vector x, N_Vector w)
{
  realtype norm;
  norm = x->ops->nvwrmsnorm(x, w);
  return(norm);
}

realtype N_VWrmsNormMask(N_Vector x, N_Vector w, N_Vector id)
{
  realtype norm;
  norm = x->ops->nvwrmsnormmask(x, w, id);
  return(norm);
}

realtype N_VMin(N_Vector x)
{
  realtype minval;
  minval = x->ops->nvmin(x);
  return(minval);
}

realtype N_VWL2Norm(N_Vector x, N_Vector w)
{
  realtype norm;
  norm = x->ops->nvwl2norm(x, w);
  return(norm);
}

realtype N_VL1Norm(N_Vector x)
{
  realtype norm;
  norm = x->ops->nvl1norm(x);
  return(norm);
}

void N_VCompare(realtype c, N_Vector x, N_Vector z)
{
  z->ops->nvcompare(c, x, z);
}

booleantype N_VInvTest(N_Vector x, N_Vector z)
{
  booleantype flag;
  flag = z->ops->nvinvtest(x, z);
  return(flag);
}

booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  booleantype flag;
  flag = x->ops->nvconstrmask(c, x, m);
  return(flag);
}

realtype N_VMinQuotient(N_Vector num, N_Vector denom)
{
  realtype quotient;
  quotient = num->ops->nvminquotient(num, denom);
  return(quotient);
}

/*
 * Additional functions exported by the generic NVECTOR
 *  N_VCloneVectorArray
 *  N_VDestroyVectorArray
 */

N_Vector *N_VCloneVectorArray(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j=0; j<count; j++) {
    vs[j] = N_VClone(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray(vs, j-1);
      return(NULL);
    }
  }
  
  return(vs);

}

void N_VDestroyVectorArray(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy(vs[j]);

  free(vs);
}
