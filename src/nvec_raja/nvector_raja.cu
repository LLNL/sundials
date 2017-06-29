/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Slaven Peles @ LLNL                               
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
 */

//#include <limits>

#include <nvector/raja/Vector.hpp>
#include <RAJA/RAJA.hxx>

//#define abs(x) ((x)<0 ? -(x) : (x))

// Need better solution than defines
#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Raja(N_Vector v)
{
  return SUNDIALS_NVEC_RAJA;
}

N_Vector N_VNewEmpty_Raja(long int length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Raja content;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Raja;
  ops->nvclone           = N_VClone_Raja;
  ops->nvcloneempty      = N_VCloneEmpty_Raja;
  ops->nvdestroy         = N_VDestroy_Raja;
  ops->nvspace           = N_VSpace_Raja;
  ops->nvgetarraypointer = NULL; //N_VGetArrayPointer_Raja;
  ops->nvsetarraypointer = NULL; //N_VSetArrayPointer_Raja;
  ops->nvlinearsum       = N_VLinearSum_Raja;
  ops->nvconst           = N_VConst_Raja;
  ops->nvprod            = N_VProd_Raja;
  ops->nvdiv             = N_VDiv_Raja;
  ops->nvscale           = N_VScale_Raja;
  ops->nvabs             = N_VAbs_Raja;
  ops->nvinv             = N_VInv_Raja;
  ops->nvaddconst        = N_VAddConst_Raja;
  ops->nvdotprod         = N_VDotProd_Raja;
  ops->nvmaxnorm         = N_VMaxNorm_Raja;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Raja;
  ops->nvwrmsnorm        = N_VWrmsNorm_Raja;
  ops->nvmin             = N_VMin_Raja;
  ops->nvwl2norm         = N_VWL2Norm_Raja;
  ops->nvl1norm          = N_VL1Norm_Raja;
  ops->nvcompare         = N_VCompare_Raja;
  ops->nvinvtest         = N_VInvTest_Raja;
  ops->nvconstrmask      = N_VConstrMask_Raja;
  ops->nvminquotient     = N_VMinQuotient_Raja;

  /* Create content */
  content = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

    
N_Vector N_VNew_Raja(long int length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja(length);
  if (v == NULL) return(NULL);

  v->content = new sunrajavec::Vector<realtype, long int>(length);

  return(v);
}


N_Vector N_VMake_Raja(N_VectorContent_Raja c)
{
  N_Vector v;
  sunrajavec::Vector<realtype, long int>* x = static_cast<sunrajavec::Vector<realtype, long int>*>(c);
  long int length = x->size();

  v = NULL;
  v = N_VNewEmpty_Raja(length);
  if (v == NULL) return(NULL);

  v->content = c;

  return(v);
}


/* ----------------------------------------------------------------------------
 * Function to create an array of new RAJA-based vectors.
 */

N_Vector *N_VCloneVectorArray_Raja(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Raja(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Raja(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new RAJA-based vectors with NULL data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_Raja(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Raja(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Raja(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Raja
 */

void N_VDestroyVectorArray_Raja(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Raja(vs[j]);

  free(vs); vs = NULL;

  return;
}


/* ----------------------------------------------------------------------------
 * Function to print the a serial vector
 */

void N_VPrint_Raja(N_Vector X)
{
  const realtype *xd = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  long int i;

  for (i = 0; i < N; ++i) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%35.32Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%19.16g\n", xd[i]);
#else
    printf("%11.8g\n", xd[i]);
#endif
  }
  printf("\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Raja(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;
  ops->nvprod            = w->ops->nvprod;
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale;
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */
  v->content = NULL;
  v->ops  = ops;

  return(v);
}

N_Vector N_VClone_Raja(N_Vector w)
{
  N_Vector v;
  sunrajavec::Vector<realtype, long int>* wdat = static_cast<sunrajavec::Vector<realtype, long int>*>(w->content);
  sunrajavec::Vector<realtype, long int>* vdat = new sunrajavec::Vector<realtype, long int>(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Raja(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Raja(N_Vector v)
{
  sunrajavec::Vector<realtype, long int>* x = static_cast<sunrajavec::Vector<realtype, long int>*>(v->content);
  if (x != NULL) {
    if (!x->isClone()) {
      delete x;
      v->content = NULL;
    }
  }

  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Raja(N_Vector X, long int *lrw, long int *liw)
{
    *lrw = sunrajavec::getSize<realtype, long int>(X);
    *liw = 1;
}

void N_VConst_Raja(realtype c, N_Vector Z)
{
  const long int N = sunrajavec::getSize<realtype, long int>(Z);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *ydata = sunrajavec::getDevData<realtype, long int>(Y);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = a*xdata[i] + b*ydata[i];
  });
}

void N_VProd_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *ydata = sunrajavec::getDevData<realtype, long int>(Y);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = xdata[i] * ydata[i];
  });
}

void N_VDiv_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *ydata = sunrajavec::getDevData<realtype, long int>(Y);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = xdata[i] / ydata[i];
  });
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = c * xdata[i];
  });
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = abs(xdata[i]);
  });
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = RCONST(1.0) / xdata[i];
  });
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = xdata[i] + b;
  });
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *ydata = sunrajavec::getDevData<realtype, long int>(Y);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result += xdata[i] * ydata[i] ;
  });

  return static_cast<realtype>(gpu_result);
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceMax<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result.max(abs(xdata[i]));
  });

  return static_cast<realtype>(gpu_result);
}

realtype N_VWrmsNorm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *wdata = sunrajavec::getDevData<realtype, long int>(W);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
  });

  return std::sqrt(static_cast<realtype>(gpu_result)/N);
}

realtype N_VWrmsNormMask_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *wdata = sunrajavec::getDevData<realtype, long int>(W);
  const realtype *iddata = sunrajavec::getDevData<realtype, long int>(ID);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i] * iddata[i]);
  });

  return std::sqrt(static_cast<realtype>(gpu_result)/N);
}

realtype N_VMin_Raja(N_Vector X)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceMin<RAJA::cuda_reduce<128>, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result.min(xdata[i]);
  });

  return static_cast<realtype>(gpu_result);
}

realtype N_VWL2Norm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const realtype *wdata = sunrajavec::getDevData<realtype, long int>(W);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
  });

  return std::sqrt(static_cast<realtype>(gpu_result));
}

realtype N_VL1Norm_Raja(N_Vector X)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(0.0);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    gpu_result += (abs(xdata[i]));
  });

  return static_cast<realtype>(gpu_result);
}

void N_VCompare_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(X);
  const long int N = sunrajavec::getSize<realtype, long int>(X);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(Z);

  RAJA::forall<RAJA::cuda_exec<256> >(0, N, [=] __device__(long int i) {
     zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
  });
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(x);
  const long int N = sunrajavec::getSize<realtype, long int>(x);
  realtype *zdata = sunrajavec::getDevData<realtype, long int>(z);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(ZERO);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    if (xdata[i] == ZERO) {
      gpu_result += ONE;
    } else {
      zdata[i] = ONE/xdata[i];
    }
  });

  return (static_cast<realtype>(gpu_result) < HALF);
}

booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m)
{
  const realtype *cdata = sunrajavec::getDevData<realtype, long int>(c);
  const realtype *xdata = sunrajavec::getDevData<realtype, long int>(x);
  const long int N = sunrajavec::getSize<realtype, long int>(x);
  realtype *mdata = sunrajavec::getDevData<realtype, long int>(m);

  RAJA::ReduceSum<RAJA::cuda_reduce<128>, realtype> gpu_result(ZERO);
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    bool test = (abs(cdata[i]) > ONEPT5 && cdata[i]*xdata[i] <= ZERO) ||
                (abs(cdata[i]) > HALF   && cdata[i]*xdata[i] <  ZERO);
    mdata[i] = test ? ONE : ZERO;
    gpu_result += mdata[i];
  });

  return (static_cast<realtype>(gpu_result) < HALF);
}

realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom)
{
  const realtype *ndata = sunrajavec::getDevData<realtype, long int>(num);
  const realtype *ddata = sunrajavec::getDevData<realtype, long int>(denom);
  const long int N = sunrajavec::getSize<realtype, long int>(num);

  RAJA::ReduceMin<RAJA::cuda_reduce<128>, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall<RAJA::cuda_exec<128> >(0, N, [=] __device__(long int i) {
    if (ddata[i] != ZERO)
      gpu_result.min(ndata[i]/ddata[i]);
  });

  return (static_cast<realtype>(gpu_result));
}


} // extern "C"
