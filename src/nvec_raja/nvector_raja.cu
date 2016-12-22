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


#include <nvector/raja/Vector.hpp>
#include <nvector/raja/VectorKernels.cuh>


/// Vector extractor
inline rvec::Vector<double, long int>* extract_raja(N_Vector v)
{ 
    return static_cast<rvec::Vector<double, long int>*>(v->content);
}



extern "C" {

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

//   ops->nvclone           = N_VClone_Raja;
//   ops->nvcloneempty      = N_VCloneEmpty_Raja;
  ops->nvdestroy         = N_VDestroy_Raja;
//   ops->nvspace           = N_VSpace_Raja;
//   ops->nvgetarraypointer = N_VGetArrayPointer_Raja;
//   ops->nvsetarraypointer = N_VSetArrayPointer_Raja;
  ops->nvlinearsum       = N_VLinearSum_Raja;
//   ops->nvconst           = N_VConst_Raja;
//   ops->nvprod            = N_VProd_Raja;
//   ops->nvdiv             = N_VDiv_Raja;
//   ops->nvscale           = N_VScale_Raja;
//   ops->nvabs             = N_VAbs_Raja;
//   ops->nvinv             = N_VInv_Raja;
//   ops->nvaddconst        = N_VAddConst_Raja;
  ops->nvdotprod         = N_VDotProd_Raja;
//   ops->nvmaxnorm         = N_VMaxNorm_Raja;
//   ops->nvwrmsnormmask    = N_VWrmsNormMask_Raja;
//   ops->nvwrmsnorm        = N_VWrmsNorm_Raja;
//   ops->nvmin             = N_VMin_Raja;
//   ops->nvwl2norm         = N_VWL2Norm_Raja;
//   ops->nvl1norm          = N_VL1Norm_Raja;
//   ops->nvcompare         = N_VCompare_Raja;
//   ops->nvinvtest         = N_VInvTest_Raja;
//   ops->nvconstrmask      = N_VConstrMask_Raja;
//   ops->nvminquotient     = N_VMinQuotient_Raja;

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

  v->content = new rvec::Vector<double, long int>(length);

  return(v);
}


N_Vector N_VMake_Raja(N_VectorContent_Raja c)
{
  N_Vector v;
  rvec::Vector<double, long int>* x = static_cast<rvec::Vector<double, long int>*>(c);
  long int length = x->size();

  v = NULL;
  v = N_VNewEmpty_Raja(length);
  if (v == NULL) return(NULL);

  v->content = c;

  return(v);
}


void N_VDestroy_Raja(N_Vector v)
{
//   if (NV_OWN_DATA_S(v) == TRUE) {
//     free(NV_DATA_S(v));
//     NV_DATA_S(v) = NULL;
//   }
//   free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VLinearSum_Raja(double a, N_Vector X, double b, N_Vector Y, N_Vector Z)
{
    rvec::linearSum(a, *extract_raja(X), b, *extract_raja(Y), *extract_raja(Z));
}

double N_VDotProd_Cuda(N_Vector X, N_Vector Y)
{
    return (rvec::dotProd(*extract_raja(X), *extract_raja(Y)));
}


} // extern "C"
