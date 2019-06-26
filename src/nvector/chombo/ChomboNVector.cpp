#include "LevelData.H"
#include "FArrayBox.H"
#include "FORT_PROTO.H"

#include "sundials/sundials_types.h"

#include "ChomboNVector.H"

#include "N_VectorOps_F.H"

#include "ChomboHelper.H"

#include <cmath>
using std::sqrt;
#include <iostream>
using std::cout;
using std::endl;

enum MPIOpType { MPIOpSum, MPIOpMax, MPIOpMin, MPIOpLast };
static realtype VAllReduce_Parallel(realtype d, MPIOpType op, MPI_Comm comm);
static int      VAllReduce_Parallel(int d, MPIOpType op, MPI_Comm comm);

///// Operations /////
// forward declaration
N_Vector N_VCreate();

// creation
N_Vector N_VNew(const DisjointBoxLayout& dp, const int nComp, const IntVect& ghost) {
  N_Vector v = N_VCreate();
  NV_CONTENT_CH(v)->Data = new LevelData<FArrayBox>;
  NV_DATA_CH(v).define(dp, nComp, ghost);
  NV_OWNDATA_CH(v) = true;
  return v;
}

N_Vector N_VAssociate(LevelData<FArrayBox>* Data) {
  N_Vector v = N_VCreate();
  NV_CONTENT_CH(v)->Data = Data;
  NV_OWNDATA_CH(v) = false;
  return v;
}

// clone
N_Vector N_VCloneEmpty_Ch(N_Vector v) {
  return N_VCreate();
}

N_Vector N_VClone_Ch(N_Vector v) {
  return N_VNew(NV_DBL_CH(v), NV_NCOMP_CH(v), NV_GHOST_CH(v));
}

// destroy
void N_VDestroy_Ch(N_Vector v) {
  if (NV_OWNDATA_CH(v)) delete NV_CONTENT_CH(v)->Data;
  delete NV_CONTENT_CH(v); v->content = NULL;
  delete v->ops; v->ops = NULL;
  delete v; v = NULL;
}

// space check
void N_VSpace_Ch(N_Vector nvSpec, long int *lrw, long int *liw) {
 // dummy
}

// linear sum
void N_VLinearSum_Ch(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    NV_DATA_CH(z)[dit()].axby(NV_DATA_CH(x)[dit()], NV_DATA_CH(y)[dit()], a, b);
  }
}

// set to constant
void N_VConst_Ch(realtype c, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    NV_DATA_CH(z)[dit()].setVal(c);
  }
}

// product
void N_VProd_Ch(N_Vector x, N_Vector y, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(z)[dit()];
    const FArrayBox& xx  = NV_DATA_CH(x)[dit()];
    const FArrayBox& yy  = NV_DATA_CH(y)[dit()];
    FArrayBox& zz  = NV_DATA_CH(z)[dit()];
    FORT_ARRAYPROD(CHF_CONST_FRA(xx), CHF_CONST_FRA(yy), CHF_FRA(zz), CHF_BOX(region));
  }
}

// division
void N_VDiv_Ch(N_Vector x, N_Vector y, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(z)[dit()];
    const FArrayBox& xx  = NV_DATA_CH(x)[dit()];
    const FArrayBox& yy  = NV_DATA_CH(y)[dit()];
    FArrayBox& zz  = NV_DATA_CH(z)[dit()];
    FORT_ARRAYDIV(CHF_CONST_FRA(xx), CHF_CONST_FRA(yy), CHF_FRA(zz), CHF_BOX(region));
  }
}

// scale
void N_VScale_Ch(realtype c, N_Vector x, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(z)[dit()];
    const FArrayBox& xx  = NV_DATA_CH(x)[dit()];
    FArrayBox& zz  = NV_DATA_CH(z)[dit()];
    FORT_ARRAYSCL(CHF_CONST_FRA(xx), CHF_FRA(zz), CHF_BOX(region), CHF_CONST_REAL(c));
  }
}

// abs
void N_VAbs_Ch(N_Vector x, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(z)[dit()];
    const FArrayBox& xx  = NV_DATA_CH(x)[dit()];
    FArrayBox& zz  = NV_DATA_CH(z)[dit()];
    FORT_ARRAYABS(CHF_CONST_FRA(xx), CHF_FRA(zz), CHF_BOX(region));
  }
}

// inv
void N_VInv_Ch(N_Vector x, N_Vector z) {
  int nonzero = 1;
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    FArrayBox& zz = NV_DATA_CH(z)[dit()];
    FORT_INVWCHK(CHF_CONST_FRA(xx), CHF_FRA(zz), CHF_BOX(region), CHF_INT(nonzero));
  }
  // in this version of the function nonzero is not returned
}

// add constant
void N_VAddConst_Ch(N_Vector x, realtype b, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    FArrayBox& zz = NV_DATA_CH(z)[dit()];
    FORT_ADDCONST(CHF_CONST_FRA(xx), CHF_FRA(zz), CHF_BOX(region), CHF_CONST_REAL(b));
  }
}

// dot product
realtype N_VDotProd_Ch(N_Vector x, N_Vector y) {
  realtype prod = 0.;
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    // don't double-count cells by using ghost cell information
    const Box& region = NV_DBL_CH(x)[dit()];
    prod += NV_DATA_CH(x)[dit()].dotProduct(NV_DATA_CH(y)[dit()], region);
  }
  return VAllReduce_Parallel(prod, MPIOpSum, NV_COMM_CH(x));
}

// max norm
realtype N_VMaxNorm_Ch(N_Vector x) {
  realtype m;
  DataIterator dit = NV_DBLITER_CH(x);
  dit.reset();
  m = NV_DATA_CH(x)[dit()](NV_DBL_CH(x)[dit()].smallEnd());
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xb = NV_DATA_CH(x)[dit()];
    FORT_ARRAYMAXNORM(CHF_CONST_FRA(xb), CHF_BOX(region), CHF_REAL(m));
  }
  return VAllReduce_Parallel(m, MPIOpMax, NV_COMM_CH(x));
}

// weighted root-mean-square norm
realtype N_VWrmsNorm_Ch(N_Vector x, N_Vector w) {
  realtype norm = 0.;
  int n = 0;
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    const FArrayBox& ww = NV_DATA_CH(w)[dit()];
    FORT_WTDSQ(CHF_CONST_FRA(xx), CHF_CONST_FRA(ww), CHF_BOX(region), CHF_REAL(norm));
    n += region.size().product();
  }
  n    = VAllReduce_Parallel(n,    MPIOpSum, NV_COMM_CH(x));
  norm = VAllReduce_Parallel(norm, MPIOpSum, NV_COMM_CH(x));
  return sqrt(norm/static_cast<Real>(n));
}

// signed weighted root-mean-square norm
realtype N_VWrmsNormMask_Ch(N_Vector x, N_Vector w, N_Vector id) {
  realtype norm = 0.;
  int n = 0;
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx  = NV_DATA_CH(x) [dit()];
    const FArrayBox& ww  = NV_DATA_CH(w) [dit()];
    const FArrayBox& idx = NV_DATA_CH(id)[dit()];
    FORT_WTDSIGNSQ(CHF_CONST_FRA(xx), CHF_CONST_FRA(ww), CHF_CONST_FRA(idx), CHF_BOX(region), CHF_REAL(norm));
    n += region.size().product();
  }
  n    = VAllReduce_Parallel(n,    MPIOpSum, NV_COMM_CH(x));
  norm = VAllReduce_Parallel(norm, MPIOpSum, NV_COMM_CH(x));
  return sqrt(norm/static_cast<Real>(n));
}

// min
realtype N_VMin_Ch(N_Vector x) {
  realtype m;
  int i = 0;
  DataIterator dit = NV_DBLITER_CH(x);
  dit.reset();
  m = NV_DATA_CH(x)[dit()](NV_DBL_CH(x)[dit()].smallEnd());
  for (dit.reset(); dit.ok(); i++, ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xb = NV_DATA_CH(x)[dit()];
    FORT_ARRAYMIN(CHF_CONST_FRA(xb), CHF_BOX(region), CHF_REAL(m));
  }
  return VAllReduce_Parallel(m, MPIOpMin, NV_COMM_CH(x));
}

// weighted Euclidean l2 norm
realtype N_VWL2Norm_Ch(N_Vector x, N_Vector w) {
  realtype sq = 0.;
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    const FArrayBox& ww = NV_DATA_CH(w)[dit()];
    FORT_WTDSQ(CHF_CONST_FRA(xx), CHF_CONST_FRA(ww), CHF_BOX(region), CHF_REAL(sq));
  }
  sq = VAllReduce_Parallel(sq, MPIOpSum, NV_COMM_CH(x));
  return sqrt(sq);
}

// l1 norm
realtype N_VL1Norm_Ch(N_Vector x) {
  const int p=1;
  realtype norm = 0.;
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    // ignore ghost cells
    const Box& region = NV_DBL_CH(x)[dit()];
    for (int comp=0; comp<NV_NCOMP_CH(x); comp++)
      norm += NV_DATA_CH(x)[dit()].norm(region, p, comp);
  }
  return VAllReduce_Parallel(norm, MPIOpSum, NV_COMM_CH(x));
}

// compare
void N_VCompare_Ch(realtype c, N_Vector x, N_Vector z) {
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    FArrayBox& zz = NV_DATA_CH(z)[dit()];
    FORT_ARRAYCOMP(CHF_CONST_FRA(xx), CHF_FRA(zz), CHF_BOX(region), CHF_CONST_REAL(c));
  }
}

// inverse with check
booleantype N_VInvTest_Ch(N_Vector x, N_Vector z) {
  int nonzero = 1;
  DataIterator dit = NV_DBLITER_CH(z);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    FArrayBox& zz = NV_DATA_CH(z)[dit()];
    FORT_INVWCHK(CHF_CONST_FRA(xx), CHF_FRA(zz), CHF_BOX(region), CHF_INT(nonzero));
  }
  return static_cast<bool>( VAllReduce_Parallel(nonzero, MPIOpMin, NV_COMM_CH(z)) );
}

// constraint check
booleantype N_VConstrMask_Ch(N_Vector c, N_Vector x, N_Vector m) {
  int allpassed = 1;
  DataIterator dit = NV_DBLITER_CH(x);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(x)[dit()];
    const FArrayBox& cc = NV_DATA_CH(c)[dit()];
    const FArrayBox& xx = NV_DATA_CH(x)[dit()];
    FArrayBox& mm = NV_DATA_CH(m)[dit()];
    FORT_CONSTRCHK(CHF_CONST_FRA(cc), CHF_CONST_FRA(xx), CHF_FRA(mm), CHF_BOX(region), CHF_INT(allpassed));
  }
  return static_cast<bool>( VAllReduce_Parallel(allpassed, MPIOpMin, NV_COMM_CH(c)) );
}

// min quotient
realtype N_VMinQuotient_Ch(N_Vector num, N_Vector denom) {
  // Assume none exists
  realtype t, q = BIG_REAL;
  int n;
  DataIterator dit = NV_DBLITER_CH(num);
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& region = NV_DBL_CH(num)[dit()];
    const FArrayBox& x = NV_DATA_CH(num)[dit()];
    const FArrayBox& y = NV_DATA_CH(denom)[dit()];
    t = BIG_REAL;
    FORT_MINQUOT(CHF_CONST_FRA(x), CHF_CONST_FRA(y), CHF_BOX(region), CHF_REAL(t), CHF_INT(n));
    if (n == 1 && t < q) q = t;
  }
  return VAllReduce_Parallel(q, MPIOpMin, NV_COMM_CH(num));
}

N_Vector N_VCreate() {
  _generic_N_Vector *v  = new _generic_N_Vector;
  v->content = (void *)new _N_VectorContent_Ch;

  // assign the communicator
  NV_CONTENT_CH(v)->comm = Chombo_MPI::comm;

  // define operations
  NV_OPS_CH(v) = new _generic_N_Vector_Ops;
  NV_OPS_CH(v)->nvclone         = N_VClone_Ch;
  NV_OPS_CH(v)->nvcloneempty    = N_VCloneEmpty_Ch;
  NV_OPS_CH(v)->nvdestroy       = N_VDestroy_Ch;
  NV_OPS_CH(v)->nvspace         = N_VSpace_Ch;
  // Skipped some optional kernel functions here
  NV_OPS_CH(v)->nvlinearsum     = N_VLinearSum_Ch;
  NV_OPS_CH(v)->nvconst         = N_VConst_Ch;
  NV_OPS_CH(v)->nvprod          = N_VProd_Ch;
  NV_OPS_CH(v)->nvdiv           = N_VDiv_Ch;
  NV_OPS_CH(v)->nvscale         = N_VScale_Ch;
  NV_OPS_CH(v)->nvabs           = N_VAbs_Ch;
  NV_OPS_CH(v)->nvinv           = N_VInv_Ch;
  NV_OPS_CH(v)->nvaddconst      = N_VAddConst_Ch;
  NV_OPS_CH(v)->nvdotprod       = N_VDotProd_Ch;
  NV_OPS_CH(v)->nvmaxnorm       = N_VMaxNorm_Ch;
  NV_OPS_CH(v)->nvwrmsnorm      = N_VWrmsNorm_Ch;
  NV_OPS_CH(v)->nvwrmsnormmask  = N_VWrmsNormMask_Ch;
  NV_OPS_CH(v)->nvmin           = N_VMin_Ch;
  NV_OPS_CH(v)->nvwl2norm       = N_VWL2Norm_Ch;
  NV_OPS_CH(v)->nvl1norm        = N_VL1Norm_Ch;
  NV_OPS_CH(v)->nvcompare       = N_VCompare_Ch;
  NV_OPS_CH(v)->nvinvtest       = N_VInvTest_Ch;
  NV_OPS_CH(v)->nvconstrmask    = N_VConstrMask_Ch;
  NV_OPS_CH(v)->nvminquotient   = N_VMinQuotient_Ch;

  return v;
}

void N_VDataExchange(N_Vector v) {
  // using Chombo's exchange feature
  NV_DATA_CH(v).exchange(NV_DATA_CH(v).interval());
}

int N_VEquate(N_Vector x, N_Vector z) {
  if (!(NV_DBL_CH(x) == NV_DBL_CH(z))) return CHSUN_SIZE_ERROR;
  NV_DATA_CH(x).copyTo(NV_DATA_CH(z));
  return CHSUN_SUCCESS;
}

/********************/
/* Helper functions */
/********************/
// VAllReduce_Parallel overloaded for real and int types
static realtype VAllReduce_Parallel(realtype d, MPIOpType op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction. The operation 
   * is over all processors in the communicator 
   */

  realtype out;

#ifdef CH_MPI
  switch (op) {
   case MPIOpSum:
     MPI_Allreduce(&d, &out, 1, CH_REAL_MPI_TYPE, MPI_SUM, comm);
     break;
   case MPIOpMax:
     MPI_Allreduce(&d, &out, 1, CH_REAL_MPI_TYPE, MPI_MAX, comm);
     break;
   case MPIOpMin:
     MPI_Allreduce(&d, &out, 1, CH_REAL_MPI_TYPE, MPI_MIN, comm);
     break;
   default:
     break;
  }
#else
  out = d;
#endif

  return out;
}

static int VAllReduce_Parallel(int d, MPIOpType op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction. The operation 
   * is over all processors in the communicator 
   */

  int out;

#ifdef CH_MPI
  switch (op) {
   case MPIOpSum:
     MPI_Allreduce(&d, &out, 1, CH_INTEGER_MPI_TYPE, MPI_SUM, comm);
     break;
   case MPIOpMax:
     MPI_Allreduce(&d, &out, 1, CH_INTEGER_MPI_TYPE, MPI_MAX, comm);
     break;
   case MPIOpMin:
     MPI_Allreduce(&d, &out, 1, CH_INTEGER_MPI_TYPE, MPI_MIN, comm);
     break;
   default:
     break;
  }
#else
  out = d;
#endif

  return out;
}

