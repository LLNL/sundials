/*
 * -----------------------------------------------------------------
 * $Revision: 1.20 $
 * $Date: 2005-01-24 22:29:12 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for a parallel MPI implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "nvector_parallel.h"
#include "sundialsmath.h"
#include "sundialstypes.h"

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error Message */

#define BAD_N1 "N_VNew_Parallel -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/* Private function prototypes */

/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_Parallel(realtype d, int op, MPI_Comm comm);
/* z=x */
static void VCopy_Parallel(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_Parallel(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_Parallel(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_Parallel(realtype a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Function to create a new parallel vector with empty data array
 */

N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, 
                              long int local_length,
                              long int global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel content;
  long int n, Nsum;

  /* Compute global length as sum of local lengths */
  n = local_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  if (Nsum != global_length) {
    printf(BAD_N);
    return(NULL);
  } 

  /* Create vector */
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvclone           = N_VClone_Parallel;
  ops->nvdestroy         = N_VDestroy_Parallel;
  ops->nvspace           = N_VSpace_Parallel;
  ops->nvgetarraypointer = N_VGetArrayPointer_Parallel;
  ops->nvsetarraypointer = N_VSetArrayPointer_Parallel;
  ops->nvlinearsum       = N_VLinearSum_Parallel;
  ops->nvconst           = N_VConst_Parallel;
  ops->nvprod            = N_VProd_Parallel;
  ops->nvdiv             = N_VDiv_Parallel;
  ops->nvscale           = N_VScale_Parallel;
  ops->nvabs             = N_VAbs_Parallel;
  ops->nvinv             = N_VInv_Parallel;
  ops->nvaddconst        = N_VAddConst_Parallel;
  ops->nvdotprod         = N_VDotProd_Parallel;
  ops->nvmaxnorm         = N_VMaxNorm_Parallel;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Parallel;
  ops->nvwrmsnorm        = N_VWrmsNorm_Parallel;
  ops->nvmin             = N_VMin_Parallel;
  ops->nvwl2norm         = N_VWL2Norm_Parallel;
  ops->nvl1norm          = N_VL1Norm_Parallel;
  ops->nvcompare         = N_VCompare_Parallel;
  ops->nvinvtest         = N_VInvTest_Parallel;
  ops->nvconstrmask      = N_VConstrMask_Parallel;
  ops->nvminquotient     = N_VMinQuotient_Parallel;

  /* Create content */
  content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length = local_length;
  content->global_length = global_length;
  content->comm = comm;
  content->own_data = FALSE;
  content->data = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops = ops;

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a new parallel vector
 */

N_Vector N_VNew_Parallel(MPI_Comm comm, 
                         long int local_length,
                         long int global_length)
{
  N_Vector v;
  realtype *data;

  v = N_VNewEmpty_Parallel(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) {N_VDestroy_Parallel(v);return(NULL);}

    /* Attach data */
    NV_OWN_DATA_P(v) = TRUE;
    NV_DATA_P(v) = data; 

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to clone from a template a new vector with empty (NULL) data array
 */

N_Vector N_VCloneEmpty_Parallel(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }
  
  ops->nvclone           = w->ops->nvclone;
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
  content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = NV_LOCLENGTH_P(w);
  content->global_length = NV_GLOBLENGTH_P(w);
  content->comm = NV_COMM_P(w);
  content->own_data = FALSE;
  content->data = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops = ops;

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parallel N_Vector with user data component 
 */

N_Vector N_VMake_Parallel(MPI_Comm comm, 
                          long int local_length,
                          long int global_length,
                          realtype *v_data)
{
  N_Vector v;

  v = N_VNewEmpty_Parallel(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  if (local_length > 0) {
    /* Attach data */
    NV_OWN_DATA_P(v) = FALSE;
    NV_DATA_P(v) = v_data;
  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parallel vectors. 
 */

N_Vector *N_VNewVectorArray_Parallel(int count, 
                                     MPI_Comm comm, 
                                     long int local_length,
                                     long int global_length)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j=0; j<count; j++) {
    vs[j] = N_VNew_Parallel(comm, local_length, global_length);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parallel vectors with empty
 * (NULL) data array.
 */

N_Vector *N_VNewVectorArrayEmpty_Parallel(int count, 
                                          MPI_Comm comm, 
                                          long int local_length,
                                          long int global_length)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j=0; j<count; j++) {
    vs[j] = N_VNewEmpty_Parallel(comm, local_length, global_length);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VNewVectorArray_Parallel
 */

void N_VDestroyVectorArray_Parallel(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Parallel(vs[j]);

  free(vs);
}

/* ---------------------------------------------------------------- 
 * Function to print a parallel vector 
 */

void N_VPrint_Parallel(N_Vector x)
{
  long int i, N;
  realtype *xd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i=0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%Lg\n", *xd++);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%lg\n", *xd++);
#else
    printf("%g\n", *xd++);
#endif
  }
  printf("\n");
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VClone_Parallel(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int local_length;

  v = N_VCloneEmpty_Parallel(w);
  if (v == NULL) return(NULL);

  local_length  = NV_LOCLENGTH_P(w);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) {N_VDestroy_Parallel(v);return(NULL);}

    /* Attach data */
    NV_OWN_DATA_P(v) = TRUE;
    NV_DATA_P(v) = data;
  }

  return(v);
}

void N_VDestroy_Parallel(N_Vector v)
{
  if ( (NV_OWN_DATA_P(v) == TRUE) && (NV_DATA_P(v) != NULL) ) 
    free(NV_DATA_P(v));
  free(v->content);
  free(v->ops);
  free(v);
}

void N_VSpace_Parallel(N_Vector v, long int *lrw, long int *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_COMM_P(v);
  MPI_Comm_size(comm, &npes);
  
  *lrw = NV_GLOBLENGTH_P(v);
  *liw = 2 * npes;
}

realtype *N_VGetArrayPointer_Parallel(N_Vector v)
{
  realtype *v_data;

  v_data = NV_DATA_P(v);

  return(v_data);
}

void N_VSetArrayPointer_Parallel(realtype *v_data, N_Vector v)
{
  if (NV_LOCLENGTH_P(v) > 0) NV_DATA_P(v) = v_data;
}

void N_VLinearSum_Parallel(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Parallel(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_Parallel(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_Parallel(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Parallel(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Parallel(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Parallel(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_Parallel(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_Parallel(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++) 
    *zd++ = a * (*xd++) + b * (*yd++);
}

void N_VConst_Parallel(realtype c, N_Vector z)
{
  long int i, N;
  realtype *zd;

  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++) 
    *zd++ = c;
}

void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) * (*yd++);
}

void N_VDiv_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) / (*yd++);
}

void N_VScale_Parallel(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy_Parallel(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Parallel(x, z);
  } else if (c == -ONE) {
    VNeg_Parallel(x, z);
  } else {
    N  = NV_LOCLENGTH_P(x);
    xd = NV_DATA_P(x);
    zd = NV_DATA_P(z);
    for (i=0; i < N; i++) 
      *zd++ = c * (*xd++);
  }
}

void N_VAbs_Parallel(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++, xd++, zd++)
    *zd = ABS(*xd);
}

void N_VInv_Parallel(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = ONE / (*xd++);
}

void N_VAddConst_Parallel(N_Vector x, realtype b, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);
  
  for (i=0; i < N; i++) *zd++ = (*xd++) + b; 
}

realtype N_VDotProd_Parallel(N_Vector x, N_Vector y)
{
  long int i, N;
  realtype sum = ZERO, *xd, *yd, gsum;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  comm = NV_COMM_P(x);

  for (i=0; i < N; i++) sum += xd[i] * yd[i];

  gsum = VAllReduce_Parallel(sum, 1, comm);
  return(gsum);
}

realtype N_VMaxNorm_Parallel(N_Vector x)
{
  long int i, N;
  realtype max, *xd, gmax;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  max = ZERO;

  for (i=0; i < N; i++, xd++) {
    if (ABS(*xd) > max) max = ABS(*xd);
  }
   
  gmax = VAllReduce_Parallel(max, 2, comm);
  return(gmax);
}

realtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w)
{
  long int i, N, N_global;
  realtype sum = ZERO, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  comm = NV_COMM_P(x);

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  gsum = VAllReduce_Parallel(sum, 1, comm);
  return(RSqrt(gsum / N_global));
}

realtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id)
{
  long int i, N, N_global;
  realtype sum = ZERO, prodi, *xd, *wd, *idd, gsum;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  idd = NV_DATA_P(id);
  comm = NV_COMM_P(x);

  for (i=0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i] * wd[i];
      sum += prodi * prodi;
    }
  }

  gsum = VAllReduce_Parallel(sum, 1, comm);
  return(RSqrt(gsum / N_global));
}

realtype N_VMin_Parallel(N_Vector x)
{
  long int i, N;
  realtype min, *xd, gmin;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  comm = NV_COMM_P(x);

  min = BIG_REAL;

  if (N > 0) {

    xd = NV_DATA_P(x);
    
    min = xd[0];
    
    xd++;
    for (i=1; i < N; i++, xd++) {
      if ((*xd) < min) min = *xd;
    }

  }
    
  gmin = VAllReduce_Parallel(min, 3, comm);
  return(gmin);
}

realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w)
{
  long int i, N;
  realtype sum = ZERO, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  comm = NV_COMM_P(x);

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  gsum = VAllReduce_Parallel(sum, 1, comm);
  return(RSqrt(gsum));
}

realtype N_VL1Norm_Parallel(N_Vector x)
{
  long int i, N;
  realtype sum = ZERO, gsum, *xd;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  for (i=0; i<N; i++,xd++) 
    sum += ABS(*xd);

  gsum = VAllReduce_Parallel(sum, 1, comm);
  return(gsum);
}

void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++, xd++, zd++) {
    *zd = (ABS(*xd) >= c) ? ONE : ZERO;
  }
}

booleantype N_VInvTest_Parallel(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);
  comm = NV_COMM_P(x);

  val = ONE;
  for (i=0; i < N; i++) {
    if (*xd == ZERO) 
      val = ZERO;
    else
      *zd++ = ONE / (*xd++);
  }

  gval = VAllReduce_Parallel(val, 3, comm);
  if (gval == ZERO)
    return(FALSE);
  else
    return(TRUE);
}

booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m)
{
  long int i, N;
  booleantype test;
  realtype *cd, *xd, *md;
  MPI_Comm comm;
 
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  cd = NV_DATA_P(c);
  md = NV_DATA_P(m);
  comm = NV_COMM_P(x);

  test = TRUE;

  for (i=0; i<N; i++, cd++, xd++, md++) {
    *md = ZERO;
    if (*cd == ZERO) continue;
    if (*cd > ONEPT5 || (*cd) < -ONEPT5) {
      if ( (*xd)*(*cd) <= ZERO) {test = FALSE; *md = ONE; }
      continue;
    }
    if ( (*cd) > HALF || (*cd) < -HALF) {
      if ( (*xd)*(*cd) < ZERO ) {test = FALSE; *md = ONE; }
    }
  }

  return((booleantype)VAllReduce_Parallel((realtype)test, 3, comm));
}

realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  long int i, N;
  realtype *nd, *dd, min;
  MPI_Comm comm;

  N  = NV_LOCLENGTH_P(num);
  nd = NV_DATA_P(num);
  dd = NV_DATA_P(denom);
  comm = NV_COMM_P(num);

  notEvenOnce = TRUE;

  for (i = 0; i < N; i++, nd++, dd++) {
    if (*dd == ZERO) continue;
    else {
      if (notEvenOnce) {
        min = *nd / *dd ;
        notEvenOnce = FALSE;
      }
      else min = MIN(min, (*nd) / (*dd));
    }
  }

  if (notEvenOnce || (N == 0)) min = BIG_REAL;

  return(VAllReduce_Parallel(min, 3, comm));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static realtype VAllReduce_Parallel(realtype d, int op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator 
   */

  realtype out;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}

static void VCopy_Parallel(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = *xd++; 
}

static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) + (*yd++);
}

static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
 
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) - (*yd++);
}

static void VNeg_Parallel(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = -(*xd++);
}

static void VScaleSum_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = c * ((*xd++) + (*yd++));
}

static void VScaleDiff_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = c * ((*xd++) - (*yd++));
}

static void VLin1_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = a * (*xd++) + (*yd++);
}

static void VLin2_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = a * (*xd++) - (*yd++);
}

static void Vaxpy_Parallel(realtype a, N_Vector x, N_Vector y)
{
  long int i, N;
  realtype *xd, *yd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);

  if (a == ONE) {
    for (i=0; i < N; i++)
      *yd++ += (*xd++);
    return;
  }
  
  if (a == -ONE) {
    for (i=0; i < N; i++)
      *yd++ -= (*xd++);
    return;
  }    
  
  for (i=0; i < N; i++)
    *yd++ += a * (*xd++);
}

static void VScaleBy_Parallel(realtype a, N_Vector x)
{
  long int i, N;
  realtype *xd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i=0; i < N; i++)
    *xd++ *= a;
}
