/*******************************************************************
 *                                                                 *
 * File          : nvector_parallel.c                              *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh,              *
 *                 Radu Serban, and Allan G. Taylor, LLNL          *
 * Version of    : 06 June 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for a parallel implementation   *
 * of the NVECTOR package. It contains the implementation of       *
 * the parallel vector specification intialization and free        *
 * routines (and of the Fortran callable interfaces to them)       *
 * and of the N_Vector kernels listed in nvector_parallel.h.       *
 * It uses MPI for message-passing.                                *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nvector_parallel.h"
#include "sundialstypes.h"
#include "sundialsmath.h" 

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error message */

#define BAD_N1  "NV_SpecInit_Parallel -- Sum of local vector lengths differs from "
#define BAD_N2  "input global length. \n\n"
#define BAD_N    BAD_N1 BAD_N2

/* Private Helper Prototypes */

/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_Parallel(realtype d, int op, NV_Spec nvspec);
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


/********************* Exported Functions ************************/

/* Parallel implementation of the vector specification
   initialization routine */

NV_Spec NV_SpecInit_Parallel(MPI_Comm comm,  integertype local_vec_length, 
                             integertype global_vec_length, int *argc, char ***argv)
{
  NV_Spec nvspec;
  int initflag, initerr;
  integertype n, Nsum;

  /* Create vector specification structure */
  nvspec = (NV_Spec) malloc (sizeof *nvspec);
  if (nvspec == NULL) return(NULL);
  
  /* Create parallel content of vector specification structure */
  nvspec->content = (NV_SpecContent_Parallel) malloc(sizeof(struct _NV_SpecContent_Parallel));
  if (nvspec->content == NULL) {
    free(nvspec);
    return(NULL);
  }

  /* Load parallel content of vector specification structure */
  NS_CONTENT_P(nvspec)->local_vec_length = local_vec_length;
  NS_CONTENT_P(nvspec)->global_vec_length = global_vec_length;
  
  MPI_Initialized(&initflag);
  if (!initflag) {
    initerr = MPI_Init(argc,argv);
    if (initerr != MPI_SUCCESS) return(NULL);
  }
  NS_CONTENT_P(nvspec)->init_by_user = initflag;
  
  NS_CONTENT_P(nvspec)->comm = comm;
  
  /* Attach vector operations */
  nvspec->ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (nvspec->ops == NULL) {
    free(nvspec->content);
    free(nvspec);
    return(NULL);
  }

  nvspec->ops->nvnew           = N_VNew_Parallel;
  nvspec->ops->nvfree          = N_VFree_Parallel;
  nvspec->ops->nvspace         = N_VSpace_Parallel;
  nvspec->ops->nvmake          = N_VMake_Parallel;
  nvspec->ops->nvdispose       = N_VDispose_Parallel;
  nvspec->ops->nvgetdata       = N_VGetData_Parallel;
  nvspec->ops->nvsetdata       = N_VSetData_Parallel;
  nvspec->ops->nvlinearsum     = N_VLinearSum_Parallel;
  nvspec->ops->nvconst         = N_VConst_Parallel;
  nvspec->ops->nvprod          = N_VProd_Parallel;
  nvspec->ops->nvdiv           = N_VDiv_Parallel;
  nvspec->ops->nvscale         = N_VScale_Parallel;
  nvspec->ops->nvabs           = N_VAbs_Parallel;
  nvspec->ops->nvinv           = N_VInv_Parallel;
  nvspec->ops->nvaddconst      = N_VAddConst_Parallel;
  nvspec->ops->nvdotprod       = N_VDotProd_Parallel;
  nvspec->ops->nvmaxnorm       = N_VMaxNorm_Parallel;
  nvspec->ops->nvwrmsnorm      = N_VWrmsNorm_Parallel;
  nvspec->ops->nvwrmsnormmask  = N_VWrmsNormMask_Parallel;
  nvspec->ops->nvmin           = N_VMin_Parallel;
  nvspec->ops->nvwl2norm       = N_VWL2Norm_Parallel;
  nvspec->ops->nvl1norm        = N_VL1Norm_Parallel;
  nvspec->ops->nvcompare       = N_VCompare_Parallel;
  nvspec->ops->nvinvtest       = N_VInvTest_Parallel;
  nvspec->ops->nvconstrprodpos = N_VConstrProdPos_Parallel;
  nvspec->ops->nvconstrmask    = N_VConstrMask_Parallel;
  nvspec->ops->nvminquotient   = N_VMinQuotient_Parallel;
  nvspec->ops->nvprint         = N_VPrint_Parallel;

  /* If local length is negative, return now */
  if (local_vec_length < 0) return(nvspec);
  
  /* Compute global length as sum of local lengths */
  n = local_vec_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  NS_CONTENT_P(nvspec)->global_vec_length = Nsum;
  
  /* Check input global length against computed value */
  if (Nsum != global_vec_length) {
    printf(BAD_N);
    NV_SpecFree_Parallel(nvspec);
    return(NULL);
  } 

  /* Attach ID tag */
  nvspec->tag = ID_TAG_P;
  
  /* Return the vector specification */
  return(nvspec);
}

/* Parallel implementation of the vector specification
   free routine */

void NV_SpecFree_Parallel(NV_Spec nvspec)
{
  if (nvspec == NULL) return;

  if (!(NS_CONTENT_P(nvspec)->init_by_user)) MPI_Finalize();

  free(nvspec->content);
  free(nvspec->ops);
  free(nvspec);
}
 
/***************************************************************************/

/* BEGIN implementation of vector operations */

N_Vector N_VNew_Parallel(NV_Spec nvspec)
{
  N_Vector v;
  int N_local, N_global;

  if (nvspec == NULL) return(NULL);

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  v->content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (v->content == NULL) {
    free(v);
    return(NULL);
  }

  N_local  = NS_CONTENT_P(nvspec)->local_vec_length;
  N_global = NS_CONTENT_P(nvspec)->global_vec_length;

  NV_CONTENT_P(v)->data = (realtype *) malloc(N_local * sizeof(realtype));
  if (NV_CONTENT_P(v)->data == NULL) {
    free(v->content);
    free(v);
    return(NULL);
  }

  NV_CONTENT_P(v)->local_length  = N_local;
  NV_CONTENT_P(v)->global_length = N_global;

  v->nvspec = nvspec;
  
  return(v);
}

void N_VSpace_Parallel(NV_Spec nvspec, long int *lrw, long int *liw)
{
  *lrw = NS_CONTENT_P(nvspec)->global_vec_length;
  *liw = 2;
}

void N_VFree_Parallel(N_Vector v)
{
  free(NV_DATA_P(v));
  free(NV_CONTENT_P(v));
  free(v);
}

N_Vector N_VMake_Parallel(realtype *v_data, NV_Spec nvspec)
{
  N_Vector v;
  int N_local, N_global;

  if (nvspec == NULL) return(NULL);

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  v->content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (v->content == NULL) {
    free(v);
    return(NULL);
  }

  N_local  = NS_CONTENT_P(nvspec)->local_vec_length;
  N_global = NS_CONTENT_P(nvspec)->global_vec_length;

  NV_CONTENT_P(v)->data = v_data;

  NV_CONTENT_P(v)->local_length  = N_local;
  NV_CONTENT_P(v)->global_length = N_global;

  v->nvspec = nvspec;
  
  return(v);
}

void N_VDispose_Parallel(N_Vector v)
{
  free(NV_CONTENT_P(v));
  free(v);
}

realtype *N_VGetData_Parallel(N_Vector v)
{
  realtype *v_data;
  v_data = NV_CONTENT_P(v)->data;
  return(v_data);
}

void N_VSetData_Parallel(realtype *v_data, N_Vector v)
{
  NV_CONTENT_P(v)->data = v_data;
}

void N_VLinearSum_Parallel(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  integertype i, N;
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
  integertype i, N;
  realtype *zd;

  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++) 
    *zd++ = c;
}


void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++, xd++, zd++)
    *zd = ABS(*xd);
}


void N_VInv_Parallel(N_Vector x, N_Vector z)
{
  integertype i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = ONE / (*xd++);
}


void N_VAddConst_Parallel(N_Vector x, realtype b, N_Vector z)
{
  integertype i, N;
  realtype *xd, *zd;
  
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);
  
  for (i=0; i < N; i++) *zd++ = (*xd++) + b; 
}


realtype N_VDotProd_Parallel(N_Vector x, N_Vector y)
{
  integertype i, N;
  realtype sum = ZERO, *xd, *yd, gsum;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  nvspec = x->nvspec; 

  for (i=0; i < N; i++) sum += xd[i] * yd[i];

  gsum = VAllReduce_Parallel(sum, 1, nvspec);
  return(gsum);
}


realtype N_VMaxNorm_Parallel(N_Vector x)
{
  integertype i, N;
  realtype max = ZERO, *xd, gmax;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  nvspec = x->nvspec;

  for (i=0; i < N; i++, xd++) {
    if (ABS(*xd) > max) max = ABS(*xd);
  }
   
  gmax = VAllReduce_Parallel(max, 2, nvspec);
  return(gmax);
}


realtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w)
{
  integertype i, N, N_global;
  realtype sum = ZERO, prodi, *xd, *wd, gsum;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  nvspec = x->nvspec;

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  gsum = VAllReduce_Parallel(sum, 1, nvspec);
  return(RSqrt(gsum / N_global));
}

realtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id)
{
  integertype i, N, N_global;
  realtype sum = ZERO, prodi, *xd, *wd, *idd, gsum;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  idd = NV_DATA_P(id);
  nvspec = x->nvspec;

  for (i=0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i] * wd[i];
      sum += prodi * prodi;
    }
  }

  gsum = VAllReduce_Parallel(sum, 1, nvspec);
  return(RSqrt(gsum / N_global));
}

realtype N_VMin_Parallel(N_Vector x)
{
  integertype i, N;
  realtype min, *xd, gmin;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  nvspec = x->nvspec;

  min = xd[0];

  xd++;
  for (i=1; i < N; i++, xd++) {
    if ((*xd) < min) min = *xd;
  }

  gmin = VAllReduce_Parallel(min, 3, nvspec);
  return(gmin);
}


realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w)
{
  integertype i, N;
  realtype sum = ZERO, prodi, *xd, *wd, gsum;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  nvspec = x->nvspec;

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  gsum = VAllReduce_Parallel(sum, 1, nvspec);
  return(RSqrt(gsum));
}


realtype N_VL1Norm_Parallel(N_Vector x)
{
  integertype i, N;
  realtype sum = ZERO, gsum, *xd;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  nvspec = x->nvspec;

  for (i=0; i<N; i++,xd++) 
    sum += ABS(*xd);

  gsum = VAllReduce_Parallel(sum, 1, nvspec);
  return(gsum);
}


void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z)
{
  integertype i, N;
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
  integertype i, N;
  realtype *xd, *zd, val, gval;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);
  nvspec = x->nvspec; 

  val = ONE;
  for (i=0; i < N; i++) {
    if (*xd == ZERO) 
      val = ZERO;
    else
      *zd++ = ONE / (*xd++);
  }

  gval = VAllReduce_Parallel(val, 3, nvspec);
  if (gval == ZERO)
    return(FALSE);
  else
    return(TRUE);
}


booleantype N_VConstrProdPos_Parallel(N_Vector c, N_Vector x)
{
  integertype i, N;
  booleantype test;
  realtype  *xd, *cd;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  cd = NV_DATA_P(c);
  nvspec = x->nvspec;

  test = TRUE;

  for (i=0; i < N; i++, xd++,cd++) {
    if (*cd != ZERO) {
      if ( (*xd)*(*cd) <= ZERO) {
        test = FALSE;
        break;
      }
    }
  }

  return((booleantype)VAllReduce_Parallel((realtype)test, 3, nvspec));
}


booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m)
{
  integertype i, N;
  booleantype test;
  realtype *cd, *xd, *md;
  NV_Spec nvspec;
 
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  cd = NV_DATA_P(c);
  md = NV_DATA_P(m);
  nvspec = x->nvspec;

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

  return((booleantype)VAllReduce_Parallel((realtype)test, 3, nvspec));
}


realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  integertype i, N;
  realtype *nd, *dd, min;
  NV_Spec nvspec;

  N  = NV_LOCLENGTH_P(num);
  nd = NV_DATA_P(num);
  dd = NV_DATA_P(denom);
  nvspec = num->nvspec;

  notEvenOnce = TRUE;

  min = 0.0;
  for (i=0; i<N; i++, nd++, dd++) {
    if (*dd == ZERO) continue;
    else {
      if (notEvenOnce) {
        min = *nd / *dd ;
        notEvenOnce = FALSE;
      }
      else 
        min = MIN(min, (*nd)/(*dd));
    }
  }
  if (notEvenOnce) min = 1.e99;
  
  return(VAllReduce_Parallel(min, 3, nvspec));
}

 
void N_VPrint_Parallel(N_Vector x)
{
  integertype i, N;
  realtype *xd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i=0; i < N; i++) printf("%g\n", *xd++);

  printf("\n");
}


/***************** Private Helper Functions **********************/

static realtype VAllReduce_Parallel(realtype d, int op, NV_Spec nvspec)
{
  /* This function does a global reduction.  The operation is
       sum if op = 1,
       max if op = 2,
       min if op = 3.
     The operation is over all processors in the group defined by
     the parameters within nvspec. */

  MPI_Comm comm;
  realtype out;

  comm = NS_CONTENT_P(nvspec)->comm;

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
  integertype i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = *xd++; 
}


static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
  realtype *xd, *zd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i=0; i < N; i++)
    *zd++ = -(*xd++);
}


static void VScaleSum_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
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
  integertype i, N;
  realtype *xd;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i=0; i < N; i++)
    *xd++ *= a;
}
