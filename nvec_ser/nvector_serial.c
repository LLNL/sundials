/*******************************************************************
 *                                                                 *
 * File          : nvector_serial.c                                *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh,              *
 *                 Radu Serban, and Allan G. Taylor, LLNL          *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for a serial implementation     *
 * of the NVECTOR package. It contains the implementation of       *
 * the serial vector specification intialization and free          *
 * routines (and of the Fortran callable interfaces to them)       *
 * and of the N_Vector kernels listed in nvector_serial.h.         *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nvector_serial.h"
#include "sundialstypes.h"
#include "sundialsmath.h" 

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define ONEPT5 RCONST(1.5)


/* Private Helper Prototypes */
/* z=x */
static void VCopy_Serial(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_Serial(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_Serial(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_Serial(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_Serial(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_Serial(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_Serial(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_Serial(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_Serial(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_Serial(realtype a, N_Vector x);

/********************* Exported Functions ************************/

/* Serial implementation of the vector specification 
   initialization routine */

NV_Spec NV_SpecInit_Serial(long int vec_length)
{
  NV_Spec nvspec;

  /* Create vector specification structure */
  nvspec = (NV_Spec) malloc(sizeof *nvspec);
  if (nvspec == NULL) return(NULL);

  /* Create serial content of vector specification structure */
  nvspec->content = (NV_SpecContent_Serial) malloc(sizeof(struct _NV_SpecContent_Serial));
  if (nvspec->content == NULL) {
    free(nvspec);
    return(NULL);
  }

  /* Load serial content of vector specification structure */
  NS_CONTENT_S(nvspec)->length = vec_length;

  /* Attach vector operations */
  nvspec->ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (nvspec->ops == NULL) {
    free(nvspec->content);
    free(nvspec);
    return(NULL);
  }

  nvspec->ops->nvnew           = N_VNew_Serial;
  nvspec->ops->nvfree          = N_VFree_Serial;
  nvspec->ops->nvspace         = N_VSpace_Serial;
  nvspec->ops->nvmake          = N_VMake_Serial;
  nvspec->ops->nvdispose       = N_VDispose_Serial;
  nvspec->ops->nvgetdata       = N_VGetData_Serial;
  nvspec->ops->nvsetdata       = N_VSetData_Serial;
  nvspec->ops->nvlinearsum     = N_VLinearSum_Serial;
  nvspec->ops->nvconst         = N_VConst_Serial;
  nvspec->ops->nvprod          = N_VProd_Serial;
  nvspec->ops->nvdiv           = N_VDiv_Serial;
  nvspec->ops->nvscale         = N_VScale_Serial;
  nvspec->ops->nvabs           = N_VAbs_Serial;
  nvspec->ops->nvinv           = N_VInv_Serial;
  nvspec->ops->nvaddconst      = N_VAddConst_Serial;
  nvspec->ops->nvdotprod       = N_VDotProd_Serial;
  nvspec->ops->nvmaxnorm       = N_VMaxNorm_Serial;
  nvspec->ops->nvwrmsnormmask  = N_VWrmsNormMask_Serial;
  nvspec->ops->nvwrmsnorm      = N_VWrmsNorm_Serial;
  nvspec->ops->nvmin           = N_VMin_Serial;
  nvspec->ops->nvwl2norm       = N_VWL2Norm_Serial;
  nvspec->ops->nvl1norm        = N_VL1Norm_Serial;
  nvspec->ops->nvcompare       = N_VCompare_Serial;
  nvspec->ops->nvinvtest       = N_VInvTest_Serial;
  nvspec->ops->nvconstrprodpos = N_VConstrProdPos_Serial;
  nvspec->ops->nvconstrmask    = N_VConstrMask_Serial;
  nvspec->ops->nvminquotient   = N_VMinQuotient_Serial;
  nvspec->ops->nvprint         = N_VPrint_Serial;

  /* Attach ID tag */
  nvspec->tag = ID_TAG_S;

  return(nvspec);

}
 
/* Serial implementation of the vector specification
   free routine */

void NV_SpecFree_Serial(NV_Spec nvspec)
{
  if (nvspec == NULL) return;

  free(nvspec->content);
  free(nvspec->ops);
  free(nvspec);
}

/***************************************************************************/

/* BEGIN implementation of vector operations */

N_Vector N_VNew_Serial(NV_Spec nvspec)
{
  N_Vector v;
  long int length;

  if (nvspec == NULL) return(NULL);

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  v->content = (N_VectorContent_Serial) malloc(sizeof(struct _N_VectorContent_Serial));
  if (v->content == NULL) {
    free(v);
    return(NULL);
  }

  length = NS_CONTENT_S(nvspec)->length;

  NV_CONTENT_S(v)->data = (realtype *) malloc(length * sizeof(realtype));
  if(NV_CONTENT_S(v)->data == NULL) {
    free(v->content);
    free(v);
    return(NULL);
  }

  NV_CONTENT_S(v)->length = length;

  v->nvspec = nvspec;

  return(v);
}

void N_VSpace_Serial(NV_Spec nvspec, long int *lrw, long int *liw)
{
  *lrw = NS_CONTENT_S(nvspec)->length;
  *liw = 1;
}

void N_VFree_Serial(N_Vector v)
{
  free(NV_DATA_S(v));
  free(NV_CONTENT_S(v));
  free(v);
}


N_Vector N_VMake_Serial(realtype *v_data, NV_Spec nvspec)
{
  N_Vector v;
  long int length;

  if (nvspec == NULL) return(NULL);

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  v->content = (N_VectorContent_Serial) malloc(sizeof(struct _N_VectorContent_Serial));
  if (v->content == NULL) {
    free(v);
    return(NULL);
  }

  length = NS_CONTENT_S(nvspec)->length;

  NV_CONTENT_S(v)->data = v_data;

  NV_CONTENT_S(v)->length = length;

  v->nvspec = nvspec;

  return(v);
}

void N_VDispose_Serial(N_Vector v)
{
  free(NV_CONTENT_S(v));
  free(v);
}

realtype *N_VGetData_Serial(N_Vector v)
{
  realtype *v_data;
  v_data = NV_CONTENT_S(v)->data;
  return(v_data);
}

void N_VSetData_Serial(realtype *v_data, N_Vector v)
{
  NV_CONTENT_S(v)->data = v_data;
}

void N_VLinearSum_Serial(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Serial(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_Serial(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_Serial(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Serial(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Serial(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Serial(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_Serial(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_Serial(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++) 
    *zd++ = a * (*xd++) + b * (*yd++);
}


void N_VConst_Serial(realtype c, N_Vector z)
{
  long int i, N;
  realtype *zd;

  N  = NV_LENGTH_S(z);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++) 
    *zd++ = c;
}


void N_VProd_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) * (*yd++);
}


void N_VDiv_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) / (*yd++);
}


void N_VScale_Serial(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy_Serial(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Serial(x, z);
  } else if (c == -ONE) {
    VNeg_Serial(x, z);
  } else {
    N  = NV_LENGTH_S(x);
    xd = NV_DATA_S(x);
    zd = NV_DATA_S(z);
    for (i=0; i < N; i++) 
      *zd++ = c * (*xd++);
  }
}


void N_VAbs_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++, xd++, zd++)
    *zd = ABS(*xd);
}


void N_VInv_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = ONE / (*xd++);
}


void N_VAddConst_Serial(N_Vector x, realtype b, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++) 
    *zd++ = (*xd++) + b; 
}


realtype N_VDotProd_Serial(N_Vector x, N_Vector y)
{
  long int i, N;
  realtype sum = ZERO, *xd, *yd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);

  for (i=0; i < N; i++)
    sum += (*xd++) * (*yd++);
  
  return(sum);
}


realtype N_VMaxNorm_Serial(N_Vector x)
{
  long int i, N;
  realtype max = ZERO, *xd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i=0; i < N; i++, xd++) {
    if (ABS(*xd) > max) max = ABS(*xd);
  }
   
  return(max);
}


realtype N_VWrmsNorm_Serial(N_Vector x, N_Vector w)
{
  long int i, N;
  realtype sum = ZERO, prodi, *xd, *wd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  wd = NV_DATA_S(w);

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  return(RSqrt(sum / N));
}

realtype N_VWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id)
{
  long int i, N;
  realtype sum = ZERO, prodi, *xd, *wd, *idd;

  N  = NV_LENGTH_S(x);
  xd  = NV_DATA_S(x);
  wd  = NV_DATA_S(w);
  idd = NV_DATA_S(id);

  for (i=0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i] * wd[i];
      sum += prodi * prodi;
    }
  }

  return(RSqrt(sum / N));
}


realtype N_VMin_Serial(N_Vector x)
{
  long int i, N;
  realtype min, *xd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  min = xd[0];

  xd++;
  for (i=1; i < N; i++, xd++) {
    if ((*xd) < min) min = *xd;
  }

  return(min);
}


realtype N_VWL2Norm_Serial(N_Vector x, N_Vector w)
{
  long int i, N;
  realtype sum = ZERO, prodi, *xd, *wd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  wd = NV_DATA_S(w);

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  return(RSqrt(sum));
}


realtype N_VL1Norm_Serial(N_Vector x)
{
  long int i, N;
  realtype sum = ZERO, *xd;
  
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  
  for (i=0; i<N; i++)  
    sum += ABS(xd[i]);

  return(sum);
}


void N_VOneMask_Serial(N_Vector x)
{
  long int i, N;
  realtype *xd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i=0; i<N; i++,xd++) {
    if (*xd != ZERO) *xd = ONE;
  }
}


void N_VCompare_Serial(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++, xd++, zd++) {
    *zd = (ABS(*xd) >= c) ? ONE : ZERO;
  }
}


booleantype N_VInvTest_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++) {
    if (*xd == ZERO) return(FALSE);
    *zd++ = ONE / (*xd++);
  }

  return(TRUE);
}


booleantype N_VConstrProdPos_Serial(N_Vector c, N_Vector x)
{
  long int i, N;
  realtype  *xd, *cd;
  booleantype test;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  cd = NV_DATA_S(c);

  test = TRUE;

  for (i=0; i < N; i++, xd++,cd++) {
    if (*cd != ZERO) {
      if ((*xd)*(*cd) <= ZERO) {
        test = FALSE;
        break;
      }
    }
  }
  return(test);
}


booleantype N_VConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m)
{
  long int i, N;
  booleantype test;
  realtype *cd, *xd, *md;
  
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  cd = NV_DATA_S(c);
  md = NV_DATA_S(m);

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
  return(test);
}


realtype N_VMinQuotient_Serial(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  long int i, N;
  realtype *nd, *dd, min;

  N  = NV_LENGTH_S(num);
  nd = NV_DATA_S(num);
  dd = NV_DATA_S(denom);

  notEvenOnce = TRUE;

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
  
  return(min);
}

 
void N_VPrint_Serial(N_Vector x)
{
  long int i, N;
  realtype *xd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i=0; i < N; i++) printf("%11.8g\n", *xd++);

  printf("\n");
}


/***************** Private Helper Functions **********************/


static void VCopy_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = *xd++; 
}


static void VSum_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) + (*yd++);
}


static void VDiff_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
 
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = (*xd++) - (*yd++);
}


static void VNeg_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);
  
  for (i=0; i < N; i++)
    *zd++ = -(*xd++);
}


static void VScaleSum_Serial(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = c * ((*xd++) + (*yd++));
}


static void VScaleDiff_Serial(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = c * ((*xd++) - (*yd++));
}


static void VLin1_Serial(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = a * (*xd++) + (*yd++);
}


static void VLin2_Serial(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i=0; i < N; i++)
    *zd++ = a * (*xd++) - (*yd++);
}

static void Vaxpy_Serial(realtype a, N_Vector x, N_Vector y)
{
  long int i, N;
  realtype *xd, *yd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);

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

static void VScaleBy_Serial(realtype a, N_Vector x)
{
  long int i, N;
  realtype *xd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i=0; i < N; i++)
    *xd++ *= a;
}
