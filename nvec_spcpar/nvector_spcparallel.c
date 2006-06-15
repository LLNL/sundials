/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2006-06-15 15:39:50 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds and Radu Serban @LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a spc_parallel MPI 
 * implementation of the NVECTOR package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "nvector_spcparallel.h"
#include "sundials_math.h" 

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private function prototypes */

/* attach data and set pointers in data */
static void VAttach_Data(N_Vector v, realtype *data);
/* z=x */
static void VCopy_SpcParallel(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_SpcParallel(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_SpcParallel(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_SpcParallel(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_SpcParallel(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_SpcParallel(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_SpcParallel(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_SpcParallel(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_SpcParallel(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_SpcParallel(realtype a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* -----------------------------------------------------------------
 * Create a new spc_parallel vector with empty data array
 */

N_Vector N_VNewEmpty_SpcParallel(MPI_Comm comm, int Ngrp, int *Nspc,
                                 long int Nx,  long int Ny,  long int Nz, 
                                 long int NGx, long int NGy, long int NGz) 
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_SpcParallel content;
  int ig;
  long int tmp, n, N;

  /* Create the new vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) {
    free(v); 
    return(NULL); 
  }

  /* Attach custom vector routines to N_Vector_Ops structure */
  ops->nvclone           = N_VClone_SpcParallel;
  ops->nvcloneempty      = N_VCloneEmpty_SpcParallel;
  ops->nvdestroy         = N_VDestroy_SpcParallel;
  ops->nvspace           = N_VSpace_SpcParallel;
  ops->nvgetarraypointer = N_VGetArrayPointer_SpcParallel;
  ops->nvsetarraypointer = N_VSetArrayPointer_SpcParallel;
  ops->nvlinearsum       = N_VLinearSum_SpcParallel;
  ops->nvconst           = N_VConst_SpcParallel;
  ops->nvprod            = N_VProd_SpcParallel;
  ops->nvdiv             = N_VDiv_SpcParallel;
  ops->nvscale           = N_VScale_SpcParallel;
  ops->nvabs             = N_VAbs_SpcParallel;
  ops->nvinv             = N_VInv_SpcParallel;
  ops->nvaddconst        = N_VAddConst_SpcParallel;
  ops->nvdotprod         = N_VDotProd_SpcParallel;
  ops->nvmaxnorm         = N_VMaxNorm_SpcParallel;
  ops->nvwrmsnorm        = N_VWrmsNorm_SpcParallel;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_SpcParallel;
  ops->nvmin             = N_VMin_SpcParallel;
  ops->nvwl2norm         = N_VWL2Norm_SpcParallel;
  ops->nvl1norm          = N_VL1Norm_SpcParallel;
  ops->nvcompare         = N_VCompare_SpcParallel;
  ops->nvinvtest         = N_VInvTest_SpcParallel;
  ops->nvconstrmask      = N_VConstrMask_SpcParallel;
  ops->nvminquotient     = N_VMinQuotient_SpcParallel;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_SpcParallel) malloc(sizeof(struct _N_VectorContent_SpcParallel));
  if (content == NULL) {free(ops); free(v); return(NULL);}

  /* Allocate space for arrays in content */
  content->Nspc  = (int *) malloc(Ngrp*sizeof(int));
  content->n1    = (long int *) malloc(Ngrp*sizeof(long int));
  content->gdata = (realtype **) malloc(Ngrp*sizeof(realtype *));

  /* Attach lengths and communicator to N_Vector_Content_SpcParallel structure */
  content->Ngrp = Ngrp;
  tmp = (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
  n = 0;
  for(ig=0; ig<Ngrp; ig++) {
    content->Nspc[ig] = Nspc[ig];
    content->n1[ig] = Nspc[ig] * tmp;
    n += Nspc[ig] * tmp;
  }
  MPI_Allreduce(&n, &N, 1, SPVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);

  content->n    = n;
  content->N    = N;
  content->Nx   = Nx;
  content->Ny   = Ny;
  content->Nz   = Nz;
  content->NGx  = NGx;
  content->NGy  = NGy;
  content->NGz  = NGz;
  content->comm = comm;

  /* Attach empty data array to content structure */
  content->data     = NULL;
  content->own_data = FALSE;

  /* Attach content and ops to generic N_Vector */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* -----------------------------------------------------------------
 * Create a new spc_parallel vector
 */

N_Vector N_VNew_SpcParallel(MPI_Comm comm, int Ngrp, int *Nspc, 
                            long int Nx,  long int Ny,  long int Nz, 
                            long int NGx, long int NGy, long int NGz)
{
  N_Vector v;
  realtype *data;
  long int n;

  /* Create the new vector */
  v = NULL;
  v = N_VNewEmpty_SpcParallel(comm, Ngrp, Nspc, Nx, Ny, Nz, NGx, NGy, NGz);
  if (v == NULL) return(NULL);

  /* Get local length */
  n = SPV_LOCLENGTH(v);

  /* Create data */
  if ( n > 0 ) {
    /* Allocate memory */
    data = NULL;
    data = (realtype *) calloc(n, sizeof(realtype));
    if(data == NULL) {
      N_VDestroy_SpcParallel(v); 
      return(NULL); 
    }
    /* Attach data */
    VAttach_Data(v,data);
    SPV_OWN_DATA(v) = TRUE;
  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parallel N_Vector and attach to it the
 * user data, which must be a contiguous array and include the ghost
 * boundaries. In other words, data must have length
 *    sum(Nspc(igrp)) * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
 */

N_Vector N_VMake_SpcParallel(MPI_Comm comm, int Ngrp, int *Nspc, 
                             long int Nx,  long int Ny,  long int Nz, 
                             long int NGx, long int NGy, long int NGz,
                             realtype *data)
{
  N_Vector v;

  /* Create the new N_Vector */
  v = NULL;
  v = N_VNewEmpty_SpcParallel(comm, Ngrp, Nspc, Nx, Ny, Nz, NGx, NGy, NGz);
  if (v == NULL) return(NULL);

  /* Attach data if it has nonzero size*/
  if ( SPV_LOCLENGTH(v) > 0 ) {
    VAttach_Data(v,data);
    SPV_OWN_DATA(v) = FALSE;
  }
  
  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to copy into a given SPCPARALLEL N_Vector the data for
 * one species. The user must indicate the group index and the
 * species index within that group. The data must be a realtype 
 * array of length Nx * Ny * Nz
 */

void N_VLoad_SpcParallel(N_Vector v, int igrp, int ispc, realtype *data)
{
  realtype *vd;
  int Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock;

  vd = SPV_GDATA(v,igrp);
  Ns = SPV_NSPECIES(v,igrp);
  Nx   = SPV_XLENGTH(v);
  Ny   = SPV_YLENGTH(v);
  Nz   = SPV_ZLENGTH(v);
  NGx  = SPV_XGHOST(v);
  NGy  = SPV_YGHOST(v);
  NGz  = SPV_ZGHOST(v);

  for (k=NGz; k<Nz+NGz; k++) {
    Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
    for (j=NGy; j<Ny+NGy; j++) { 
      Yblock = j * (Nx + 2*NGx) * Ns;
      for (i=NGx; i<Nx+NGx; i++) {
        Xblock = i * Ns;
        vd[Zblock+Yblock+Xblock+ispc] = *(data++);
      }
    }
  }
  
}

/* ---------------------------------------------------------------- 
 * Function to create an array of SPCPARALLEL vectors with empty
 * (NULL) data array by cloning from a given vector w.
 */

N_Vector *N_VCloneVectorArrayEmpty_SpcParallel(int count, N_Vector w)
{
  N_Vector *vs;
  int j;


  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_SpcParallel(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_SpcParallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of SPCPARALLEL vectors by cloning
 * from a given vector w
 */

N_Vector *N_VCloneVectorArray_SpcParallel(int count, N_Vector w) 
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_SpcParallel(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_SpcParallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array of SPCPARALLEL vectors created with 
 * N_VCloneVectorArray_SpcParallel
 */

void N_VDestroyVectorArray_SpcParallel(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_SpcParallel(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ---------------------------------------------------------------- 
 * Function to print a parallel vector to stdout or to a file
 */

void N_VPrint_SpcParallel(N_Vector v, char *fname)
{
  FILE *fp;
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *vd;

  vd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(v);
  Nx   = SPV_XLENGTH(v);
  Ny   = SPV_YLENGTH(v);
  Nz   = SPV_ZLENGTH(v);
  NGx  = SPV_XGHOST(v);
  NGy  = SPV_YGHOST(v);
  NGz  = SPV_ZGHOST(v);

  /* open output file */
  if (fname[0] == '\0')
    fp = stdout;
  else
    fp = fopen(fname, "w");

  for (ig=0; ig<Ngrp; ig++) {
    vd = SPV_GDATA(v,ig);
    Ns = SPV_NSPECIES(v,ig);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            fprintf(fp, "(%d,%ld,%ld,%ld,%d) = %e\n", 
                    ig, i-NGx+1, j-NGy+1, k-NGz+1, is, vd[loc]);
          }
        }
      }
    }
  }

  if (fname[0] != '\0') fclose(fp);
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

/* 
 * N_VCloneEmpty_SpcParallel returns a new N_Vector of the same form 
 * as the input N_Vector, but with empty data container 
 */
 
N_Vector N_VCloneEmpty_SpcParallel(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_SpcParallel content;
  int ig, Ngrp;

  /* Check that w has been created */
  if (w == NULL) return(NULL);

  /* Create the new vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *w);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  /* Attach operations */
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
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_SpcParallel) malloc(sizeof(struct _N_VectorContent_SpcParallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  Ngrp = SPV_NGROUPS(w);

  /* Allocate space for arrays in content */  
  content->Nspc  = (int *) malloc(Ngrp*sizeof(int));
  content->n1    = (long int *) malloc(Ngrp*sizeof(long int));
  content->gdata = (realtype **) malloc(Ngrp*sizeof(realtype *));

  /* Attach lengths and communicator to content structure */
  content->Ngrp = SPV_NGROUPS(w);
  content->Nx   = SPV_XLENGTH(w);
  content->Ny   = SPV_YLENGTH(w);
  content->Nz   = SPV_ZLENGTH(w);
  content->NGx  = SPV_XGHOST(w);
  content->NGy  = SPV_YGHOST(w);
  content->NGz  = SPV_ZGHOST(w);
  for(ig=0; ig<Ngrp; ig++) {
    content->Nspc[ig] = SPV_NSPECIES(w,ig);
    content->n1[ig]   = SPV_LOCLENGTH1(w,ig);
  }
  content->n    = SPV_LOCLENGTH(w);
  content->N    = SPV_GLOBLENGTH(w);
  content->comm = SPV_COMM(w);

  /* Attach empty data array to content structure */
  content->data     = NULL;
  content->own_data = FALSE;

  /* Attach content and ops to generic N_Vector */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/*
 * N_VClone_SpcParallel returns a new N_Vector of the same 
 * form as the input N_Vector.
 */ 

N_Vector N_VClone_SpcParallel(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int n;

  /* Create vector */
  v = NULL;
  v = N_VCloneEmpty_SpcParallel(w);
  if (v == NULL) return(NULL);

  /* Get local length */
  n  = SPV_LOCLENGTH(w);

  /* Create data */
  if ( n > 0 ) {
    /* Allocate memory */
    data = NULL;
    data = (realtype *) calloc(n, sizeof(realtype));
    if(data == NULL) { 
      N_VDestroy_SpcParallel(v); 
      return(NULL); 
    }
    /* Attach data */
    VAttach_Data(v,data);
    SPV_OWN_DATA(v) = TRUE;
  }

  return(v);
}

/* 
 * N_VDestroy_SpcParallel frees the data storage for an N_Vector 
*/ 

void N_VDestroy_SpcParallel(N_Vector v)
{
  N_VectorContent_SpcParallel content;

  if ( SPV_OWN_DATA(v) == TRUE ) 
    if (SPV_DATA(v) != NULL) {
      free(SPV_DATA(v));
      SPV_DATA(v) = NULL;
    }

  content = SPV_CONTENT(v);

  free(content->Nspc); content->Nspc = NULL;
  free(content->n1); content->n1 = NULL;
  free(content->gdata); content->gdata = NULL;
  free(content); content = NULL;

  free(v->ops); v->ops = NULL;
 
  free(v); v = NULL;
}

/* 
 * N_VSpace_SpcParallel returns the local space requirements for one N_Vector. 
 * The amount of realtype data is given in lrw, and long int data in liw.  
 * Note: this includes ghost cell data storage as well 
 */

void N_VSpace_SpcParallel(N_Vector v, long int *lrw, long int *liw)
{
  *lrw = SPV_LOCLENGTH(v);
  *liw = 11 + 3*SPV_NGROUPS(v);
}

/* 
 * N_VGetArrayPointer_SpcParallel extracts the data component array 
 * from the N_Vector v
 */

realtype *N_VGetArrayPointer_SpcParallel(N_Vector v)
{
  return((realtype *) SPV_DATA(v));
}

/* 
 * N_VSetArrayPointer_SpcParallel attaches the data component array 
 * v_data to the N_Vector v
 */

void N_VSetArrayPointer_SpcParallel(realtype *v_data, N_Vector v)
{
  SPV_DATA(v) = v_data;
}

/* 
 * N_VLinearSum_SpcParallel calculates z = a*x + b*y 
*/

void N_VLinearSum_SpcParallel(realtype a, N_Vector x, realtype b, 
                              N_Vector y, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (z == y)) {    
    Vaxpy_SpcParallel(a,x,y);
    return;
  }
  
  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (z == x)) {    
    Vaxpy_SpcParallel(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE)) {
    VSum_SpcParallel(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */
  if (( test = ((a == ONE) && (b == -ONE))) 
            || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_SpcParallel(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, */
  /*        (2) a == other or 0.0, b == 1.0  */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_SpcParallel(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */
  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_SpcParallel(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have 
     called N_VConst */
  if (a == b) {
    VScaleSum_SpcParallel(a, x, y, z);
    return;
  }

  /* Case: a == -b */
  if (a == -b) {
    VScaleDiff_SpcParallel(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */  

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = a * xd[loc] + b * yd[loc];
          }
        }
      }
    }
  }

}

/* 
 * N_VConst_SpcParallel calculates z[i] = c for all i 
 */

void N_VConst_SpcParallel(realtype c, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *zd;

  zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(z);
  Nx   = SPV_XLENGTH(z);
  Ny   = SPV_YLENGTH(z);
  Nz   = SPV_ZLENGTH(z);
  NGx  = SPV_XGHOST(z);
  NGy  = SPV_YGHOST(z);
  NGz  = SPV_ZGHOST(z);

  for(ig=0; ig<Ngrp; ig++) {
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(z,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = c;
          }
        }
      }
    }
  }

}

/* 
 * N_VProd_SpcParallel calculates z[i] = x[i]*y[i] 
 */

void N_VProd_SpcParallel(N_Vector x, N_Vector y, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = xd[loc] * yd[loc];
          }
        }
      }
    }
  }

}

/* 
 * N_VDiv_SpcParallel calculates z[i] = x[i]/y[i] 
 */

void N_VDiv_SpcParallel(N_Vector x, N_Vector y, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);  

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = xd[loc] / yd[loc];
          }
        }
      }
    }
  }

}

/* 
 * N_VScale_SpcParallel calculates z = c*x 
 */

void N_VScale_SpcParallel(realtype c, N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* BLAS usage: scale x <- cx */
  if (z == x) {
    VScaleBy_SpcParallel(c, x);
    return;
  }

  if (c == ONE) {

    VCopy_SpcParallel(x, z);

  } else if (c == -ONE) {

    VNeg_SpcParallel(x, z);

  } else {

    /* Get mesh info & data from vector */
    Ngrp = SPV_NGROUPS(x);
    Nx   = SPV_XLENGTH(x);
    Ny   = SPV_YLENGTH(x);
    Nz   = SPV_ZLENGTH(x);
    NGx  = SPV_XGHOST(x);
    NGy  = SPV_YGHOST(x);
    NGz  = SPV_ZGHOST(x);

    for(ig=0; ig<Ngrp; ig++) {
      xd = SPV_GDATA(x,ig);
      zd = SPV_GDATA(z,ig);
      Ns  = SPV_NSPECIES(x,ig);   
      for (k=NGz; k<Nz+NGz; k++) {
        Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
        for (j=NGy; j<Ny+NGy; j++) { 
          Yblock = j * (Nx + 2*NGx) * Ns;
          for (i=NGx; i<Nx+NGx; i++) {
            Xblock = i * Ns;
            for (is=0; is<Ns; is++) {
              loc = Zblock + Yblock + Xblock + is;
              zd[loc] = c * xd[loc];
            }
          }
        }
      }
    }

  }

}

/* 
 * N_VAbs_SpcParallel or (nvabs) calculates z[i] = |x[i]| 
 */

void N_VAbs_SpcParallel(N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);  

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = ABS(xd[loc]);
          }
        }
      }
    }
  }

}

/* 
 * N_VInv_SpcParallel calculates z[i] = 1/x[i].  
 *  Note: it does not check for division by 0.  It should be called only 
 *  with an N_Vector x which is guaranteed to have all non-zero components
 */

void N_VInv_SpcParallel(N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = ONE / xd[loc];
          }
        }
      }
    }
  }

}

/* 
 * N_VAddConst_SpcParallel calculates z[i] = x[i] + b 
 */

void N_VAddConst_SpcParallel(N_Vector x, realtype b, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;
  
  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = xd[loc] + b;
          }
        }
      }
    }
  }

}

/* 
 * N_VDotProd_SpcParallel returns the value of the ordinary dot product 
 *  of x and y, i.e. sum (i=0 to N-1) {x[i] * y[i]} 
 */

realtype N_VDotProd_SpcParallel(N_Vector x, N_Vector y)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype sum, gsum, *xd, *yd;
  MPI_Comm comm;

  xd = yd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  sum = ZERO;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            sum += xd[loc] * yd[loc];
          }
        }
      }
    }
  }

  /* obtain global sum from local sums */
  gsum = ZERO;
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  return(gsum);
}

/* 
 * N_VMaxNorm_SpcParallel returns the maximum norm of x, 
 *  i.e. max(i=1 to N-1) |x[i]|  
 */

realtype N_VMaxNorm_SpcParallel(N_Vector x)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype max, *xd, gmax;
  MPI_Comm comm;

  xd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  max = ZERO;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            if (ABS(xd[loc]) > max) max = ABS(xd[loc]);
          }
        }
      }
    }
  }

  /* Obtain global max from local maxima */
  gmax = ZERO;  
  MPI_Allreduce(&max, &gmax, 1, SPVEC_REAL_MPI_TYPE, MPI_MAX, comm);

  return(gmax);
}

/*
 * N_VWrmsNorm_SpcParallel returns the weighted root mean square norm 
 * of x with weight factor w, i.e. 
 * sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N] 
 */

realtype N_VWrmsNorm_SpcParallel(N_Vector x, N_Vector w)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz, Nglobal;
  long int Yblock, Zblock, Xblock, loc;
  realtype sum, gsum, prodi, *xd, *wd;
  MPI_Comm comm;

  xd = wd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  Nglobal = SPV_GLOBLENGTH(x);
  comm = SPV_COMM(x);

  sum = ZERO;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    wd = SPV_GDATA(w,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            prodi = xd[loc] * wd[loc];
            sum += prodi * prodi;
          }
        }
      }
    }
  }

  /* Obtain global sum from local sums */
  gsum = ZERO;
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  /* scale correctly and return */
  return(RSqrt(gsum / Nglobal));
}

/* 
 * N_VWrmsNormMask_SpcParallel or (nvwrmsnormmask) returns ??? 
*/

realtype N_VWrmsNormMask_SpcParallel(N_Vector x, N_Vector w, N_Vector id)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz, Nglobal;
  long int Yblock, Zblock, Xblock, loc;
  realtype sum, gsum, prodi, *xd, *wd, *idd;
  MPI_Comm comm;

  xd = wd = idd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x); 
  Nglobal = SPV_GLOBLENGTH(x);
  comm = SPV_COMM(x);

  sum = ZERO;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    wd = SPV_GDATA(w,ig);
    idd = SPV_GDATA(id,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            if (idd[loc] > ZERO) {
              prodi = xd[loc] * wd[loc];
              sum += prodi * prodi;
            }
          }
        }
      }
    }
  }

  /* Obtain global sum from local sums */
  gsum = ZERO;
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  /* scale result and return */
  return(RSqrt(gsum / Nglobal));
}

/* 
 * N_VMin_SpcParallel returns the smallest element of x 
 */

realtype N_VMin_SpcParallel(N_Vector x)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype min, gmin, *xd; 
  MPI_Comm comm;

  xd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);  
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  min = BIG_REAL;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            if (xd[loc] < min) min = xd[loc];
          }
        }
      }
    }
  }

  /* Obtain global min from local mins */
  gmin = BIG_REAL;
  MPI_Allreduce(&min, &gmin, 1, SPVEC_REAL_MPI_TYPE, MPI_MIN, comm);

  return(gmin);
}

/* 
 * N_VWL2Norm_SpcParallel returns the weighted Euclidean L2 norm of x 
 * with weight factor w, 
 * i.e. sqrt [(sum (i=0 to N-1) {(x[i]*w[i])^2}) ] 
 */

realtype N_VWL2Norm_SpcParallel(N_Vector x, N_Vector w)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  xd = wd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);  
  comm = SPV_COMM(x);

  sum = ZERO;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    wd = SPV_GDATA(w,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            prodi = xd[loc] * wd[loc];
            sum += prodi * prodi;
          }
        }
      }
    }
  }

  /* Obtain global sum from local sums */
  gsum = ZERO;
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  return(RSqrt(gsum));
}

/* 
 * N_VL1Norm_SpcParallel returns the L1 norm of x, 
 * i.e. sum (i=0 to N-1) {|x[i]|} 
 */

realtype N_VL1Norm_SpcParallel(N_Vector x)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype sum, gsum, *xd;
  MPI_Comm comm;

  xd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  sum = ZERO;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            sum += ABS(xd[loc]);
          }
        }
      }
    }
  }

  /* Obtain global sum from local sums */
  gsum = ZERO;
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  return(gsum);
}

/* 
 * N_VCompare_SpcParallel calculates z[i] = 1 if |x[i]| > c, z[i] = 0 otherwise 
 */

void N_VCompare_SpcParallel(realtype c, N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;
  
  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = (ABS(xd[loc]) >= c) ? ONE : ZERO;
          }
        }
      }
    }
  }

}

/*
 * N_VInvTest_SpcParallel computes z[i] = 1/x[i] with a test for x[i] == 0 
 * before inverting x[i].  This routine returns TRUE if all components 
 * of x are nonzero (successful inversion) and returns FALSE otherwise. 
 */

booleantype N_VInvTest_SpcParallel(N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  /* Initialize return value */
  val = ONE;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            if (xd[loc] == ZERO) val = ZERO;
            else                 zd[loc] = ONE / xd[loc];
          }
        }
      }
    }
  }

  /* Obtain global return value from local values */
  MPI_Allreduce(&val, &gval, 1, SPVEC_REAL_MPI_TYPE, MPI_MIN, comm);

  if (gval == ZERO)  return(FALSE);
  else               return(TRUE);
}

/* 
 * N_VConstrMask_SpcParallel returns a boolean FALSE if any element fails 
 * the constraint test, and TRUE if all passed.  The constraint test 
 * is as follows: 
 *       if c[i] =  2.0, then x[i] must be >  0.0
 *       if c[i] =  1.0, then x[i] must be >= 0.0
 *       if c[i] = -1.0, then x[i] must be <= 0.0
 *       if c[i] = -2.0, then x[i] must be <  0.0
 * It also sets a mask vector m, with elements equal to 1.0 where the 
 * corresponding constraint test failed, and equal to 0.0 where the 
 * constraint test passed.  This routine is specialized in that it is 
 * used only for constraint checking. 
 */

booleantype N_VConstrMask_SpcParallel(N_Vector c, N_Vector x, N_Vector m)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  booleantype test, alltest;
  realtype *cd, *xd, *md;
  MPI_Comm comm;
 
  cd = xd = md = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  /* Initialize output variable */
  test = TRUE;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    cd = SPV_GDATA(c,ig);
    md = SPV_GDATA(m,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            md[loc] = ZERO;
            if (cd[loc] == ZERO) continue;
            if (cd[loc] > ONEPT5 || cd[loc] < -ONEPT5) {
              if (xd[loc]*cd[loc] <= ZERO) {
                test = FALSE; 
                md[loc] = ONE; 
              }
              continue;
            }
            if (cd[loc] > HALF || cd[loc] < -HALF) {
              if (xd[loc]*cd[loc] < ZERO) {
                test = FALSE;
                md[loc] = ONE;
              }
            }
          }
        }
      }
    }
  }

  /* Obtain global return value from local return values */
  MPI_Allreduce(&test, &alltest, 1, MPI_INT, MPI_MIN, comm);

  return(alltest);
}

/* 
 * N_VMinQuotient_SpcParallel returns min(x[i]/y[i]) over all i 
 * such that denom[i] != 0. 
 */

realtype N_VMinQuotient_SpcParallel(N_Vector x, N_Vector y)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  booleantype notEvenOnce;
  realtype *xd, *yd, min, gmin;
  MPI_Comm comm;

  xd = yd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);
  comm = SPV_COMM(x);

  /* Initialize output value, minimum */
  notEvenOnce = TRUE;
  min = BIG_REAL;
  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            if (yd[loc] == ZERO) continue;
            else {
              if (!notEvenOnce) min = MIN(min, xd[loc] / yd[loc]);
              else {
                min = xd[loc] / yd[loc] ;
                notEvenOnce = FALSE;
              }
            }
          }
        }
      }
    }
  }

  /* Obtain global min from local minima */
  gmin = BIG_REAL;
  MPI_Allreduce(&min, &gmin, 1, SPVEC_REAL_MPI_TYPE, MPI_MIN, comm);

  return(gmin);
}
 
/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

/*
 * Attach data to a vector and set pointers to start points of 
 * group data
 */

static void VAttach_Data(N_Vector v, realtype *data)
{
  int ig, Ngrp;

  SPV_DATA(v) = data;
  Ngrp = SPV_NGROUPS(v);

  SPV_GDATA(v,0) = SPV_DATA(v);
  for (ig=1; ig<Ngrp; ig++)
    SPV_GDATA(v,ig) = SPV_GDATA(v,ig-1) + SPV_LOCLENGTH1(v,ig-1);

}


/*
 * VCopy_SpcParallel is a private helper function that copies the 
 * local parts of the vector x into the local parts of the vector z. 
 * Note: the mesh dimensions must match. 
 */

static void VCopy_SpcParallel(N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = xd[loc];
          }
        }
      }
    }
  }

}

/*
 * VSum_SpcParallel is a private helper function that calculates the 
 * sum z[i] = x[i] + y[i] over the local processor domain. 
 */

static void VSum_SpcParallel(N_Vector x, N_Vector y, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = xd[loc] + yd[loc];
          }
        }
      }
    }
  }

}

/*
 * VDiff_SpcParallel is a private helper function that calculates the 
 * difference z[i] = x[i] - y[i] over the local processor domain. 
 */

static void VDiff_SpcParallel(N_Vector x, N_Vector y, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;
 
  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = xd[loc] - yd[loc];
          }
        }
      }
    }
  }

}

/*
 * VNeg_SpcParallel is a private helper function that calculates the 
 * negation z[i] = -x[i] over the local processor domain. 
 */

static void VNeg_SpcParallel(N_Vector x, N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = -xd[loc];
          }
        }
      }
    }
  }

}

/*
 * VScaleSum_SpcParallel is a private helper function that calculates 
 * the scaled sum z[i] = c*(x[i] + y[i]) over the local processor domain 
 */

static void VScaleSum_SpcParallel(realtype c, N_Vector x, N_Vector y, 
                                  N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = c * (xd[loc] + yd[loc]);
          }
        }
      }
    }
  }

}

/*
 * VScaleDiff_SpcParallel is a private helper function that calculates
 * the scaled difference z[i] = c*(x[i]-y[i]) over the local processor 
 * domain 
 */

static void VScaleDiff_SpcParallel(realtype c, N_Vector x, N_Vector y, 
				    N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = c * (xd[loc] - yd[loc]);
          }
        }
      }
    }
  }

}

/*
 * VLin1_SpcParallel is a private helper function that calculates the 
 * linear operation z[i] = a*x[i] + y[i] over the local processor domain. 
 */

static void VLin1_SpcParallel(realtype a, N_Vector x, N_Vector y, 
			       N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = a*xd[loc] + yd[loc];
          }
        }
      }
    }
  }

}

/*
 * VLin2_SpcParallel is a private helper function that calculates the
 * linear operation z[i] = a*x[i] - y[i] over the local processor domain. 
 */

static void VLin2_SpcParallel(realtype a, N_Vector x, N_Vector y, 
			       N_Vector z)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    zd = SPV_GDATA(z,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            zd[loc] = a*xd[loc] - yd[loc];
          }
        }
      }
    }
  }

}

/*
 * Vaxpy_SpcParallel is a private helper function that calculates the
 * axpy operation y[i] = a*x[i] + y[i] over the local processor domain. 
 */

static void Vaxpy_SpcParallel(realtype a, N_Vector x, N_Vector y)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd, *yd;

  xd = yd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  if (a == ONE) {

    for(ig=0; ig<Ngrp; ig++) {
      xd = SPV_GDATA(x,ig);
      yd = SPV_GDATA(y,ig);
      Ns  = SPV_NSPECIES(x,ig);   
      for (k=NGz; k<Nz+NGz; k++) {
        Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
        for (j=NGy; j<Ny+NGy; j++) { 
          Yblock = j * (Nx + 2*NGx) * Ns;
          for (i=NGx; i<Nx+NGx; i++) {
            Xblock = i * Ns;
            for (is=0; is<Ns; is++) {
              loc = Zblock + Yblock + Xblock + is;
              yd[loc] += xd[loc];
            }
          }
        }
      }
    }

    return;

  }

  if (a == -ONE) {

    for(ig=0; ig<Ngrp; ig++) {
      xd = SPV_GDATA(x,ig);
      yd = SPV_GDATA(y,ig);
      Ns  = SPV_NSPECIES(x,ig);   
      for (k=NGz; k<Nz+NGz; k++) {
        Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
        for (j=NGy; j<Ny+NGy; j++) { 
          Yblock = j * (Nx + 2*NGx) * Ns;
          for (i=NGx; i<Nx+NGx; i++) {
            Xblock = i * Ns;
            for (is=0; is<Ns; is++) {
              loc = Zblock + Yblock + Xblock + is;
              yd[loc] -= xd[loc];
            }
          }
        }
      }
    }

    return;

  }

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    yd = SPV_GDATA(y,ig);
    Ns  = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            yd[loc] += a * xd[loc];
          }
        }
      }
    }
  }

}

/*
 * VScaleBy_SpcParallel is a private helper function that calculates the 
 * scaled product x[i] = a*x[i] over the local processor domain. 
 */

static void VScaleBy_SpcParallel(realtype a, N_Vector x)
{
  int ig, Ngrp, is, Ns;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Xblock, loc;
  realtype *xd;

  xd = NULL;

  /* Get mesh info & data from vector */
  Ngrp = SPV_NGROUPS(x);
  Nx   = SPV_XLENGTH(x);
  Ny   = SPV_YLENGTH(x);
  Nz   = SPV_ZLENGTH(x);
  NGx  = SPV_XGHOST(x);
  NGy  = SPV_YGHOST(x);
  NGz  = SPV_ZGHOST(x);

  for(ig=0; ig<Ngrp; ig++) {
    xd = SPV_GDATA(x,ig);
    Ns = SPV_NSPECIES(x,ig);   
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy) * Ns;
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx) * Ns;
	for (i=NGx; i<Nx+NGx; i++) {
          Xblock = i * Ns;
          for (is=0; is<Ns; is++) {
            loc = Zblock + Yblock + Xblock + is;
            xd[loc] *= a;
          }
        }
      }
    }
  }

}
