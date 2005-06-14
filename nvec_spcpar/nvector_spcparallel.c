/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-06-14 19:00:57 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds and Radu Serban @LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for a spc_parallel MPI 
 * implementation of the NVECTOR package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "nvector_spcparallel.h"
#include "sundialstypes.h"
#include "sundialsmath.h" 

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private function prototypes */

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

N_Vector N_VNewEmpty_SpcParallel(MPI_Comm comm, int nspc,
                                 long int Nx,  long int Ny,  long int Nz, 
                                 long int NGx, long int NGy, long int NGz) 
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_SpcParallel content;

  long int n1, n1g, n, ng, Ng;

  v = NULL;
  ops = NULL;
  content = NULL;
  
  /* Compute various local lengths */

  n1  = Nx * Ny * Nz;
  n1g = (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
  n   = nspc * n1;
  ng  = nspc * n1g;

  /* Compute global length */

  MPI_Allreduce(&ng, &Ng, 1, SPVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);

  /* Create vector */

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */

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

  content = (N_VectorContent_SpcParallel) malloc(sizeof(struct _N_VectorContent_SpcParallel));
  if (content == NULL) {free(ops); free(v); return(NULL);}

  /* Attach lengths and communicator to N_Vector_Content_SpcParallel structure */

  content->nspc   = nspc;

  content->Nx   = Nx;
  content->Ny   = Ny;
  content->Nz   = Nz;
  content->NGx  = NGx;
  content->NGy  = NGy;
  content->NGz  = NGz;

  content->n1   = n1;
  content->n1g  = n1g;
  content->n    = n;
  content->ng   = ng;

  content->Ng   = Ng;

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

N_Vector N_VNew_SpcParallel(MPI_Comm comm, int nspc, 
                            long int Nx,  long int Ny,  long int Nz, 
                            long int NGx, long int NGy, long int NGz)
{
  N_Vector v;
  realtype *data;
  long int ng;

  v = NULL;
  data = NULL;

  /* Create the new N_Vector */

  v = N_VNewEmpty_SpcParallel(comm, nspc, Nx, Ny, Nz, NGx, NGy, NGz);
  if (v == NULL) return(NULL);

  /* Compute local length */

  ng = nspc * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);

  /* Create data */

  if ( ng > 0 ) {

    /* Allocate memory */
    data = (realtype *) malloc(ng * sizeof(realtype *));
    if(data == NULL) { 
      N_VDestroy_SpcParallel(v); 
      return(NULL); 
    }

    /* Attach data */
    SPV_OWN_DATA(v) = TRUE;
    SPV_DATA(v)     = data;
  
  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parallel N_Vector and attach to it the
 * user data, which must be a contiguous array and include the ghost
 * boundaries. In other words, data must have length
 *    nspc * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
 */

N_Vector N_VAttach_SpcParallel(MPI_Comm comm, int nspc, 
                               long int Nx,  long int Ny,  long int Nz, 
                               long int NGx, long int NGy, long int NGz,
                               realtype *data)
{
  N_Vector v;

  v = NULL;

  /* Create the new N_Vector */
  v = N_VNewEmpty_SpcParallel(comm, nspc, Nx, Ny, Nz, NGx, NGy, NGz);
  if (v == NULL) return(NULL);

  /* Attach data if it has nonzero size*/
  if ( SPV_LOCLENGTH(v) > 0 ) {
    SPV_OWN_DATA(v) = FALSE;
    SPV_DATA(v)     = data;
  }
  
  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parallel N_Vector and copy into it some
 * user data (which must be an array of nspc vectors, each of length
 * Nx*Ny*Nz)
 */

N_Vector N_VLoad_SpcParallel(MPI_Comm comm, int nspc, 
                             long int Nx,  long int Ny,  long int Nz, 
                             long int NGx, long int NGy, long int NGz,
                             realtype **data)
{
  N_Vector v;
  realtype *vd, *sdata;
  int s;
  long int i, j, k, Yblock, Zblock, loc;

  v = NULL;

  /* Create the new N_Vector */
  v = N_VNew_SpcParallel(comm, nspc, Nx, Ny, Nz, NGx, NGy, NGz);
  if (v == NULL) return(NULL);

  /* If data has nonzero size, copy user data */

  if ( SPV_DATA(v) != NULL ) {

    for(s=0; s<nspc; s++) {

      vd = SPV_SDATA(v,s);
      sdata = data[s];

      for (k=NGz; k<Nz+NGz; k++) {
        Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
        for (j=NGy; j<Ny+NGy; j++) { 
          Yblock = j * (Nx + 2*NGx);
          for (i=NGx; i<Nx+NGx; i++) {
            loc = Zblock + Yblock + i;
            vd[loc] = *(sdata++);
          }
        }
      }

    }
      
  }
  
  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new spcparallel vectors. 
 */

N_Vector *N_VNewVectorArray_SpcParallel(int count, 
                                        MPI_Comm comm, int nspc,
                                        long int Nx,  long int Ny,  long int Nz, 
                                        long int NGx, long int NGy, long int NGz)
{
  N_Vector *vs;
  int j;

  vs = NULL;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = N_VNew_SpcParallel(comm, nspc, Nx, Ny, Nz, NGx, NGy, NGz);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new spcparallel vectors with empty
 * (NULL) data array.
 */

N_Vector *N_VNewVectorArrayEmpty_SpcParallel(int count, 
                                             MPI_Comm comm, int nspc, 
                                             long int Nx,  long int Ny,  long int Nz, 
                                             long int NGx, long int NGy, long int NGz) 
{
  N_Vector *vs;
  int j;

  vs = NULL;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = N_VNewEmpty_SpcParallel(comm, nspc, Nx, Ny, Nz, NGx, NGy, NGz);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VNewVectorArray_SpcParallel
 */

void N_VDestroyVectorArray_SpcParallel(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_SpcParallel(vs[j]);

  free(vs);

  return;
}

/* ---------------------------------------------------------------- 
 * Function to print a parallel vector to stdout
 */

void N_VPrint_SpcParallel(N_Vector v)
{
  int nspc, s;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, loc;
  realtype *vd;

  vd = NULL;

  /* Get mesh info & data from vector */
  nspc  = SPV_NSPECIES(v);
  Nx  = SPV_XLENGTH(v);
  Ny  = SPV_YLENGTH(v);
  Nz  = SPV_ZLENGTH(v);
  NGx = SPV_XGHOST(v);
  NGy = SPV_YGHOST(v);
  NGz = SPV_ZGHOST(v);

  for (s=0; s<nspc; s++) {
    vd = SPV_SDATA(v,s);
    printf("\nSpecies %d\n\n",s+1);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Zblock + Yblock + i;
	  printf("(%ld,%ld,%ld,%d) = %e\n", 
                 i-NGx+1, j-NGy+1, k-NGz+1, s+1, vd[loc]);
	}
      }
    }
  }

  printf("\n");
}

/* ---------------------------------------------------------------- 
 * Function to print a parallel vector to a file
 */

void N_VPrintFile_SpcParallel(char *fname, N_Vector v)
{
  FILE *fp;
  int nspc, s;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, loc;
  realtype *vd;

  vd = NULL;

  /* Get mesh info & data from vector */
  nspc  = SPV_NSPECIES(v);
  Nx  = SPV_XLENGTH(v);
  Ny  = SPV_YLENGTH(v);
  Nz  = SPV_ZLENGTH(v);
  NGx = SPV_XGHOST(v);
  NGy = SPV_YGHOST(v);
  NGz = SPV_ZGHOST(v);

  /* open output file */
  fp = fopen(fname, "w");

  for (s=0; s<nspc; s++) {
    vd = SPV_SDATA(v,s);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Zblock + Yblock + i;
	  fprintf(fp, "(%ld,%ld,%ld,%d) = %e\n", 
		  i-NGx+1, j-NGy+1, k-NGz+1, s+1, vd[loc]);
	}
      }
    }
  }

  fclose(fp);
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

  int nspc;
  long int Nx, Ny, Nz, NGx, NGy, NGz;
  long int n1, n1g, n, ng;

  v = NULL;
  ops = NULL;
  content = NULL;

  /* Check that w has been created */
  if (w == NULL) return(NULL);

  /* Create vector */
  v = (N_Vector) malloc(sizeof *w);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
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

  content = (N_VectorContent_SpcParallel) malloc(sizeof(struct _N_VectorContent_SpcParallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Obtain various dimensions from w */

  nspc = SPV_NSPECIES(w);

  Nx = SPV_XLENGTH(w);
  Ny = SPV_YLENGTH(w);
  Nz = SPV_ZLENGTH(w);

  NGx = SPV_XGHOST(w);
  NGy = SPV_YGHOST(w);
  NGz = SPV_ZGHOST(w);

  n1  = Nx * Ny * Nz;
  n1g = (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
  n   = nspc * n1;
  ng  = nspc * n1g;

  /* Attach lengths and communicator to content structure */

  content->nspc   = nspc;

  content->Nx   = Nx;
  content->Ny   = Ny;
  content->Nz   = Nz;
  content->NGx  = NGx;
  content->NGy  = NGy;
  content->NGz  = NGz;

  content->n1   = n1;
  content->n1g  = n1g;
  content->n    = n;
  content->ng   = ng;

  content->Ng   = SPV_GLOBLENGTH(w);

  content->comm = SPV_COMM(w);

  /* set data arrays to null */
  content->own_data = FALSE;
  content->data     = NULL;

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
  long int ng;

  v = NULL;
  data = NULL;

  /* Create vector */

  v = N_VCloneEmpty_SpcParallel(w);
  if (v == NULL) return(NULL);

  /* Get local length */

  ng  = SPV_LOCLENGTH(w);

  /* Create data */

  if ( ng > 0 ) {

    /* Allocate memory */
    data = (realtype *) malloc(ng * sizeof(realtype *));
    if(data == NULL) { 
      N_VDestroy_SpcParallel(v); 
      return(NULL); 
    }

    /* Attach data */
    SPV_OWN_DATA(v) = TRUE;
    SPV_DATA(v)     = data;
  }

  return(v);
}

/* 
 * N_VDestroy_SpcParallel frees the data storage for an N_Vector 
*/ 

void N_VDestroy_SpcParallel(N_Vector v)
{
  if ( SPV_OWN_DATA(v) == TRUE ) 
    if (SPV_DATA(v) != NULL) free(SPV_DATA(v));

  free(v->content);
  free(v->ops);
  free(v);
}

/* 
 * N_VSpace_SpcParallel returns the space requirements for one N_Vector. 
 * The amount of realtype data is given in lrw, and long int data in liw.  
 * Note: this includes ghost cell data storage as well 
 */

void N_VSpace_SpcParallel(N_Vector v, long int *lrw, long int *liw)
{
  MPI_Comm comm;
  int npes;

  int NGlobal, NxLoc, NyLoc, NzLoc, NSpecies;
  int NProcs, NLoc;

  comm = SPV_COMM(v);
  MPI_Comm_size(comm, &npes);

  *lrw = SPV_GLOBLENGTH(v);
  *liw = 14*npes;
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
  if (SPV_LOCLENGTH(v) > 0) 
    SPV_DATA(v) = v_data;
}


/* 
 * N_VLinearSum_SpcParallel calculates z = a*x + b*y 
*/

void N_VLinearSum_SpcParallel(realtype a, N_Vector x, realtype b, 
		      N_Vector y, N_Vector z)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
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
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  /* Calculates:  for (i=0; i < N; i++),  *zd++ = a*(*xd++) + b*(*yd++); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = a * xd[loc] + b * yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *zd;

  zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(z);
  Ny  = SPV_YLENGTH(z);
  Nz  = SPV_ZLENGTH(z);
  nspc  = SPV_NSPECIES(z);
  NGx = SPV_XGHOST(z);
  NGy = SPV_YGHOST(z);
  NGz = SPV_ZGHOST(z);
  zd  = SPV_DATA(z);

  /*  Calculates:  for (i=0; i < N; i++),  *zd++ = c; */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = c;
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  /*  Calculates:  for (i=0; i < N; i++),  *zd++ = (*xd++) * (*yd++); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = xd[loc] * yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  /*  for (i=0; i < N; i++),  *zd++ = (*xd++) / (*yd++); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = xd[loc] / yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
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
    Nx  = SPV_XLENGTH(x);
    Ny  = SPV_YLENGTH(x);
    Nz  = SPV_ZLENGTH(x);
    nspc  = SPV_NSPECIES(x);
    NGx = SPV_XGHOST(x);
    NGy = SPV_YGHOST(x);
    NGz = SPV_ZGHOST(x);
    xd  = SPV_DATA(x);
    zd  = SPV_DATA(z);

    /*  for (i=0; i < N; i++),  *zd++ = c * (*xd++); */
    for (l=0; l<nspc; l++) {
      Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
      for (k=NGz; k<Nz+NGz; k++) {
	Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
	for (j=NGy; j<Ny+NGy; j++) { 
	  Yblock = j * (Nx + 2*NGx);
	  for (i=NGx; i<Nx+NGx; i++) {
	    loc = Sblock + Zblock + Yblock + i;
	    zd[loc] = c * xd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);

  /*  for (i=0; i < N; i++, xd++, zd++),  *zd = ABS(*xd); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = ABS(xd[loc]);
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);

  /*  for (i=0; i < N; i++),  *zd++ = ONE / (*xd++); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = ONE / xd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd;
  
  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);
  
  /*  for (i=0; i < N; i++) *zd++ = (*xd++) + b; */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = xd[loc] + b;
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype sum, *xd, *yd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = yd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  comm = SPV_COMM(x);

  /*  for (i=0; i < N; i++) sum += xd[i] * yd[i]; */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  sum += xd[loc] * yd[loc];
	}
      }
    }
  }

  /* obtain global sum from local sums */
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  return(gsum);
}

/* 
 * N_VMaxNorm_SpcParallel returns the maximum norm of x, 
 *  i.e. max(i=1 to N-1) |x[i]|  
 */

realtype N_VMaxNorm_SpcParallel(N_Vector x)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype max, *xd, gmax;
  MPI_Comm comm;

  max = gmax = ZERO;
  xd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  comm = SPV_COMM(x);

  /*  for (i=0; i < N; i++, xd++),  if (ABS(*xd) > max) max = ABS(*xd); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  if (ABS(xd[loc]) > max) max = ABS(xd[loc]);
	}
      }
    }
  }
   
  /* Obtain global max from local maxima */
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc, N_global;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  wd  = SPV_DATA(w);
  N_global = SPV_GLOBLENGTH(x);
  comm = SPV_COMM(x);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  prodi = xd[loc] * wd[loc];
	  sum += prodi * prodi;
	}
      }
    }
  }
  
  /* Obtain global sum from local sums */
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  /* scale correctly and return */
  return(RSqrt(gsum / N_global));
}

/* 
 * N_VWrmsNormMask_SpcParallel or (nvwrmsnormmask) returns ??? 
*/

realtype N_VWrmsNormMask_SpcParallel(N_Vector x, N_Vector w, N_Vector id)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc, N_global;
  realtype sum, prodi, *xd, *wd, *idd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = idd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  wd  = SPV_DATA(w);
  idd = SPV_DATA(id);
  N_global = SPV_GLOBLENGTH(x);
  comm = SPV_COMM(x);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  if (idd[loc] > ZERO) {
	    prodi = xd[loc] * wd[loc];
	    sum += prodi * prodi;
	  }
	}
      }
    }
  }

  /* Obtain global sum from local sums */
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  /* scale result and return */
  return(RSqrt(gsum / N_global));
}

/* 
 * N_VMin_SpcParallel returns the smallest element of x 
 */

realtype N_VMin_SpcParallel(N_Vector x)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype min, *xd, gmin; 
  MPI_Comm comm;

  gmin = min = BIG_REAL;
  xd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  comm = SPV_COMM(x);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  if (xd[loc] < min) min = xd[loc];
	}
      }
    }
  }

  /* Obtain global min from local mins */
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  wd  = SPV_DATA(w);
  comm = SPV_COMM(x);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  prodi = xd[loc] * wd[loc];
	  sum += prodi * prodi;
	}
      }
    }
  }

  /* Obtain global sum from local sums */
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  return(RSqrt(gsum));
}

/* 
 * N_VL1Norm_SpcParallel returns the L1 norm of x, 
 * i.e. sum (i=0 to N-1) {|x[i]|} 
 */

realtype N_VL1Norm_SpcParallel(N_Vector x)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype sum, gsum, *xd;
  MPI_Comm comm;

  gsum = sum = ZERO;
  xd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  comm = SPV_COMM(x);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  sum += ABS(xd[loc]);
	}
      }
    }
  }

  /* Obtain global sum from local sums */
  MPI_Allreduce(&sum, &gsum, 1, SPVEC_REAL_MPI_TYPE, MPI_SUM, comm);

  return(gsum);
}

/* 
 * N_VCompare_SpcParallel calculates z[i] = 1 if |x[i]| > c, z[i] = 0 otherwise 
 */

void N_VCompare_SpcParallel(realtype c, N_Vector x, N_Vector z)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;
  
  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = (ABS(xd[loc]) >= c) ? ONE : ZERO;
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);
  comm = SPV_COMM(x);

  /* Initialize return value */
  val = ONE;

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  if (xd[loc] == ZERO)  val = ZERO;
	  else  zd[loc] = ONE / xd[loc];
	}
      }
    }
  }

  /* Obtain global return value from local values */
  MPI_Allreduce(&val, &gval, 1, SPVEC_REAL_MPI_TYPE, MPI_MIN, comm);

  if (gval == ZERO)  return(FALSE);
  else  return(TRUE);
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  booleantype test, alltest;
  realtype *cd, *xd, *md;
  MPI_Comm comm;
 
  cd = xd = md = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  cd  = SPV_DATA(c);
  md  = SPV_DATA(m);
  comm = SPV_COMM(x);

  /* Initialize output variable */
  test = TRUE;

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
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
  
  /* Obtain global return value from local return values */
  MPI_Allreduce(&test, &alltest, 1, MPI_INT, MPI_MIN, comm);

  return(alltest);
}

/* 
 * N_VMinQuotient_SpcParallel returns min(num[i]/denom[i]) over all i 
 * such that denom[i] != 0. 
 */

realtype N_VMinQuotient_SpcParallel(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *nd, *dd, min, gmin;
  MPI_Comm comm;

  nd = dd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(num);
  Ny  = SPV_YLENGTH(num);
  Nz  = SPV_ZLENGTH(num);
  nspc  = SPV_NSPECIES(num);
  NGx = SPV_XGHOST(num);
  NGy = SPV_YGHOST(num);
  NGz = SPV_ZGHOST(num);
  nd  = SPV_DATA(num);
  dd  = SPV_DATA(denom);
  comm = SPV_COMM(num);

  /* Initialize output value, minimum */
  notEvenOnce = TRUE;
  min = BIG_REAL;

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  if (dd[loc] == ZERO) continue;
	  else {
	    if (!notEvenOnce) min = MIN(min, nd[loc] / dd[loc]);
	    else {
	      min = nd[loc] / dd[loc] ;
	      notEvenOnce = FALSE;
	    }
	  }
	}
      }
    }
  }

  /* Obtain global min from local minima */
  MPI_Allreduce(&min, &gmin, 1, SPVEC_REAL_MPI_TYPE, MPI_MIN, comm);

  return(gmin);
}
 
/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

/*
 * VCopy_SpcParallel is a private helper function that copies the 
 * local parts of the vector x into the local parts of the vector z. 
 * Note: the mesh dimensions must match. 
 */

static void VCopy_SpcParallel(N_Vector x, N_Vector z)
{
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = xd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = xd[loc] + yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;
 
  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = xd[loc] - yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *zd;

  xd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = -xd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  /*  for (i=0; i < N; i++),  *zd++ = c * ((*xd++) + (*yd++)); */
  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = c * (xd[loc] + yd[loc]);
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = c * (xd[loc] - yd[loc]);
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = a*xd[loc] + yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);
  zd  = SPV_DATA(z);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  zd[loc] = a*xd[loc] - yd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd, *yd;

  xd = yd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);
  yd  = SPV_DATA(y);

  if (a == ONE) {
    for (l=0; l<nspc; l++) {
      Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
      for (k=NGz; k<Nz+NGz; k++) {
	Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
	for (j=NGy; j<Ny+NGy; j++) { 
	  Yblock = j * (Nx + 2*NGx);
	  for (i=NGx; i<Nx+NGx; i++) {
	    loc = Sblock + Zblock + Yblock + i;
	    yd[loc] += xd[loc];
	  }
	}
      }
    }
    return;
  }

  if (a == -ONE) {
    for (l=0; l<nspc; l++) {
      Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
      for (k=NGz; k<Nz+NGz; k++) {
	Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
	for (j=NGy; j<Ny+NGy; j++) { 
	  Yblock = j * (Nx + 2*NGx);
	  for (i=NGx; i<Nx+NGx; i++) {
	    loc = Sblock + Zblock + Yblock + i;
	    yd[loc] -= xd[loc];
	  }
	}
      }
    }
    return;
  }

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  yd[loc] += a * xd[loc];
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
  int l, nspc;
  long int i, j, k, Nx, Ny, Nz, NGx, NGy, NGz;
  long int Yblock, Zblock, Sblock, loc;
  realtype *xd;

  xd = NULL;

  /* Get mesh info & data from vector */
  Nx  = SPV_XLENGTH(x);
  Ny  = SPV_YLENGTH(x);
  Nz  = SPV_ZLENGTH(x);
  nspc  = SPV_NSPECIES(x);
  NGx = SPV_XGHOST(x);
  NGy = SPV_YGHOST(x);
  NGz = SPV_ZGHOST(x);
  xd  = SPV_DATA(x);

  for (l=0; l<nspc; l++) {
    Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
    for (k=NGz; k<Nz+NGz; k++) {
      Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
      for (j=NGy; j<Ny+NGy; j++) { 
	Yblock = j * (Nx + 2*NGx);
	for (i=NGx; i<Nx+NGx; i++) {
	  loc = Sblock + Zblock + Yblock + i;
	  xd[loc] *= a;
	} 
      } 
    } 
  }
}

