/*----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 -----------------------------------------------------------------
 This is the implementation file for a grid parallel MPI 
 implementation of the NVECTOR package.
 ---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_parallel_grid.h>
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error Message */

#define BAD_N1 "N_VNew_Parallel_Grid -- Sum of local vector lengths "
#define BAD_N2 "differs from input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/* Domain looping macros */
#define NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) \
for (i5=NV_OFFSET_PG(x,5); i5<NV_OFFSET_PG(x,5)+NV_ACTIVELEN_PG(x,5); i5++) {	\
  for (i4=NV_OFFSET_PG(x,4); i4<NV_OFFSET_PG(x,4)+NV_ACTIVELEN_PG(x,4); i4++) { \
    for (i3=NV_OFFSET_PG(x,3); i3<NV_OFFSET_PG(x,3)+NV_ACTIVELEN_PG(x,3); i3++) { \
      for (i2=NV_OFFSET_PG(x,2); i2<NV_OFFSET_PG(x,2)+NV_ACTIVELEN_PG(x,2); i2++) { \
        for (i1=NV_OFFSET_PG(x,1); i1<NV_OFFSET_PG(x,1)+NV_ACTIVELEN_PG(x,1); i1++) { \
          for (i0=NV_OFFSET_PG(x,0); i0<NV_OFFSET_PG(x,0)+NV_ACTIVELEN_PG(x,0); i0++) { \
  	    i = ((((i5*NV_ARRAYLEN_PG(x,4) + i4)*NV_ARRAYLEN_PG(x,3) + i3)*NV_ARRAYLEN_PG(x,2) + i2)*NV_ARRAYLEN_PG(x,1) + i1)*NV_ARRAYLEN_PG(x,0) + i0;

#define NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) \
for (i0=NV_OFFSET_PG(x,0); i0<NV_OFFSET_PG(x,0)+NV_ACTIVELEN_PG(x,0); i0++) { \
  for (i1=NV_OFFSET_PG(x,1); i1<NV_OFFSET_PG(x,1)+NV_ACTIVELEN_PG(x,1); i1++) { \
    for (i2=NV_OFFSET_PG(x,2); i2<NV_OFFSET_PG(x,2)+NV_ACTIVELEN_PG(x,2); i2++) { \
      for (i3=NV_OFFSET_PG(x,3); i3<NV_OFFSET_PG(x,3)+NV_ACTIVELEN_PG(x,3); i3++) { \
        for (i4=NV_OFFSET_PG(x,4); i4<NV_OFFSET_PG(x,4)+NV_ACTIVELEN_PG(x,4); i4++) { \
	  for (i5=NV_OFFSET_PG(x,5); i5<NV_OFFSET_PG(x,5)+NV_ACTIVELEN_PG(x,5); i5++) { \
	    i = ((((i0*NV_ARRAYLEN_PG(x,1) + i1)*NV_ARRAYLEN_PG(x,2) + i2)*NV_ARRAYLEN_PG(x,3) + i3)*NV_ARRAYLEN_PG(x,4) + i4)*NV_ARRAYLEN_PG(x,5) + i5;


#define NV_LOOPEND_PG()	\
          } \
        } \
      } \
    } \
  } \
}

/* Private function prototypes */

/* vector compatability check */
static booleantype VCheck_Compatible(N_Vector x, N_Vector z);
/* vector data length */
static long int NV_DATALEN_PG(N_Vector x);
/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_Parallel_Grid(realtype d, int op, MPI_Comm comm);
/* z=x */
static void VCopy_Parallel_Grid(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_Parallel_Grid(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_Parallel_Grid(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_Parallel_Grid(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_Parallel_Grid(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_Parallel_Grid(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_Parallel_Grid(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_Parallel_Grid(realtype a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Function to create a new parallel grid vector with empty data array
 */

N_Vector N_VNewEmpty_Parallel_Grid(MPI_Comm comm, 
				   long int dims,
				   long int *dim_length,
				   long int *dim_alength,
				   long int *dim_offset,
				   long int F_ordering,
				   long int global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel_Grid content;
  long int n, Nsum;
  int i;

  /* ensure that 0 < dims <= MAX_DIMS */
  if (dims > MAX_DIMS) {
    printf("N_VNew_Parallel_Grid error -- dims must be at most %i\n\n",
	   MAX_DIMS);
    return(NULL);
  }
  if (dims == 0) {
    printf("N_VNew_Parallel_Grid error -- dims must be at least 1\n\n");
    return(NULL);
  }

  /* Compute total local length */
  n = 1;
  for (i=0; i<dims; i++)  n *= dim_alength[i];

  /* Compute global length as sum of local lengths */
  MPI_Allreduce(&n, &Nsum, 1, PGVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  if (Nsum != global_length) {
    printf(BAD_N);
    return(NULL);
  } 

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvclone           = N_VClone_Parallel_Grid;
  ops->nvcloneempty      = N_VCloneEmpty_Parallel_Grid;
  ops->nvdestroy         = N_VDestroy_Parallel_Grid;
  ops->nvspace           = N_VSpace_Parallel_Grid;
  ops->nvgetarraypointer = N_VGetArrayPointer_Parallel_Grid;
  ops->nvsetarraypointer = N_VSetArrayPointer_Parallel_Grid;
  ops->nvlinearsum       = N_VLinearSum_Parallel_Grid;
  ops->nvconst           = N_VConst_Parallel_Grid;
  ops->nvprod            = N_VProd_Parallel_Grid;
  ops->nvdiv             = N_VDiv_Parallel_Grid;
  ops->nvscale           = N_VScale_Parallel_Grid;
  ops->nvabs             = N_VAbs_Parallel_Grid;
  ops->nvinv             = N_VInv_Parallel_Grid;
  ops->nvaddconst        = N_VAddConst_Parallel_Grid;
  ops->nvdotprod         = N_VDotProd_Parallel_Grid;
  ops->nvmaxnorm         = N_VMaxNorm_Parallel_Grid;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Parallel_Grid;
  ops->nvwrmsnorm        = N_VWrmsNorm_Parallel_Grid;
  ops->nvmin             = N_VMin_Parallel_Grid;
  ops->nvwl2norm         = N_VWL2Norm_Parallel_Grid;
  ops->nvl1norm          = N_VL1Norm_Parallel_Grid;
  ops->nvcompare         = N_VCompare_Parallel_Grid;
  ops->nvinvtest         = N_VInvTest_Parallel_Grid;
  ops->nvconstrmask      = N_VConstrMask_Parallel_Grid;
  ops->nvminquotient     = N_VMinQuotient_Parallel_Grid;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Parallel_Grid) 
    malloc(sizeof(struct _N_VectorContent_Parallel_Grid));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach vector components */
  content->dims          = dims;
  for (i=0; i<dims; i++)  content->dim_length[i]  = dim_length[i];
  for (i=0; i<dims; i++)  content->dim_alength[i] = dim_alength[i];
  for (i=0; i<dims; i++)  content->dim_offset[i]  = dim_offset[i];
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = FALSE;
  content->F_ordering    = (F_ordering == 1) ? TRUE : FALSE;
  content->data          = NULL;

  /* set additional dimensions to have length 1, offset 0 */
  for (i=dims; i<MAX_DIMS; i++)  content->dim_length[i]  = 1;
  for (i=dims; i<MAX_DIMS; i++)  content->dim_alength[i] = 1;
  for (i=dims; i<MAX_DIMS; i++)  content->dim_offset[i]  = 0;

  /* verify that array lengths are legal */
  for (i=0; i<MAX_DIMS; i++) 
    if (content->dim_alength[i]+content->dim_offset[i] > content->dim_length[i])
      {
	printf("N_VNewEmpty_Parallel_Grid: illegal inputs for dimension %i:\n",i);
	printf("  total length (%li) must exceed active length (%li) plus offset (%li)\n",
	       content->dim_length[i],content->dim_alength[i],content->dim_offset[i]);
	return NULL;
      }

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a new parallel grid vector
 */

N_Vector N_VNew_Parallel_Grid(MPI_Comm comm, 
			      long int dims,
			      long int *dim_length,
			      long int *dim_alength,
			      long int *dim_offset,
			      long int F_ordering,
			      long int global_length)
{
  N_Vector v;
  realtype *data;
  long int i, local_length;

  v = NULL;
  v = N_VNewEmpty_Parallel_Grid(comm, dims, dim_length, dim_alength, 
				dim_offset, F_ordering, global_length);
  if (v == NULL) return(NULL);

  /* Compute total local length */
  local_length = (dims == 0) ? 0 : 1;
  for (i=0; i<dims; i++)  local_length *= dim_length[i];

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory, initialize to zero */
    data = NULL;
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Parallel_Grid(v); return(NULL); }
    for (i=0; i<local_length; i++)  data[i] = 0.0;

    /* Attach data */
    NV_OWN_DATA_PG(v) = TRUE;
    NV_DATA_PG(v)     = data; 

  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parallel grid N_Vector with user data component 
 */

N_Vector N_VMake_Parallel_Grid(MPI_Comm comm, 
			       long int dims,
			       long int *dim_length,
			       long int *dim_alength,
			       long int *dim_offset,
			       long int F_ordering,
			       long int global_length,
			       realtype *v_data)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Parallel_Grid(comm, dims, dim_length, dim_alength, 
				dim_offset, F_ordering, global_length);
  if (v == NULL) return(NULL);

  if (NV_DATALEN_PG(v) > 0) {
    /* Attach data */
    NV_OWN_DATA_PG(v) = FALSE;
    NV_DATA_PG(v)     = v_data;
  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parallel grid vectors. 
 */

N_Vector *N_VCloneVectorArray_Parallel_Grid(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Parallel_Grid(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel_Grid(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parallel grid vectors with empty
 * (NULL) data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_Parallel_Grid(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Parallel_Grid(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel_Grid(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Parallel_Grid
 */

void N_VDestroyVectorArray_Parallel_Grid(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Parallel_Grid(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ---------------------------------------------------------------- 
 * Function to print the active portion of a parallel grid vector 
 */

void N_VPrint_Parallel_Grid(N_Vector x)
{
  long int i, i0, i1, i2, i3, i4, i5;
  long int dims, N[MAX_DIMS], n[MAX_DIMS], o[MAX_DIMS];
  booleantype Forder;
  realtype *xd = NULL;

  /* get array dimensions */
  Forder = NV_FORDER_PG(x);
  dims = NV_DIMS_PG(x);
  for (i=0; i<MAX_DIMS; i++)  N[i] = NV_ARRAYLEN_PG(x,i);
  for (i=0; i<MAX_DIMS; i++)  n[i] = NV_ACTIVELEN_PG(x,i);
  for (i=0; i<MAX_DIMS; i++)  o[i] = NV_OFFSET_PG(x,i);

  /* access array data */
  xd = NV_DATA_PG(x);

  /* iterate over the vector and output, depending on dimensionality 
     and data ordering */
  if (dims == 1) 
    for (i0=o[0]; i0<o[0]+n[0]; i0++) 
      printf("%li : %"DOUT"\n", i0-o[0], xd[i0]);
  if (dims == 2) {
    if (Forder) {
      for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	for (i0=o[0]; i0<o[0]+n[0]; i0++) 
	  printf("%li, %li : %"DOUT"\n", i0-o[0], i1-o[1], xd[i1*N[0]+i0]);
    } else {
      for (i0=o[0]; i0<o[0]+n[0]; i0++) 
	for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	  printf("%li, %li : %"DOUT"\n", i0-o[0], i1-o[1], xd[i0*N[1]+i1]);
    }
  }
  if (dims == 3) {
    if (Forder) {
      for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	  for (i0=o[0]; i0<o[0]+n[0]; i0++) {
	    i = (i2*N[1] + i1)*N[0] + i0;
	    printf("%li, %li, %li : %"DOUT"\n", 
		   i0-o[0], i1-o[1], i2-o[2], xd[i]);
	  } 
    } else {
      for (i0=o[0]; i0<o[0]+n[0]; i0++) 
	for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	  for (i2=o[2]; i2<o[2]+n[2]; i2++) {
	    i = (i0*N[1] + i1)*N[2] + i2;
	    printf("%li, %li, %li : %"DOUT"\n", 
		   i0-o[0], i1-o[1], i2-o[2], xd[i]);
	  }
    }
  }
  if (dims == 4) {
    if (Forder) {
      for (i3=o[3]; i3<o[3]+n[3]; i3++) 
	for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	  for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	    for (i0=o[0]; i0<o[0]+n[0]; i0++) {
	      i = ((i3*N[2] + i2)*N[1] + i1)*N[0] + i0;
	      printf("%li, %li, %li, %li : %"DOUT"\n", 
		     i0-o[0], i1-o[1], i2-o[2], i3-o[3], xd[i]);
	    }
    } else {
      for (i0=o[0]; i0<o[0]+n[0]; i0++) 
	for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	  for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	    for (i3=o[3]; i3<o[3]+n[3]; i3++) {
	      i = ((i0*N[1] + i1)*N[2] + i2)*N[3] + i3;
	      printf("%li, %li, %li, %li : %"DOUT"\n", 
		     i0-o[0], i1-o[1], i2-o[2], i3-o[3], xd[i]);
	    }
    }
  }
  if (dims == 5) {
    if (Forder) {
      for (i4=o[4]; i4<o[4]+n[4]; i4++) 
	for (i3=o[3]; i3<o[3]+n[3]; i3++) 
	  for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	    for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	      for (i0=o[0]; i0<o[0]+n[0]; i0++) {
		i = (((i4*N[3] + i3)*N[2] + i2)*N[1] + i1)*N[0] + i0;
		printf("%li, %li, %li, %li, %li : %"DOUT"\n", 
		       i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], xd[i]);
	      }
    } else {
      for (i0=o[0]; i0<o[0]+n[0]; i0++) 
	for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	  for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	    for (i3=o[3]; i3<o[3]+n[3]; i3++) 
	      for (i4=o[4]; i4<o[4]+n[4]; i4++) {
		i = (((i0*N[1] + i1)*N[2] + i2)*N[3] + i3)*N[4] + i4;
		printf("%li, %li, %li, %li, %li : %"DOUT"\n", 
		       i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], xd[i]);
	      }
    }
  }
  if (dims == 6) {
    if (Forder) {
      for (i5=o[5]; i5<o[5]+n[5]; i5++) 
	for (i4=o[4]; i4<o[4]+n[4]; i4++) 
	  for (i3=o[3]; i3<o[3]+n[3]; i3++) 
	    for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	      for (i1=o[1]; i1<o[1]+n[1]; i1++) 
		for (i0=o[0]; i0<o[0]+n[0]; i0++) {
		  i = ((((i5*N[4] + i4)*N[3] + i3)*N[2] + i2)*N[1] + i1)*N[0] + i0;
		  printf("%li, %li, %li, %li, %li, %li : %"DOUT"\n", 
			 i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], i5-o[5], xd[i]);
		}
    } else {
    for (i0=o[0]; i0<o[0]+n[0]; i0++) 
      for (i1=o[1]; i1<o[1]+n[1]; i1++) 
	for (i2=o[2]; i2<o[2]+n[2]; i2++) 
	  for (i3=o[3]; i3<o[3]+n[3]; i3++) 
	    for (i4=o[4]; i4<o[4]+n[4]; i4++) 
	      for (i5=o[5]; i5<o[5]+n[5]; i5++) {
		i = ((((i0*N[1] + i1)*N[2] + i2)*N[3] + i3)*N[4] + i4)*N[5] + i5;
		printf("%li, %li, %li, %li, %li, %li : %"DOUT"\n", 
		       i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], i5-o[5], xd[i]);
	      }
    }
  }
  
  printf("\n");
  return;
}

/* ---------------------------------------------------------------- 
 * Function to print all data in a parallel grid vector 
 */

void N_VPrintAll_Parallel_Grid(N_Vector x)
{
  long int i, i0, i1, i2, i3, i4, i5;
  long int dims, N[MAX_DIMS], o[MAX_DIMS];
  booleantype Forder;
  realtype *xd = NULL;

  /* get array dimensions */
  dims = NV_DIMS_PG(x);
  Forder = NV_FORDER_PG(x);
  for (i=0; i<MAX_DIMS; i++)  N[i] = NV_ARRAYLEN_PG(x,i);
  for (i=0; i<MAX_DIMS; i++)  o[i] = NV_OFFSET_PG(x,i);

  /* access array data */
  xd = NV_DATA_PG(x);

  /* iterate over the vector and output, depending on dimensionality 
     and data ordering */
  if (dims == 1) 
    for (i0=0; i0<N[0]; i0++) 
      printf("%li : %"DOUT"\n", i0-o[0], xd[i0]);
  if (dims == 2) {
    if (Forder) {
      for (i1=0; i1<N[1]; i1++) 
	for (i0=0; i0<N[0]; i0++) 
	  printf("%li, %li : %"DOUT"\n", i0-o[0], i1-o[1], xd[i1*N[0]+i0]);
    } else {
      for (i0=0; i0<N[0]; i0++) 
	for (i1=0; i1<N[1]; i1++) 
	  printf("%li, %li : %"DOUT"\n", i0-o[0], i1-o[1], xd[i0*N[1]+i1]);
    }
  }
  if (dims == 3) {
    if (Forder) {
      for (i2=0; i2<N[2]; i2++) 
	for (i1=0; i1<N[1]; i1++) 
	  for (i0=0; i0<N[0]; i0++) {
	    i = (i2*N[1] + i1)*N[0] + i0;
	    printf("%li, %li, %li : %"DOUT"\n", 
		   i0-o[0], i1-o[1], i2-o[2], xd[i]);
	  } 
    } else {
      for (i0=0; i0<N[0]; i0++) 
	for (i1=0; i1<N[1]; i1++) 
	  for (i2=0; i2<N[2]; i2++) {
	    i = (i0*N[1] + i1)*N[2] + i2;
	    printf("%li, %li, %li : %"DOUT"\n", 
		   i0-o[0], i1-o[1], i2-o[2], xd[i]);
	  }
    }
  }
  if (dims == 4) {
    if (Forder) {
      for (i3=0; i3<N[3]; i3++) 
	for (i2=0; i2<N[2]; i2++) 
	  for (i1=0; i1<N[1]; i1++) 
	    for (i0=0; i0<N[0]; i0++) {
	      i = ((i3*N[2] + i2)*N[1] + i1)*N[0] + i0;
	      printf("%li, %li, %li, %li : %"DOUT"\n", 
		     i0-o[0], i1-o[1], i2-o[2], i3-o[3], xd[i]);
	    }
    } else {
      for (i0=0; i0<N[0]; i0++) 
	for (i1=0; i1<N[1]; i1++) 
	  for (i2=0; i2<N[2]; i2++) 
	    for (i3=0; i3<N[3]; i3++) {
	      i = ((i0*N[1] + i1)*N[2] + i2)*N[3] + i3;
	      printf("%li, %li, %li, %li : %"DOUT"\n", 
		     i0-o[0], i1-o[1], i2-o[2], i3-o[3], xd[i]);
	    }
    }
  }
  if (dims == 5) {
    if (Forder) {
      for (i4=0; i4<N[4]; i4++) 
	for (i3=0; i3<N[3]; i3++) 
	  for (i2=0; i2<N[2]; i2++) 
	    for (i1=0; i1<N[1]; i1++) 
	      for (i0=0; i0<N[0]; i0++) {
		i = (((i4*N[3] + i3)*N[2] + i2)*N[1] + i1)*N[0] + i0;
		printf("%li, %li, %li, %li, %li : %"DOUT"\n", 
		       i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], xd[i]);
	      }
    } else {
      for (i0=0; i0<N[0]; i0++) 
	for (i1=0; i1<N[1]; i1++) 
	  for (i2=0; i2<N[2]; i2++) 
	    for (i3=0; i3<N[3]; i3++) 
	      for (i4=0; i4<N[4]; i4++) {
		i = (((i0*N[1] + i1)*N[2] + i2)*N[3] + i3)*N[4] + i4;
		printf("%li, %li, %li, %li, %li : %"DOUT"\n", 
		       i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], xd[i]);
	      }
    }
  }
  if (dims == 6) {
    if (Forder) {
      for (i5=0; i5<N[5]; i5++) 
	for (i4=0; i4<N[4]; i4++) 
	  for (i3=0; i3<N[3]; i3++) 
	    for (i2=0; i2<N[2]; i2++) 
	      for (i1=0; i1<N[1]; i1++) 
		for (i0=0; i0<N[0]; i0++) {
		  i = ((((i5*N[4] + i4)*N[3] + i3)*N[2] + i2)*N[1] + i1)*N[0] + i0;
		  printf("%li, %li, %li, %li, %li, %li : %"DOUT"\n", 
			 i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], i5-o[5], xd[i]);
		}
    } else {
    for (i0=0; i0<N[0]; i0++) 
      for (i1=0; i1<N[1]; i1++) 
	for (i2=0; i2<N[2]; i2++) 
	  for (i3=0; i3<N[3]; i3++) 
	    for (i4=0; i4<N[4]; i4++) 
	      for (i5=0; i5<N[5]; i5++) {
		i = ((((i0*N[1] + i1)*N[2] + i2)*N[3] + i3)*N[4] + i4)*N[5] + i5;
		printf("%li, %li, %li, %li, %li, %li : %"DOUT"\n", 
		       i0-o[0], i1-o[1], i2-o[2], i3-o[3], i4-o[4], i5-o[5], xd[i]);
	      }
    }
  }
  
  printf("\n");
  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Parallel_Grid(N_Vector w)
{
  int i;
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel_Grid content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }
  
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
  content = NULL;
  content = (N_VectorContent_Parallel_Grid) 
    malloc(sizeof(struct _N_VectorContent_Parallel_Grid));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->dims          = NV_DIMS_PG(w);
  content->global_length = NV_GLOBLENGTH_PG(w);
  content->comm          = NV_COMM_PG(w);
  content->own_data      = FALSE;
  content->data          = NULL;
  content->F_ordering    = NV_FORDER_PG(w);
  for (i=0; i<MAX_DIMS; i++)  content->dim_length[i]  = NV_ARRAYLEN_PG(w,i);
  for (i=0; i<MAX_DIMS; i++)  content->dim_alength[i] = NV_ACTIVELEN_PG(w,i);
  for (i=0; i<MAX_DIMS; i++)  content->dim_offset[i]  = NV_OFFSET_PG(w,i);

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

N_Vector N_VClone_Parallel_Grid(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int local_length, i;

  v = NULL;
  v = N_VCloneEmpty_Parallel_Grid(w);
  if (v == NULL) return(NULL);

  local_length = NV_DATALEN_PG(w);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Parallel_Grid(v); return(NULL); }
    for (i=0; i<local_length; i++)  data[i] = 0.0;

    /* Attach data */
    NV_OWN_DATA_PG(v) = TRUE;
    NV_DATA_PG(v)     = data;
  }

  return(v);
}

void N_VDestroy_Parallel_Grid(N_Vector v)
{
  if ((NV_OWN_DATA_PG(v) == TRUE) && (NV_DATA_PG(v) != NULL)) {
    free(NV_DATA_PG(v));
    NV_DATA_PG(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Parallel_Grid(N_Vector v, long int *lrw, long int *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_COMM_PG(v);
  MPI_Comm_size(comm, &npes);
  
  *lrw = NV_GLOBLENGTH_PG(v);
  *liw = (2+3*MAX_DIMS)*npes;

  return;
}

realtype *N_VGetArrayPointer_Parallel_Grid(N_Vector v)
{
  return((realtype *) NV_DATA_PG(v));
}

void N_VSetArrayPointer_Parallel_Grid(realtype *v_data, N_Vector v)
{
  if (NV_DATALEN_PG(v) > 0) NV_DATA_PG(v) = v_data;

  return;
}

void N_VLinearSum_Parallel_Grid(realtype a, N_Vector x, realtype b, 
				N_Vector y, N_Vector z)
{
  long int i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;
  xd = yd = zd = NULL;

  /* first check that N_Vectors are compatible */
  if (!VCheck_Compatible(x,y)) {
    fprintf(stderr,"N_VLinearSum_Parallel_Grid error: x,y incompatible\n");
    return;
  }
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VLinearSum_Parallel_Grid error: x,z incompatible\n");
    return;
  }


  /* usage scenarios */

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Parallel_Grid(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_Parallel_Grid(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_Parallel_Grid(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Parallel_Grid(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Parallel_Grid(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Parallel_Grid(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_Parallel_Grid(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_Parallel_Grid(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  /* access array data */
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* set total local vector length */
  N = NV_DATALEN_PG(x);

  /* iterate over entire vector data, including grid cells */
  for (i=0; i<N; i++)  zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}

void N_VConst_Parallel_Grid(realtype c, N_Vector z)
{
  long int i, N;
  realtype *zd = NV_DATA_PG(z);

  /* set all entries of z to the constant (including grid zones) */
  N = NV_DATALEN_PG(z);
  for (i=0; i<N; i++)  zd[i] = c;

  return;
}

void N_VProd_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z)
{
  realtype *xd, *yd, *zd;
  long int i, N;
;
  xd = yd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,y)) {
    fprintf(stderr,"N_VProd_Parallel_Grid error: x,y incompatible\n");
    return;
  }
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VProd_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  /* access data arrays, length */
  N = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform product (even on grid zones) */
  for (i=0; i<N; i++)  zd[i] = xd[i]*yd[i];

  return;
}

void N_VDiv_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z)
{
  realtype *xd, *yd, *zd;
  long int i0, i1, i2, i3, i4, i5, i;
  xd = yd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,y)) {
    fprintf(stderr,"N_VDiv_Parallel_Grid error: x,y incompatible\n");
    return;
  }
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VDiv_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform division on domain interior */
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      zd[i] = xd[i]/yd[i];
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      zd[i] = xd[i]/yd[i];
    }
    NV_LOOPEND_PG();
  }

  return;
}

void N_VScale_Parallel_Grid(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  xd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VScale_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy_Parallel_Grid(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Parallel_Grid(x, z);
  } else if (c == -ONE) {
    VNeg_Parallel_Grid(x, z);
  } else {
    N  = NV_DATALEN_PG(x);
    xd = NV_DATA_PG(x);
    zd = NV_DATA_PG(z);
    for (i=0; i<N; i++)   zd[i] = c*xd[i];
  }

  return;
}

void N_VAbs_Parallel_Grid(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  xd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VAbs_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);
  for (i=0; i<N; i++)   zd[i] = SUN_ABS(xd[i]);

  return;
}

void N_VInv_Parallel_Grid(N_Vector x, N_Vector z)
{
  realtype *xd, *zd;
  long int i0, i1, i2, i3, i4, i5, i;
  xd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VAbs_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);

  /* perform inversion on domain interior */
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      zd[i] = ONE/xd[i];
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      zd[i] = ONE/xd[i];
    }
    NV_LOOPEND_PG();
  }


  return;
}

void N_VAddConst_Parallel_Grid(N_Vector x, realtype b, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  xd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VAddConst_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);
  for (i=0; i<N; i++)  zd[i] = xd[i]+b;

  return;
}

realtype N_VDotProd_Parallel_Grid(N_Vector x, N_Vector y)
{
  realtype sum, *xd, *yd, gsum;
  MPI_Comm comm;
  long int i0, i1, i2, i3, i4, i5, i;
  xd = yd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,y)) {
    fprintf(stderr,"N_VAddConst_Parallel_Grid error: x,y incompatible\n");
    return(-1.0);
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);

  /* perform dot product on domain interior */
  sum = ZERO;
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      sum += xd[i]*yd[i];
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      sum += xd[i]*yd[i];
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global sum */
  comm = NV_COMM_PG(x);
  gsum = VAllReduce_Parallel_Grid(sum, 1, comm);

  return(gsum);
}

realtype N_VMaxNorm_Parallel_Grid(N_Vector x)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype max, *xd, gmax;
  MPI_Comm comm;
  xd = NULL;

  /* access data array */
  xd = NV_DATA_PG(x);

  /* perform max norm on domain interior */
  max = ZERO;
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      if (SUN_ABS(xd[i]) > max) max = SUN_ABS(xd[i]);
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      if (SUN_ABS(xd[i]) > max) max = SUN_ABS(xd[i]);
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global max */
  comm = NV_COMM_PG(x);
  gmax = VAllReduce_Parallel_Grid(max, 2, comm);

  return(gmax);
}

realtype N_VWrmsNorm_Parallel_Grid(N_Vector x, N_Vector w)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype sum, *xd, *wd, gsum, prodi;
  MPI_Comm comm;
  xd = wd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,w)) {
    fprintf(stderr,"N_VWrmsNorm_Parallel_Grid error: x,w incompatible\n");
    return(-1.0);
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  wd = NV_DATA_PG(w);

  /* perform wrms norm on domain interior */
  sum = ZERO;
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      prodi = xd[i]*wd[i];
      sum += SUN_SQR(prodi);
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      prodi = xd[i]*wd[i];
      sum += SUN_SQR(prodi);
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global sum */
  comm = NV_COMM_PG(x);
  gsum = VAllReduce_Parallel_Grid(sum, 1, comm);

  return(SUN_SQRT(gsum/NV_GLOBLENGTH_PG(x)));
}

realtype N_VWrmsNormMask_Parallel_Grid(N_Vector x, N_Vector w, N_Vector id)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype sum, prodi, *xd, *wd, *idd, gsum;
  MPI_Comm comm;
  xd = wd = idd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,w)) {
    fprintf(stderr,"N_VWrmsNormMask_Parallel_Grid error: x,w incompatible\n");
    return(-1.0);
  }
  if (!VCheck_Compatible(x,id)) {
    fprintf(stderr,"N_VWrmsNormMask_Parallel_Grid error: x,id incompatible\n");
    return(-1.0);
  }

  /* access data arrays */
  xd  = NV_DATA_PG(x);
  wd  = NV_DATA_PG(w);
  idd = NV_DATA_PG(id);

  /* perform masked wrms norm on domain interior */
  sum = ZERO;
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      if (idd[i] > ZERO) {
	prodi = xd[i]*wd[i];
	sum += SUN_SQR(prodi);
      }
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      if (idd[i] > ZERO) {
	prodi = xd[i]*wd[i];
	sum += SUN_SQR(prodi);
      }
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global sum */
  comm = NV_COMM_PG(x);
  gsum = VAllReduce_Parallel_Grid(sum, 1, comm);

  return(SUN_SQRT(gsum/NV_GLOBLENGTH_PG(x)));
}

realtype N_VMin_Parallel_Grid(N_Vector x)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype gmin;
  MPI_Comm comm;
  realtype *xd = NULL;
  realtype min = BIG_REAL;

  if (NV_DATALEN_PG(x) > 0) {

    /* access data array */
    xd = NV_DATA_PG(x);

    /* initialize result */
    min = xd[0];

    /* iterate over active values to compute minimum */
    if (NV_FORDER_PG(x)) {
      NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
	if (xd[i] < min)  min = xd[i];
      }
      NV_LOOPEND_PG();
    } else {
      NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
	if (xd[i] < min)  min = xd[i];
      }
      NV_LOOPEND_PG();
    }

  }

  /* communicate to obtain global min */
  comm = NV_COMM_PG(x);
  gmin = VAllReduce_Parallel_Grid(min, 3, comm);

  return(gmin);
}

realtype N_VWL2Norm_Parallel_Grid(N_Vector x, N_Vector w)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype prodi, *xd, *wd, gsum;
  MPI_Comm comm;
  realtype sum = ZERO;
  xd = wd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,w)) {
    fprintf(stderr,"N_VWl2Norm_Parallel_Grid error: x,w incompatible\n");
    return(-1.0);
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  wd = NV_DATA_PG(w);
  
  /* iterate over active values to compute norm */
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      prodi = xd[i]*wd[i];
      sum += SUN_SQR(prodi);
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      prodi = xd[i]*wd[i];
      sum += SUN_SQR(prodi);
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global sum */
  comm = NV_COMM_PG(x);
  gsum = VAllReduce_Parallel_Grid(sum, 1, comm);

  return(SUN_SQRT(gsum));
}

realtype N_VL1Norm_Parallel_Grid(N_Vector x)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype gsum;
  MPI_Comm comm;
  realtype sum = ZERO;
  realtype *xd = NULL;

  /* access data array */
  xd = NV_DATA_PG(x);

  /* iterate over active values to compute norm */
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      sum += SUN_ABS(xd[i]);
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      sum += SUN_ABS(xd[i]);
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global sum */
  comm = NV_COMM_PG(x);
  gsum = VAllReduce_Parallel_Grid(sum, 1, comm);

  return(gsum);
}

void N_VCompare_Parallel_Grid(realtype c, N_Vector x, N_Vector z)
{
  long int N, i;
  realtype *xd, *zd;
  xd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VCompare_Parallel_Grid error: x,z incompatible\n");
    return;
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);

  /* perform operation on full domain (including grids) */
  N = NV_DATALEN_PG(x);
  for (i=0; i<N; i++) 
    zd[i] = (SUN_ABS(xd[i]) >= c) ? ONE : ZERO;

  return;
}

booleantype N_VInvTest_Parallel_Grid(N_Vector x, N_Vector z)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;
  xd = zd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,z)) {
    fprintf(stderr,"N_VInvTest_Parallel_Grid error: x,z incompatible\n");
    return(TRUE);
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);

  /* perform test on domain interior */
  val = ONE;
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      if (xd[i] == ZERO) 
	val = ZERO;
      else
	zd[i] = ONE/xd[i];
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      if (xd[i] == ZERO) 
	val = ZERO;
      else
	zd[i] = ONE/xd[i];
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global min */
  comm = NV_COMM_PG(x);
  gval = VAllReduce_Parallel_Grid(val, 3, comm);

  if (gval == ZERO)
    return(FALSE);
  else
    return(TRUE);
}

booleantype N_VConstrMask_Parallel_Grid(N_Vector c, N_Vector x, N_Vector m)
{
  long int i0, i1, i2, i3, i4, i5, i;
  realtype temp;
  realtype *cd, *xd, *md;
  MPI_Comm comm;
  cd = xd = md = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(x,c)) {
    fprintf(stderr,"N_VConstrMask_Parallel_Grid error: x,c incompatible\n");
    return(TRUE);
  }
  if (!VCheck_Compatible(x,m)) {
    fprintf(stderr,"N_VConstrMask_Parallel_Grid error: x,m incompatible\n");
    return(TRUE);
  }

  /* access data arrays */
  xd = NV_DATA_PG(x);
  cd = NV_DATA_PG(c);
  md = NV_DATA_PG(m);


  /* perform operation on domain interior */
  temp = ONE;
  if (NV_FORDER_PG(x)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      md[i] = ZERO;
      if (cd[i] == ZERO) continue;
      if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
	if (xd[i]*cd[i] <= ZERO) { temp = ZERO; md[i] = ONE; }
	continue;
      }
      if (cd[i] > HALF || cd[i] < -HALF) {
	if (xd[i]*cd[i] < ZERO ) { temp = ZERO; md[i] = ONE; }
      }
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, x) {
      md[i] = ZERO;
      if (cd[i] == ZERO) continue;
      if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
	if (xd[i]*cd[i] <= ZERO) { temp = ZERO; md[i] = ONE; }
	continue;
      }
      if (cd[i] > HALF || cd[i] < -HALF) {
	if (xd[i]*cd[i] < ZERO ) { temp = ZERO; md[i] = ONE; }
      }
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global min */
  comm = NV_COMM_PG(x);
  temp = VAllReduce_Parallel_Grid(temp, 3, comm);

  if (temp == ONE) return(TRUE);
  else return(FALSE);
}

realtype N_VMinQuotient_Parallel_Grid(N_Vector num, N_Vector denom)
{
  long int i0, i1, i2, i3, i4, i5, i;
  booleantype notEvenOnce;
  realtype *nd, *dd, min;
  MPI_Comm comm;
  nd = dd = NULL;

  /* check for compatibility */
  if (!VCheck_Compatible(num,denom)) {
    fprintf(stderr,"N_VMinQuotient_Parallel_Grid error: num,denom incompatible\n");
    return(-1.0);
  }

  /* access data arrays */
  nd = NV_DATA_PG(num);
  dd = NV_DATA_PG(denom);

  /* perform operation on domain interior */
  notEvenOnce = TRUE;
  min = BIG_REAL;
  if (NV_FORDER_PG(num)) {
    NV_FLOOP_PG(i0, i1, i2, i3, i4, i5, i, num) {
      if (dd[i] == ZERO) continue;
      else {
	if (!notEvenOnce) min = SUN_MIN(min, nd[i]/dd[i]);
	else {
	  min = nd[i]/dd[i];
	  notEvenOnce = FALSE;
	}
      }
    }
    NV_LOOPEND_PG();
  } else {
    NV_CLOOP_PG(i0, i1, i2, i3, i4, i5, i, num) {
      if (dd[i] == ZERO) continue;
      else {
	if (!notEvenOnce) min = SUN_MIN(min, nd[i]/dd[i]);
	else {
	  min = nd[i]/dd[i];
	  notEvenOnce = FALSE;
	}
      }
    }
    NV_LOOPEND_PG();
  }

  /* communicate to obtain global min */
  comm = NV_COMM_PG(num);
  return(VAllReduce_Parallel_Grid(min, 3, comm));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static booleantype VCheck_Compatible(N_Vector x, N_Vector y)
{
  /* This function checks that the two input vector layouts match */

  long int N, M, i;

  /* check for matching dims */
  N = NV_DIMS_PG(x);
  M = NV_DIMS_PG(y);
  if (N != M) {
    fprintf(stderr,"VCheck_Compatible: x,y dims mismatch (%li vs %li)\n",N,M);
    return FALSE;
  }

  /* check for matching data size */
  for (i=0; i<MAX_DIMS; i++) {
    N = NV_ARRAYLEN_PG(x,i);
    M = NV_ARRAYLEN_PG(y,i);
    if (N != M) {
      fprintf(stderr,"VCheck_Compatible: x,y size mismatch (dim %li: %li vs %li)\n",i,N,M);
      return FALSE;
    }
  }

  /* check for matching active data size */
  for (i=0; i<MAX_DIMS; i++) {
    N = NV_ACTIVELEN_PG(x,i);
    M = NV_ACTIVELEN_PG(y,i);
    if (N != M) {
      fprintf(stderr,"VCheck_Compatible: x,y active size mismatch (dim %li: %li vs %li)\n",i,N,M);
      return FALSE;
    }
  }

  /* check for matching offsets */
  for (i=0; i<MAX_DIMS; i++) {
    N = NV_OFFSET_PG(x,i);
    M = NV_OFFSET_PG(y,i);
    if (N != M) {
      fprintf(stderr,"VCheck_Compatible: x,y offset mismatch (dim %li: %li vs %li)\n",i,N,M);
      return FALSE;
    }
  }

  /* check for matching order */
  if (NV_FORDER_PG(x) != NV_FORDER_PG(y)) {
    fprintf(stderr,"VCheck_Compatible: x,y ordering mismatch\n");
    return FALSE;
  }

  /* if we made it here, the arrays are compatible */
  return TRUE;
}

static long int NV_DATALEN_PG(N_Vector x) {

  /* simple routine to output the local vector data length */
  long int N, i;
  N = 1;
  for (i=0; i<MAX_DIMS; i++)  N *= NV_ARRAYLEN_PG(x,i);

  return N;
}

static realtype VAllReduce_Parallel_Grid(realtype d, int op, MPI_Comm comm)
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
   case 1: MPI_Allreduce(&d, &out, 1, PGVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PGVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PGVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}

static void VCopy_Parallel_Grid(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  xd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = xd[i]; 

  return;
}

static void VSum_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = xd[i]-yd[i];

  return;
}

static void VNeg_Parallel_Grid(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;
  xd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = -xd[i];

  return;
}

static void VScaleSum_Parallel_Grid(realtype c, N_Vector x, 
				    N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)   zd[i] = c*(xd[i]+yd[i]);

  return;
}

static void VScaleDiff_Parallel_Grid(realtype c, N_Vector x, 
				     N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_Parallel_Grid(realtype a, N_Vector x, 
				N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = (a*xd[i])+yd[i];

  return;
}

static void VLin2_Parallel_Grid(realtype a, N_Vector x, 
				N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);
  zd = NV_DATA_PG(z);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  zd[i] = (a*xd[i])-yd[i];

  return;
}

static void Vaxpy_Parallel_Grid(realtype a, N_Vector x, N_Vector y)
{
  long int i, N;
  realtype *xd, *yd;
  xd = yd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);
  yd = NV_DATA_PG(y);

  /* perform operation on all local data */ 
  if (a == ONE) {
    for (i=0; i<N; i++)  yd[i] += xd[i];
    return;
  }
  
  if (a == -ONE) {
    for (i=0; i<N; i++)  yd[i] -= xd[i];
    return;
  }    
  
  for (i=0; i<N; i++)  yd[i] += a*xd[i];

  return;
}

static void VScaleBy_Parallel_Grid(realtype a, N_Vector x)
{
  long int i, N;
  realtype *xd;
  xd = NULL;

  /* access data arrays */
  N  = NV_DATALEN_PG(x);
  xd = NV_DATA_PG(x);

  /* perform operation on all local data */ 
  for (i=0; i<N; i++)  xd[i] *= a;

  return;
}

/* ----------------------------------------------------------------- */
