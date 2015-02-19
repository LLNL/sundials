/*----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 -----------------------------------------------------------------
 This is the main header file for the MPI-enabled implementation
 of the grid NVECTOR module.
 
  Part I contains declarations specific to the parallel grid
  implementation of the supplied NVECTOR module.
 
 Part II defines accessor macros that allow the user to efficiently
 use the type N_Vector without making explicit references to the
 underlying data structure.
 
 Part III contains the prototype for the constructor
 N_VNew_Parallel_Grid as well as implementation-specific prototypes
 for various useful vector operations.
 
 Notes:
 
   - The definition of the generic N_Vector structure can be
     found in the header file sundials_nvector.h.
 
   - The definition of the type realtype can be found in the
     header file sundials_types.h, and it may be changed (at the 
     configuration stage) according to the user's needs. 
     The sundials_types.h file also contains the definition
     for the type booleantype.
 
   - N_Vector arguments to arithmetic vector operations need not
     be distinct. For example, the following call:
 
        N_VLinearSum_Parallel_Grid(a,x,b,y,y);
 
     (which stores the result of the operation a*x+b*y in y)
     is legal.
 ---------------------------------------------------------------*/

#ifndef _NVECTOR_PARALLEL_GRID_H
#define _NVECTOR_PARALLEL_GRID_H

#include <mpi.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*
 * -----------------------------------------------------------------
 * PART I: PARALLEL grid implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* define maximum dimensionality of multi-dimensional arrays */
#define MAX_DIMS 6


/* define MPI data types and output formatting specifiers */

#if defined(SUNDIALS_SINGLE_PRECISION)

#define PGVEC_REAL_MPI_TYPE MPI_FLOAT
#define DOUT "g"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define PGVEC_REAL_MPI_TYPE MPI_DOUBLE
#define DOUT "lg"

#elif defined(SUNDIALS_EXTENDED_PRECISION)

#define PGVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE
#define DOUT "Lg"

#endif

#define PGVEC_INTEGER_MPI_TYPE MPI_LONG

/* parallel grid implementation of the N_Vector 'content' structure 
   contains the vector dimensionality, the total allocated array length 
   in each direction, the active array length in each direction, the
   offset to the beginning of the active data in each direction, the 
   global length of the vector (all MPI tasks), a flag indicating 
   ownership of the data, a pointer to a contiguous array of 'realtype' 
   components, the MPI communicator, and a flag indicating whether the 
   data is stored in Fortran (column-major) or C (row-major) ordering.

   Note: the total length of the allocated array will be the product of 
   dim_length[0] through dim_length[dims-1].  We assume that the 
   multidimensional data held in the actual "data" array adheres to the 
   choice of C or Fortran major ordering held in the 'ordering' flag. */

struct _N_VectorContent_Parallel_Grid {
  long int dims;                  /* vector dimensionality                 */
  long int dim_length[MAX_DIMS];  /* total array length in each dimension  */
  long int dim_alength[MAX_DIMS]; /* active array length in each dimension */
  long int dim_offset[MAX_DIMS];  /* offset to active data in each dim     */
  long int global_length;         /* global vector length                  */
  booleantype own_data;           /* ownership of data                     */
  realtype *data;                 /* local data array                      */
  MPI_Comm comm;                  /* pointer to MPI communicator           */
  booleantype F_ordering;         /* Fortran or C ordering                 */
};

typedef struct _N_VectorContent_Parallel_Grid *N_VectorContent_Parallel_Grid;

/*
 * -----------------------------------------------------------------
 * PART II: macros NV_CONTENT_PG, NV_DATA_PG, NV_OWN_DATA_PG,
 *          NV_ARRAYLEN_PG, NV_ACTIVELEN_PG, NV_OFFSET_PG, 
 *          NV_GLOBLENGTH_PG, NV_COMM_PG, NV_FORDER_PG
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * long int i[] = {i0, i1, ..., id};  [has length dims]
 *
 * (1) NV_CONTENT_PG
 *
 *     This routines gives access to the contents of the parallel
 *     grid vector N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_P(v) sets v_cont to be
 *     a pointer to the parallel grid N_Vector content structure.
 *
 * (2) NV_DATA_PG, NV_DIMS_PG, NV_ORDER_PG, NV_OWN_DATA_PG, NV_ARRAYLEN_PG, 
 *     NV_ACTIVELEN_PG, NV_OFFSET_PG, NV_GLOBLENGTH_PG and NV_COMM_PG
 *
 *     These routines give access to the individual parts of
 *     the content structure of a parallel grid N_Vector.
 *
 *     The assignment v_data = NV_DATA_PG(v) sets v_data to be
 *     a pointer to the first component of the local data for
 *     the vector v. The assignment NV_DATA_PG(v) = data_v sets
 *     the component array of v to be data_V by storing the
 *     pointer data_v.
 *
 *     The assignment v_dims = NV_DIMS_PG(v) sets v_dims to be
 *     the dimensionality of the vector v. The assignment 
 *     NV_DIMS_PG(v) = dims_v sets the dimensionality of v to be dims_v.
 *
 *     The assignment v_order = NV_FORDER_PG(v) sets v_order to be 
 *     either 0 or 1, where 0 corresponds to C ordering and 1 corresponds 
 *     to Fortran ordering of the multi-dimensional array.  The assignment
 *     NV_FORDER_PG(v) = order_v sets the ordering of v to be order_v.
 *
 *     The assignment v_owndata = NV_OWN_DATA_PG(v) sets v_owndata to be 
 *     either 0 or 1, where 0 signifies that v does not own the data.  
 *     The assignment NV_OWN_DATA_PG(v) = owndata_v sets the ownership of 
 *     v's data to be owndata_v.
 *
 *     The assignment v_llen = NV_ARRAYLEN_PG(v,d) sets v_llen to
 *     be the length of the local part of the vector v in the dimension 
 *     d. The call NV_ARRAYLEN_PG(v,d) = llen_v sets the local length
 *     of v in the dimension d to be llen_v.
 *
 *     The assignment v_alen = NV_ACTIVELEN_PG(v,d) sets v_alen to
 *     be the length of the active portion of the local part of the vector 
 *     v in the dimension d. The call NV_ACTIVELEN_PG(v,d) = alen_v sets 
 *     the local length of v in the dimension d to be alen_v.
 *
 *     The assignment v_off = NV_OFFSET_PG(v,d) sets v_off to
 *     be the length of the offset to the active portion of the local
 *     part of the vector v in the dimension d. The call 
 *     NV_OFFSET_PG(v,d) = off_v sets the offset of v in the dimension d 
 *     to be off_v.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PG(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PG(v) = glen_v sets the global length of v to
 *     be glen_v.
 *
 *     The assignment v_comm = NV_COMM_PG(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_PG(v) = comm_v sets the MPI communicator of v to be
 *     comm_v.
 *
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PG(v)     ( (N_VectorContent_Parallel_Grid)(v->content) )

#define NV_DATA_PG(v)        ( NV_CONTENT_PG(v)->data )

#define NV_DIMS_PG(v)        ( NV_CONTENT_PG(v)->dims )

#define NV_FORDER_PG(v)      ( NV_CONTENT_PG(v)->F_ordering )

#define NV_OWN_DATA_PG(v)    ( NV_CONTENT_PG(v)->own_data )

#define NV_ARRAYLEN_PG(v,d)  ( (NV_CONTENT_PG(v)->dim_length)[d] )

#define NV_ACTIVELEN_PG(v,d) ( (NV_CONTENT_PG(v)->dim_alength)[d] )

#define NV_OFFSET_PG(v,d)    ( (NV_CONTENT_PG(v)->dim_offset)[d] )

#define NV_GLOBLENGTH_PG(v)  ( NV_CONTENT_PG(v)->global_length )

#define NV_COMM_PG(v)        ( NV_CONTENT_PG(v)->comm )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by nvector_parallel_grid
 * 
 * CONSTRUCTORS:
 *    N_VNew_Parallel_Grid
 *    N_VNewEmpty_Parallel_Grid
 *    N_VMake_Parallel_Grid
 *    N_VCloneVectorArray_Parallel_Grid
 *    N_VCloneVectorArrayEmpty_Parallel_Grid
 * DESTRUCTORS:
 *    N_VDestroy_Parallel_Grid
 *    N_VDestroyVectorArray_Parallel_Grid
 * OTHER:
 *    N_VPrint_Parallel_Grid
 * OTHER:
 *    N_VPrintAll_Parallel_Grid
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_Parallel_Grid
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a parallel 
 * grid vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_Parallel_Grid(MPI_Comm comm, 
					      long int dims,
					      long int *dim_length,
					      long int *dim_alength,
					      long int *dim_offset,
					      long int F_ordering,
					      long int global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_Parallel_Grid
 * -----------------------------------------------------------------
 * This function creates a new parallel grid N_Vector with an 
 * empty (NULL) data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Parallel_Grid(MPI_Comm comm, 
						   long int dims,
						   long int *dim_length,
						   long int *dim_alength,
						   long int *dim_offset,
						   long int F_ordering,
						   long int global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_Parallel_Grid
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a parallel grid 
 * vector with a user-supplied data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_Parallel_Grid(MPI_Comm comm, 
					       long int dims,
					       long int *dim_length,
					       long int *dim_alength,
					       long int *dim_offset,
					       long int F_ordering,
					       long int global_length,
					       realtype *v_data);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_Parallel_Grid
 * -----------------------------------------------------------------
 * This function creates an array of 'count' PARALLEL grid 
 * vectors by cloning a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_Parallel_Grid(int count, 
							    N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_Parallel_Grid
 * -----------------------------------------------------------------
 * This function creates an array of 'count' PARALLEL grid 
 * vectors each with an empty (NULL) data array by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_Parallel_Grid(int count, 
								 N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_Parallel_Grid
 * -----------------------------------------------------------------
 * This function frees an array of N_Vector created with 
 * N_VCloneVectorArray_Parallel_Grid or 
 * N_VCloneVectorArrayEmpty_Parallel_Grid.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_Parallel_Grid(N_Vector *vs, 
							 int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_Parallel_Grid
 * -----------------------------------------------------------------
 * This function prints the content of a parallel grid vector to 
 * stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_Parallel_Grid(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintAll_Parallel_Grid
 * -----------------------------------------------------------------
 * This function prints the content of a parallel grid vector to 
 * stdout, including all grid zone values.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrintAll_Parallel_Grid(N_Vector v);

/*
 * -----------------------------------------------------------------
 * grid parallel implementations of the vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Parallel_Grid(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Parallel_Grid(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Parallel_Grid(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Parallel_Grid(N_Vector v, long int *lrw, long int *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_Parallel_Grid(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Parallel_Grid(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_Parallel_Grid(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Parallel_Grid(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Parallel_Grid(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Parallel_Grid(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Parallel_Grid(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Parallel_Grid(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Parallel_Grid(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Parallel_Grid(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Parallel_Grid(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Parallel_Grid(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Parallel_Grid(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Parallel_Grid(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Parallel_Grid(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Parallel_Grid(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Parallel_Grid(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Parallel_Grid(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Parallel_Grid(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Parallel_Grid(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
