/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the main header file for the MPI-enabled implementation
 * of the NVECTOR module.
 *
 * Part I contains declarations specific to the parallel
 * implementation of the supplied NVECTOR module.
 *
 * Part II defines accessor macros that allow the user to efficiently
 * use the type N_Vector without making explicit references to the
 * underlying data structure.
 *
 * Part III contains the prototype for the constructor
 * N_VNew_petsc as well as implementation-specific prototypes
 * for various useful vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type realtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type booleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_petsc(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------
 */

#ifndef _NVECTOR_PARALLEL_H
#define _NVECTOR_PARALLEL_H

#include <mpi.h>
#include <petscvec.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*
 * -----------------------------------------------------------------
 * PART I: PARALLEL implementation of N_Vector               
 * -----------------------------------------------------------------
 */

/* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)

#define PVEC_REAL_MPI_TYPE MPI_FLOAT

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

#elif defined(SUNDIALS_EXTENDED_PRECISION)

#define PVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE

#endif

#define PVEC_INTEGER_MPI_TYPE MPI_LONG

/* parallel implementation of the N_Vector 'content' structure
   contains the global and local lengths of the vector, a pointer
   to an array of 'realtype components', the MPI communicator,
   and a flag indicating ownership of the data */

struct _N_VectorContent_petsc {
  long int local_length;   /* local vector length         */
  long int global_length;  /* global vector length        */
  booleantype own_data;    /* ownership of data           */
  Vec *pvec;               /* PETSc vector                */
  MPI_Comm comm;           /* pointer to MPI communicator */
};

typedef struct _N_VectorContent_petsc *N_VectorContent_petsc;

/*
 * -----------------------------------------------------------------
 * PART II: macros NV_CONTENT_PTC, NV_DATA_PTC, NV_OWN_DATA_PTC,
 *          NV_LOCLENGTH_PTC, NV_GLOBLENGTH_PTC,NV_COMM_PTC, and NV_Ith_PTC
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * long int v_len, s_len, i;
 *
 * (1) NV_CONTENT_PTC
 *
 *     This routines gives access to the contents of the parallel
 *     vector N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_PTC(v) sets v_cont to be
 *     a pointer to the parallel N_Vector content structure.
 *
 * (2) NV_DATA_PTC, NV_OWN_DATA_PTC, NV_LOCLENGTH_PTC, NV_GLOBLENGTH_PTC,
 *     and NV_COMM_PTC
 *
 *     These routines give access to the individual parts of
 *     the content structure of a parallel N_Vector.
 *
 *     The assignment v_data = NV_DATA_PTC(v) sets v_data to be
 *     a pointer to the first component of the local data for
 *     the vector v. The assignment NV_DATA_PTC(v) = data_v sets
 *     the component array of v to be data_V by storing the
 *     pointer data_v.
 *
 *     The assignment v_llen = NV_LOCLENGTH_PTC(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_LOCLENGTH_PTC(v) = llen_v sets the local length
 *     of v to be llen_v.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PTC(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PTC(v) = glen_v sets the global length of v to
 *     be glen_v.
 *
 *     The assignment v_comm = NV_COMM_PTC(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_C(v) = comm_v sets the MPI communicator of v to be
 *     comm_v.
 *
 * (3) NV_Ith_PTC
 *
 *     In the following description, the components of the
 *     local part of an N_Vector are numbered 0..n-1, where n
 *     is the local length of (the local part of) v.
 *
 *     The assignment r = NV_Ith_PTC(v,i) sets r to be the value
 *     of the ith component of the local part of the vector v.
 *     The assignment NV_Ith_PTC(v,i) = r sets the value of the
 *     ith local component of v to be r.
 *
 * Note: When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_PTC(v) and then access v_data[i] within the
 * loop than it is to use NV_Ith_PTC(v,i) within the loop.
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PTC(v)    ( (N_VectorContent_petsc)(v->content) )

#define NV_LOCLENGTH_PTC(v)  ( NV_CONTENT_PTC(v)->local_length )

#define NV_GLOBLENGTH_PTC(v) ( NV_CONTENT_PTC(v)->global_length )

#define NV_OWN_DATA_PTC(v)   ( NV_CONTENT_PTC(v)->own_data )

#define NV_PVEC_PTC(v)       ( NV_CONTENT_PTC(v)->pvec )

#define NV_COMM_PTC(v)       ( NV_CONTENT_PTC(v)->comm )

// This needs to be reworked 
//#define NV_Ith_PTC(v,i)      ( NV_DATA_PTC(v)[i] )
// Above needs to be reworked

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by nvector_parallel
 * 
 * CONSTRUCTORS:
 *    N_VNewEmpty_petsc
 *    N_VNew_petsc
 *    N_VMake_petsc
 *    N_VCloneVectorArray_petsc
 *    N_VCloneVectorArrayEmpty_petsc
 * DESTRUCTORS:
 *    N_VDestroy_petsc
 *    N_VDestroyVectorArray_petsc
 * OTHER:
 *    N_VPrint_petsc
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_petsc
 * -----------------------------------------------------------------
 * This function creates a new parallel N_Vector with an empty
 * (NULL) data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_petsc(MPI_Comm comm, 
                                           long int local_length,
                                           long int global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_petsc
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a parallel vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_petsc(MPI_Comm comm, 
                                      long int local_length,
                                      long int global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_petsc
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a parallel vector
 * with a user-supplied data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_petsc(MPI_Comm comm, 
                                       long int local_length,
                                       long int global_length,
                                       realtype *v_data);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_petsc
 * -----------------------------------------------------------------
 * This function creates an array of 'count' PARALLEL vectors by
 * cloning a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_petsc(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_petsc
 * -----------------------------------------------------------------
 * This function creates an array of 'count' PARALLEL vectors each 
 * with an empty (NULL) data array by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_petsc(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_petsc
 * -----------------------------------------------------------------
 * This function frees an array of N_Vector created with 
 * N_VCloneVectorArray_petsc or N_VCloneVectorArrayEmpty_petsc.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_petsc(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_petsc
 * -----------------------------------------------------------------
 * This function prints the content of a parallel vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_petsc(N_Vector v);

/*
 * -----------------------------------------------------------------
 * parallel implementations of the vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_petsc(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_petsc(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_petsc(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_petsc(N_Vector v, long int *lrw, long int *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_petsc(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_petsc(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_petsc(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_petsc(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_petsc(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_petsc(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_petsc(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_petsc(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_petsc(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_petsc(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_petsc(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_petsc(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_petsc(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_petsc(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_petsc(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
