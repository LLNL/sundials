/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ SMU
 *                Slaven Peles @ LLNL
 * ----------------------------------------------------------------- 
 * Based on work by: Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                   and Aaron Collier @ LLNL
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
 * N_VNew_ParHyp as well as implementation-specific prototypes
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
 *        N_VLinearSum_ParHyp(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------
 */

#ifndef _NVECTOR_PARHYP_H
#define _NVECTOR_PARHYP_H

#include <mpi.h>
#include <sundials/sundials_nvector.h>
/* hypre header files */
#include <seq_mv.h>
#include <_hypre_parcsr_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <_hypre_parcsr_ls.h>

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

/* 
 * Parallel implementation of the N_Vector 'content' structure
 * contains the global and local lengths of the vector, a pointer
 * to an array of 'realtype components', the MPI communicator,
 * and a flag indicating ownership of the data. 
 */
struct _N_VectorContent_ParHyp {
  long int local_length;         /* local vector length         */
  long int global_length;        /* global vector length        */
  booleantype own_data;          /* ownership of data           */
  booleantype own_parvector;     /* ownership of vector pointer */
  MPI_Comm comm;                 /* pointer to MPI communicator */

  hypre_ParVector *x; /* The actual hypre_ParVector object the local data and    */
                      /* other relevant info (such as partitioning) is stored in */
};

typedef struct _N_VectorContent_ParHyp *N_VectorContent_ParHyp;

/*
 * -----------------------------------------------------------------
 * PART II: macros NV_CONTENT_PH, NV_DATA_PH, NV_OWN_DATA_PH,
 *          NV_LOCLENGTH_PH, NV_GLOBLENGTH_PH,NV_COMM_PH, and 
 *          NV_Ith_PH
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * long int v_len, s_len, i;
 *
 * (1) NV_CONTENT_PH
 *
 *     This routines gives access to the contents of the hypre
 *     vector wrapper (the N_Vector).
 *
 *     The assignment v_cont = NV_CONTENT_PH(v) sets v_cont to be
 *     a pointer to the N_Vector content structure.
 *
 * (2) NV_DATA_PH, NV_OWN_DATA_PH, NV_LOCLENGTH_PH, NV_GLOBLENGTH_PH,
 *     and NV_COMM_PH
 *
 *     These routines give access to the individual parts of
 *     the content structure of a parhyp N_Vector.
 *
 *     The assignment v_llen = NV_LOCLENGTH_PH(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_LOCLENGTH_PH(v) = llen_v generally should NOT be used! It 
 *     will change locally stored value with the hypre local vector 
 *     length, but it will NOT change the length of the actual hypre
 *     local vector.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PH(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PH(v) = glen_v generally should NOT be used! It 
 *     will change locally stored value with the hypre parallel vector 
 *     length, but it will NOT change the length of the actual hypre
 *     parallel vector.
 *
 *     The assignment v_comm = NV_COMM_PH(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_C(v) = comm_v sets the MPI communicator of v to be
 *     NV_COMM_PH(v) = comm_v generally should NOT be used! It 
 *     will change locally stored value with the hypre parallel vector 
 *     communicator, but it will NOT change the communicator of the 
 *     actual hypre parallel vector.
 * 
 * (3) NV_DATA_PH, NV_HYPRE_PARVEC_PH
 *     
 *     The assignment v_data = NV_DATA_PH(v) sets v_data to be
 *     a pointer to the first component of the data inside the 
 *     local vector of the hypre_parhyp vector for
 *     the vector v. The assignment NV_DATA_PH(v) = data_v sets
 *     the component array of the local vector of the hypre_parhyp
 *     vector for the vector v to be data_v by storing the
 *     pointer data_v. This is dangerous operation. Use only if
 *     you know what you are doing! Use NV_HYPRE_PARVEC_PH instead.
 *
 *     The assignment v_parhyp = NV_HYPRE_PARVEC_PH(v) sets v_parhyp
 *     to be a pointer to hypre_ParVector of vector v. The assignment
 *     NV_HYPRE_PARVEC_PH(v) = parhyp_v sets pointer to 
 *     hypre_ParVector of vector v to be parhyp_v.
 *
 * (4) NV_Ith_PH
 *
 *     In the following description, the components of the
 *     local part of an N_Vector are numbered 0..n-1, where n
 *     is the local length of (the local part of) v.
 *
 *     The assignment r = NV_Ith_PH(v,i) sets r to be the value
 *     of the ith component of the local part of the vector v.
 *     The assignment NV_Ith_PH(v,i) = r sets the value of the
 *     ith local component of v to be r.
 *
 * Note: When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_PH(v) and then access v_data[i] within the
 * loop than it is to use NV_Ith_PH(v,i) within the loop.
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PH(v)    ( (N_VectorContent_ParHyp)(v->content) )

#define NV_LOCLENGTH_PH(v)  ( NV_CONTENT_PH(v)->local_length )

#define NV_GLOBLENGTH_PH(v) ( NV_CONTENT_PH(v)->global_length )

#define NV_OWN_DATA_PH(v)   ( NV_CONTENT_PH(v)->own_data )

#define NV_OWN_PARVEC_PH(v) ( NV_CONTENT_PH(v)->own_parvector )

#define NV_HYPRE_PARVEC_PH(v) ( NV_CONTENT_PH(v)->x )

#define NV_DATA_PH(v)       ( NV_HYPRE_PARVEC_PH(v) == NULL ? NULL : hypre_VectorData(hypre_ParVectorLocalVector(NV_HYPRE_PARVEC_PH(v))) )

#define NV_COMM_PH(v)       ( NV_CONTENT_PH(v)->comm )

#define NV_Ith_PH(v,i)      ( NV_DATA_PH(v)[i] )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by nvector_ParHyp
 * 
 * CONSTRUCTORS:
 *    N_VNewEmpty_ParHyp
 *    N_VMake_ParHyp
 *    N_VCloneVectorArray_ParHyp
 *    N_VCloneVectorArrayEmpty_ParHyp
 * DESTRUCTORS:
 *    N_VDestroy_ParHyp
 *    N_VDestroyVectorArray_ParHyp
 * OTHER:
 *    N_VPrint_ParHyp
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_ParHyp
 * -----------------------------------------------------------------
 * This function creates a new hypre vector wrapper without the 
 * hypre vector itself.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm, 
                                            long int local_length,
                                            long int global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_ParHyp
 * -----------------------------------------------------------------
 * This function creates a hypre vector wrapper around user-supplied 
 * hypre vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_ParHyp(hypre_ParVector *x);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_ParHyp
 * -----------------------------------------------------------------
 * This function creates an array of 'count' N_Vectors by cloning a 
 * given vector w. Both, the wrapper and the underlying hypre vector
 * are cloned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_ParHyp(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_ParHyp
 * -----------------------------------------------------------------
 * This function creates an array of 'count' empty hypre vector
 * wrappers by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_ParHyp(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_ParHyp
 * -----------------------------------------------------------------
 * This function frees an array of N_Vector created with 
 * N_VCloneVectorArray_ParHyp or N_VCloneVectorArrayEmpty_ParHyp.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_ParHyp(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_ParHyp
 * -----------------------------------------------------------------
 * This function prints the content of a parallel vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_ParHyp(N_Vector v);

/*
 * -----------------------------------------------------------------
 * parallel implementations of the vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_ParHyp(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_ParHyp(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_ParHyp(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_ParHyp(N_Vector v, long int *lrw, long int *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_ParHyp(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_ParHyp(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_ParHyp(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_ParHyp(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_ParHyp(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_ParHyp(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_ParHyp(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_ParHyp(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_ParHyp(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_ParHyp(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_ParHyp(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
