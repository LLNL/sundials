/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2013/08/31 12:43:00 $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner and Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR 
 *                   Serial module by Scott D. Cohen, Alan C. 
 *                   Hindmarsh, Radu Serban, and Aaron Collier 
 *                   @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the openMP implementation of the
 * NVECTOR module.
 *
 * Part I contains declarations specific to the openMP
 * implementation of the supplied NVECTOR module.
 *
 * Part II defines accessor macros that allow the user to
 * efficiently use the type N_Vector without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor N_VNew_openMP
 * as well as implementation-specific prototypes for various useful
 * vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_openMP(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------
 */

#ifndef _NVECTOR_OPENMP_H
#define _NVECTOR_OPENMP_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_nvector.h>

/*
 * -----------------------------------------------------------------
 * PART I: OPENMP implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* openMP implementation of the N_Vector 'content' structure
   contains the length of the vector, a pointer to an array
   of 'realtype' components, and a flag indicating ownership of
   the data */

struct _N_VectorContent_openMP {
  long int length;
  booleantype own_data;
  realtype *data;
};

typedef struct _N_VectorContent_openMP *N_VectorContent_openMP;

/*
 * -----------------------------------------------------------------
 * PART II: macros NV_CONTENT_OMP, NV_DATA_OMP, NV_OWN_DATA_OMP,
 *          NV_LENGTH_OMP, and NV_Ith_OMP
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * long int i;
 *
 * (1) NV_CONTENT_OMP
 *
 *     This routines gives access to the contents of the openMP
 *     vector N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_OMP(v) sets v_cont to be
 *     a pointer to the openMP N_Vector content structure.
 *
 * (2) NV_DATA_OMP NV_OWN_DATA_OMP and NV_LENGTH_OMP
 *
 *     These routines give access to the individual parts of
 *     the content structure of a openMP N_Vector.
 *
 *     The assignment v_data = NV_DATA_OMP(v) sets v_data to be
 *     a pointer to the first component of v. The assignment
 *     NV_DATA_OMP(v) = data_V sets the component array of v to
 *     be data_v by storing the pointer data_v.
 *
 *     The assignment v_len = NV_LENGTH_OMP(v) sets v_len to be
 *     the length of v. The call NV_LENGTH_OMP(v) = len_v sets
 *     the length of v to be len_v.
 *
 * (3) NV_Ith_OMP
 *
 *     In the following description, the components of an
 *     N_Vector are numbered 0..n-1, where n is the length of v.
 *
 *     The assignment r = NV_Ith_OMP(v,i) sets r to be the value of
 *     the ith component of v. The assignment NV_Ith_OMP(v,i) = r
 *     sets the value of the ith component of v to be r.
 *
 * Note: When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_OMP(v) and then access v_data[i] within the
 * loop than it is to use NV_Ith_OMP(v,i) within the loop.
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_OMP(v)  ( (N_VectorContent_openMP)(v->content) )

#define NV_LENGTH_OMP(v)   ( NV_CONTENT_OMP(v)->length )

#define NV_OWN_DATA_OMP(v) ( NV_CONTENT_OMP(v)->own_data )

#define NV_DATA_OMP(v)     ( NV_CONTENT_OMP(v)->data )

#define NV_Ith_OMP(v,i)    ( NV_DATA_OMP(v)[i] )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by nvector_openMP
 * 
 * CONSTRUCTORS:
 *    N_VNew_openMP
 *    N_VNewEmpty_openMP
 *    N_VMake_openMP
 *    N_VCloneVectorArray_openMP
 *    N_VCloneVectorArrayEmpty_openMP
 * DESTRUCTORS:
 *    N_VDestroy_openMP
 *    N_VDestroyVectorArray_openMP
 * OTHER:
 *    N_VPrint_openMP
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_openMP
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a openMP vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_openMP(long int vec_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_openMP
 * -----------------------------------------------------------------
 * This function creates a new openMP N_Vector with an empty (NULL)
 * data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_openMP(long int vec_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_openMP
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a openMP vector
 * with a user-supplied data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_openMP(long int vec_length, realtype *v_data);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_openMP
 * -----------------------------------------------------------------
 * This function creates an array of 'count' OPENMP vectors by
 * cloning a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_openMP(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_openMP
 * -----------------------------------------------------------------
 * This function creates an array of 'count' OPENMP vectors each
 * with an empty (NULL) data array by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_openMP(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_openMP
 * -----------------------------------------------------------------
 * This function frees an array of OPENMP vectors created with 
 * N_VCloneVectorArray_openMP or N_VCloneVectorArrayEmpty_openMP.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_openMP(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_openMP
 * -----------------------------------------------------------------
 * This function prints the content of a openMP vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_openMP(N_Vector v);

/*
 * -----------------------------------------------------------------
 * openMP implementations of various useful vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_openMP(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_openMP(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_openMP(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_openMP(N_Vector v, long int *lrw, long int *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_openMP(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_openMP(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_openMP(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_openMP(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_openMP(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_openMP(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_openMP(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_openMP(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_openMP(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_openMP(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_openMP(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_openMP(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_openMP(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_openMP(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_openMP(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_openMP(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_openMP(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_openMP(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_openMP(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_openMP(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_openMP(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
