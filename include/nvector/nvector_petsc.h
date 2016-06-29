/*
 * -----------------------------------------------------------------
 * $Revision:$
 * $Date:$
 * ----------------------------------------------------------------- 
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the main header file for the PETSc vector wrapper 
 * for NVECTOR module.
 *
 * Part I contains declarations specific to the PETSc vector wrapper
 * implementation.
 *
 * Part II contains the prototype for the constructor
 * N_VMake_petsc as well as PETSc-specific prototypes
 * for various useful vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type realtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     build configuration stage) according to the user's needs. 
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

#ifndef _NVECTOR_PETSC_H
#define _NVECTOR_PETSC_H

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
   to the underlying PETSc vector, the MPI communicator,
   and a flag indicating ownership of the PETSc vector */

struct _N_VectorContent_petsc {
  long int local_length;   /* copy of local vector length  */
  long int global_length;  /* copy of global vector length */
  booleantype own_data;    /* ownership of data            */
  Vec *pvec;               /* pointer to PETSc vector      */
  MPI_Comm comm;           /* copy of MPI communicator     */
};

typedef struct _N_VectorContent_petsc *N_VectorContent_petsc;

/*
 * -----------------------------------------------------------------
 * PART II: functions exported by nvector_petsc
 * 
 * CONSTRUCTORS:
 *    N_VNewEmpty_petsc
 *    N_VMake_petsc
 *    N_VCloneVectorArray_petsc
 *    N_VCloneVectorArrayEmpty_petsc
 * DESTRUCTORS:
 *    N_VDestroyVectorArray_petsc
 * OTHER:
 *    N_VGetVector_petsc
 *    N_VPrint_petsc
 * -----------------------------------------------------------------
 */


/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_petsc
 * -----------------------------------------------------------------
 * This function creates a new N_Vector wrapper around an empty
 * (NULL) PETSc vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_petsc(MPI_Comm comm, 
                                           long int local_length,
                                           long int global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_petsc
 * -----------------------------------------------------------------
 * This function is not supported for PETSc vector wrapper.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_petsc(Vec *v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetArrayPointer_petsc
 * -----------------------------------------------------------------
 * This function is not supported for PETSc vector wrapper.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype *N_VGetArrayPointer_petsc(N_Vector v);

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
 * Function : N_VGetVector_petsc
 * -----------------------------------------------------------------
 * Extracts PETSc vector from N_Vector wrapper.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT Vec *N_VGetVector_petsc(N_Vector v);

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
 * PETSc implementations of the vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_petsc(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_petsc(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_petsc(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_petsc(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_petsc(N_Vector v, long int *lrw, long int *liw);
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

#endif /* _NVECTOR_PETSC_H */
