/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Header file for the generic ARKSLS linear solver module.
 *---------------------------------------------------------------*/

#ifndef _ARKSPARSE_H
#define _ARKSPARSE_H

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  ARKSPARSE CONSTANTS
===============================================================*/

/* ARKSLS return values */
#define ARKSLS_SUCCESS           0
#define ARKSLS_MEM_NULL         -1
#define ARKSLS_LMEM_NULL        -2
#define ARKSLS_ILL_INPUT        -3
#define ARKSLS_MEM_FAIL         -4
#define ARKSLS_JAC_NOSET        -5
#define ARKSLS_MASS_NOSET       -6
#define ARKSLS_PACKAGE_FAIL     -7
#define ARKSLS_MASSMEM_NULL     -8

/* Additional last_flag values */
#define ARKSLS_JACFUNC_UNRECVR  -9
#define ARKSLS_JACFUNC_RECVR    -10
#define ARKSLS_MASSFUNC_UNRECVR -11
#define ARKSLS_MASSFUNC_RECVR   -12


/*===============================================================
  FUNCTION TYPES
===============================================================*/

/*-----------------------------------------------------------------
 Types : ARKSlsSparseJacFn

 A sparse Jacobian approximation function Jac must be of type
 ARKSlsSparseJacFn.  Its parameters are:

 t   is the current value of the independent variable t.
                                                                
 y   is the current value of the dependent variable vector,
     namely the predicted value of y(t)
                                                                
 fy  is the vector f(t,y).
     namely the predicted value of y'(t)
                                                                
 JacMat is the compressed sparse column matrix (of type SlsMat)
     to be loaded by an ARKSlsSparseJacFn routine with an 
     approximation to the system Jacobian matrix
            J = (df_i/dy_j) at the point (t,y). 
     Note that JacMat is NOT preset to zero!
     Matrix data is for the nonzero entries of the Jacobian stored 
     in compressed column format.  Row indices of entries in 
     column j are stored in JacMat->rowvals[colptrs[j]] 
     through JacMat->rowvals[colptrs[j+i]-1]
     and corresponding numerical values of the Jacobian are stored 
     in the same entries of JacMat->data.

 J_data is a pointer to user Jacobian data - the same as the
     user_data parameter passed to ARKodeSetFdata.
                                                                
 tmp1, tmp2, tmp3 are pointers to memory allocated for
     N_Vectors which can be used by an ARKSparseJacFn routine
     as temporary storage or work space.
                                                                
 A ARKSlsSparseJacFn should return 0 if successful, a positive 
 value if a recoverable error occurred, and a negative value if 
 a nonrecoverable error occurred. 

 NOTE: If the user's Jacobian routine needs other quantities,
     they are accessible as follows: hcur (the current stepsize)
     and ewt (the error weight vector) are accessible through
     ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively
     (see arkode.h). The unit roundoff is available as
     UNIT_ROUNDOFF defined in sundials_types.h.

---------------------------------------------------------------*/
typedef int (*ARKSlsSparseJacFn)(realtype t,
				 N_Vector y, N_Vector fy, 
				 SlsMat JacMat, void *user_data,
				 N_Vector tmp1, N_Vector tmp2, 
				 N_Vector tmp3);

/*-----------------------------------------------------------------
 Types : ARKSlsSparseMassFn

 A sparse mass matrix approximation function Mass must be of type
 ARKSlsSparseMassFn.  Its parameters are:

 t   is the current value of the independent variable t.
                                                                
 MassMat is the compressed sparse column matrix (of type SlsMat)
     to be loaded by an ARKSlsSparseJacFn routine with an 
     approximation to the system mass matrix.  Note that MassMat 
     is NOT preset to zero! Matrix data is for the nonzero entries 
     of the mass matrix stored in compressed column format.  Row 
     indices of entries in column j are stored in 
     MassMat->rowvals[colptrs[j]] through 
     MassMat->rowvals[colptrs[j+i]-1]
     and corresponding numerical values of the mass matrix are 
     stored in the same entries of MassMat->data.

 M_data is a pointer to user mass matrix data - the same as the
     user_data parameter passed to ARKodeSetFdata.
                                                                
 tmp1, tmp2, tmp3 are pointers to memory allocated for
     N_Vectors which can be used by an ARKSparseJacFn routine
     as temporary storage or work space.
                                                                
 A ARKSlsSparseMassFn should return 0 if successful, a positive 
 value if a recoverable error occurred, and a negative int if a
 nonrecoverable error occurred. 

---------------------------------------------------------------*/
typedef int (*ARKSlsSparseMassFn)(realtype t, SlsMat MassMat, 
				  void *user_data, N_Vector tmp1, 
				  N_Vector tmp2, N_Vector tmp3);


/*===============================================================
  EXPORTED FUNCTIONS
===============================================================*/

/*---------------------------------------------------------------
 Optional inputs to the ARKSLS linear solver:

 ARKSlsSetSparseJacFn specifies the Jacobian approximation
 routine to be used for a sparse direct linear solver.

 The return value is one of:
    ARKSLS_SUCCESS   if successful
    ARKSLS_MEM_NULL  if the ARKODE memory was NULL
    ARKSLS_LMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSlsSetSparseJacFn(void *arkode_mem, 
					 ARKSlsSparseJacFn jac);

/*---------------------------------------------------------------
 Optional inputs to the ARKSLS linear solver:

 ARKSlsSetSparseMassFn specifies the sparse mass matrix 
 approximation routine to be used for a sparse direct linear 
 solver.

 The return value is one of:
    ARKSLS_SUCCESS   if successful
    ARKSLS_MEM_NULL  if the ARKODE memory was NULL
    ARKSLS_MASSMEM_NULL if the mass matrix solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSlsSetSparseMassFn(void *arkode_mem, 
					  ARKSlsSparseMassFn mass);

/*---------------------------------------------------------------
 Optional outputs from the ARKSLS linear solver:

 ARKSlsGetNumJacEvals returns the number of calls made to the
                      Jacobian evaluation routine jac.
 ARKSlsGetNumMassEvals returns the number of calls made to the
                      mass matrix evaluation routine Mass.
 ARKSlsGetLastFlag    returns the last error flag set by any of
                      the ARKSLS interface functions.
 ARKSlsGetLastMassFlag returns the last error flag set by any of 
                      the ARKSLS interface mass matrix functions.

 The return value of ARKSlsGet* is one of:
    ARKSLS_SUCCESS   if successful
    ARKSLS_MEM_NULL  if the IDA memory was NULL
    ARKSLS_LMEM_NULL if the linear solver memory was NULL
    ARKSLS_MASSMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSlsGetNumJacEvals(void *arkode_mem, 
					 long int *njevals);
SUNDIALS_EXPORT int ARKSlsGetNumMassEvals(void *arkode_mem, 
					  long int *nmevals);
SUNDIALS_EXPORT int ARKSlsGetLastFlag(void *arkode_mem, 
				      long int *flag);
SUNDIALS_EXPORT int ARKSlsGetLastMassFlag(void *arkode_mem, 
					  long int *flag);

/*---------------------------------------------------------------
 The following function returns the name of the constant 
 associated with a ARKSLS return flag
---------------------------------------------------------------*/
SUNDIALS_EXPORT char *ARKSlsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
