/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * Common header file for the direct linear solver interface in 
 * ARKODE.
 *--------------------------------------------------------------*/

#ifndef _ARKDLS_H
#define _ARKDLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  ARKDLS Constants
===============================================================*/

/* ARKDLS return values */
#define ARKDLS_SUCCESS          0
#define ARKDLS_MEM_NULL        -1
#define ARKDLS_LMEM_NULL       -2
#define ARKDLS_ILL_INPUT       -3
#define ARKDLS_MEM_FAIL        -4
#define ARKDLS_MASSMEM_NULL    -5

/* Additional last_flag values */
#define ARKDLS_JACFUNC_UNRECVR  -6
#define ARKDLS_JACFUNC_RECVR    -7
#define ARKDLS_MASSFUNC_UNRECVR -8
#define ARKDLS_MASSFUNC_RECVR   -9

#define ARKDLS_SUNMAT_FAIL      -10


  
/*===============================================================
 ARKDLS user-supplied function prototypes
===============================================================*/

/*---------------------------------------------------------------
 Type: ARKDlsJacFn

 A Jacobian approximation function Jac must be of type 
 ARKDlsJacFn. Its parameters are:

 Jac is the SUNMatrix that will be loaded by a ARKDlsJacFn with
     an approximation to the Jacobian matrix 
     J = (df_i/dy_j) at the point (t,y). 

 t   is the current value of the independent variable.

 y   is the current value of the dependent variable vector,
     namely the predicted value of y(t).

 fy  is the vector f(t,y).

 user_data is a pointer to user data - the same as the user_data
     parameter passed to ARKodeSetUserData.

 tmp1, tmp2, and tmp3 are pointers to memory allocated for
     vectors of length N which can be used by a ARKDlsJacFn
     as temporary storage or work space.

 A ARKDlsJacFn should return 0 if successful, a positive 
 value if a recoverable error occurred, and a negative value if 
 an unrecoverable error occurred.

 NOTE: See the relevant SUNMatrix implementation header files
     and documentation for mechanisms to inquire about matrix 
     dimensions, and for efficient ways to set matrix entries.
                                                                
 NOTE: If the user's Jacobian routine needs other quantities,   
     they are accessible as follows: hcur (the current stepsize)
     and ewt (the error weight vector) are accessible through   
     ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively 
     (see arkode.h). The unit roundoff is available as 
     UNIT_ROUNDOFF defined in sundials_types.h.

---------------------------------------------------------------*/
typedef int (*ARKDlsJacFn)(realtype t, N_Vector y, N_Vector fy, 
                           SUNMatrix Jac, void *user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  

/*---------------------------------------------------------------
 Type: ARKDlsMassFn

 A mass matrix approximation function Mass must be of type 
 ARKDlsMassFn. Its parameters are:

 t   is the current value of the independent variable.

 M   is the SUNMatrix that will be loaded by a ARKDlsMassFn 
     with an approximation to the mass matrix.

 user_data is a pointer to user data - the same as the user_data
     parameter passed to ARKodeSetUserData.

 tmp1, tmp2, and tmp3 are pointers to memory allocated for
 vectors of length N which can be used by a ARKDlsMassFn
 as temporary storage or work space.

 A ARKDlsMassFn should return 0 if successful, and a 
 negative value if an unrecoverable error occurred.

---------------------------------------------------------------*/
typedef int (*ARKDlsMassFn)(realtype t, SUNMatrix M, void *user_data, 
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  

/*===============================================================
  ARKDLS Exported functions
===============================================================*/

/*---------------------------------------------------------------
 Required inputs for the ARKDLS linear solver interface:

 ARKDlsSetLinearSolver specifies the direct SUNLinearSolver 
 object that ARKode should use.  This is required if ARKode is 
 solving an implicit or IMEX IVP, and using a modified Newton 
 solver.

 ARKDlsSetMassLinearSolver specifies the direct SUNLinearSolver 
 object (and user-provided function to fill the mass matrix) that
 ARKode should use when solving mass-matrix linear systems.  This 
 is required if ARKode is solving a problem with non-identity 
 mass matrix and the user wishes to use a direct solver for these 
 systems.

 NOTE: when solving an implicit or IMEX IVP with non-identity mass
 matrix and direct linear solver, both the system and mass matrices 
 must have the same type (i.e. you cannot combine a direct system 
 solver with an iterative mass matrix solver, or a dense and banded 
 matrices, etc.).

 The return value is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the ARKODE memory was NULL
    ARKDLS_ILL_INPUT if the arguments are incompatible
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDlsSetLinearSolver(void *arkode_mem, 
                                          SUNLinearSolver LS,
                                          SUNMatrix A);

SUNDIALS_EXPORT int ARKDlsSetMassLinearSolver(void *arkode_mem, 
                                              SUNLinearSolver LS,
                                              SUNMatrix M,
                                              booleantype time_dep);

/*---------------------------------------------------------------
 Optional inputs to the ARKDLS linear solver interface:

 ARKDlsSetJacFn specifies the dense/band/sparse Jacobian
 approximation routine to be used for a direct linear solver.

 By default, a difference quotient approximation is used for 
 dense/band; no default exists for sparse (so this must be 
 user-supplied).

 ARKDlsSetMassFn specifies the mass matrix approximation routine 
 to be used for a direct linear solver.  As this must have already
 been supplied to attach the direct mass linear solver to ARKode, 
 this routine is provided to allow changing this function after it 
 is initially set.

 The return value is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the ARKODE memory was NULL
    ARKDLS_LMEM_NULL if the linear solver memory was NULL
    ARKDLS_MASSMEM_NULL if the mass matrix solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDlsSetJacFn(void *arkode_mem, ARKDlsJacFn jac);
SUNDIALS_EXPORT int ARKDlsSetMassFn(void *arkode_mem, ARKDlsMassFn mass);


/*---------------------------------------------------------------
 Optional outputs from the ARKDLS linear solver:

 ARKDlsGetWorkSpace   returns the real and integer workspace used
                      by the direct linear solver.
 ARKDlsGetNumJacEvals returns the number of calls made to the
                      Jacobian evaluation routine jac.
 ARKDlsGetNumRhsEvals returns the number of calls to the user
                      f routine due to finite difference Jacobian
                      evaluation.
 ARKDlsGetLastFlag    returns the last error flag set by any of
                      the ARKDLS interface functions.

 ARKDlsGetMassWorkSpace returns the real/integer workspace used
                        by the mass matrix direct linear solver.
 ARKDlsGetNumMassSetups returns the number of calls made to the
                        mass matrix solver setup routine
 ARKDlsGetNumMassSolves returns the number of calls made to the
                        mass matrix solver 'solve' routine
 ARKDlsGetNumMassMult   returns the number of calls made to the
                        mass matrix-times-vector routine
 ARKDlsGetLastMassFlag  returns the last error flag set by any of
                        the ARKDLS mass functions

 The return value of ARKDlsGet* is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the ARKODE memory was NULL
    ARKDLS_LMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDlsGetWorkSpace(void *arkode_mem, 
				       long int *lenrwLS, 
				       long int *leniwLS);
SUNDIALS_EXPORT int ARKDlsGetNumJacEvals(void *arkode_mem, 
					 long int *njevals);
SUNDIALS_EXPORT int ARKDlsGetNumRhsEvals(void *arkode_mem, 
					 long int *nfevalsLS);
SUNDIALS_EXPORT int ARKDlsGetLastFlag(void *arkode_mem, 
				      long int *flag);

SUNDIALS_EXPORT int ARKDlsGetMassWorkSpace(void *arkode_mem, 
					   long int *lenrwMLS, 
					   long int *leniwMLS);
SUNDIALS_EXPORT int ARKDlsGetNumMassSetups(void *arkode_mem, 
                                           long int *nmsetups);
SUNDIALS_EXPORT int ARKDlsGetNumMassSolves(void *arkode_mem, 
                                           long int *nmsolves);
SUNDIALS_EXPORT int ARKDlsGetNumMassMult(void *arkode_mem, 
                                         long int *nmmults);
SUNDIALS_EXPORT int ARKDlsGetLastMassFlag(void *arkode_mem, 
					  long int *flag);


/*---------------------------------------------------------------
 The following function returns the name of the constant 
 associated with a ARKDLS return flag
---------------------------------------------------------------*/
SUNDIALS_EXPORT char *ARKDlsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
