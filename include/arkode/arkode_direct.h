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
 * Common header file for the direct linear solvers in ARKODE.
 *--------------------------------------------------------------*/

#ifndef _ARKDLS_H
#define _ARKDLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  ARKDIRECT CONSTANTS
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


/*===============================================================
  FUNCTION TYPES
===============================================================*/

/*---------------------------------------------------------------
 Type: ARKDlsDenseJacFn

 A dense Jacobian approximation function Jac must be of type 
 ARKDlsDenseJacFn. Its parameters are:

 N   is the problem size.

 Jac is the dense matrix (of type DlsMat) that will be loaded
     by a ARKDlsDenseJacFn with an approximation to the Jacobian 
     matrix J = (df_i/dy_j) at the point (t,y). 

 t   is the current value of the independent variable.

 y   is the current value of the dependent variable vector,
     namely the predicted value of y(t).

 fy  is the vector f(t,y).

 user_data is a pointer to user data - the same as the user_data
     parameter passed to ARKodeSetFdata.

 tmp1, tmp2, and tmp3 are pointers to memory allocated for
 vectors of length N which can be used by a ARKDlsDenseJacFn
 as temporary storage or work space.

 A ARKDlsDenseJacFn should return 0 if successful, a positive 
 value if a recoverable error occurred, and a negative value if 
 an unrecoverable error occurred.

 NOTE: The following are two efficient ways to load a dense Jac:         
 (1) (with macros - no explicit data structure references)      
     for (j=0; j < Neq; j++) {                                  
       col_j = DENSE_COL(Jac,j);                                 
       for (i=0; i < Neq; i++) {                                
         generate J_ij = the (i,j)th Jacobian element           
         col_j[i] = J_ij;                                       
       }                                                        
     }                                                          
 (2) (without macros - explicit data structure references)      
     for (j=0; j < Neq; j++) {                                  
       col_j = (Jac->data)[j];                                   
       for (i=0; i < Neq; i++) {                                
         generate J_ij = the (i,j)th Jacobian element           
         col_j[i] = J_ij;                                       
       }                                                        
     }                                                          
 A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
 efficient in general.  It is only appropriate for use in small 
 problems in which efficiency of access is NOT a major concern. 
                                                                
 NOTE: If the user's Jacobian routine needs other quantities,   
     they are accessible as follows: hcur (the current stepsize)
     and ewt (the error weight vector) are accessible through   
     ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively 
     (see arkode.h). The unit roundoff is available as 
     UNIT_ROUNDOFF defined in sundials_types.h.

---------------------------------------------------------------*/
typedef int (*ARKDlsDenseJacFn)(long int N, realtype t,
				N_Vector y, N_Vector fy, 
				DlsMat Jac, void *user_data,
				N_Vector tmp1, N_Vector tmp2, 
				N_Vector tmp3);
  

/*---------------------------------------------------------------
 Type: ARKDlsDenseMassFn

 A dense mass matrix approximation function Mass must be of type 
 ARKDlsDenseMassFn. Its parameters are:

 N   is the problem size.

 t   is the current value of the independent variable.

 M   is the dense matrix (of type DlsMat) that will be loaded
     by a ARKDlsDenseMassFn with an approximation to the mass matrix.

 user_data is a pointer to user data - the same as the user_data
     parameter passed to ARKodeSetFdata.

 tmp1, tmp2, and tmp3 are pointers to memory allocated for
 vectors of length N which can be used by a ARKDlsDenseMassFn
 as temporary storage or work space.

 A ARKDlsDenseMassFn should return 0 if successful, and a 
 negative value if an unrecoverable error occurred.

---------------------------------------------------------------*/
typedef int (*ARKDlsDenseMassFn)(long int N, realtype t, DlsMat M, 
				 void *user_data, N_Vector tmp1, 
				 N_Vector tmp2, N_Vector tmp3);
  

/*---------------------------------------------------------------
 Type: ARKDlsBandJacFn

 A band Jacobian approximation function Jac must have the
 prototype given below. Its parameters are:

 N is the length of all vector arguments.

 mupper is the upper half-bandwidth of the approximate banded
 Jacobian. This parameter is the same as the mupper parameter
 passed by the user to the linear solver initialization function.

 mlower is the lower half-bandwidth of the approximate banded
 Jacobian. This parameter is the same as the mlower parameter
 passed by the user to the linear solver initialization function.

 t is the current value of the independent variable.

 y is the current value of the dependent variable vector,
      namely the predicted value of y(t).

 fy is the vector f(t,y).

 Jac is the band matrix (of type DlsMat) that will be loaded
 by a ARKDlsBandJacFn with an approximation to the Jacobian matrix
 Jac = (df_i/dy_j) at the point (t,y).
 Three efficient ways to load J are:

 (1) (with macros - no explicit data structure references)
    for (j=0; j < n; j++) {
       col_j = BAND_COL(Jac,j);
       for (i=j-mupper; i <= j+mlower; i++) {
         generate J_ij = the (i,j)th Jacobian element
         BAND_COL_ELEM(col_j,i,j) = J_ij;
       }
     }

 (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)
    for (j=0; j < n; j++) {
       col_j = BAND_COL(Jac,j);
       for (k=-mupper; k <= mlower; k++) {
         generate J_ij = the (i,j)th Jacobian element, i=j+k
         col_j[k] = J_ij;
       }
     }

 (3) (without macros - explicit data structure references)
     offset = Jac->smu;
     for (j=0; j < n; j++) {
       col_j = ((Jac->data)[j])+offset;
       for (k=-mupper; k <= mlower; k++) {
         generate J_ij = the (i,j)th Jacobian element, i=j+k
         col_j[k] = J_ij;
       }
     }
 Caution: Jac->smu is generally NOT the same as mupper.

 The BAND_ELEM(A,i,j) macro is appropriate for use in small
 problems in which efficiency of access is NOT a major concern.

 user_data is a pointer to user data - the same as the user_data
          parameter passed to ARKodeSetFdata.

 NOTE: If the user's Jacobian routine needs other quantities,
     they are accessible as follows: hcur (the current stepsize)
     and ewt (the error weight vector) are accessible through
     ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively
     (see arkode.h). The unit roundoff is available as
     UNIT_ROUNDOFF defined in sundials_types.h

 tmp1, tmp2, and tmp3 are pointers to memory allocated for
 vectors of length N which can be used by a ARKDlsBandJacFn
 as temporary storage or work space.

 A ARKDlsBandJacFn should return 0 if successful, a positive value
 if a recoverable error occurred, and a negative value if an 
 unrecoverable error occurred.
---------------------------------------------------------------*/
typedef int (*ARKDlsBandJacFn)(long int N, long int mupper, 
			       long int mlower, realtype t, 
			       N_Vector y, N_Vector fy, 
			       DlsMat Jac, void *user_data,
			       N_Vector tmp1, N_Vector tmp2, 
			       N_Vector tmp3);

/*---------------------------------------------------------------
 Type: ARKDlsBandMassFn

 A band mass matrix approximation function Mass must have the
 prototype given below. Its parameters are:

 N is the length of all vector arguments.

 mupper is the upper half-bandwidth of the approximate banded
 Jacobian. This parameter is the same as the mupper parameter
 passed by the user to the linear solver initialization function.

 mlower is the lower half-bandwidth of the approximate banded
 Jacobian. This parameter is the same as the mlower parameter
 passed by the user to the linear solver initialization function.

 t is the current value of the independent variable.

 M is the band matrix (of type DlsMat) that will be loaded
 by a ARKDlsBandMassFn with an approximation to the mass matrix

 user_data is a pointer to user data - the same as the user_data
     parameter passed to ARKodeSetFdata.

 tmp1, tmp2, and tmp3 are pointers to memory allocated for
 vectors of length N which can be used by a ARKDlsBandMassFn
 as temporary storage or work space.

 A ARKDlsBandMassFn should return 0 if successful, and a negative 
 value if an unrecoverable error occurred.
---------------------------------------------------------------*/
typedef int (*ARKDlsBandMassFn)(long int N, long int mupper, 
				long int mlower, realtype t, 
				DlsMat M, void *user_data, 
				N_Vector tmp1, N_Vector tmp2, 
				N_Vector tmp3);


/*===============================================================
  EXPORTED FUNCTIONS
===============================================================*/

/*---------------------------------------------------------------
 Optional inputs to the ARKDLS linear solver:

 ARKDlsSetDenseJacFn specifies the dense Jacobian approximation
 routine to be used for a direct dense linear solver.

 ARKDlsSetBandJacFn specifies the band Jacobian approximation
 routine to be used for a direct band linear solver.

 By default, a difference quotient approximation, supplied with
 the solver is used.

 The return value is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the ARKODE memory was NULL
    ARKDLS_LMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDlsSetDenseJacFn(void *arkode_mem, 
					ARKDlsDenseJacFn jac);
SUNDIALS_EXPORT int ARKDlsSetBandJacFn(void *arkode_mem, 
				       ARKDlsBandJacFn jac);


/*---------------------------------------------------------------
 Optional inputs to the ARKDLS linear solver:

 ARKDlsSetDenseMassFn specifies the dense mass matrix 
 approximation routine to be used for a direct dense solver.

 ARKDlsSetBandMassFn specifies the band mass matrix approximation
 routine to be used for a direct band solver.

 The return value is one of:
    ARKDLS_SUCCESS      if successful
    ARKDLS_MEM_NULL     if the ARKODE memory was NULL
    ARKDLS_MASSMEM_NULL if the mass matrix solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDlsSetDenseMassFn(void *arkode_mem, 
					 ARKDlsDenseMassFn mass);
SUNDIALS_EXPORT int ARKDlsSetBandMassFn(void *arkode_mem, 
					ARKDlsBandMassFn mass);


/*---------------------------------------------------------------
 Optional outputs from the ARKDLS linear solver:

 ARKDlsGetWorkSpace   returns the real and integer workspace used
                     by the direct linear solver.
 ARKDlsGetMassWorkSpace   returns the real and integer workspace used
                     by the mass matrix direct linear solver.
 ARKDlsGetNumJacEvals returns the number of calls made to the
                     Jacobian evaluation routine jac.
 ARKDlsGetNumMassEvals returns the number of calls made to the
                     mass matrix evaluation routine Mass.
 ARKDlsGetNumRhsEvals returns the number of calls to the user
                     f routine due to finite difference Jacobian
                     evaluation.
 ARKDlsGetLastFlag    returns the last error flag set by any of
                     the ARKDLS interface functions.
 ARKDlsGetLastMassFlag  returns the last error flag set by any of
                     the ARKDLS interface mass matrix functions.

 The return value of ARKDlsGet* is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the ARKODE memory was NULL
    ARKDLS_LMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDlsGetWorkSpace(void *arkode_mem, 
				       long int *lenrwLS, 
				       long int *leniwLS);
SUNDIALS_EXPORT int ARKDlsGetMassWorkSpace(void *arkode_mem, 
					   long int *lenrwMLS, 
					   long int *leniwMLS);
SUNDIALS_EXPORT int ARKDlsGetNumJacEvals(void *arkode_mem, 
					 long int *njevals);
SUNDIALS_EXPORT int ARKDlsGetNumMassEvals(void *arkode_mem, 
					  long int *nmevals);
SUNDIALS_EXPORT int ARKDlsGetNumRhsEvals(void *arkode_mem, 
					 long int *nfevalsLS);
SUNDIALS_EXPORT int ARKDlsGetLastFlag(void *arkode_mem, 
				      long int *flag);
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
