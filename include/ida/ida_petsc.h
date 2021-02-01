/*
 * ----------------------------------------------------------------- 
 * Programmers: Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the PETSc Scaled and Preconditioned
 * Iterative Linear Solvers for use with IDA.
 * -----------------------------------------------------------------
 */

#ifndef _IDA_PETSC_H
#define _IDA_PETSC_H

#include <sundials/sundials_nvector.h>
#include <petscmat.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* 
 * -----------------------------------------------------------------
 * IDAPETScKSP return values 
 * -----------------------------------------------------------------
 */

#define PETSC_KSP_SUCCESS     0
#define PETSC_KSP_MEM_NULL   -1 
#define PETSC_KSP_LMEM_NULL  -2 
#define PETSC_KSP_ILL_INPUT  -3
#define PETSC_KSP_MEM_FAIL   -4
#define PETSC_KSP_PMEM_NULL  -5


/*
 * -----------------------------------------------------------------
 *                                                                
 * Function : IDAPETScKSP                                            
 * -----------------------------------------------------------------
 * A call to the IDAPETScKSP function links the main integrator with 
 * the PETSc linear solver module.  Its parameters are as      
 * follows:                                                       
 *                                                                
 * ida_mem   is the pointer to memory block returned by IDACreate.
 *                                                                
 * comm      is the MPI communicator for the computation.
 *
 * JacMat    is pointer to PETSc matrix containing Jacobian.     
 *                                                                
 * The return values of IDASpgmr are:                             
 *    IDASPILS_SUCCESS    if successful                            
 *    IDASPILS_MEM_NULL   if the ida memory was NULL
 *    IDASPILS_MEM_FAIL   if there was a memory allocation failure 
 *    IDASPILS_ILL_INPUT  if there was illegal input.              
 * The above constants are defined in ida_spils.h
 *                                                                
 * -----------------------------------------------------------------
 */                                                                

SUNDIALS_EXPORT int IDAPETScKSP(void *ida_mem, MPI_Comm comm, Mat *JacMat);

/*
 * -----------------------------------------------------------------
 * Types : IDAPETScJacFn
 * -----------------------------------------------------------------
 *
 * A sparse Jacobian approximation function jaceval must be of type 
 * IDAPETScJacFn.
 * Its parameters are:                     
 *                                                                
 * t   is the current value of the independent variable t.        
 *                                                                
 * c_j is the scalar in the system Jacobian, proportional to 
 *     the inverse of the step size h.
 *                                                                
 * y   is the current value of the dependent variable vector,     
 *     namely the predicted value of y(t).                     
 *                                                                
 * yp  is the current value of the derivative vector y',          
 *     namely the predicted value of y'(t).                    
 *                                                                
 * r   is the residual vector F(tt,yy,yp).                     
 *                                                                
 * JacMat is the sparse matrix (of PETSc Mat type) to be loaded 
 *     by an IDAPETScJacFn routine with an approximation to the 
 *     system Jacobian matrix
 *            J = dF/dy' + c_j*dF/dy                            
 *     at the given point (t,y,y'), where the DAE system is    
 *     given by F(t,y,y') = 0.
 *     Note that JacMat is NOT preset to zero!
 * 
 * user_data is a pointer to user Jacobian data - the same as the    
 *     user_data parameter passed to IDASetRdata.                     
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an IDASparseJacFn routine 
 *     as temporary storage or work space.                     
 *                                                                
 * A IDAPETScJacFn should return                                
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 * In the case of a recoverable error return, the integrator will 
 * attempt to recover by reducing the stepsize (which changes cj).
 *
 * -----------------------------------------------------------------
 *
  * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     IDAGetCurrentStep and IDAGetErrWeights, respectively 
 *     (see ida.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*IDAPETScJacFn)(realtype t, realtype c_j,
                             N_Vector y, N_Vector yp, N_Vector r, 
                             Mat JacMat, void *user_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);



/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDAPETScKSP linear solver                  
 * -----------------------------------------------------------------
 *                                                                
 * IDAPETScSetJacFn specifies the Jacobian evaluation function.        
 *           This function must be supplied when using PETSc linear 
 *           solver.                           
 * IDAPETScSetGSType specifies the type of Gram-Schmidt           
 *           orthogonalization to be used. Currently disabled.
 * IDAPETScSetMaxRestarts specifies the maximum number of restarts
 *           to be used in the GMRES algorithm. Currently disabled.
 * IDASpbcgSetMaxl specifies the maximum Krylov subspace size. 
 *           Default is 5. Currently disabled.
 * IDAPETScSetEpsLin specifies the factor in the linear iteration 
 *           convergence test constant. Currently disabled.                     
 * IDAPETScSetIncrementFactor specifies a factor in the increments
 *           to yy used in the difference quotient approximations 
 *           to matrix-vector products Jv. Currently disabled.                       
 *                                                                
 * The return value of IDAPETScSet* is one of:
 *    PETSC_KSP_SUCCESS   if successful
 *    PETSC_KSP_MEM_NULL  if the ida memory was NULL
 *    PETSC_KSP_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAPETScSetJacFn(void *ida_mem, IDAPETScJacFn jac);

SUNDIALS_EXPORT int IDAPETScSetGSType(void *ida_mem, int gstype);
SUNDIALS_EXPORT int IDAPETScSetMaxRestarts(void *ida_mem, int maxrs);
SUNDIALS_EXPORT int IDAPETScSetMaxl(void *ida_mem, int maxl);
SUNDIALS_EXPORT int IDAPETScSetEpsLin(void *ida_mem, realtype eplifac);
SUNDIALS_EXPORT int IDAPETScSetIncrementFactor(void *ida_mem, realtype dqincfac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDAPETScKSP linear solver               
 *----------------------------------------------------------------
 *                                                                
 * IDAPETScGetWorkSpace returns the real and integer workspace used 
 *     by IDAPETScKSP. TODO: Currently disabled.                                                 
 * IDAPETScGetNumPrecEvals returns the number of preconditioner   
 *     evaluations, i.e. the number of calls made to PrecSetup    
 *     with jok==SUNFALSE. TODO: Currently disabled.                                     
 * IDAPETScGetNumPrecSolves returns the number of calls made to   
 *     PrecSolve. TODO: Currently disabled.                 
 * IDAPETScGetNumLinIters returns the number of linear iterations.
 * IDAPETScGetNumConvFails returns the number of linear           
 *     convergence failures.                                      
 * IDAPETScGetNumJacEvals returns the number of Jacobian evaluations.
 * IDAPETScGetNumResEvals returns the number of calls to the user 
 *     res routine due to finite difference Jacobian times vector 
 *     evaluation. TODO: Currently disabled, will be removed.                                               
 * IDAPETScGetLastFlag returns the last error flag set by any of
 *     the IDAPETScKSP interface functions. TODO: Currently disabled.
 *                                                                
 * The return value of IDAPETScGet* is one of:
 *    PETSC_KSP_SUCCESS   if successful
 *    PETSC_KSP_MEM_NULL  if the ida memory was NULL
 *    PETSC_KSP_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */                                                                

SUNDIALS_EXPORT int IDAPETScGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int IDAPETScGetNumPrecEvals(void *ida_mem, long int *npevals);
SUNDIALS_EXPORT int IDAPETScGetNumPrecSolves(void *ida_mem, long int *npsolves);
SUNDIALS_EXPORT int IDAPETScGetNumLinIters(void *ida_mem, long int *nliters);
SUNDIALS_EXPORT int IDAPETScGetNumConvFails(void *ida_mem, long int *nlcfails);
SUNDIALS_EXPORT int IDAPETScGetNumJacEvals(void *ida_mem, long int *njacevals);
SUNDIALS_EXPORT int IDAPETScGetNumResEvals(void *ida_mem, long int *nrevalsLS); 
SUNDIALS_EXPORT int IDAPETScGetLastFlag(void *ida_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a PETSC_KSP return flag. TODO: Currently disabled.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *IDAPETScGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
