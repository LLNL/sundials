/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:56 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for the IDAS Scaled Preconditioned GMRES     
 * linear solver module, IDASPGMR.                                 
 * -----------------------------------------------------------------
 */

#ifndef _IDASPGMR_H
#define _IDASPGMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idas_spils.h"
#include "sundials_spgmr.h"

/*
 * -----------------------------------------------------------------
 *                                                                
 * Function : IDASpgmr                                            
 *----------------------------------------------------------------
 * A call to the IDASpgmr function links the main integrator with 
 * the IDASPGMR linear solver module.  Its parameters are as      
 * follows:                                                       
 *                                                                
 * IDA_mem   is the pointer to memory block returned by IDACreate.
 *                                                                
 * maxl      is the maximum Krylov subspace dimension, an         
 *           optional input.  Pass 0 to use the default value,    
 *           MIN(Neq, 5).  Otherwise pass a positive integer.     
 *                                                                
 * The return values of IDASpgmr are:                             
 *    IDASPGMR_SUCCESS    if successful                            
 *    IDASPGMR_MEM_NULL   if the ida memory was NULL
 *    IDASPGMR_MEM_FAIL   if there was a memory allocation failure 
 *    IDASPGMR_ILL_INPUT  if there was illegal input.              
 *                                                                
 * -----------------------------------------------------------------
 */                                                                

int IDASpgmr(void *ida_mem, int maxl);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDASPGMR linear solver                  
 *----------------------------------------------------------------
 *                                                                
 * IDASpgmrSetPreconditioner specifies the PrecSetup and PrecSolve
 *                           functions, and a pointer to user
 *                           preconditioner data (passed to PrecSetup
 *                           and PrecSolve every time these routines
 *                           are called)
 * IDASpgmrSetJacTimesVecFn specifies the jtimes function.        
 *                          Default is to use an internal finite
 *                          difference approximation routine.
 *                          Also takes as an argument a pointer
 *                          to user Jacobian data (passed to jtimes
 *                          every time it is called).   Default is
 *                          NULL.
 * IDASpgmrSetGSType specifies the type of Gram-Schmidt           
 *           orthogonalization to be used. This must be one of    
 *           the two enumeration constants MODIFIED_GS or         
 *           CLASSICAL_GS defined in iterativ.h. These correspond 
 *           to using modified Gram-Schmidt and classical         
 *           Gram-Schmidt, respectively.                          
 *           Default value is MODIFIED_GS.                        
 * IDASpgmrSetMaxRestarts specifies the maximum number of restarts
 *           to be used in the GMRES algorithm.  maxrs must be a  
 *           non-negative integer.  Pass 0 to specify no restarts.
 *           Default is 5.                                        
 * IDASpgmrSetEpsLin specifies the factor in the linear iteration 
 *           convergence test constant.                           
 *           Default is 0.05                                      
 * IDASpgmrSetIncrementFactor specifies a factor in the increments
 *           to yy used in the difference quotient approximations 
 *           to matrix-vector products Jv.                        
 *           Default is 1.0                                       
 *                                                                
 * The return value of IDASpgmrSet* is one of:
 *    IDASPGMR_SUCCESS   if successful
 *    IDASPGMR_MEM_NULL  if the ida memory was NULL
 *    IDASPGMR_LMEM_NULL if the idaspgmr memory was NULL
 * -----------------------------------------------------------------
 */

int IDASpgmrSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset,
			      IDASpilsPrecSolveFn psolve, void *prec_data);
int IDASpgmrSetJacTimesVecFn(void *ida_mem, IDASpilsJacTimesVecFn jtimes, void *jac_data);
int IDASpgmrSetGSType(void *ida_mem, int gstype);
int IDASpgmrSetMaxRestarts(void *ida_mem, int maxrs);
int IDASpgmrSetEpsLin(void *ida_mem, realtype eplifac);
int IDASpgmrSetIncrementFactor(void *ida_mem, realtype dqincfac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDASPGMR linear solver               
 *----------------------------------------------------------------
 *                                                                
 * IDASpgmrGetWorkSpace returns the real and integer workspace used 
 *     by IDASPGMR.                                                  
 * IDASpgmrGetNumPrecEvals returns the number of preconditioner   
 *     evaluations, i.e. the number of calls made to PrecSetup    
 *     with jok==FALSE.                                           
 * IDASpgmrGetNumPrecSolves returns the number of calls made to   
 *     PrecSolve.                                                 
 * IDASpgmrGetNumLinIters returns the number of linear iterations.
 * IDASpgmrGetNumConvFails returns the number of linear           
 *     convergence failures.                                      
 * IDASpgmrGetNumJtimesEvals returns the number of calls to jtimes
 * IDASpgmrGetNumResEvals returns the number of calls to the user 
 *     res routine due to finite difference Jacobian times vector 
 *     evaluation.                                                
 * IDASpgmrGetLastFlag returns the last error flag set by any of
 *     the IDASPGMR interface functions.
 *                                                                
 * The return value of IDASpgmrGet* is one of:
 *    IDASPGMR_SUCCESS   if successful
 *    IDASPGMR_MEM_NULL  if the ida memory was NULL
 *    IDASPGMR_LMEM_NULL if the idaspgmr memory was NULL
 * -----------------------------------------------------------------
 */                                                                

int IDASpgmrGetWorkSpace(void *ida_mem, long int *lenrwSG, long int *leniwSG);
int IDASpgmrGetNumPrecEvals(void *ida_mem, long int *npevals);
int IDASpgmrGetNumPrecSolves(void *ida_mem, long int *npsolves);
int IDASpgmrGetNumLinIters(void *ida_mem, long int *nliters);
int IDASpgmrGetNumConvFails(void *ida_mem, long int *nlcfails);
int IDASpgmrGetNumJtimesEvals(void *ida_mem, long int *njvevals);
int IDASpgmrGetNumResEvals(void *ida_mem, long int *nrevalsSG); 
int IDASpgmrGetLastFlag(void *ida_mem, int *flag);

/* IDASPGMR return values */

#define IDASPGMR_SUCCESS     0
#define IDASPGMR_MEM_NULL   -1 
#define IDASPGMR_LMEM_NULL  -2 
#define IDASPGMR_ILL_INPUT  -3
#define IDASPGMR_MEM_FAIL   -4

#ifdef __cplusplus
}
#endif

#endif
