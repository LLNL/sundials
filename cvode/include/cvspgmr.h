/*
 * -----------------------------------------------------------------
 * $Revision: 1.16 $
 * $Date: 2004-10-08 15:11:08 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES scaled,
 * preconditioned GMRES linear solver, CVSPGMR.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvspgmr_h
#define _cvspgmr_h

#include "spgmr.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * CVSPGMR solver constants                                       
 * -----------------------------------------------------------------
 * CVSPGMR_MAXL    : default value for the maximum Krylov         
 *                   dimension.                                   
 * 
 * CVSPGMR_MSBPRE  : maximum number of steps between              
 *                   preconditioner evaluations                   
 *                                                                
 * CVSPGMR_DGMAX   : maximum change in gamma between              
 *                   preconditioner evaluations                   
 *                                                                
 * CVSPGMR_DELT    : default value for factor by which the        
 *                   tolerance on the nonlinear iteration is      
 *                   multiplied to get a tolerance on the linear  
 *                   iteration                                    
 * -----------------------------------------------------------------
 */

#define CVSPGMR_MAXL    5

#define CVSPGMR_MSBPRE  50 

#define CVSPGMR_DGMAX   RCONST(0.2)  

#define CVSPGMR_DELT    RCONST(0.05) 
 
/*
 * -----------------------------------------------------------------
 * Type : CVSpgmrPrecSetupFn                                      
 * -----------------------------------------------------------------
 * The user-supplied preconditioner setup function PrecSetup and  
 * the user-supplied preconditioner solve function PrecSolve      
 * together must define left and right preconditoner matrices     
 * P1 and P2 (either of which may be trivial), such that the      
 * product P1*P2 is an approximation to the Newton matrix         
 * M = I - gamma*J.  Here J is the system Jacobian J = df/dy,     
 * and gamma is a scalar proportional to the integration step     
 * size h.  The solution of systems P z = r, with P = P1 or P2,   
 * is to be carried out by the PrecSolve function, and PrecSetup  
 * is to do any necessary setup operations.                       
 *                                                                
 * The user-supplied preconditioner setup function PrecSetup      
 * is to evaluate and preprocess any Jacobian-related data        
 * needed by the preconditioner solve function PrecSolve.         
 * This might include forming a crude approximate Jacobian,       
 * and performing an LU factorization on the resulting            
 * approximation to M.  This function will not be called in       
 * advance of every call to PrecSolve, but instead will be called 
 * only as often as necessary to achieve convergence within the   
 * Newton iteration.  If the PrecSolve function needs no          
 * preparation, the PrecSetup function can be NULL.               
 *                                                                
 * For greater efficiency, the PrecSetup function may save        
 * Jacobian-related data and reuse it, rather than generating it  
 * from scratch.  In this case, it should use the input flag jok  
 * to decide whether to recompute the data, and set the output    
 * flag *jcurPtr accordingly.                                     
 *                                                                
 * Each call to the PrecSetup function is preceded by a call to   
 * the RhsFn f with the same (t,y) arguments.  Thus the PrecSetup 
 * function can use any auxiliary data that is computed and       
 * saved by the f function and made accessible to PrecSetup.      
 *                                                                
 * A function PrecSetup must have the prototype given below.      
 * Its parameters are as follows:                                 
 *                                                                
 * t       is the current value of the independent variable.      
 *                                                                
 * y       is the current value of the dependent variable vector, 
 *           namely the predicted value of y(t).                  
 *                                                                
 * fy      is the vector f(t,y).                                  
 *                                                                
 * jok     is an input flag indicating whether Jacobian-related   
 *         data needs to be recomputed, as follows:               
 *           jok == FALSE means recompute Jacobian-related data   
 *                  from scratch.                                 
 *           jok == TRUE  means that Jacobian data, if saved from 
 *                  the previous PrecSetup call, can be reused    
 *                  (with the current value of gamma).            
 *         A Precset call with jok == TRUE can only occur after   
 *         a call with jok == FALSE.                              
 *                                                                
 * jcurPtr is a pointer to an output integer flag which is        
 *         to be set by PrecSetup as follows:                     
 *         Set *jcurPtr = TRUE if Jacobian data was recomputed.   
 *         Set *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                        but saved data was reused.  
 *                                                                
 * gamma   is the scalar appearing in the Newton matrix.          
 *                                                                
 * P_data is a pointer to user data - the same as the P_data      
 *         parameter passed to CVSpgmr.                           
 *                                                                
 * tmp1, tmp2, and tmp3 are pointers to memory allocated          
 *         for N_Vectors which can be used by CVSpgmrPrecSetupFn  
 *         as temporary storage or work space.                    
 *                                                                
 * NOTE: If the user's preconditioner needs other quantities,     
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively). 
 *     The unit roundoff is available as UNIT_ROUNDOFF defined in 
 *     sundialstypes.h                                            
 *                                                                
 * Returned value:                                                
 * The value to be returned by the PrecSetup function is a flag   
 * indicating whether it was successful.  This value should be    
 *   0   if successful,                                           
 *   > 0 for a recoverable error (step will be retried),          
 *   < 0 for an unrecoverable error (integration is halted).      
 * -----------------------------------------------------------------
 */
  
typedef int (*CVSpgmrPrecSetupFn)(realtype t, N_Vector y, N_Vector fy, 
                                  booleantype jok, booleantype *jcurPtr,
                                  realtype gamma, void *P_data,
                                  N_Vector tmp1, N_Vector tmp2,
                                  N_Vector tmp3);
 
/*
 * -----------------------------------------------------------------
 * Type : CVSpgmrPrecSolveFn                                      
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve function PrecSolve      
 * is to solve a linear system P z = r in which the matrix P is   
 * one of the preconditioner matrices P1 or P2, depending on the  
 * type of preconditioning chosen.                                
 *                                                                
 * A function PrecSolve must have the prototype given below.      
 * Its parameters are as follows:                                 
 *                                                                
 * t      is the current value of the independent variable.       
 *                                                                
 * y      is the current value of the dependent variable vector.  
 *                                                                
 * fy     is the vector f(t,y).                                   
 *                                                                
 * r      is the right-hand side vector of the linear system.     
 *                                                                
 * z      is the output vector computed by PrecSolve.             
 *                                                                
 * gamma  is the scalar appearing in the Newton matrix.           
 *                                                                
 * delta  is an input tolerance for use by PSolve if it uses      
 *        an iterative method in its solution.  In that case,     
 *        the residual vector Res = r - P z of the system         
 *        should be made less than delta in weighted L2 norm,     
 *        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta .          
 *        Note: the error weight vector ewt can be obtained       
 *        through a call to the routine CVodeGetErrWeights.       
 *                                                                
 * lr     is an input flag indicating whether PrecSolve is to use 
 *        the left preconditioner P1 or right preconditioner      
 *        P2: lr = 1 means use P1, and lr = 2 means use P2.       
 *                                                                
 * P_data is a pointer to user data - the same as the P_data      
 *        parameter passed to CVSpgmr.                            
 *                                                                
 * tmp    is a pointer to memory allocated for an N_Vector        
 *        which can be used by PSolve for work space.             
 *                                                                
 * Returned value:                                                
 * The value to be returned by the PrecSolve function is a flag   
 * indicating whether it was successful.  This value should be    
 *   0 if successful,                                             
 *   positive for a recoverable error (step will be retried),     
 *   negative for an unrecoverable error (integration is halted). 
 * -----------------------------------------------------------------
 */
  
typedef int (*CVSpgmrPrecSolveFn)(realtype t, N_Vector y, N_Vector fy, 
                                  N_Vector r, N_Vector z, 
                                  realtype gamma, realtype delta,
                                  int lr, void *P_data, N_Vector tmp);
 
/*
 * -----------------------------------------------------------------
 * Type : CVSpgmrJacTimesVecFn                                    
 * -----------------------------------------------------------------
 * The user-supplied function jtimes is to generate the product   
 * J*v for given v, where J is the Jacobian df/dy, or an          
 * approximation to it, and v is a given vector. It should return 
 * 0 if successful and a nonzero int otherwise.                   
 *                                                                
 * A function jtimes must have the prototype given below. Its     
 * parameters are as follows:                                     
 *                                                                
 *   v        is the N_Vector to be multiplied by J.              
 *                                                                
 *   Jv       is the output N_Vector containing J*v.              
 *                                                                
 *   t        is the current value of the independent variable.   
 *                                                                
 *   y        is the current value of the dependent variable      
 *            vector.                                             
 *                                                                
 *   fy       is the vector f(t,y).                               
 *                                                                
 *   jac_data is a pointer to user Jacobian data, the same as the 
 *            pointer passed to CVSpgmr.                          
 *                                                                
 *   tmp      is a pointer to memory allocated for an N_Vector    
 *            which can be used by Jtimes for work space.         
 * -----------------------------------------------------------------
 */

typedef int (*CVSpgmrJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
                                    N_Vector y, N_Vector fy, 
                                    void *jac_data, N_Vector tmp);
 
/*
 * -----------------------------------------------------------------
 * Function : CVSpgmr                                             
 * -----------------------------------------------------------------
 * A call to the CVSpgmr function links the main CVODE integrator 
 * with the CVSPGMR linear solver.                                
 *                                                                
 * cvode_mem is the pointer to the integrator memory returned by  
 *           CVodeCreate.                                         
 *                                                                
 * pretype   is the type of user preconditioning to be done.      
 *           This must be one of the four enumeration constants   
 *           NONE, LEFT, RIGHT, or BOTH defined in iterative.h.    
 *           These correspond to no preconditioning,              
 *           left preconditioning only, right preconditioning     
 *           only, and both left and right preconditioning,       
 *           respectively.                                        
 *                                                                
 * maxl      is the maximum Krylov dimension. This is an          
 *           optional input to the CVSPGMR solver. Pass 0 to      
 *           use the default value CVSPGMR_MAXL=5.                
 *                                                                
 * The return value of CVSpgmr is one of:                              
 *    CVSPGMR_SUCCESS   if successful                                
 *    CVSPGMR_MEM_NULL  if the cvode memory was NULL
 *    CVSPGMR_MEM_FAIL  if there was a memory allocation failure     
 *    CVSPGMR_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

int CVSpgmr(void *cvode_mem, int pretype, int maxl);

/*
 * -----------------------------------------------------------------
 * Function: CVSpgmrResetPrecType                                 
 * -----------------------------------------------------------------
 * CVSpgmrResetPrecType resets the type of preconditioner, pretype,
 *     from the value set in a prior call to CVSpgmr.     
 *     This must be one of NONE, LEFT, RIGHT, or BOTH.            
 * -----------------------------------------------------------------
 */

int CVSpgmrResetPrecType(void *cvode_mem, int pretype);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVSPGMR linear solver                   
 * -----------------------------------------------------------------
 *                                                                
 * CVSpgmrSetGSType specifies the type of Gram-Schmidt            
 *           orthogonalization to be used. This must be one of    
 *           the two enumeration constants MODIFIED_GS or         
 *           CLASSICAL_GS defined in iterative.h. These correspond 
 *           to using modified Gram-Schmidt and classical         
 *           Gram-Schmidt, respectively.                          
 *           Default value is MODIFIED_GS.                        
 * CVSpgmrSetDelt specifies the factor by which the tolerance on  
 *           the nonlinear iteration is multiplied to get a       
 *           tolerance on the linear iteration. This is an        
 *           optional input to the CVSPGMR solver.                
 *           Default value is 0.05.                               
 * CVSpgmrSetPrecSetupFn specifies the PrecSetup function.        
 *           Default is NULL.                                     
 * CVSpgmrSetPrecSolveFn specifies the PrecSolve function.        
 *           Default is NULL.                                     
 * CVSpgmrSetPrecData specifies a pointer to user preconditioner  
 *           data. This pointer is passed to PrecSetup and        
 *           PrecSolve every time these routines are called.      
 *           Default is NULL.                                     
 * CVSpgmrSetJacTimesVecFn specifies the jtimes function.         
 *           Default is to use an internal finite difference      
 *           approximation routine.                               
 * CVSpgmrSetJacData specifies a pointer to user Jacobian data.   
 *           This pointer is passed to jtimes every time this     
 *           routine is called.                                   
 *           Default is NULL.                                     
 *
 * The return value of CVSpgmrSet* is one of:
 *    CVSPGMR_SUCCESS   if successful
 *    CVSPGMR_MEM_NULL  if the cvode memory was NULL
 *    CVSPGMR_LMEM_NULL if the cvspgmr memory was NULL
 *    CVSPGMR_ILL_INPUT if an input has an illegal value
 * -----------------------------------------------------------------
 */

int CVSpgmrSetGSType(void *cvode_mem, int gstype);
int CVSpgmrSetDelt(void *cvode_mem, realtype delt);
int CVSpgmrSetPrecSetupFn(void *cvode_mem, CVSpgmrPrecSetupFn pset);
int CVSpgmrSetPrecSolveFn(void *cvode_mem, CVSpgmrPrecSolveFn psolve);
int CVSpgmrSetPrecData(void *cvode_mem, void *P_data);
int CVSpgmrSetJacTimesVecFn(void *cvode_mem, CVSpgmrJacTimesVecFn jtimes);
int CVSpgmrSetJacData(void *cvode_mem, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSPGMR linear solver                
 * -----------------------------------------------------------------
 *                                                                
 * CVSpgmrGetWorkSpace returns the real and integer workspace used
 *     by CVSPGMR.                                                   
 * CVSpgmrGetNumPrecEvals returns the number of preconditioner    
 *     evaluations, i.e. the number of calls made to PrecSetup    
 *     with jok==FALSE.                                           
 * CVSpgmrGetNumPrecSolves returns the number of calls made to    
 *     PrecSolve.                                                 
 * CVSpgmrGetNumLinIters returns the number of linear iterations. 
 * CVSpgmrGetNumConvFails returns the number of linear            
 *     convergence failures.                                      
 * CVSpgmrGetNumJtimesEvals returns the number of calls to jtimes.
 * CVSpgmrGetNumRhsEvals returns the number of calls to the user  
 *     f routine due to finite difference Jacobian times vector   
 *     evaluation.
 * CVSpgmrGetLastFlag returns the last error flag set by any of
 *     the CVSPGMR interface functions.
 *           
 * The return value of CVSpgmrGet* is one of:
 *    CVSPGMR_SUCCESS   if successful
 *    CVSPGMR_MEM_NULL  if the cvode memory was NULL
 *    CVSPGMR_LMEM_NULL if the cvspgmr memory was NULL
 * -----------------------------------------------------------------
 */                                                                

int CVSpgmrGetWorkSpace(void *cvode_mem, long int *lenrwSG, long int *leniwSG);
int CVSpgmrGetNumPrecEvals(void *cvode_mem, long int *npevals);
int CVSpgmrGetNumPrecSolves(void *cvode_mem, long int *npsolves);
int CVSpgmrGetNumLinIters(void *cvode_mem, long int *nliters);
int CVSpgmrGetNumConvFails(void *cvode_mem, long int *nlcfails);
int CVSpgmrGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
int CVSpgmrGetNumRhsEvals(void *cvode_mem, long int *nfevalsSG); 
int CVSpgmrGetLastFlag(void *cvode_mem, int *flag);

/* CVSPGMR return values */

#define CVSPGMR_SUCCESS     0
#define CVSPGMR_MEM_NULL   -1 
#define CVSPGMR_LMEM_NULL  -2 
#define CVSPGMR_ILL_INPUT  -3
#define CVSPGMR_MEM_FAIL   -4

#define CVSPGMR_DATA_NULL  -10

#endif

#ifdef __cplusplus
}
#endif
