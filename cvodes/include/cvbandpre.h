/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2004-07-22 21:19:23 $
 * ----------------------------------------------------------------- 
 * Programmers: Michael Wittman, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for the CVBANDPRE module, which         
 * provides a banded difference quotient Jacobian-based            
 * preconditioner and solver routines for use with CVSPGMR.        
 *                                                                 
 * Summary:                                                        
 * These routines provide a band matrix preconditioner based on    
 * difference quotients of the ODE right-hand side function f.     
 * The user supplies parameters                                    
 *   mu = upper half-bandwidth (number of super-diagonals)         
 *   ml = lower half-bandwidth (number of sub-diagonals)           
 * The routines generate a band matrix of bandwidth ml + mu + 1    
 * and use this to form a preconditioner for use with the Krylov   
 * linear solver in CVSPGMR.  Although this matrix is intended     
 * to approximate the Jacobian df/dy, it may be a very crude       
 * approximation.  The true Jacobian need not be banded, or its    
 * true bandwith may be larger than ml + mu + 1, as long as the    
 * banded approximation generated here is sufficiently accurate    
 * to speed convergence as a preconditioner.                       
 *                                                                 
 * Usage:                                                          
 *   The following is a summary of the usage of this module.       
 *   Details of the calls to CVodeCreate, CVodeMalloc, CVSpgmr,    
 *   and CVode are available in the User Guide.                    
 *   To use these routines, the sequence of calls in the user      
 *   main program should be as follows:                            
 *                                                                 
 *   #include "cvbandpre.h"                                        
 *   #include "nvector_serial.h"                                   
 *   ...                                                           
 *   void *bp_data;                                                
 *   ...                                                           
 *   y0 = N_VNew_Serial(...);                             
 *   ...                                                           
 *   cvode_mem = CVodeCreate(...);                                 
 *   ier = CVodeMalloc(...);                                       
 *   ...                                                           
 *   bp_data = CVBandPrecAlloc(N, f, f_data, mu, ml, cvode_mem);   
 *   ...                                                           
 *   flag = CVBPSpgmr(cvode_mem, pretype, maxl, bp_data);          
 *   ...                                                           
 *   flag = CVode(...);                                            
 *   ...                                                           
 *   CVBandPrecFree(bp_data);                                      
 *   ...                                                           
 *   N_VDestroy_Serial(y0);                                   
 *   ...                                                           
 *   CVodeFree(cvode_mem);                                         
 *                                                                 
 * Notes:                                                          
 * (1) Include this file for the CVBandPrecData type definition.   
 * (2) In the CVBandPrecAlloc call, the arguments N, f, and f_data 
 *     are the same as in the call to CVodeMalloc.                 
 * (3) In the CVSpgmr call, the user is free to specify the inputs 
 *     pretype and gstype, and the optional inputs maxl and delt.  
 *     But the last three arguments must be as shown, with the     
 *     last argument being the pointer returned by CVBandPreAlloc. 
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvbandpre_h
#define _cvbandpre_h

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Function : CVBandPrecAlloc                                     
 * -----------------------------------------------------------------
 * CVBandPrecAlloc allocates and initializes a CVBandPrecData     
 * structure to be passed to CVSpgmr (and subsequently used by    
 * CVBandPrecSetup and CVBandPrecSolve).                          
 *                                                                
 * The parameters of CVBandPrecAlloc are as follows:              
 *                                                                
 * N       is the length of all vector arguments.                 
 *                                                                
 * mu      is the upper half bandwidth.                           
 *                                                                
 * ml      is the lower half bandwidth.                           
 *                                                                
 * CVBandPrecAlloc returns the storage pointer of type            
 * CVBandPrecData or NULL if the request for storage cannot be    
 * satisfied.                                                     
 *                                                                
 * NOTE: The band preconditioner assumes a serial implementation  
 *       of the NVECTOR package. Therefore, CVBandPrecAlloc will  
 *       first test for a compatible N_Vector internal            
 *       representation by checking (1) the vector specification  
 *       ID tag and (2) that the functions N_VGetData, and        
 *       N_VSetData are implemented.                              
 * -----------------------------------------------------------------
 */

void *CVBandPrecAlloc(void *cvode_mem, long int N, 
                      long int mu, long int ml);

/*
 * -----------------------------------------------------------------
 * Function : CVBPSpgmr                                           
 * -----------------------------------------------------------------
 * CVBPSpgmr links the CVBANDPPRE preconditioner to the CVSPGMR   
 * linear solver. It performs the following actions:              
 *  1) Calls the CVSPGMR specification routine and attaches the   
 *     CVSPGMR linear solver to the integrator memory;            
 *  2) Sets the preconditioner data structure for CVSPGMR         
 *  3) Sets the preconditioner setup routine for CVSPGMR          
 *  4) Sets the preconditioner solve routine for CVSPGMR          
 *                                                                
 * Its first 3 arguments are the same as for CVSpgmr (see         
 * cvspgmr.h). The last argument is the pointer to the CVBANDPPRE 
 * memory block returned by CVBandPrecAlloc.                      
 * Note that the user need not call CVSpgmr anymore.              
 *                                                                
 * Possible return values are:                                    
 *                   SUCCESS                                      
 *                   LIN_NO_MEM                                   
 *                   LMEM_FAIL                                    
 *                   LIN_NO_LMEM                                  
 *                   LIN_ILL_INPUT                                
 *   Additionaly, if CVBandPrecAlloc was not previously called,   
 *   CVBPSpgmr returns BP_NO_PDATA (defined below).               
 * -----------------------------------------------------------------
 */

int CVBPSpgmr(void *cvode_mem, int pretype, int maxl, void *p_data);

/*
 * -----------------------------------------------------------------
 * Function : CVBandPrecFree                                      
 * -----------------------------------------------------------------
 * CVBandPrecFree frees the memory allocated by CVBandPrecAlloc   
 * in the argument pdata.                                         
 * -----------------------------------------------------------------
 */

void CVBandPrecFree(void *bp_data);

/*
 * -----------------------------------------------------------------
 * Function : CVBandPrecGet*                                      
 * -----------------------------------------------------------------
 * CVBandPrecGetIntWorkSpace returns the integer workspace used   
 *    by CVBANDPRE.                                               
 * CVBandPrecGetRealWorkSpace returns the real workspace used     
 *    by CVBANDPRE.                                               
 * CVBandPrecGetNumRhsEvals returns the number of calls made from 
 *    CVBANDPRE to the user's right hand side routine f.          
 * -----------------------------------------------------------------
 */

int CVBandPrecGetIntWorkSpace(void *bp_data, long int *leniwBP);
int CVBandPrecGetRealWorkSpace(void *bp_data, long int *lenrwBP);
int CVBandPrecGetNumRhsEvals(void *bp_data, long int *nfevalsBP);

/* Return values for CVBandPrecGet* functions */
/* OKAY = 0 */
enum { BP_NO_PDATA = -11 };

#endif

#ifdef __cplusplus
}
#endif
