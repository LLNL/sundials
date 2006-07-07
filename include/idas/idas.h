/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-07-07 16:49:25 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see the LICENSE file
 * -----------------------------------------------------------------
 * This is the header (include) file for the main IDAS solver.
 * -----------------------------------------------------------------
 *
 * IDAS is used to solve numerically the initial value problem     
 * for the differential algebraic equation (DAE) system           
 *   F(t,y,y') = 0,                                               
 * given initial conditions                                       
 *   y(t0) = y0,   y'(t0) = yp0.                                  
 * Here y and F are vectors of length N.                          
 *
 * Additionally, IDAS can perform forward or adjoint sensitivity
 * analysis.
 * -----------------------------------------------------------------
 */

#ifndef _IDAS_H
#define _IDAS_H

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              I D A S     C O N S T A N T S
   * =================================================================
   */

  /*
   * ----------------------------------------------------------------
   * Inputs to:
   *  IDAMalloc, IDAReInit, 
   *  IDASensMalloc, IDASensReInit, 
   *  IDAQuadMalloc, IDAQuadReInit,
   *  IDACalcIC, IDASolve,
   *  IDAadjMalloc
   * ----------------------------------------------------------------
   */

  /* itol */
#define IDA_SS               1
#define IDA_SV               2
#define IDA_WF               3
#define IDA_EE               4 

  /* itask */
#define IDA_NORMAL           1
#define IDA_ONE_STEP         2
#define IDA_NORMAL_TSTOP     3 
#define IDA_ONE_STEP_TSTOP   4

  /* icopt */
#define IDA_YA_YDP_INIT      1 
#define IDA_Y_INIT           2

  /* ism */
#define IDA_SIMULTANEOUS 1
#define IDA_STAGGERED    2

  /* interp */
#define IDA_HERMITE    1
#define IDA_POLYNOMIAL 2

  /*
   * ===============================================================
   * IDAS RETURN VALUES
   * ===============================================================
   */

#define IDA_SUCCESS          0
#define IDA_TSTOP_RETURN     1
#define IDA_ROOT_RETURN      2

#define IDA_WARNING          99

#define IDA_MEM_NULL        -1
#define IDA_ILL_INPUT       -2
#define IDA_NO_MALLOC       -3
#define IDA_TOO_MUCH_WORK   -4
#define IDA_TOO_MUCH_ACC    -5
#define IDA_ERR_FAIL        -6
#define IDA_CONV_FAIL       -7
#define IDA_LINIT_FAIL      -8
#define IDA_LSETUP_FAIL     -9
#define IDA_LSOLVE_FAIL     -10
#define IDA_RES_FAIL        -11
#define IDA_CONSTR_FAIL     -12
#define IDA_REP_RES_ERR     -13

#define IDA_MEM_FAIL        -14

#define IDA_BAD_T           -15

#define IDA_BAD_EWT         -16
#define IDA_FIRST_RES_FAIL  -17
#define IDA_LINESEARCH_FAIL -18
#define IDA_NO_RECOVERY     -19

#define IDA_RTFUNC_FAIL     -20

#define IDA_NO_QUAD         -30
#define IDA_QRHS_FAIL       -31
#define IDA_FIRST_QRHS_ERR  -32
#define IDA_REP_QRHS_ERR    -33

#define IDA_BAD_IS          -40
#define IDA_NO_SENS         -41
#define IDA_SRES_FAIL       -42
#define IDA_REP_SRES_ERR    -43

#define IDA_ADJMEM_NULL -101
#define IDA_BAD_TB0     -103
#define IDA_BCKMEM_NULL -104
#define IDA_REIFWD_FAIL -105
#define IDA_FWD_FAIL    -106
#define IDA_BAD_ITASK   -107
#define IDA_BAD_TBOUT   -108
#define IDA_GETY_BADT   -109

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * ----------------------------------------------------------------
   * Type : IDAResFn                                                   
   * ----------------------------------------------------------------
   * The F function which defines the DAE system   F(t,y,y')=0      
   * must have type IDAResFn.                                          
   * Symbols are as follows: 
   *                  t  <-> t        y <-> yy               
   *                  y' <-> yp       F <-> rr
   * A IDAResFn takes as input the independent variable value t,    
   * the dependent variable vector yy, and the derivative (with     
   * respect to t) of the yy vector, yp.  It stores the result of   
   * F(t,y,y') in the vector rr. The yy, yp, and rr arguments are of 
   * type N_Vector. The res_data parameter is the pointer res_data 
   * passed by the user to the IDASetRdata routine. This user-supplied 
   * pointer is passed to the user's res function every time it is called, 
   * to provide access in res to user data.                                    
   *                                                                
   * A IDAResFn res should return a value of 0 if successful, a positive
   * value if a recoverable error occured (e.g. yy has an illegal value),
   * or a negative value if a nonrecoverable error occured. In the latter
   * case, the program halts. If a recoverable error occured, the integrator
   * will attempt to correct and retry.
   * ----------------------------------------------------------------
   */

  typedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp,
                          N_Vector rr, void *res_data);

  /*
   * -----------------------------------------------------------------
   * Type : IDARootFn
   * -----------------------------------------------------------------
   * A function g, which defines a set of functions g_i(t,y,y') whose
   * roots are sought during the integration, must have type IDARootFn.
   * The function g takes as input the independent variable value t,
   * the dependent variable vector y, and its t-derivative yp (= y').
   * It stores the nrtfn values g_i(t,y,y') in the realtype array gout.
   * (Allocation of memory for gout is handled within IDA.)
   * The g_data parameter is the same as that passed by the user
   * to the IDARootInit routine.  This user-supplied pointer is
   * passed to the user's g function every time it is called.
   *
   * An IDARootFn should return 0 if successful or a non-zero value
   * if an error occured (in which case the integration will be halted).
   * -----------------------------------------------------------------
   */

  typedef int (*IDARootFn)(realtype t, N_Vector y, N_Vector yp,
                           realtype *gout, void *g_data);

  /*
   * -----------------------------------------------------------------
   * Type : IDAEwtFn
   * -----------------------------------------------------------------
   * A function e, which sets the error weight vector ewt, must have
   * type IDAEwtFn.
   * The function e takes as input the current dependent variable y.
   * It must set the vector of error weights used in the WRMS norm:
   * 
   *   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]
   *
   * Typically, the vector ewt has components:
   * 
   *   ewt_i = 1 / (reltol * |y_i| + abstol_i)
   *
   * The e_data parameter is the same as that passed by the user
   * to the IDASetEdata routine.  This user-supplied pointer is
   * passed to the user's e function every time it is called.
   * An IDAEwtFn e must return 0 if the error weight vector has been
   * successfuly set and a non-zero value otherwise.
   * -----------------------------------------------------------------
   */

  typedef int (*IDAEwtFn)(N_Vector y, N_Vector ewt, void *e_data);

  /*
   * -----------------------------------------------------------------
   * Type : IDAErrHandlerFn
   * -----------------------------------------------------------------
   * A function eh, which handles error messages, must have type
   * IDAErrHandlerFn.
   * The function eh takes as input the error code, the name of the
   * module reporting the error, the error message, and a pointer to
   * user data, the same as that passed to CVodeSetErrHandlerFn.
   * 
   * All error codes are negative, except IDA_WARNING which indicates 
   * a warning (the solver continues).
   *
   * An IDAErrHandlerFn has no return value.
   * -----------------------------------------------------------------
   */

  typedef void (*IDAErrHandlerFn)(int error_code, 
                                  const char *module, const char *function, 
                                  char *msg, void *eh_data); 


  /*
   * -----------------------------------------------------------------
   * Type : IDAQuadRhsFn
   * -----------------------------------------------------------------
   */
  
  typedef int (*IDAQuadRhsFn)(realtype tres, 
                              N_Vector yy, N_Vector yp,
                              N_Vector ypQ,
                              void *rdataQ);

  /*
   * -----------------------------------------------------------------
   * Type : IDASensResFn
   * -----------------------------------------------------------------
   */

  typedef int (*IDASensResFn)(int Ns, realtype tres, 
                              N_Vector yy, N_Vector yp, N_Vector resval,
                              N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                              void *rdataS,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /*
   * -----------------------------------------------------------------
   * Type : IDAResFnB
   * -----------------------------------------------------------------
   */

  typedef int (*IDAResFnB)(realtype tt, 
                           N_Vector yy, N_Vector yp,
                           N_Vector yyB, N_Vector ypB, N_Vector rrB,
                           void *rdataB);

  /*
   * -----------------------------------------------------------------
   * Type : IDAQuadRhsFnB
   * -----------------------------------------------------------------
   */

  typedef void (*IDAQuadRhsFnB)(realtype tt, 
                                N_Vector yy, N_Vector yp, 
                                N_Vector yyB, N_Vector ypB,
                                N_Vector ypQB, void *rdataQB);


  /*
   * ================================================================
   *          U S E R - C A L L A B L E   R O U T I N E S           
   * ================================================================
   */

  /* 
   * ----------------------------------------------------------------
   * Function : IDACreate                                           
   * ----------------------------------------------------------------
   * IDACreate creates an internal memory block for a problem to    
   * be solved by IDA.                                              
   *                                                                
   * If successful, IDACreate returns a pointer to initialized      
   * problem memory. This pointer should be passed to IDAMalloc.    
   * If an initialization error occurs, IDACreate prints an error   
   * message to standard err and returns NULL.                      
   *                                                                
   * ----------------------------------------------------------------
   */

  void *IDACreate(void);

  /*
   * ----------------------------------------------------------------
   * Integrator optional input specification functions              
   * ----------------------------------------------------------------
   * The following functions can be called to set optional inputs   
   * to values other than the defaults given below:                 
   *                                                                
   *                      |                                         
   * Function             |  Optional input / [ default value ]     
   *                      |                                          
   * ---------------------------------------------------------------- 
   *                      |                                          
   * IDASetErrHandlerFn   | user-provided ErrHandler function.
   *                      | [internal]
   *                      |
   * IDASetErrFile        | the file pointer for an error file
   *                      | where all CVODE warning and error
   *                      | messages will be written if the default
   *                      | internal error handling function is used. 
   *                      | This parameter can be stdout (standard 
   *                      | output), stderr (standard error), or a 
   *                      | file pointer (corresponding to a user 
   *                      | error file opened for writing) returned 
   *                      | by fopen.
   *                      | If not called, then all messages will
   *                      | be written to the standard error stream.
   *                      | [stderr]
   *                      |                                          
   * IDASetRdata          | a pointer to user data that will be     
   *                      | passed to the user's res function every 
   *                      | time res is called.                     
   *                      | [NULL]                                  
   *                      |         
   * IDASetEwtFn          | user-provide EwtSet function e and 
   *                      | a pointer to user data that will be
   *                      | passed to the user's e function every
   *                      | time e is called.
   *                      | [NULL]
   *                      | [NULL]
   *                      |
   * IDASetMaxOrd         | maximum lmm order to be used by the     
   *                      | solver.                                 
   *                      | [5]                                      
   *                      |                                          
   * IDASetMaxNumSteps    | maximum number of internal steps to be  
   *                      | taken by the solver in its attempt to   
   *                      | reach tout.                             
   *                      | [500]                                   
   *                      |                                          
   * IDASetInitStep       | initial step size.                      
   *                      | [estimated by IDA]                       
   *                      |                                          
   * IDASetMaxStep        | maximum absolute value of step size     
   *                      | allowed.                                
   *                      | [infinity]                              
   *                      |                                          
   * IDASetStopTime       | the independent variable value past     
   *                      | which the solution is not to proceed.   
   *                      | [infinity]                              
   *                      |                                          
   * IDASetNonlinConvCoef | Newton convergence test  constant       
   *                      | for use during integration.             
   *                      | [0.33]                                  
   *                      |                                          
   * IDASetMaxErrTestFails| Maximum number of error test failures   
   *                      | in attempting one step.                 
   *                      | [10]                                    
   *                      |                                         
   * IDASetMaxNonlinIters | Maximum number of nonlinear solver      
   *                      | iterations at one solution.             
   *                      | [4]                                     
   *                      |                                         
   * IDASetMaxConvFails   | Maximum number of allowable conv.       
   *                      | failures in attempting one step.        
   *                      | [10]                                    
   *                      |                                         
   * IDASetSuppressAlg    | flag to indicate whether or not to      
   *                      | suppress algebraic variables in the     
   *                      | local error tests:                      
   *                      | FALSE = do not suppress;                 
   *                      | TRUE = do suppress;                     
   *                      | [FALSE]                                 
   *                      | NOTE: if suppressed algebraic variables 
   *                      | is selected, the nvector 'id' must be   
   *                      | supplied for identification of those    
   *                      | algebraic components (see IDASetId).    
   *                      |                                          
   * IDASetId             | an N_Vector, which states a given       
   *                      | element to be either algebraic or       
   *                      | differential.                           
   *                      | A value of 1.0 indicates a differential 
   *                      | variable while a 0.0 indicates an       
   *                      | algebraic variable. 'id' is required    
   *                      | if optional input SUPPRESSALG is set,   
   *                      | or if IDACalcIC is to be called with    
   *                      | icopt = IDA_YA_YDP_INIT.               
   *                      |                                         
   * IDASetConstraints    | an N_Vector defining inequality         
   *                      | constraints for each component of the   
   *                      | solution vector y. If a given element   
   *                      | of this vector has values +2 or -2,     
   *                      | then the corresponding component of y   
   *                      | will be constrained to be > 0.0 or      
   *                      | <0.0, respectively, while if it is +1   
   *                      | or -1, the y component is constrained   
   *                      | to be >= 0.0 or <= 0.0, respectively.   
   *                      | If a component of constraints is 0.0,   
   *                      | then no constraint is imposed on the    
   *                      | corresponding component of y.           
   *                      | The presence of a non-NULL constraints  
   *                      | vector that is not 0.0 (ZERO) in all    
   *                      | components will cause constraint        
   *                      | checking to be performed.               
   *                      |                                         
   * -----------------------------------------------------------------
   *                      |
   * IDASetTolerances     | Changes the integration tolerances
   *                      | between calls to IDASolve().
   *                      | [set by IDAMalloc/IDAReInit]
   *                      |
   * ---------------------------------------------------------------- 
   * Return flag:
   *   IDA_SUCCESS   if successful
   *   IDA_MEM_NULL  if the IDAS memory is NULL
   *   IDA_ILL_INPUT if an argument has an illegal value
   *
   * ----------------------------------------------------------------
   */

  int IDASetErrHandlerFn(void *ida_mem, IDAErrHandlerFn ehfun, void *eh_data);
  int IDASetErrFile(void *ida_mem, FILE *errfp);
  int IDASetRdata(void *ida_mem, void *res_data);
  int IDASetEwtFn(void *ida_mem, IDAEwtFn efun, void *edata);
  int IDASetMaxOrd(void *ida_mem, int maxord);
  int IDASetMaxNumSteps(void *ida_mem, long int mxsteps);
  int IDASetInitStep(void *ida_mem, realtype hin);
  int IDASetMaxStep(void *ida_mem, realtype hmax);
  int IDASetStopTime(void *ida_mem, realtype tstop);
  int IDASetNonlinConvCoef(void *ida_mem, realtype epcon);
  int IDASetMaxErrTestFails(void *ida_mem, int maxnef);
  int IDASetMaxNonlinIters(void *ida_mem, int maxcor);
  int IDASetMaxConvFails(void *ida_mem, int maxncf);
  int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg);
  int IDASetId(void *ida_mem, N_Vector id);
  int IDASetConstraints(void *ida_mem, N_Vector constraints);

  int IDASetTolerances(void *ida_mem, int itol, realtype rtol, void *atol);

  /*
   * ----------------------------------------------------------------
   * Function : IDAMalloc                                           
   * ----------------------------------------------------------------
   * IDAMalloc allocates and initializes memory for a problem to    
   * to be solved by IDAS.                                           
   *                                                                
   * res     is the residual function F in F(t,y,y') = 0.                     
   *                                                                
   * t0      is the initial value of t, the independent variable.   
   *                                                                
   * yy0     is the initial condition vector y(t0).                 
   *                                                                
   * yp0     is the initial condition vector y'(t0)                 
   *                                                                
   * itol    is the type of tolerances to be used.                  
   *            The legal values are:                               
   *               IDA_SS (scalar relative and absolute  tolerances),   
   *               IDA_SV (scalar relative tolerance and vector         
   *                       absolute tolerance).                         
   *               IDA_WF (user-provided weight function)                       
   *                                         
   * rtol    is the relative tolerance scalar.         
   *                                                                
   * atol    is a pointer (void) to the absolute tolerance scalar or
   *            an N_Vector tolerance or an IDAEwtFn funciton.                              
   * (ewt)                                                          
   *         Both rtol and atol are used to compute the error weight
   *         vector, ewt. The error test required of a correction   
   *         delta is that the weighted-RMS norm of delta be less   
   *         than or equal to 1.0. Other convergence tests use the  
   *         same norm. The weighting vector used in this norm is   
   *         ewt. The components of ewt are defined by              
   *         ewt[i] = 1.0/(rtol*yy[i] + atol[i]). Here, yy is the   
   *         current approximate solution.  See the routine         
   *         N_VWrmsNorm for the norm used in this error test.      
   *                                                                
   * Note: The tolerance values may be changed in between calls to  
   *       IDASolve for the same problem. These values refer to     
   *       (*rtol) and either (*atol), for a scalar absolute        
   *       tolerance, or the components of atol, for a vector       
   *       absolute tolerance.                                      
   *                                                                 
   *  IDA_SUCCESS if successful
   *  IDA_MEM_NULL if the IDAS memory was NULL
   *  IDA_MEM_FAIL if a memory allocation failed
   *  IDA_ILL_INPUT f an argument has an illegal value.
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDAMalloc(void *ida_mem, IDAResFn res,
                realtype t0, N_Vector yy0, N_Vector yp0, 
                int itol, realtype rtol, void *atol);

  /*
   * ----------------------------------------------------------------
   * Function : IDAReInit                                           
   * ----------------------------------------------------------------
   * IDAReInit re-initializes IDAS for the solution of a problem,    
   * where a prior call to IDAMalloc has been made.                 
   * IDAReInit performs the same input checking and initializations 
   * that IDAMalloc does.                                           
   * But it does no memory allocation, assuming that the existing   
   * internal memory is sufficient for the new problem.             
   *                                                                
   * The use of IDAReInit requires that the maximum method order,   
   * maxord, is no larger for the new problem than for the problem  
   * specified in the last call to IDAMalloc.  This condition is    
   * automatically fulfilled if the default value for maxord is     
   * specified.                                                     
   *                                                                
   * Following the call to IDAReInit, a call to the linear solver   
   * specification routine is necessary if a different linear solver
   * is chosen, but may not be otherwise.  If the same linear solver
   * is chosen, and there are no changes in its input parameters,   
   * then no call to that routine is needed.                        
   *                                                                
   * The first argument to IDAReInit is:                            
   *                                                                
   * ida_mem = pointer to IDA memory returned by IDACreate.         
   *                                                                
   * All the remaining arguments to IDAReInit have names and        
   * meanings identical to those of IDAMalloc.                      
   *                                                                
   * The return value of IDAReInit is equal to SUCCESS = 0 if there 
   * were no errors; otherwise it is a negative int equal to:       
   *   IDA_MEM_NULL   indicating ida_mem was NULL, or            
   *   IDA_NO_MALLOC  indicating that ida_mem was not allocated. 
   *   IDA_ILL_INPUT  indicating an input argument was illegal   
   *                  (including an attempt to increase maxord). 
   * In case of an error return, an error message is also printed.  
   * ----------------------------------------------------------------
   */                                                                

  int IDAReInit(void *ida_mem, IDAResFn res,
                realtype t0, N_Vector yy0, N_Vector yp0,
                int itol, realtype rtol, void *atol);
 
  /* ----------------------------------------------------------------
   * Initial Conditions optional input specification functions      
   * ----------------------------------------------------------------
   * The following functions can be called to set optional inputs   
   * to control the initial conditions calculations.                
   *                                                                
   *                        |                                        
   * Function               |  Optional input / [ default value ]   
   *                        |                                        
   * -------------------------------------------------------------- 
   *                        |                                        
   * IDASetNonlinConvCoefIC | positive coeficient in the Newton     
   *                        | convergence test.  This test uses a   
   *                        | weighted RMS norm (with weights       
   *                        | defined by the tolerances, as in      
   *                        | IDASolve).  For new initial value     
   *                        | vectors y and y' to be accepted, the  
   *                        | norm of J-inverse F(t0,y,y') is       
   *                        | required to be less than epiccon,     
   *                        | where J is the system Jacobian.       
   *                        | [0.01 * 0.33]                          
   *                        |                                        
   * IDASetMaxNumStepsIC    | maximum number of values of h allowed 
   *                        | when icopt = IDA_YA_YDP_INIT, where  
   *                        | h appears in the system Jacobian,     
   *                        | J = dF/dy + (1/h)dF/dy'.              
   *                        | [5]                                   
   *                        |                                        
   * IDASetMaxNumJacsIC     | maximum number of values of the       
   *                        | approximate Jacobian or preconditioner
   *                        | allowed, when the Newton iterations   
   *                        | appear to be slowly converging.       
   *                        | [4]                                    
   *                        |                                        
   * IDASetMaxNumItersIC    | maximum number of Newton iterations   
   *                        | allowed in any one attempt to solve   
   *                        | the IC problem.                       
   *                        | [10]                                  
   *                        |                                        
   * IDASetLineSearchOffIC  | a boolean flag to turn off the        
   *                        | linesearch algorithm.                 
   *                        | [FALSE]                               
   *                        |                                        
   * IDASetStepToleranceIC  | positive lower bound on the norm of   
   *                        | a Newton step.                        
   *                        | [(unit roundoff)^(2/3)                
   *                                                                
   * ---------------------------------------------------------------- 
   * Return flag:
   *   IDA_SUCCESS   if successful
   *   IDA_MEM_NULL  if the IDAS memory is NULL
   *   IDA_ILL_INPUT if an argument has an illegal value
   *
   * ----------------------------------------------------------------
   */

  int IDASetNonlinConvCoefIC(void *ida_mem, realtype epiccon);
  int IDASetMaxNumStepsIC(void *ida_mem, int maxnh);
  int IDASetMaxNumJacsIC(void *ida_mem, int maxnj);
  int IDASetMaxNumItersIC(void *ida_mem, int maxnit);
  int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff);
  int IDASetStepToleranceIC(void *ida_mem, realtype steptol);

  /*
   * -----------------------------------------------------------------
   * Function : IDARootInit
   * -----------------------------------------------------------------
   * IDARootInit initializes a rootfinding problem to be solved
   * during the integration of the DAE system.  It must be called
   * after IDACreate, and before IDASolve.  The arguments are:
   *
   * ida_mem = pointer to IDA memory returned by IDACreate.
   *
   * nrtfn   = number of functions g_i, an int >= 0.
   *
   * g       = name of user-supplied function, of type IDARootFn,
   *           defining the functions g_i whose roots are sought.
   *
   * g_data  = a pointer to user data that will be passed to the 
   *           user's g function every time g is called.
   *
   * If a new problem is to be solved with a call to IDAReInit,
   * where the new problem has no root functions but the prior one
   * did, then call IDARootInit with nrtfn = 0.
   *
   * The return value of IDARootInit is IDA_SUCCESS = 0 if there were
   * no errors; otherwise it is a negative int equal to:
   *   IDA_MEM_NULL     indicating ida_mem was NULL, or
   *   IDA_MEM_FAIL     indicating a memory allocation failed.
   *                    (including an attempt to increase maxord).
   *   IDA_ILL_INPUT    indicating nrtfn > 0 but g = NULL.
   * In case of an error return, an error message is also printed.
   * -----------------------------------------------------------------
   */

  int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g, void *g_data);

  /*
   * -----------------------------------------------------------------
   * Quadrature optional input specification functions
   * -----------------------------------------------------------------
   * The following functions can be called to set optional inputs
   * to values other than the defaults given below:
   *
   * Function             |  Optional input / [ default value ]
   * --------------------------------------------------------------
   *                      |
   * IDASetQuadRdata      | a pointer to user data that will be
   *                      | passed to the user's rhsQ function every
   *                      | time rhsQ is called.
   *                      | [NULL]
   *                      |
   * IDASetQuadErrCon     | are quadrature variables considered in
   *                      | the error control?
   *                      | If yes, set tolerances for quadrature
   *                      | integration. 
   *                      | [errconQ = FALSE]
   *                      | [ not tolerances]
   *                      |
   * -----------------------------------------------------------------
   * If successful, these functions return IDA_SUCCESS. If an argument
   * has an illegal value, they print an error message to the
   * file specified by errfp and return one of the error flags
   * defined for the CVodeSet* routines.
   * -----------------------------------------------------------------
   */

  int IDASetQuadRdata(void *ida_mem, void *rhsQ_data);
  int IDASetQuadErrCon(void *ida_mem, booleantype errconQ, 
                       int itolQ, realtype reltolQ, void *abstolQ);

  /*
   * ----------------------------------------------------------------
   * Function : IDAQuadMalloc                                       
   * ----------------------------------------------------------------
   * IDAQuadMalloc allocates and initializes memory related to      
   * quadrature integration.                                        
   *                                                                
   * ida_mem is a pointer to IDAS memory returned by IDACreate      
   *                                                                
   * rhsQ  is the user-provided integrand routine.                  
   *                                                                
   * yQ0   is a pointer to a vector specification structure         
   *       for N_Vectors containing quadrature variables.           
   *                                                                
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDAQuadMalloc(void *ida_mem, IDAQuadRhsFn rhsQ, N_Vector yQ0);

  /* 
   * ----------------------------------------------------------------
   * Function : IDAQuadReInit                                       
   * ----------------------------------------------------------------
   * IDAQuadReInit re-initializes IDAS's quadrature related         
   * memory for a problem, assuming it has already been allocated   
   * in prior calls to IDAMalloc and IDAQuadMalloc.                 
   *                                                                
   * All problem specification inputs are checked for errors.       
   * The number of quadratures Nq is assumed to be unchanged        
   * since the previous call to IDAQuadMalloc.                      
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDAQuadReInit(void *ida_mem, IDAQuadRhsFn rhsQ, N_Vector yQ0);


  /* 
   * ----------------------------------------------------------------
   * Forward sensitivity optional input specification functions     
   * ----------------------------------------------------------------
   * The following functions can be called to set optional inputs   
   * to other values than the defaults given below:                 
   *                                                                
   * Function                 |  Optional input / [ default value ]     
   *                          |                                          
   * -------------------------------------------------------------- 
   *                          |                                         
   * IDASetSensResFn          | sensitivity residual function and
   *                          | user data pointer.
   *                          | This function must compute residuals    
   *                          | for all sensitivity equations.          
   *                          | [IDAS difference quotient approx.]      
   *                          | [internal]
   *                          |                                         
   * IDASetSensRho            | controls the selection of finite        
   *                          | difference schemes used in evaluating   
   *                          | the sensitivity right hand sides.       
   *                          | [0.0]                                   
   *                          |                                         
   * IDASetSensParams         |   parameter information:
   *                          | p: pointer to problem parameters
   *                          | plist: list of parameters with respect
   *                          |        to which sensitivities are to be
   *                          |        computed.
   *                          | pbar: order of magnitude info. 
   *                          |       Typically, if p[plist[i]] is nonzero, 
   *                          |       pbar[i]=p[plist[i]].
   *                          | [p=NULL]
   *                          | [plist=NULL]
   *                          | [pbar=NULL]                               
   *                          |                                         
   * IDASetSensErrCon         | are sensitivity variables considered in 
   *                          | the error control?                      
   *                          | [TRUE]                                  
   *                          |                                         
   * IDASetSensTolerances     | a pointer to the sensitivity relative   
   *                          | tolerance scalar and one for the        
   *                          | absolute tolerance                      
   *                          | [estimated by IDAS]                     
   *                          |                                         
   * IDASetSensMaxNonlinIters | Maximum number of nonlinear solver  
   *                          | iterations for sensitivity systems  
   *                          | (staggered)                         
   *                          | [4]                                 
   *                          |                                     
   * -------------------------------------------------------------- 
   * If successful, these functions return IDA_SUCCESS. If an argument  
   * has an illegal value, they return one of the error flags      
   * defined for the IDASet* routines.                              
   * ----------------------------------------------------------------
   */

  int IDASetSensResFn(void *ida_mem, IDASensResFn resS, void *rdataS);
  int IDASetSensRho(void *ida_mem, realtype rho);
  int IDASetSensParams(void *ida_mem, realtype *p, realtype *pbar, int *plist);
  int IDASetSensErrCon(void *ida_mem, booleantype errconS);
  int IDASetSensTolerances(void *ida_mem, int itolS, 
                           realtype reltolS, void *abstolS);
  int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS);
  
  /*
   * ----------------------------------------------------------------
   * Function : IDASensMalloc                                       
   * ----------------------------------------------------------------
   * IDASensMalloc allocates and initializes memory related to      
   * sensitivity computations.                                      
   *                                                                
   * ida_mem is a pointer to IDAS memory returned by IDACreate.     
   *                                                                
   * Ns        is the number of sensitivities to be computed.       
   *                                                                
   * ism       is the type of corrector used in sensitivity         
   *           analysis. The legal values are: SIMULTANEOUS
   *           and STAGGERED (see previous description) 
   *                                                                
   * yS0       is the array of initial condition vectors for        
   *           sensitivity variables.                                
   *                                                                
   * ypS0      is the array of initial condition vectors for        
   *           sensitivity derivatives.                              
   *                                                                
   * If successful, IDASensMalloc returns SUCCESS. If an            
   * initialization error occurs, IDASensMalloc returns one of      
   * the error flags defined above.                                 
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDASensMalloc(void *ida_mem, int Ns, int ism, N_Vector *yS0, N_Vector *ypS0);
  
  
  /*
   * ----------------------------------------------------------------
   * Function : IDASensReInit                                       
   * ----------------------------------------------------------------
   * IDASensReInit re-initializes the IDAS sensitivity related      
   * memory for a problem, assuming it has already been allocated   
   * in prior calls to IDAMalloc and IDASensMalloc.                 
   *                                                                
   * All problem specification inputs are checked for errors.       
   * The number of sensitivities Ns is assumed to be unchanged      
   * since the previous call to IDASensMalloc.                      
   * If any error occurs during initialization, it is reported to   
   * the file whose file pointer is errfp.                          
   *                                                                
   * IDASensReInit potentially does some minimal memory allocation  
   * (for the sensitivity absolute tolerance).
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDASensReInit(void *ida_mem, int ism, N_Vector *yS0, N_Vector *ypS0);

  /*
   * -----------------------------------------------------------------
   * Function : IDASensToggle
   * -----------------------------------------------------------------
   * IDASensToggle activates or deactivates sensitivity calculations.
   * It does NOT deallocate sensitivity-related memory.
   * It is allowed to set sensi=TRUE only if IDASensMalloc has been
   * previously called.
   * 
   * The return value is equal to IDA_SUCCESS = 0 if there were no
   * errors; otherwise it is a negative int equal to:
   *   IDa_MEM_NULL  indicating ida_mem was NULL
   *   IDA_NO_SENS   indicating there was not a prior call to
   *                 IDASensMalloc.   
   * -----------------------------------------------------------------
   */

  int IDASensToggle(void *ida_mem, booleantype sensi);


  /*
   * ----------------------------------------------------------------
   * Function : IDACalcIC                                           
   * ----------------------------------------------------------------
   * IDACalcIC calculates corrected initial conditions for the DAE  
   * system for a class of index-one problems of semi-implicit form.
   * It uses Newton iteration combined with a Linesearch algorithm. 
   * Calling IDACalcIC is optional. It is only necessary when the   
   * initial conditions do not solve the given system.  I.e., if    
   * y0 and yp0 are known to satisfy F(t0, y0, yp0) = 0, then       
   * a call to IDACalcIC is NOT necessary (for index-one problems). 
   *                                                                
   * A call to IDACalcIC must be preceded by a successful call to   
   * IDAMalloc or IDAReInit for the given DAE problem, and by a     
   * successful call to the linear system solver specification      
   * routine.                                                       
   * In addition, IDACalcIC assumes that the vectors y0, yp0, and   
   * (if relevant) id and constraints that were set through         
   * IDASetConstraints remain unaltered since that call.            
   *                                                                
   * The call to IDACalcIC should precede the call(s) to IDASolve   
   * for the given problem.                                         
   *                                                                
   * The arguments to IDACalcIC are as follows.  The first three -- 
   * ida_mem, icopt, tout1 -- are required; the others are optional.
   * A zero value passed for any optional input specifies that the  
   * default value is to be used.                                   
   *                                                                
   * IDA_mem is the pointer to IDA memory returned by IDACreate.    
   *                                                                
   * t0      is the initial value of t, the independent variable.   
   *                                                                
   * yy0     is the initial condition vector y(t0).                 
   *                                                                
   * yp0     is the initial condition vector y'(t0)                 
   *                                                                
   * icopt  is the option of IDACalcIC to be used.                  
   *        icopt = IDA_YA_YDP_INIT   directs IDACalcIC to compute 
   *                the algebraic components of y and differential  
   *                components of y', given the differential        
   *                components of y.  This option requires that the 
   *                N_Vector id was set through a call to IDASetId  
   *                specifying the differential and algebraic       
   *                components.                                     
   *        icopt = IDA_Y_INIT   directs IDACalcIC to compute all  
   *                components of y, given y'.  id is not required. 
   *                                                                
   * tout1  is the first value of t at which a soluton will be      
   *        requested (from IDASolve).  (This is needed here to     
   *        determine the direction of integration and rough scale  
   *        in the independent variable t.)                          
   *                                                                
   *                                                                
   * IDACalcIC returns an int flag.  Its symbolic values and their  
   * meanings are as follows.  (The numerical return values are set 
   * above in this file.)  All unsuccessful returns give a negative 
   * return value.  If IFACalcIC failed, y0 and yp0 contain         
   * (possibly) altered values, computed during the attempt.        
   *                                                                
   * SUCCESS             IDACalcIC was successful.  The corrected   
   *                     initial value vectors are in y0 and yp0.    
   *                                                                
   * IDA_MEM_NULL        The argument ida_mem was NULL.             
   *                                                                
   * IDA_ILL_INPUT       One of the input arguments was illegal.    
   *                     See printed message.                       
   *                                                                
   * IDA_LINIT_FAIL      The linear solver's init routine failed.   
   *                                                                
   * IDA_BAD_EWT         Some component of the error weight vector  
   *                     is zero (illegal), either for the input    
   *                     value of y0 or a corrected value.          
   *                                                                
   * IDA_RES_FAIL        The user's residual routine returned 
   *                     a non-recoverable error flag.              
   *                                                                
   * IDA_FIRST_RES_FAIL  The user's residual routine returned 
   *                     a recoverable error flag on the first call,
   *                     but IDACalcIC was unable to recover.       
   *                                                                
   * IDA_LSETUP_FAIL     The linear solver's setup routine had a    
   *                     non-recoverable error.                     
   *                                                                
   * IDA_LSOLVE_FAIL     The linear solver's solve routine had a    
   *                     non-recoverable error.                     
   *                                                                
   * IDA_NO_RECOVERY     The user's residual routine, or the linear 
   *                     solver's setup or solve routine had a      
   *                     recoverable error, but IDACalcIC was       
   *                     unable to recover.                         
   *                                                                
   * IDA_CONSTR_FAIL     IDACalcIC was unable to find a solution    
   *                     satisfying the inequality constraints.     
   *                                                                
   * IDA_LINESEARCH_FAIL The Linesearch algorithm failed to find a  
   *                     solution with a step larger than steptol   
   *                     in weighted RMS norm.                      
   *                                                                
   * IDA_CONV_FAIL       IDACalcIC failed to get convergence of the 
   *                     Newton iterations.                         
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDACalcIC (void *ida_mem, realtype t0, N_Vector yy0, N_Vector yp0, 
                 int icopt, realtype tout1); 

  /*
   * ----------------------------------------------------------------
   * Function : IDASolve                                            
   * ----------------------------------------------------------------
   * IDASolve integrates the DAE over an interval in t, the         
   * independent variable. If itask is IDA_NORMAL, then the solver      
   * integrates from its current internal t value to a point at or  
   * beyond tout, then interpolates to t = tout and returns y(tret) 
   * in the user-allocated vector yret. In general, tret = tout.    
   * If itask is IDA_ONE_STEP, then the solver takes one internal
   * step of the independent variable and returns in yret the value
   * of y at the new internal independent variable value. In this
   * case, tout is used only during the first call to IDASolve to         
   * determine the direction of integration and the rough scale of  
   * the problem. If itask is IDA_NORMAL_TSTOP or IDA_ONE_STEP_TSTOP,
   * then IDA returns the solution at tstop if that comes sooner than
   * tout or the end of the next internal step, respectively.
   * In any case, the independent variable value reached by the solver
   * is placed in (*tret).  The user is responsible for allocating the
   * memory for this value.          
   *                                                                
   * ida_mem is the pointer (void) to IDA memory returned by        
   *         IDACreate.
   *                                                                
   * tout    is the next independent variable value at which a      
   *         computed solution is desired.                          
   *                                                                
   * tret    is a pointer to a real location.  IDASolve sets (*tret)
   *         to the actual t value reached, corresponding to the
   *         solution vector yret.  In IDA_NORMAL mode, with no
   *         errors and no roots found, (*tret) = tout.
   *
   * yret    is the computed solution vector.  With no errors,
   *         yret = y(tret).                                        
   *                                                                
   * ypret   is the derivative of the computed solution at t = tret.
   *                                                                
   * Note: yret and ypret may be the same N_Vectors as y0 and yp0   
   * in the call to IDAMalloc or IDAReInit.                         
   *                                                                
   * itask   is IDA_NORMAL, IDA_NORMAL_TSTOP, IDA_ONE_STEP, or
   *         IDA_ONE_STEP_TSTOP.   These modes are described above.
   *
   *
   * The return values for IDASolve are described below.            
   * (The numerical return values are defined above in this file.)  
   * All unsuccessful returns give a negative return value.         
   *                                                                
   * IDA_SUCCESS
   *   IDASolve succeeded and no roots were found.                       
   *
   * IDA_ROOT_RETURN:  IDASolve succeeded, and found one or more roots.
   *   If nrtfn > 1, call IDAGetRootInfo to see which g_i were found
   *   to have a root at (*tret).
   *
   * IDA_TSTOP_RETURN: 
   *   IDASolve returns computed results for the independent variable 
   *   value tstop. That is, tstop was reached.                            
   *                                                                
   * IDA_MEM_NULL: 
   *   The ida_mem argument was NULL.            
   *                                                                
   * IDA_ILL_INPUT: 
   *   One of the inputs to IDASolve is illegal. This includes the 
   *   situation when a component of the error weight vectors 
   *   becomes < 0 during internal stepping.  It also includes the
   *   situation where a root of one of the root functions was found
   *   both at t0 and very near t0.  The ILL_INPUT flag          
   *   will also be returned if the linear solver function IDA---
   *   (called by the user after calling IDACreate) failed to set one 
   *   of the linear solver-related fields in ida_mem or if the linear 
   *   solver's init routine failed. In any case, the user should see 
   *   the printed error message for more details.                
   *                                                                
   * IDA_TOO_MUCH_WORK: 
   *   The solver took mxstep internal steps but could not reach tout. 
   *   The default value for mxstep is MXSTEP_DEFAULT = 500.                
   *                                                                
   * IDA_TOO_MUCH_ACC: 
   *   The solver could not satisfy the accuracy demanded by the user 
   *   for some internal step.   
   *                                                                
   * IDA_ERR_FAIL:
   *   Error test failures occurred too many times (=MXETF = 10) during 
   *   one internal step.  
   *                                                                
   * IDA_CONV_FAIL: 
   *   Convergence test failures occurred too many times (= MXNCF = 10) 
   *   during one internal step.                                          
   *                                                                
   * IDA_LSETUP_FAIL: 
   *   The linear solver's setup routine failed  
   *   in an unrecoverable manner.                    
   *                                                                
   * IDA_LSOLVE_FAIL: 
   *   The linear solver's solve routine failed  
   *   in an unrecoverable manner.                    
   *                                                                
   * IDA_CONSTR_FAIL:
   *    The inequality constraints were violated, 
   *    and the solver was unable to recover.         
   *                                                                
   * IDA_REP_RES_ERR: 
   *    The user's residual function repeatedly returned a recoverable 
   *    error flag, but the solver was unable to recover.                 
   *                                                                
   * IDA_RES_FAIL:
   *    The user's residual function returned a nonrecoverable error 
   *    flag.
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDASolve(void *ida_mem, realtype tout, realtype *tret,
               N_Vector yret, N_Vector ypret, int itask);

  /*
   * ----------------------------------------------------------------
   * Function: IDAGetSolution                                       
   * ----------------------------------------------------------------
   *                                                                
   * This routine evaluates y(t) and y'(t) as the value and         
   * derivative of the interpolating polynomial at the independent  
   * variable t, and stores the results in the vectors yret and     
   * ypret.  It uses the current independent variable value, tn,    
   * and the method order last used, kused. This function is        
   * called by IDASolve with t = tout, t = tn, or t = tstop.        
   *                                                                
   * If kused = 0 (no step has been taken), or if t = tn, then the  
   * order used here is taken to be 1, giving yret = phi[0],        
   * ypret = phi[1]/psi[0].                                         
   *                                                                
   * The return values are:                                         
   *   IDA_SUCCESS:  succeess.                                  
   *   IDA_BAD_T:    t is not in the interval [tn-hu,tn].                   
   *   IDA_MEM_NULL: The ida_mem argument was NULL.     
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDAGetSolution(void *ida_mem, realtype t, 
                     N_Vector yret, N_Vector ypret);

  /* ----------------------------------------------------------------
   * Integrator optional output extraction functions                
   * ----------------------------------------------------------------
   *                                                                
   * The following functions can be called to get optional outputs  
   * and statistics related to the main integrator.                 
   * ---------------------------------------------------------------- 
   *                                                                
   * IDAGetWorkSpace returns the IDA real and integer workspace sizes      
   * IDAGetNumSteps returns the cumulative number of internal       
   *       steps taken by the solver                                
   * IDAGetNumRhsEvals returns the number of calls to the user's    
   *       res function                                             
   * IDAGetNumLinSolvSetups returns the number of calls made to     
   *       the linear solver's setup routine                        
   * IDAGetNumErrTestFails returns the number of local error test   
   *       failures that have occured                               
   * IDAGetNumBacktrackOps returns the number of backtrack          
   *       operations done in the linesearch algorithm in IDACalcIC 
   * IDAGetLastOrder returns the order used during the last         
   *       internal step                                            
   * IDAGetCurentOrder returns the order to be used on the next     
   *       internal step                                            
   * IDAGetActualInitStep returns the actual initial step size      
   *       used by IDA                                              
   * IDAGetLAstStep returns the step size for the last internal     
   *       step (if from IDASolve), or the last value of the        
   *       artificial step size h (if from IDACalcIC)               
   * IDAGetCurrentStep returns the step size to be attempted on the 
   *       next internal step                                       
   * IDAGetCurrentTime returns the current internal time reached    
   *       by the solver                                            
   * IDAGetTolScaleFactor returns a suggested factor by which the   
   *       user's tolerances should be scaled when too much         
   *       accuracy has been requested for some internal step       
   * IDAGetErrWeights returns the current state error weight vector.        
   *       The user must allocate space for eweight.
   * IDAGetEstLocalErrors returns the estimated local errors. The user
   *       must allocate space for the vector ele.
   * IDAGetNumGEvals returns the number of calls to the user's
   *       g function (for rootfinding)
   * IDAGetRootInfo returns the indices for which g_i was found to 
   *       have a root. The user must allocate space for rootsfound.
   *       For i = 0 ... nrtfn-1, rootsfound[i] = 1 if g_i has a root,
   *       and rootsfound[i]= 0 if not.
   *                                                                
   * IDAGet* return values:
   *   IDA_SUCCESS   if succesful
   *   IDA_MEM_NULL  if the IDAS memory was NULL
   *
   * ----------------------------------------------------------------
   */

  int IDAGetWorkSpace(void *ida_mem, long int *lenrw, long int *leniw);
  int IDAGetNumSteps(void *ida_mem, long int *nsteps);
  int IDAGetNumResEvals(void *ida_mem, long int *nrevals);
  int IDAGetNumLinSolvSetups(void *ida_mem, long int *nlinsetups);
  int IDAGetNumErrTestFails(void *ida_mem, long int *netfails);
  int IDAGetNumBacktrackOps(void *ida_mem, long int *nbacktr);
  int IDAGetLastOrder(void *ida_mem, int *klast);
  int IDAGetCurrentOrder(void *ida_mem, int *kcur);
  int IDAGetActualInitStep(void *ida_mem, realtype *hinused);
  int IDAGetLastStep(void *ida_mem, realtype *hlast);
  int IDAGetCurrentStep(void *ida_mem, realtype *hcur);
  int IDAGetCurrentTime(void *ida_mem, realtype *tcur);
  int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact);
  int IDAGetErrWeights(void *ida_mem, N_Vector eweight);
  int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele);
  int IDAGetNumGEvals(void *ida_mem, long int *ngevals);
  int IDAGetRootInfo(void *ida_mem, int *rootsfound);

  /*
   * ----------------------------------------------------------------
   * As a convenience, the following function provides the          
   * optional outputs in a group.                                   
   * ----------------------------------------------------------------
   */

  int IDAGetIntegratorStats(void *ida_mem, long int *nsteps, 
                            long int *nrevals, long int *nlinsetups, 
                            long int *netfails, int *qlast, int *qcur, 
                            realtype *hlast, realtype *hcur, 
                            realtype *tcur);

  /*
   * ----------------------------------------------------------------
   * Nonlinear solver optional output extraction functions          
   * ----------------------------------------------------------------
   *                                                                
   * The following functions can be called to get optional outputs  
   * and statistics related to the nonlinear solver.                
   * -------------------------------------------------------------- 
   *                                                                
   * IDAGetNumNonlinSolvIters returns the number of nonlinear       
   *       solver iterations performed.                             
   * IDAGetNumNonlinSolvConvFails returns the number of nonlinear   
   *       convergence failures.                                    
   *                                                                
   * ----------------------------------------------------------------
   */

  int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters);
  int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails);

  /*
   * ----------------------------------------------------------------
   * As a convenience, the following function provides the          
   * nonlinear solver optional outputs in a group.                                   
   * ----------------------------------------------------------------
   */

  int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters, 
                            long int *nncfails);

  /*
   * -----------------------------------------------------------------
   * Quadrature integration solution extraction routine
   * -----------------------------------------------------------------
   * The following function can be called to obtain the quadrature
   * variables after a successful integration step.
   * If quadratures were not computed, it returns IDA_NO_QUAD.
   * -----------------------------------------------------------------
   */

  int IDAGetQuad(void *ida_mem, realtype t, N_Vector yQout);

  /*
   * -----------------------------------------------------------------
   * Quadrature integration optional output extraction routines
   * -----------------------------------------------------------------
   * The following functions can be called to get optional outputs
   * and statistics related to the integration of quadratures.
   * -----------------------------------------------------------------
   * IDAGetQuadNumRhsEvals returns the number of calls to the
   *                         user function fQ defining the right hand
   *                         side of the quadrature variables.
   * IDAGetQuadNumErrTestFails returns the number of local error
   *                             test failures for quadrature variables.
   * IDAGetQuadErrWeights returns the vector of error weights for
   *                        the quadrature variables. The user must
   *                        allocate space for ewtQ.
   * -----------------------------------------------------------------
   */

  int IDAGetQuadNumRhsEvals(void *ida_mem, long int *nrhsQevals);
  int IDAGetQuadNumErrTestFails(void *ida_mem, long int *nQetfails);
  int IDAGetQuadErrWeights(void *ida_mem, N_Vector eQweight);

  /*
   * -----------------------------------------------------------------
   * As a convenience, the following function provides the
   * optional outputs in a group.
   * -----------------------------------------------------------------
   */

  int IDAGetQuadStats(void *ida_mem, long int *nrhsQevals, long int *nQetfails);

  /*
   * -----------------------------------------------------------------
   * Sensitivity solution extraction routine
   * -----------------------------------------------------------------
   * The following functions can be called to obtain the sensitivity
   * variables after a successful integration step.
   * IDAGetSens returns all sensitivity vectors.
   * IDAGetSens1 returns only the 'is'-th sensitivity vectors.
   * If sensitivities were not computed, they returns IDA_NO_SENS.
   * -----------------------------------------------------------------
   */

  int IDAGetSens(void *ida_mem, realtype t, N_Vector *yySout, N_Vector *ypSout);
  int IDAGetSens1(void *ida_mem, realtype t, int is, N_Vector yySret, N_Vector ypSret);

  /*
   * -----------------------------------------------------------------
   * Forward sensitivity optional output extraction routines
   * -----------------------------------------------------------------
   * The following functions can be called to get optional outputs
   * and statistics related to the integration of sensitivities.
   * -----------------------------------------------------------------
   * IDAGetNumSensResEvals returns the number of calls to the
   *                       sensitivity residual function.
   * IDAGetNumResEvalsSens returns the number of calls to the
   *                       user res routine due to finite difference
   *                       evaluations of the sensitivity equations.
   * IDAGetNumSensErrTestFails returns the number of local error
   *                           test failures for sensitivity variables.
   * IDAGetNumSensLinSolvSetups returns the number of calls made
   *                            to the linear solver's setup routine
   *                            due to sensitivity computations.
   * IDAGetSensErrWeights returns the sensitivity error weight
   *                      vectors. The user need not allocate space
   *                      for ewtS.
   * -----------------------------------------------------------------
   */

  int IDAGetSensNumResEvals(void *ida_mem, long int *nresSevals);
  int IDAGetNumResEvalsSens(void *ida_mem, long int *nresevalsS);
  int IDAGetSensNumErrTestFails(void *ida_mem, long int *nSetfails);
  int IDAGetSensNumLinSolvSetups(void *ida_mem, long int *nlinsetupsS);
  int IDAGetSensErrWeights(void *ida_mem, N_Vector_S eSweight);

  /*
   * -----------------------------------------------------------------
   * As a convenience, the following function provides the
   * optional outputs in a group.
   * -----------------------------------------------------------------
   */

  int IDAGetSensStats(void *ida_mem, long int *nresSevals, long int *nresevalsS, 
                      long int *nSetfails, long int *nlinsetupsS);

  /*
   * ----------------------------------------------------------------
   * Sensitivity nonlinear solver optional output extraction functions          
   * ----------------------------------------------------------------
   */

  int IDAGetNumSensNonlinSolvIters(void *ida_mem, long int *nSniters);
  int IDAGetNumSensNonlinSolvConvFails(void *ida_mem, long int *nSncfails);
  int IDAGetSensNonlinSolvStats(void *ida_mem, long int *nSniters, long int *nSncfails);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with an IDAS return flag
   * -----------------------------------------------------------------
   */

  char *IDAGetReturnFlagName(int flag);

  /*
   * ----------------------------------------------------------------
   * Function : IDAFree                                             
   * ----------------------------------------------------------------
   * IDAFree frees the problem memory IDA_mem allocated by          
   * IDAMalloc.  Its only argument is the pointer idamem            
   * returned by IDAMalloc.                                         
   * ----------------------------------------------------------------
   */

  void IDAFree(void **ida_mem);

  /*
   * -----------------------------------------------------------------
   * Function : IDAQuadFree
   * -----------------------------------------------------------------
   * IDAQuadFree frees the problem memory in ida_mem allocated
   * for quadrature integration. Its only argument is the pointer
   * ida_mem returned by IDACreate.
   * -----------------------------------------------------------------
   */

  void IDAQuadFree(void *ida_mem);

  /*
   * -----------------------------------------------------------------
   * Function : IDASensFree
   * -----------------------------------------------------------------
   * IDASensFree frees the problem memory in ida_mem allocated
   * for sensitivity analysis. Its only argument is the pointer
   * ida_mem returned by IDACreate.
   * -----------------------------------------------------------------
   */

  void IDASensFree(void *ida_mem);

  /* 
   * -----------------------------------------------------------------
   * Initialization and optional input for ADJOINT module 
   * -----------------------------------------------------------------
   */

  void *IDAadjMalloc(void *ida_mem, long int steps, int interp);

  int IDAadjSetInterpType(void *idaadj_mem, int interp);

  /*
   * -----------------------------------------------------------------
   * Forward solution function 
   * -----------------------------------------------------------------
   */

  int IDASolveF(void *idaadj_mem, realtype tout, realtype *tret,
                N_Vector yret, N_Vector ypret, int itask, int *ncheckPtr);

  /*
   * -----------------------------------------------------------------
   * Initialization and optional input for backward integration 
   * -----------------------------------------------------------------
   */

  int IDACreateB(void *idaadj_mem);
  int IDAMallocB(void *idaadj_mem, IDAResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0, 
                 int itolB, realtype *reltolB, void *abstolB);
  int IDAReInitB(void *idaadj_mem, IDAResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0,
                 int itolB, realtype *reltolB, void *abstolB);

  int IDASetErrHandlerFnB(void *idaadj_mem, IDAErrHandlerFn ehfunB, void *eh_dataB);
  int IDASetErrFileB(void *idaadj_mem, FILE *errfpB);
  int IDASetRdataB(void *idaadj_mem, void *res_dataB);
  int IDASetMaxOrdB(void *idaadj_mem, int maxordB);
  int IDASetMaxNumStepsB(void *idaadj_mem, long int mxstepsB);
  int IDASetInitStepB(void *idaadj_mem, realtype hinB);
  int IDASetMaxStepB(void *idaadj_mem, realtype hmaxB);
  int IDASetSuppressAlgB(void *idaadj_mem, booleantype suppressalgB);
  int IDASetIdB(void *idaadj_mem, N_Vector idB);
  int IDASetConstraintsB(void *idaadj_mem, N_Vector constraintsB);

  int IDASetQuadFdataB(void *idaadj_mem, void *rhsQ_dataB);
  int IDASetQuadErrConB(void *idaadj_mem, booleantype errconQB, 
                        int itolQB, realtype reltolQB, void *abstolQB);
  int IDAQuadMallocB(void *idaadj_mem, IDAQuadRhsFnB rhsQB, N_Vector yQB0);
  int IDAQuadReInitB(void *idaadj_mem, IDAQuadRhsFnB rhsQB, N_Vector yQB0);

  /*
   * -----------------------------------------------------------------
   * Backward solution function 
   * -----------------------------------------------------------------
   */

  int IDASolveB(void *idaadj_mem, realtype tBout, realtype *tBret,
                N_Vector yBret, N_Vector ypBret, int itaskB);

  /*
   * -----------------------------------------------------------------
   * Optional output from backward integration 
   * -----------------------------------------------------------------
   */

  int IDAGetQuadB(void *idaadj_mem, N_Vector qB);

  /* 
   * -----------------------------------------------------------------
   * Deallocation of ADJOINT module 
   * -----------------------------------------------------------------
   */

  void IDAadjFree(void **idaadj_mem);

  /*
   * -----------------------------------------------------------------
   * Optional output from ADJOINT module 
   * -----------------------------------------------------------------
   */

  void *IDAadjGetIDABmem(void *idaadj_mem);

  char *IDAadjGetReturnFlagName(int flag);

  /*
   * -----------------------------------------------------------------
   * IDAadjGetY
   *    Returns the interpolated forward solution at time t. This
   *    function is a wrapper around the interpType-dependent internal
   *    function.
   *    The calling function must allocate space for y.
   * -----------------------------------------------------------------
   */

  int IDAadjGetY(void *idaadj_mem, realtype t, N_Vector y, N_Vector yp);

  /* 
   * ===============================================================
   * DEVELOPMENT USER-CALLABLE FUNCTIONS
   * ===============================================================
   */

  /*
  typedef struct {
    void *my_addr;
    void *next_addr;
    realtype t0;
    realtype t1;
    long int nstep;
    int order;
    realtype step;
  } IDAadjCheckPointRec;

  int IDAadjGetCheckPointsInfo(void *idaadj_mem, IDAadjCheckPointRec *ckpnt);
  int IDAadjGetCurrentCheckPoint(void *idaadj_mem, void **addr);

  int IDAadjGetDataPointHermite(void *idaadj_mem, long int which,
                                realtype *t, N_Vector y, N_Vector yd);
  
  int CVadjGetDataPointPolynomial(void *idaadj_mem, long int which,
                                  realtype *t, int *order, N_Vector y);
  
  */

#ifdef __cplusplus
}
#endif

#endif
