/******************************************************************
 *                                                                *
 * File          : ida.h                                          *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and        *
 *                 Radu Serban @ LLNL                             *
 * Version of    : 5 MArch 2002                                   *
 *----------------------------------------------------------------*
 * This is the header (include) file for the main IDA solver.     *
 *                                                                *
 ******************************************************************/

/*......................................................................

                            LEGAL NOTICES

This work was performed at the University of California, Lawrence
Livermore National Laboratory (UC LLNL) under contract no.
W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy
(DOE) and The Regents of the University of California (the University)
for the operation of UC LLNL.  The rights of the Federal Government are
reserved under Contract 48 subject to the restrictions agreed upon by the
DOE and University as allowed under DOE Acquisition Letter 97-1.

This work was prepared as an account of work sponsored by an agency of
the United States Government.  Neither the United States Government
nor the University of California nor any of their empolyees makes any
warranty, express or implied, or assumes any liability or
responsibility for the accuracy, completeness, or usefulness of any
information, apparatus, product, or process disclosed, or represents
that its use would not infringe privately owned rights.  Reference
herein to any specific commercial products, process, or service by
trade name, trademark, manufacturer, or otherwise, does not
necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or the University of
California.  The views and opinions of authors expressed herein do not
necessarily state or reflect those of the United States Government or
the University of California, and shall not be used for advertising or
product endorsement purposes.

......................................................................*/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _ida_h
#define _ida_h

#include <stdio.h>
#include "llnltyps.h"
#include "nvector.h"

/******************************************************************
 *                                                                *
 * IDA is used to solve numerically the initial value problem     *
 * for the differential algebraic equation (DAE) system           *
 *   F(t,y,y') = 0,                                               *
 * given initial conditions                                       *
 *   y(t0) = y0,   y'(t0) = yp0.                                  *
 * Here y and F are vectors of length N.                          *
 *                                                                *
 ******************************************************************/

/******************************************************************
 *                                                                *
 * Type : ResFn                                                   *
 *----------------------------------------------------------------*        
 * The F function which defines the DAE system   F(t,y,y')=0      *
 * must have type ResFn.                                          *
 * Symbols are as follows: t  <-> tres     y <-> yy               *
 *                         y' <-> yp       F <-> res (type ResFn) *
 * A ResFn takes as input the problem size Neq, the independent   *
 * variable value tres, the dependent variable vector yy, and the *
 * derivative (with respect to t) of the yy vector, yp.  It       *
 * stores the result of F(t,y,y') in the vector resval.  The      *
 * yy, yp, and resval arguments are of type N_Vector.             *
 * The rdata parameter is to be of the same type as the rdata     *
 * parameter passed by the user to the IDAMalloc routine. This    *
 * user-supplied pointer is passed to the user's res function     *
 * every time it is called, to provide access in res to user data.*
 *                                                                *
 * A ResFn res will return the value ires, which has possible     *
 * values RES_ERROR_RECVR = 1, RES_ERROR_NONRECVR = -1,           *
 * and SUCCESS = 0. The file ida.h may be used to obtain these    *
 * values but is not required; returning 0, +1, or -1 suffices.   *
 * RES_ERROR_NONRECVR will ensure that the program halts.         *
 * RES_ERROR_RECVR should be returned if, say, a yy or other input*
 * value is illegal. IDA will attempt to correct and retry.       *
 *                                                                *
 ******************************************************************/

typedef int (*ResFn)(integer Neq, real tres, N_Vector yy, 
              N_Vector yp, N_Vector resval, void *rdata);
 

/******************************************************************
 *                                                                *
 * Enumerations for input parameters to IDAMalloc and IDASolve.   *
 *----------------------------------------------------------------*/
/* iopt indices */

enum { MAXORD, MXSTEP, SUPPRESSALG ,
       NST, NRE, NNI, NCFN, NETF, NSETUPS, KUSED, KNEXT,
       LENRW, LENIW, NBACKTR };

/* ropt indices */

enum { HINIT, HMAX, NCONFAC, HMIN, 
       HUSED, HNEXT, TNOW, TOLSF };


enum { SS, SV };               /* itol */

enum { NORMAL, ONE_STEP, NORMAL_TSTOP, ONE_STEP_TSTOP }; /* itask */



/* Enumerations for res, lsetup and lsolve return values */

enum {SUCCESS               = 0,
      RES_ERROR_RECVR       = 1, RES_ERROR_NONRECVR      = -1,
      LSETUP_ERROR_RECVR    = 2, LSETUP_ERROR_NONRECVR   = -2,
      LSOLVE_ERROR_RECVR    = 3, LSOLVE_ERROR_NONRECVR   = -3};


/* IDACalcIC icopt values  */

enum { CALC_YA_YDP_INIT = 1 , CALC_Y_INIT = 2 };


/* IDACalcIC return values */

/* The following three values are IDASolve return values, repeated
   here for convenience. 
       SETUP_FAILURE=-7,  SOLVE_FAILURE=-8,  RES_NONRECOV_ERR=-11  */

enum { IC_IDA_NO_MEM =     -20,   IC_ILL_INPUT =       -21,
       IC_LINIT_FAIL =     -22,   IC_BAD_EWT =         -23,
       IC_FIRST_RES_FAIL = -24,   IC_NO_RECOVERY =     -25,
       IC_FAILED_CONSTR =  -26,   IC_FAILED_LINESRCH = -27,
       IC_CONV_FAILURE =   -28 };


/* IDASolve return values */

enum { NORMAL_RETURN=0,     INTERMEDIATE_RETURN=1, TSTOP_RETURN=2,
       IDA_NO_MEM=-1,       ILL_INPUT=-2,          TOO_MUCH_WORK=-3,      
       TOO_MUCH_ACC=-4,     ERR_FAILURE=-5,        CONV_FAILURE=-6,       
       SETUP_FAILURE=-7,    SOLVE_FAILURE=-8,      CONSTR_FAILURE=-9,     
       REP_RES_REC_ERR=-10, RES_NONRECOV_ERR=-11 };

 
/******************************************************************
 *                                                                *
 * Function : IDAMalloc                                           *
 *----------------------------------------------------------------*
 * IDAMalloc allocates and initializes memory for a problem to    *
 * to be solved by IDA.                                           *
 *                                                                *
 * Neq     is the number of equations in the DAE system.          *
 *         (In the parallel case, Neq is the global system size,  *
 *         not the local size.)                                   *
 *                                                                *
 * res     is the residual function F in F(t,y,y') = 0.           *          
 *                                                                *
 * rdata   is the data memory block (supplied by user) for res.   *
 *         It is passed as a void pointer and is to be cast before*
 *         use in res.                                            *
 *                                                                *
 * t0      is the initial value of t, the independent variable.   *
 *                                                                *
 * y0      is the initial condition vector y(t0).                 *
 *                                                                *
 * yp0     is the initial condition vector y'(t0)                 *
 *                                                                *
 * itol    is the type of tolerances to be used.                  *
 *            The legal values are:                               *
 *               SS (scalar relative and absolute  tolerances),   *
 *               SV (scalar relative tolerance and vector         *
 *                   absolute tolerance).                         *
 *                                                                *
 * rtol    is a pointer to the relative tolerance scalar.         *
 *                                                                *
 * atol    is a pointer (void) to the absolute tolerance scalar or*
 *            an N_Vector tolerance.                              *
 * (ewt)                                                          *
 *         Both rtol and atol are used to compute the error weight*
 *         vector, ewt. The error test required of a correction   *
 *         delta is that the weighted-RMS norm of delta be less   *
 *         than or equal to 1.0. Other convergence tests use the  *
 *         same norm. The weighting vector used in this norm is   *
 *         ewt. The components of ewt are defined by              *
 *         ewt[i] = 1.0/(rtol*yy[i] + atol[i]). Here, yy is the   *
 *         current approximate solution.  See the routine         *
 *         N_VWrmsNorm for the norm used in this error test.      *
 *                                                                *
 * id      is an N_Vector, required conditionally, which states a *
 *         given element to be either algebraic or differential.  *
 *         A value of 1.0 indicates a differential variable while *
 *         a 0.0 indicates an algebraic variable. 'id' is required*
 *         if optional input SUPPRESSALG is set, or if IDACalcIC  *
 *         is to be called with icopt = CALC_YA_YDP_INIT.         *
 *         Otherwise, 'id' may be NULL.                           *
 *                                                                *
 * constraints  is an N_Vector defining inequality constraints    *
 *         for each component of the solution vector y. If a given*
 *         element of this vector has values +2 or -2, then the   *
 *         corresponding component of y will be constrained to be *
 *         > 0.0 or < 0.0, respectively, while if it is +1 or -1, *
 *         the y component is constrained to be >= 0.0 or <= 0.0, *
 *         respectively. If a component of constraints is 0.0,    *
 *         then no constraint is imposed on the corresponding     *
 *         component of y. The presence of a non-NULL constraints *
 *         vector that is not 0.0 (ZERO) in all components will   *
 *         cause constraint checking to be performed.             *
 *                                                                *
 * errfp   is the file pointer for an error file where all IDA    *
 *            warning and error messages will be written. This    *
 *            parameter can be stdout (standard output), stderr   *
 *            (standard error), a file pointer (corresponding to  *
 *            a user error file opened for writing) returned by   *
 *            fopen, or NULL. If the user passes NULL, then all   *
 *            messages will be written to standard output.        *
 *                                                                *
 * optIn   is a flag (boole) indicating whether there are any     *
 *            optional inputs from the user in the arrays         *
 *            iopt and iopt.                                      *
 *            Pass FALSE to indicate no optional inputs and TRUE  *
 *            to indicate that optional inputs are present.       *
 *                                                                *
 * iopt    is the user-allocated array (of size OPT_SIZE given    *
 *            later) that will hold optional integer inputs and   *
 *            outputs.  The user can pass NULL if he/she does not *
 *            wish to use optional integer inputs or outputs.     *
 *            If optIn is TRUE, the user should preset to 0 those *
 *            locations for which default values are to be used.  *
 *                                                                *
 * ropt    is the user-allocated array (of size OPT_SIZE given    *
 *            later) that will hold optional real inputs and      *
 *            outputs.  The user can pass NULL if he/she does not *
 *            wish to use optional real inputs or outputs.        *
 *            If optIn is TRUE, the user should preset to 0.0 the *
 *            locations for which default values are to be used.  *
 *                                                                *
 * machEnv is a pointer to a block of machine environment-specific*
 *            information.  You may pass NULL for this argument   *
 *            in a serial computer environment.                   *
 *                                                                *
 * Note: The tolerance values may be changed in between calls to  *
 *       IDASolve for the same problem. These values refer to     *
 *       (*rtol) and either (*atol), for a scalar absolute        *
 *       tolerance, or the components of atol, for a vector       *
 *       absolute tolerance.                                      *
 *                                                                * 
 * If successful, IDAMalloc returns a pointer to initialized      *
 * problem memory. This pointer should be passed to IDA. If       *
 * an initialization error occurs, IDAMalloc prints an error      *
 * message to the file specified by errfp and returns NULL.       *
 *                                                                *
 ******************************************************************/

void *IDAMalloc(integer Neq, ResFn res, void *rdata, real t0,
                N_Vector y0, N_Vector yp0, int itol, real *rtol, void *atol, 
                N_Vector id, N_Vector constraints, FILE *errfp, boole optIn, 
                long int iopt[], real ropt[], M_Env machEnv);

 
/******************************************************************
 *                                                                *
 * Optional Inputs and Outputs                                    *
 *----------------------------------------------------------------*
 * The user should declare two arrays for optional inputs to      *
 * IDAMalloc and optional outputs from IDACalcIC and IDASolve:    *
 * an long int array iopt for optional integer input/output       *
 * and a real array ropt for optional real input/output.          *
 * The size of both these arrays should be OPT_SIZE.              *
 * So the user's declaration should look like:                    *
 *   long int iopt[OPT_SIZE];                                     *
 *   real     ropt[OPT_SIZE];                                     *
 *                                                                *
 * The enumerations listed earlier are indices into the           * 
 * iopt and ropt arrays. Here is a brief description of the       *
 * contents of these positions:                                   *
 *                                                                *
 * iopt[MAXORD] : maximum order to be used by the solver.         *
 *                Optional input. (Default = 5)                   *
 *                                                                *
 * iopt[MXSTEP] : maximum number of internal steps to be taken by *
 *                the solver in its attempt to reach tout.        *
 *                Optional input. (Default = 500).                *
 *                                                                *
 * iopt[SUPPRESSALG]: flag to indicate whether or not to suppress *
 *                algebraic variables in the local error tests:   *
 *                0 = do not suppress ; 1 = do suppress;          *
 *                the default is 0.  Optional input.              *      
 *                NOTE: if suppressed algebraic variables is      * 
 *                selected, the nvector 'id' must be supplied for *
 *                identification of those algebraic components.   *
 *                                                                *
 * iopt[NST]    : cumulative number of internal steps taken by    *
 *                the solver (total so far).  Optional output.    *
 *                                                                *
 * iopt[NRE]    : number of calls to the user's residual function.*
 *                Optional output.                                *
 *                                                                *
 * iopt[NNI]     : number of Newton iterations performed.         *
 *                 Optional output.                               *
 *                                                                *
 * iopt[NCFN]    : number of nonlinear convergence failures       *
 *                 that have occurred. Optional output.           *
 *                                                                *
 * iopt[NETF]    : number of local error test failures that       *
 *                 have occurred. Optional output.                *
 *                                                                *
 * iopt[NSETUPS] : number of calls to lsetup routine.             *
 *                                                                *
 * iopt[NBACKTR] : number of backtrack operations done in the     *
 *                 linesearch algorithm in IDACalcIC.             *
 *                                                                *
 * iopt[KUSED]   : order used during the last internal step.      *
 *                 Optional output.                               *
 *                                                                *
 * iopt[KNEXT]   : order to be used on the next internal step.    *
 *                 Optional output.                               *
 *                                                                *
 * iopt[LENRW]   : size of required IDA internal real work        *
 *                 space, in real words.  Optional output.        *
 *                                                                *
 * iopt[LENIW]   : size of required IDA internal integer work     *
 *                 space, in integer words.  Optional output.     *
 *                                                                *
 * ropt[HINIT]   : initial step size. Optional input.             *
 *                                                                *
 * ropt[HMAX]    : maximum absolute value of step size allowed.   *
 *                 Optional input. (Default is infinity).         *
 *                                                                *
 * ropt[NCONFAC] : factor in nonlinear convergence test for use   *
 *                 during integration.  Optional input.           *
 *                 The default value is 1.                        *
 *                                                                *
 * ropt[HUSED]   : step size for the last internal integration    *
 *                 step (if from IDASolve), or the last value of  *
 *                 the artificial step size h (if from IDACalcIC).*
 *                 Optional output.                               *
 *                                                                *
 * ropt[HNEXT]   : step size to be attempted on the next internal *
 *                 step. Optional output.                         *
 *                                                                *
 * ropt[TNOW]    : current internal independent variable value    *
 *                 reached by the solver.  Optional output.       *
 *                                                                *
 * ropt[TOLSF]   : a suggested factor by which the user's         *
 *                 tolerances should be scaled when too much      *
 *                 accuracy has been requested for some internal  *
 *                 step. Optional output.                         *
 *                                                                *
 ******************************************************************/


/******************************************************************
 *                                                                *
 * Function : IDACalcIC                                           *
 *                                                                *
 *----------------------------------------------------------------*
 * IDACalcIC calculates corrected initial conditions for the DAE  *
 * system for a class of index-one problems of semi-implicit form.*
 * It uses Newton iteration combined with a Linesearch algorithm. *
 * Calling IDACalcIC is optional. It is only necessary when the   *
 * initial conditions do not solve the given system.  I.e., if    *
 * y0 and yp0 are known to satisfy F(t0, y0, yp0) = 0, then       *
 * a call to IDACalcIC is NOT necessary (for index-one problems). *
 *                                                                *
 * A call to IDACalcIC must be preceded by a successful call to   *
 * IDAMalloc, and by a successful call to the linear system       *
 * solver specification routine.  In addition, IDACalcIC assumes  *
 * that the vectors y0, yp0 and (if relevant) id and constraints  *
 * that were passed to IDAMalloc remain unaltered since that call.*
 *                                                                *
 * The call to IDACalcIC should precede the call(s) to IDASolve   *
 * for the given problem.                                         *  
 *                                                                *
 * The arguments to IDACalcIC are as follows.  The first three -- *
 * ida_mem, icopt, tout1 -- are required; the others are optional.*
 * A zero value passed for any optional input specifies that the  *
 * default value is to be used.                                   *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by IDAMalloc.    *
 *                                                                *
 * icopt  is the option of IDACalcIC to be used.                  *
 *        icopt = CALC_YA_YDP_INIT   directs IDACalcIC to compute *
 *                the algebraic components of y and differential  *
 *                components of y', given the differential        *
 *                components of y.  This option requires that the *
 *                N_Vector id was input to IDAMalloc, specifying  *
 *                the differential and algebraic components.      *
 *        icopt = CALC_Y_INIT   directs IDACalcIC to compute all  *
 *                components of y, given y'.  id is not required. *
 *                                                                *
 * tout1  is the first value of t at which a soluton will be      *
 *        requested (from IDASolve).  (This is needed here to     *
 *        determine the direction of integration and rough scale  *
 *        in the independent variable t.                          *
 *                                                                *
 * epicfac is a positive scalar factor in the Newton convergence  *
 *        test.  This test uses a weighted RMS norm (with weights *
 *        defined by the tolerances, as in IDASolve).  For new    *
 *        initial value vectors y and y' to be accepted, the norm *
 *        of J-inverse F(t0,y,y') is required to be less than     *
 *        epicfac*0.33, where J is the system Jacobian.           *
 *        The default is epicfac = .01, specified by passing 0.   *
 *                                                                *
 * maxnh  is the maximum number of values of h allowed in the     *
 *        algorithm for icopt = CALC_YA_YDP_INIT, where h appears *
 *        in the system Jacobian, J = dF/dy + (1/h)dF/dy'.        *
 *        The default is maxnh = 5, specified by passing 0.       *
 *                                                                *
 * maxnj  is the maximum number of values of the approximate      *
 *        Jacobian or preconditioner allowed, when the Newton     *
 *        iterations appear to be slowly converging.              *
 *        The default is maxnj = 4, specified by passing 0.       *
 *                                                                *
 * maxnit is the maximum number of Newton iterations allowed in   *
 *        any one attempt to solve the IC problem.                *
 *        The default is maxnit = 10, specified by passing 0.     *
 *                                                                *
 * lsoff  is an integer flag to turn off the linesearch algorithm *
 *        (lsoff = 1). The default is lsoff = 0 (linesearch done).*
 *                                                                *
 * steptol is a positive lower bound on the norm of a Newton step.*
 *        The default value is steptol = (unit roundoff)^(2/3),   *
 *        specified by passing 0.                                 *
 *                                                                *
 *                                                                *
 * IDACalcIC returns an int flag.  Its symbolic values and their  *
 * meanings are as follows.  (The numerical return values are set *
 * above in this file.)  All unsuccessful returns give a negative *
 * return value.  If IFACalcIC failed, y0 and yp0 contain         *
 * (possibly) altered values, computed during the attempt.        *
 *                                                                *
 * SUCCESS             IDACalcIC was successful.  The corrected   *
 *                     initial value vectors are in y0 and yp0.   * 
 *                                                                *
 * IC_IDA_NO_MEM       The argument ida_mem was NULL.             *
 *                                                                *
 * IC_ILL_INPUT        One of the input arguments was illegal.    *
 *                     See printed message.                       *
 *                                                                *
 * IC_LINIT_FAIL       The linear solver's init routine failed.   *
 *                                                                *
 * IC_BAD_EWT          Some component of the error weight vector  *
 *                     is zero (illegal), either for the input    *
 *                     value of y0 or a corrected value.          *
 *                                                                *
 * RES_NONRECOV_ERR    The user's ResFn residual routine returned *
 *                     a non-recoverable error flag.              *
 *                                                                *
 * IC_FIRST_RES_FAIL   The user's ResFn residual routine returned *
 *                     a recoverable error flag on the first call,*
 *                     but IDACalcIC was unable to recover.       *
 *                                                                *
 * SETUP_FAILURE       The linear solver's setup routine had a    *
 *                     non-recoverable error.                     *
 *                                                                *
 * SOLVE_FAILURE       The linear solver's solve routine had a    *
 *                     non-recoverable error.                     *
 *                                                                *
 * IC_NO_RECOVERY      The user's residual routine, or the linear *
 *                     solver's setup or solve routine had a      *
 *                     recoverable error, but IDACalcIC was       *
 *                     unable to recover.                         *
 *                                                                *
 * IC_FAILED_CONSTR    IDACalcIC was unable to find a solution    *
 *                     satisfying the inequality constraints.     *
 *                                                                *
 * IC_FAILED_LINESRCH  The Linesearch algorithm failed to find a  *
 *                     solution with a step larger than steptol   *
 *                     in weighted RMS norm.                      *
 *                                                                *
 * IC_CONV_FAILURE     IDACalcIC failed to get convergence of the *
 *                     Newton iterations.                         *
 *                                                                *
 *                                                                *
 * The following optional outputs provided by IDACalcIC are       *
 * available in the iopt and ropt arrays.                         *
 *                                                                *
 * iopt[NRE]     = number of calls to the user residual function. *
 *                                                                *
 * iopt[NNI]     = number of Newton iterations performed.         *
 *                                                                *
 * iopt[NCFN]    = number of nonlinear convergence failures.      *
 *                                                                *
 * iopt[NSETUPS] = number of calls to linear solver setup routine.*
 *                                                                *
 * iopt[NBACKTR] = number of backtracks in Linesearch algorithm.  *
 *                                                                *
 * ropt[HUSED]   = value of h last used in the system Jacobian J  *
 *                 if icopt = CALC_YA_YDP_INIT.                   *
 *                                                                *
 ******************************************************************/

 int IDACalcIC (void *ida_mem, int icopt, real tout1, real epicfac, 
               int maxnh, int maxnj, int maxnit, int lsoff, real steptol);


/******************************************************************
 *                                                                *
 * Function : IDASolve                                            *
 *----------------------------------------------------------------*
 * IDASolve integrates the DAE over an interval in t, the         *
 * independent variable. If itask is NORMAL, then the solver      *
 * integrates from its current internal t value to a point at or  *
 * beyond tout, then interpolates to t = tout and returns y(tret) *
 * in the user-allocated vector yret. In general, tret = tout.    *
 * If itask is ONE_STEP, then the solver takes one internal step  *
 * of the independent variable and returns in yret the value of y *
 * at the new internal independent variable value. In this case,  *
 * tout is used only during the first call to IDASolve to         *
 * determine the direction of integration and the rough scale of  *
 * the problem. In either case, the independent variable value    *
 * reached by the solver is placed in (*tret). The user is        * 
 * responsible for allocating the memory for this value. There are*
 * two other itask options: NORMAL_TSTOP and ONE_STEP_TSTOP. Each *
 * option acts as described above with the exception that the     *
 * solution does not go past the independent variable value tstop.*
 *                                                                *
 * IDA_mem is the pointer (void) to IDA memory returned by        *
 * IDAMalloc.                                                     *
 *                                                                *
 * tout  is the next independent variable value  at which a       *
 * computed solution is desired                                   *
 *                                                                *
 * tstop is the (optional) independent variable value past which  *
 * the solution is not to proceed (see itask options)             *
 *                                                                *
 * *tret  is the actual independent variable value corresponding  *
 * to the solution vector yret.                                   *
 *                                                                *
 * yret  is the computed solution vector. With no errors,         *
 *            yret=y(tret).                                       *
 *                                                                *
 * ypret is the derivative of the computed solution vector at tret*
 *                                                                *
 * Note: yret and ypret may be the same N_Vectors as y0 and yp0   *
 * in the call to IDAMalloc.                                      *
 *                                                                *
 * itask is NORMAL, NORMAL_TSTOP, ONE_STEP, or ONE_STEP_TSTOP.    *
 *        These modes are described above.                        *
 *                                                                *
 *                                                                *
 * The return values for IDASolve are described below.            *
 * (The numerical return values are defined above in this file.)  *
 * All unsuccessful returns give a negative return value.         *
 *                                                                *
 * NORMAL_RETURN      : IDASolve succeeded.                       *
 *                                                                *
 * INTERMEDIATE_RETURN: IDASolve returns computed results for the *
 *                  last single step (itask = ONE_STEP).          *
 *                                                                *
 * TSTOP_RETURN       : IDASolve returns computed results for the *
 *                  independent variable value tstop. That is,    *
 *                  tstop was reached.                            *
 *                                                                *
 * IDA_NO_MEM         : The IDA_mem argument was NULL.            *
 *                                                                *
 * ILL_INPUT          : One of the inputs to IDASolve is illegal. *
 *                 This includes the situation when a component   *
 *                 of the error weight vectors becomes < 0 during *
 *                 internal stepping. The ILL_INPUT flag          *
 *                 will also be returned if the linear solver     *
 *                 routine IDA--- (called by the user after       *
 *                 calling IDAMalloc) failed to set one of the    *
 *                 linear solver-related fields in IDA_mem or     *
 *                 if the linear solver's init routine failed. In *
 *                 any case, the user should see the printed      *
 *                 error message for more details.                *
 *                                                                *
 * TOO_MUCH_WORK      : The solver took mxstep internal steps but *
 *                 could not reach tout. The default value for    *
 *                 mxstep is MXSTEP_DEFAULT = 500.                *
 *                                                                *
 * TOO_MUCH_ACC       : The solver could not satisfy the accuracy *
 *                 demanded by the user for some internal step.   *
 *                                                                *
 * ERR_FAILURE        : Error test failures occurred too many     *
 *                 times (=MXETF = 10) during one internal step.  *
 *                                                                *
 * CONV_FAILURE       : Convergence test failures occurred too    *
 *                 many times (= MXNCF = 10) during one internal  * 
 *                 step.                                          *
 *                                                                *
 * SETUP_FAILURE      : The linear solver's setup routine failed  *
 *                 in an unrecoverable manner.                    *
 *                                                                *
 * SOLVE_FAILURE      : The linear solver's solve routine failed  *
 *                 in an unrecoverable manner.                    *
 *                                                                *
 * CONSTR_FAILURE     : The inequality constraints were violated, *
 *                  and the solver was unable to recover.         *
 *                                                                *
 * REP_RES_REC_ERR    : The user's residual function repeatedly   *
 *                  returned a recoverable error flag, but the    *
 *                  solver was unable to recover.                 *
 *                                                                *
 * RES_NONRECOV_ERR   : The user's residual function returned a   *
 *                  nonrecoverable error flag.                    *
 *                                                                *
 ******************************************************************/


int IDASolve(void *ida_mem, real tout, real tstop, real *tret,
             N_Vector yret, N_Vector ypret, int itask);



/******************************************************************
 *                                                                *
 * Function : IDAFree                                             *
 *----------------------------------------------------------------*
 * IDAFree frees the problem memory IDA_mem allocated by          *
 * IDAMalloc.  Its only argument is the pointer idamem            *
 * returned by IDAMalloc.                                         *
 *                                                                *
 ******************************************************************/

void IDAFree(void *ida_mem);


/* iopt, ropt array sizes */

#define OPT_SIZE 40
 

/* iopt and ropt offsets                                          *
 * The constants IDA_IOPT_SIZE and IDA_ROPT_SIZE are equal to     *
 * the number of integer and real optional inputs and outputs     *
 * actually accessed in ida.c.  The locations beyond these        *
 * values are used by the linear solvers.                         */

#define IDA_IOPT_SIZE 14
#define IDA_ROPT_SIZE 8


/* Basic IDA constants */

#define MXORDP1       6     /* max value 'small' dimension for 
			  phi array size   (other dimension is Neq)  */

/******************************************************************
 *                                                                *
 * Types : struct IDAMemRec, IDAMem                               *
 *----------------------------------------------------------------*
 * The type IDAMem is type pointer to struct IDAMemRec. This      *
 * structure contains fields to keep track of problem state.      *
 *                                                                *
 ******************************************************************/

typedef struct IDAMemRec {

  real ida_uround;    /* machine unit roundoff */

  /* Problem Specification Data */

  integer  ida_Neq;            /* DAE system size                    */
  ResFn    ida_res;            /* F(t,y(t),y'(t))=0; the function F  */
  void    *ida_rdata;          /* user pointer passed to res         */
  int      ida_itol;           /* itol = SS or SV                    */
  real    *ida_rtol;           /* ptr to relative tolerance          */
  void    *ida_atol;           /* ptr to absolute tolerance          */  
  boole    ida_setupNonNull;   /* Does setup do something?           */
  boole    ida_constraintsSet; /* constraints vector present: 
                                 do constraints calc                 */
  boole    ida_suppressalg;    /* true means suppress algebraic vars
                                    in local error tests             */

  /* Divided differences array and associated minor arrays */

  N_Vector ida_phi[MXORDP1];
           /* phi = Neq x (maxord+1) array of divided differences */

  real ida_psi[MXORDP1]; /* differences in t (sums of successive step sizes) */
  real ida_alpha[MXORDP1]; /* ratios of current stepsize to psi values       */
  real ida_beta[MXORDP1]; 
                  /* ratios of current to previous product of psi values     */
  real ida_sigma[MXORDP1];/* product successive alpha values and factorial   */
  real ida_gamma[MXORDP1]; /* sum of reciprocals of psi values               */

  /* Vectors of length Neq */

  N_Vector ida_ewt;         /* error weight vector                        */
  N_Vector ida_mskewt;      /* masked error weight vector (uses ID)       */
  N_Vector ida_y0;          /* initial y vector (user-supplied)           */
  N_Vector ida_yp0;         /* initial y' vector (user-supplied)          */
  N_Vector ida_yy;          /* work space for y vector (= user's yret)    */
  N_Vector ida_yp;          /* work space for y' vector (= user's ypret)  */
  N_Vector ida_delta;       /* residual vector                            */
  N_Vector ida_id;          /* bit vector for diff./algebraic components  */
  N_Vector ida_constraints; /* vector of inequality constraint options    */
  N_Vector ida_savres;      /* saved residual vector (= tempv1)           */
  N_Vector ida_ee;          /* accumulated corrections to y               */
  N_Vector ida_mm;          /* mask vector in constraints tests (= tempv2)*/
  N_Vector ida_tempv1;      /* work space vector                          */
  N_Vector ida_tempv2;      /* work space vector                          */
  N_Vector ida_ynew;     /* work vector for y in IDACalcIC (= tempv2)     */
  N_Vector ida_ypnew;    /* work vector for yp in IDACalcIC (= ee)        */
  N_Vector ida_delnew;   /* work vector for delta in IDACalcIC (= phi[2]) */
  N_Vector ida_dtemp;    /* work vector in IDACalcIC (= phi[3])           */

/* Scalars for use by IDACalcIC*/

  int ida_icopt;        /* IC calculation user option                   */
  int ida_lsoff;        /* IC calculation linesearch turnoff option     */
  int ida_maxnj;        /* max. number of J tries in IC calculation     */
  int ida_maxnit;       /* max. number of Netwon iterations in IC calc. */
  int ida_nbacktr;      /* number of IC linesearch backtrack operations */
  int ida_sysindex;     /* computed system index (0 or 1)               */
  real ida_epsic;       /* IC nonlinear convergence test constant       */
  real ida_steptol;     /* minimum Newton step size in IC calculation   */
  real ida_tscale;      /* time scale factor = abs(tout1 - t0)          */


  /* Step Data */

  int ida_kk;        /* current BDF method order                          */
  int ida_kused;     /* method order used on last successful step         */
  int ida_knew;      /* order for next step from order decrease decision  */
  int ida_phase;     /* flag to trigger step doubling in first few steps  */
  int ida_ns;        /* counts steps at fixed stepsize and order          */

  real ida_hh;       /* current step size h                               */
  real ida_hused;    /* step size used on last successful step            */
  real ida_rr;       /* rr = hnext / hused                                */
  real ida_tn;       /* current internal value of t                       */
  real ida_tretp;    /* value of tret previously returned by IDASolve     */
  real ida_cj;       /* current value of scalar (-alphas/hh) in Jacobian  */
  real ida_cjlast;   /* cj value saved from last successful step          */
  real ida_cjold;    /* cj value saved from last call to lsetup           */
  real ida_cjratio;  /* ratio of cj values: cj/cjold                      */
  real ida_ss;       /* scalar used in Newton iteration convergence test  */
  real ida_epsNewt;  /* test constant in Newton convergence test          */
  real ida_nconfac;  /* optional factor in Newton covergence test constant*/
  real ida_toldel;   /* tolerance in direct test on Newton corrections    */


 /* Limits */

  int ida_maxord;    /* max value of method order k:                    */
  int ida_mxstep;    /* max number of internal steps for one user call  */
  real ida_hmax_inv; /* inverse of max. step size hmax (default = 0.0)  */

  /* Counters */

  long int ida_nst;     /* number of internal steps taken             */
  long int ida_nre;     /* number of function (res) calls             */
  long int ida_ncfn;    /* number of corrector convergence failures   */
  long int ida_netf;    /* number of error test failures              */
  long int ida_nni;     /* number of Newton iterations performed      */
  long int ida_lrw;     /* number of real words in IDA work vectors   */
  long int ida_liw;     /* no. of integer words in IDA work vectors   */
  long int ida_nsetups; /* number of lsetup calls                     */

  /* Saved Values */

  real ida_tolsf;        /* tolerance scale factor         */

  /* Arrays for Optional Input and Optional Output */

  long int *ida_iopt;  /* long int optional input, output */
  real     *ida_ropt;  /* real optional input, output     */


  /* Linear Solver Data */

  /* Linear Solver functions to be called */

  int (*ida_linit)(struct IDAMemRec *idamem, boole *setupNonNull);

  int (*ida_lsetup)(struct IDAMemRec *idamem, N_Vector yyp, 
                    N_Vector ypp, N_Vector resp, 
                    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3); 

  int (*ida_lsolve)(struct IDAMemRec *idamem, N_Vector b, N_Vector ycur,
		   N_Vector ypcur, N_Vector rescur);

  int (*ida_lperf)(struct IDAMemRec *idamem, int perftask);

  int (*ida_lfree)(struct IDAMemRec *idamem);

  /* Linear Solver specific memory */

  void *ida_lmem;           

  /* Flag to indicate successful ida_linit call */

  boole ida_linitOK;


  /* Error File */

  FILE *ida_errfp;      /* IDA error messages are sent to errfp */

  /* Pointer to Machine Environment-Specific Information   */

  M_Env ida_machenv;

} *IDAMem;


/******************************************************************
 *                                                                *
 * Communication between ida.c and an IDA Linear Solver Module    *
 *----------------------------------------------------------------*
 * (1) ida_linit return values                                    *
 *                                                                *
 * LINIT_OK    : The ida_linit routine succeeded.                 *
 *                                                                *
 * LINIT_ERR   : The ida_linit routine failed. Each linear solver *
 *               init routine should print an appropriate error   *
 *               message to (idamem->errfp).                      *
 *                                                                *
 * (2) Parameter documentation, as well as a brief description    *
 *     of purpose, for each IDA linear solver routine to be       *
 *     called in IDA.c is given below the constant declarations   *
 *     that follow.                                               *
 *                                                                *
 ******************************************************************/

/* ida_linit return values */

#define LINIT_OK        0
#define LINIT_ERR      -1

/*******************************************************************
 *                                                                 *
 * int (*ida_linit)(IDAMem IDA_mem, boole *setupNonNull);          *
 *-----------------------------------------------------------------*
 * The purpose of ida_linit is to allocate memory for the          *
 * solver-specific fields in the structure *(idamem->ida_lmem) and *
 * perform any needed initializations of solver-specific memory,   *
 * such as counters/statistics. The ida_linit routine should set   *
 * *setupNonNull to be TRUE if the setup operation for the linear  *
 * solver is non-empty and FALSE if the setup operation does       *
 * nothing. An (*ida_linit) should return LINIT_OK (== 0) if it has*
 * successfully initialized the IDA linear solver and LINIT_ERR    *
 * (== -1) otherwise. These constants are defined above. If an     *
 * error does occur, an appropriate message should be sent to      *
 * (idamem->errfp).                                                *
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * int (*ida_lsetup)(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,   *
 *                  N_Vector resp,                                 *
 *            N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);  *
 *-----------------------------------------------------------------*
 * The job of ida_lsetup is to prepare the linear solver for       *
 * subsequent calls to ida_lsolve. Its parameters are as follows:  *
 *                                                                 *
 * idamem - problem memory pointer of type IDAMem. See the big     *
 *          typedef earlier in this file.                          *
 *                                                                 *
 *                                                                 *
 * yyp   - the predicted y vector for the current IDA internal     *
 *         step.                                                   *
 *                                                                 *
 * ypp   - the predicted y' vector for the current IDA internal    *
 *         step.                                                   *
 *                                                                 *
 * resp  - F(tn, yyp, ypp).                                        *
 *                                                                 *
 * tempv1, tempv2, tempv3 - temporary N_Vectors provided for use   *
 *         by ida_lsetup.                                          *
 *                                                                 *
 * The ida_lsetup routine should return SUCCESS (=0) if successful,*
 * the positive value LSETUP_ERROR_RECVR for a recoverable error,  *
 * and the negative value LSETUP_ERROR_NONRECVR for an             *
 * unrecoverable error.  The code should include the file ida.h .  *
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * int (*ida_lsolve)(IDAMem IDA_mem, N_Vector b, N_Vector ycur,    *
 *                 NPVector ypcur, N_Vector rescur);               *
 *-----------------------------------------------------------------*
 * ida_lsolve must solve the linear equation P x = b, where        *
 * P is some approximation to the system Jacobian                  *
 *                  J = (dF/dy) + cj (dF/dy')                      *
 * evaluated at (tn,ycur,ypcur) and the RHS vector b is input.     *
 * The N-vector ycur contains the solver's current approximation   *
 * to y(tn), ypcur contains that for y'(tn), and the vector rescur *
 * contains the N-vector residual F(tn,ycur,ypcur).                *
 * The solution is to be returned in the vector b. ida_lsolve      *
 * returns the positive value LSOLVE_ERROR_RECVR for a             *
 * recoverable error and the negative value LSOLVE_ERROR_NONRECVR  *
 * for an unrecoverable error. Success is indicated by a return    *
 * value SUCCESS = 0. The code should include the file ida.h .     *
 *                                                                 *
 *******************************************************************/
/*******************************************************************
 *                                                                 *
 * int (*ida_lperf)(IDAMem IDA_mem, int perftask);                 *
 *                                                                 *
 *-----------------------------------------------------------------*
 * ida_lperf is called two places in IDA where linear solver       *
 * performance data is required by IDA. For perftask = 0, an       *
 * initialization of performance variables is performed, while for *
 * perftask = 1, the performance is evaluated.                     *
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * int (*ida_lfree)(IDAMem IDA_mem);                               *
 *-----------------------------------------------------------------*
 * ida_lfree should free up any memory allocated by the linear     *
 * solver. This routine is called once a problem has been          *
 * completed and the linear solver is no longer needed.            *
 *                                                                 *
 *******************************************************************/

#endif

#ifdef __cplusplus
}
#endif
