/*
 * -----------------------------------------------------------------
 * $Revision: 1.16 $
 * $Date: 2004-04-29 19:16:28 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, Radu Serban
 *              and Dan Shumaker @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * This is the interface file for the main CVODE integrator.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvode_h
#define _cvode_h

#include <stdio.h>
#include "sundialstypes.h"
#include "nvector.h"

/******************************************************************
 *                                                                *
 * CVODE is used to solve numerically the ordinary initial value  *
 * problem :                                                      *
 *                                                                *
 *                 y' = f(t,y),                                   *
 *                 y(t0) = y0,                                    *
 *                                                                *
 *  where t0, y0 in R^N, and f: R x R^N -> R^N are given.         *
 *                                                                *
 ******************************************************************/

 
/******************************************************************
 *                                                                *
 * Enumerations for inputs to CVodeMalloc, CVodeReInit, and CVode.*
 *----------------------------------------------------------------*
 * Symbolic constants for the lmm, iter, and itol input           *
 * parameters to CVodeMalloc and CVodeReInit, as well as the input*
 * parameter itask to CVode, are given below.                     *
 *                                                                *
 * lmm  : The user of the CVODE package specifies whether to use  *
 *        the ADAMS or BDF (backward differentiation formula)     *
 *        linear multistep method. The BDF method is recommended  *
 *        for stiff problems, and the ADAMS method is recommended *
 *        for nonstiff problems.                                  *
 *                                                                *
 * iter : At each internal time step, a nonlinear equation must   *
 *        be solved. The user can specify either FUNCTIONAL       *
 *        iteration, which does not require linear algebra, or a  *
 *        NEWTON iteration, which requires the solution of linear *
 *        systems. In the NEWTON case, the user also specifies a  *
 *        CVODE linear solver. NEWTON is recommended in case of   *
 *        stiff problems.                                         *
 *                                                                *
 * itol : This parameter specifies the relative and absolute      *
 *        tolerance types to be used. The SS tolerance type means *
 *        a scalar relative and absolute tolerance, while the SV  *
 *        tolerance type means a scalar relative tolerance and a  *
 *        vector absolute tolerance (a potentially different      *
 *        absolute tolerance for each vector component).          *
 *                                                                *
 * itask : The itask input parameter to CVode indicates the job   *
 *         of the solver for the next user step. The NORMAL       *
 *         itask is to have the solver take internal steps until  *
 *         it has reached or just passed the user specified tout  *
 *         parameter. The solver then interpolates in order to    *
 *         return an approximate value of y(tout). The ONE_STEP   *
 *         option tells the solver to just take one internal step *
 *         and return the solution at the point reached by that   *
 *         step. The NORMAL_TSTOP and ONE_STEP_TSTOP modes are    *
 *         similar to NORMAL and ONE_STEP, respectively, except   *
 *         that the integration never proceeds past the value     *
 *         tstop (specified through the routine CVodeSetStopTime).*
 *                                                                *
 ******************************************************************/

enum { ADAMS=1, BDF=2 };                                   /* lmm */

enum { FUNCTIONAL=1, NEWTON=2 };                          /* iter */

enum { SS=1, SV=2 };                                      /* itol */

enum { NORMAL=1, ONE_STEP=2, 
       NORMAL_TSTOP=3, ONE_STEP_TSTOP=4 };               /* itask */

/******************************************************************
 *                                                                *
 * Type : RhsFn                                                   *
 *----------------------------------------------------------------*        
 * The f function which defines the right hand side of the ODE    *
 * system y' = f(t,y) must have type RhsFn.                       *
 * f takes as input the independent variable value t, and the     *
 * dependent variable vector y.  It stores the result of f(t,y)   *
 * in the vector ydot.  The y and ydot arguments are of type      *
 * N_Vector.(Allocation of memory for ydot is handled within      *
 * CVODE). The f_data parameter is the same as the f_data         *
 * parameter passed by the user to the CVodeMalloc routine. This  *
 * user-supplied pointer is passed to the user's f function       *
 * every time it is called.                                       *
 * A RhsFn f does not have a return value.                        *
 *                                                                *
 ******************************************************************/

typedef void (*RhsFn)(realtype t, N_Vector y, N_Vector ydot, 
                      void *f_data);
 
/******************************************************************
 *                                                                *
 * Type : RootFn                                                  *
 *----------------------------------------------------------------*
 * A function g, which defines a set of functions g_i(t,y) whose  *
 * roots are sought during the integration, must have type RootFn.*
 * The function g takes as input the independent variable value   *
 * t, and the dependent variable vector y.  It stores the nrtfn   *
 * values g_i(t,y) in the realtype array gout.                    *
 * (Allocation of memory for gout is handled within CVODE.)       *
 * The g_data parameter is the same as that passed by the user    *
 * to the CVodeSetGdata routine.  This user-supplied pointer is   *
 * passed to the user's g function every time it is called.       *
 * A RootFn g does not have a return value.                       *
 *                                                                *
 ******************************************************************/

typedef void (*RootFn)(realtype t, N_Vector y, realtype *gout, 
                       void *g_data);
 
/*================================================================*
 *                                                                *
 *          U S E R - C A L L A B L E   R O U T I N E S           *
 *                                                                *
 *================================================================*/

/*----------------------------------------------------------------*
 * Function : CVodeCreate                                         *
 *----------------------------------------------------------------*
 * CVodeCreate creates an internal memory block for a problem to  *
 * be solved by CVODE.                                            *
 *                                                                *
 * lmm     is the type of linear multistep method to be used.     *
 *            The legal values are ADAMS and BDF (see previous    *
 *            description).                                       *
 *                                                                *
 * iter    is the type of iteration used to solve the nonlinear   *
 *            system that arises during each internal time step.  *
 *            The legal values are FUNCTIONAL and NEWTON.         *
 *                                                                *
 * If successful, CVodeCreate returns a pointer to initialized    *
 * problem memory. This pointer should be passed to CVodeMalloc.  *
 * If an initialization error occurs, CVodeCreate prints an error *
 * message to standard err and returns NULL.                      *
 *----------------------------------------------------------------*/

void *CVodeCreate(int lmm, int iter);

/*----------------------------------------------------------------*
 * Function : CVodeResetIterType                                  *
 *----------------------------------------------------------------*
 * CVodeResetIterType changes the cuurent nonlinear iteration     *
 * type. The legal values for iter are FUNCTIONAL or NEWTON.      *
 *                                                                *
 * If successful, CVodeResetIterType returns SUCCESS. Otherwise,  *
 * it returns one of the error values defined below for the       *
 * optional input specification routines.                         *
 *----------------------------------------------------------------*/

int CVodeResetIterType(void *cvode_mem, int iter);

/*----------------------------------------------------------------*
 * Integrator optional input specification functions              *
 *----------------------------------------------------------------*
 * The following functions can be called to set optional inputs   *
 * to values other than the defaults given below:                 *
 *                                                                *
 * Function             |  Optional input / [ default value ]     *
 * -------------------------------------------------------------- *
 *                      |                                         * 
 * CVodeSetErrFile      | the file pointer for an error file      *
 *                      | where all CVODE warning and error       *
 *                      | messages will be written. This parameter*
 *                      | can be stdout (standard output), stderr *
 *                      | (standard error), a file pointer        *
 *                      | (corresponding to a user error file     *
 *                      | opened for writing) returned by fopen.  *
 *                      | If not called, then all messages will   *
 *                      | be written to standard output.          *
 *                      | [stderr]                                *
 *                      |                                         * 
 * CVodeSetFdata        | a pointer to user data that will be     *
 *                      | passed to the user's f function every   *
 *                      | time f is called.                       *
 *                      | [NULL]                                  *
 *                      |                                         * 
 * CVodeSetGdata        | a pointer to user data that will be     *
 *                      | passed to the user's g function every   *
 *                      | time g is called.                       *
 *                      | [NULL]                                  *
 *                      |                                         * 
 * CVodeSetMaxOrd       | maximum lmm order to be used by the     *
 *                      | solver.                                 *
 *                      | [12 for Adams , 5 for BDF]              * 
 *                      |                                         * 
 * CVodeSetMaxNumSteps  | maximum number of internal steps to be  *
 *                      | taken by the solver in its attempt to   *
 *                      | reach tout.                             *
 *                      | [500]                                   *
 *                      |                                         * 
 * CVodeSetMaxHnilWarns | maximum number of warning messages      *
 *                      | issued by the solver that t+h==t on the *
 *                      | next internal step. A value of -1 means *
 *                      | no such messages are issued.            *
 *                      | [10]                                    * 
 *                      |                                         * 
 * CVodeSetStabLimDet   | flag to turn on/off stability limit     *
 *                      | detection (TRUE = on, FALSE = off).     *
 *                      | When BDF is used and order is 3 or      *
 *                      | greater, CVsldet is called to detect    *
 *                      | stability limit.  If limit is detected, *
 *                      | the order is reduced.                   *
 *                      | [FALSE]                                 *
 *                      |                                         * 
 * CVodeSetInitStep     | initial step size.                      *
 *                      | [estimated by CVODE]                    * 
 *                      |                                         * 
 * CVodeSetMinStep      | minimum absolute value of step size     *
 *                      | allowed.                                *
 *                      | [0.0]                                   *
 *                      |                                         * 
 * CVodeSetMaxStep      | maximum absolute value of step size     *
 *                      | allowed.                                *
 *                      | [infinity]                              *
 *                      |                                         * 
 * CVodeSetStopTime     | the independent variable value past     *
 *                      | which the solution is not to proceed.   *
 *                      | [infinity]                              *
 *                       \                                        * 
 *                        \                                       * 
 * CVodeSetMaxErrTestFails | Maximum number of error test failures*
 *                         | in attempting one step.              *
 *                         | [7]                                  *
 *                         |                                      *
 * CVodeSetMaxNonlinIters  | Maximum number of nonlinear solver   *
 *                         | iterations at one solution.          *
 *                         | [3]                                  *
 *                         |                                      *
 * CVodeSetMaxConvFails    | Maximum number of allowable conv.    *
 *                         | failures in attempting one step.     *
 *                         | [10]                                 *
 *                         |                                      *
 * CVodeSetNonlinConvCoef  | Coeficient in the nonlinear conv.    *
 *                         | test.                                *
 *                         | [0.1]                                *
 *                         |                                      *
 * -------------------------------------------------------------- *
 * If successful, these functions return SUCCESS. If an argument  *
 * has an illegal value, they print an error message to the       *
 * file specified by errfp and return one of the error flags      *  
 * defined below.                                                 *
 *----------------------------------------------------------------*/

int CVodeSetErrFile(void *cvode_mem, FILE *errfp);
int CVodeSetFdata(void *cvode_mem, void *f_data);
int CVodeSetGdata(void *cvode_mem, void *g_data);
int CVodeSetMaxOrd(void *cvode_mem, int maxord);
int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps);
int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil);
int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet);
int CVodeSetInitStep(void *cvode_mem, realtype hin);
int CVodeSetMinStep(void *cvode_mem, realtype hmin);
int CVodeSetMaxStep(void *cvode_mem, realtype hmax);
int CVodeSetStopTime(void *cvode_mem, realtype tstop);
int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef);
int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor);
int CVodeSetMaxConvFails(void *cvode_mem, int maxncf);
int CVodeSetNonlinConvCoef(void *cvode_mem, realtype nlscoef);

/* Error return values for CVodeSet* functions */
/* SUCCESS = 0*/
enum {CVS_NO_MEM = -1, CVS_ILL_INPUT = -2};

/*----------------------------------------------------------------*
 * Function : CVodeMalloc                                         *
 *----------------------------------------------------------------*
 * CVodeMalloc allocates and initializes memory for a problem to  *
 * to be solved by CVODE.                                         *
 *                                                                *
 * cvode_mem is pointer to CVODE memory returned by CVodeCreate.  *
 *                                                                *
 * f       is the right hand side function in y' = f(t,y).        *          
 *                                                                *
 * t0      is the initial value of t.                             *
 *                                                                *
 * y0      is the initial condition vector y(t0).                 *
 *                                                                *
 * itol    is the type of tolerances to be used.                  *
 *            The legal values are:                               *
 *               SS (scalar relative and absolute  tolerances),   *
 *               SV (scalar relative tolerance and vector         *
 *                   absolute tolerance).                         *
 *                                                                *
 * reltol  is a pointer to the relative tolerance scalar.         *
 *                                                                *
 * abstol  is a pointer to the absolute tolerance scalar or       *
 *            an N_Vector of absolute tolerances.                 *
 *                                                                *
 * nvspec  is a pointer to a vector specification structure       *
 *                                                                *
 * The parameters itol, reltol, and abstol define a vector of     *
 * error weights, ewt, with components                            *
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)   (if itol = SS), or  *
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (if itol = SV).  *
 * This vector is used in all error and convergence tests, which  *
 * use a weighted RMS norm on all error-like vectors v:           *
 *    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),    *
 * where N is the problem dimension.                              *
 *                                                                *
 * Note: The tolerance values may be changed in between calls to  *
 *       CVode for the same problem. These values refer to        *
 *       (*reltol) and either (*abstol), for a scalar absolute    *
 *       tolerance, or the components of abstol, for a vector     *
 *       absolute tolerance.                                      *
 *                                                                * 
 * If successful, CVodeMalloc returns SUCCESS. If an argument has *
 * an illegal value, CVodeMalloc prints an error message to the   *
 * file specified by errfp and returns one of the error flags     *  
 * defined below.                                                 *
 *----------------------------------------------------------------*/

int CVodeMalloc(void *cvode_mem, RhsFn f,
                realtype t0, N_Vector y0, 
                int itol, realtype *reltol, void *abstol, 
                NV_Spec nvspec);

/* Error return values for CVodeMalloc */
/* SUCCESS = 0 */
enum {CVM_NO_MEM = -1, CVM_MEM_FAIL=-2, CVM_ILL_INPUT = -3};

/*----------------------------------------------------------------*
 * Function : CVodeReInit                                         *
 *----------------------------------------------------------------*
 * CVodeReInit re-initializes CVode for the solution of a problem,*
 * where a prior call to CVodeMalloc has been made with the same  *
 * problem size N. CVodeReInit performs the same input checking   *
 * and initializations that CVodeMalloc does.                     *
 * But it does no memory allocation, assuming that the existing   *
 * internal memory is sufficient for the new problem.             *
 *                                                                *
 * The use of CVodeReInit requires that the maximum method order, *
 * maxord, is no larger for the new problem than for the problem  *
 * specified in the last call to CVodeMalloc.  This condition is  *
 * automatically fulfilled if the multistep method parameter lmm  *
 * is unchanged (or changed from ADAMS to BDF) and the default    *
 * value for maxord is specified.                                 *
 *                                                                *
 * If iter=NEWTON, then following the call to CVodeReInit, a call *
 * to the linear solver specification routine is necessary if a   *
 * different linear solver is chosen, but may not be otherwise.   *
 * If the same linear solver is chosen, and there are no changes  *
 * in the input parameters to the specification routine, then no  *
 * call to that routine is needed. Similarly for the optional     *
 * inputs to the linear solver.                                   *
 * If there are changes in parameters, but they do not increase   *
 * the linear solver memory size, then a call to the corresponding*
 * CVReInit<linsol> routine must made to communicate the new      *
 * parameters; in that case the linear solver memory is reused.   *
 * If the parameter changes do increase the linear solver memory  *
 * size, then the main linear solver specification routine must be*
 * called.  See the linear solver documentation for full details. *
 *                                                                *
 * The first argument to CVodeReInit is:                          *
 *                                                                *
 * cvode_mem = pointer to CVODE memory returned by CVodeCreate.   *
 *                                                                *
 * All the remaining arguments to CVodeReInit have names and      *
 * meanings identical to those of CVodeMalloc.                    *
 *                                                                *
 * The return value of CVodeReInit is equal to SUCCESS = 0 if     *
 * there were no errors; otherwise it is a negative int equal to: *
 *   CVREI_NO_MEM     indicating cvode_mem was NULL, or           *
 *   CVREI_ILL_INPUT  indicating an input argument was illegal    *
 *                    (including an attempt to increase maxord).  *
 * In case of an error return, an error message is also printed.  *
 *----------------------------------------------------------------*/

int CVodeReInit(void *cvode_mem, RhsFn f,
                realtype t0, N_Vector y0, 
                int itol, realtype *reltol, void *abstol);

/* CVodeReInit return values: */
/* SUCCESS = 0 */ 
enum {CVREI_NO_MEM = -1, CVREI_NO_MALLOC = -2, CVREI_ILL_INPUT = -3};

/*----------------------------------------------------------------*
 * Function : CVodeRootInit                                       *
 *----------------------------------------------------------------*
 * CVodeRootInit initializes a rootfinding problem to be solved   *
 * during the integration of the ODE system.  It must be called   *
 * after CVodeCreate, and before CVode.  The arguments are:       *
 *                                                                *
 * cvode_mem = pointer to CVODE memory returned by CVodeCreate.   *
 *                                                                *
 * g         = name of user-supplied function, of type RootFn,    *
 *             defining the functions g_i whose roots are sought. *
 *                                                                *
 * nrtfn     = number of functions g_i, an int >= 0.              *
 *                                                                *
 * If a new problem is to be solved with a call to CVodeReInit,   *
 * where the new problme has no root functions but the prior one  *
 * did, then call CVodeRootInit with nrtfn = 0.                   *
 *                                                                *
 * The return value of CVodeRootInit is SUCCESS = 0 if there were *
 * no errors; otherwise it is a negative int equal to:            *
 *   CVRT_NO_MEM     indicating cvode_mem was NULL, or            *
 *   CVRT_MEM_FAIL   indicating a memory allocation failed.       *
 *                    (including an attempt to increase maxord).  *
 * In case of an error return, an error message is also printed.  *
 *----------------------------------------------------------------*/

int CVodeRootInit(void *cvode_mem, RootFn g, int nrtfn);
/* CVodeRootInit return values: */
/* SUCCESS = 0 */ 
enum {CVRT_NO_MEM = -1, CVRT_MEM_FAIL = -2};

/******************************************************************
 * Function : CVode                                               *
 *----------------------------------------------------------------*
 * CVode integrates the ODE over an interval in t.                *
 * If itask is NORMAL, then the solver integrates from its        *
 * current internal t value to a point at or beyond tout, then    *
 * interpolates to t = tout and returns y(tout) in the user-      *
 * allocated vector yout. If itask is ONE_STEP, then the solver   *
 * takes one internal time step and returns in yout the value of  *
 * y at the new internal time. In this case, tout is used only    *
 * during the first call to CVode to determine the direction of   *
 * integration and the rough scale of the problem. In either      *
 * case, the time reached by the solver is placed in (*t). The    *
 * user is responsible for allocating the memory for this value.  *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *              CVodeCreate.                                      *
 *                                                                *
 * tout  is the next time at which a computed solution is desired *
 *                                                                *
 * yout  is the computed solution vector. In NORMAL mode with no  *
 *          errors and no roots found, yout = y(tout).            *
 *                                                                *
 * tret  is a pointer to a realtype location. CVode sets (*tret)  *
 *          to the time reached by the solver and returns         *
 *          yout = y(*tret).                                      *
 *                                                                *
 * itask is NORMAL, ONE_STEP, NORMAL_TSTOP, or ONE_STEP_TSTOP.    *
 *          These four modes are described above.                 *
 *                                                                *
 * The return values for CVode are defined later in this file.    *
 * Here is a brief description of each return value:              *
 *                                                                *
 * SUCCESS       : CVode succeeded, and no roots were found.      *
 *                                                                *
 * ROOT_RETURN   : CVode succeeded, and found one or more roots.  *
 *                 If nrtfn > 1, call CVodeGetRootInfo to see     *
 *                 which g_i were found to have a root at (*tret).*
 *                                                                *
 * TSTOP_RETURN  : CVode succeeded and returned at tstop.         *
 *                                                                *
 * CVODE_NO_MEM  : The cvode_mem argument was NULL.               *
 *                                                                *
 * CVODE_NO_MALLOC: cvode_mem was not allocated.                  *
 *                                                                *
 * ILL_INPUT     : One of the inputs to CVode is illegal. This    *
 *                 includes the situation when a component of the *
 *                 error weight vectors becomes < 0 during        *
 *                 internal time-stepping. The ILL_INPUT flag     *
 *                 will also be returned if the linear solver     *
 *                 routine CV--- (called by the user after        *
 *                 calling CVodeMalloc) failed to set one of the  *
 *                 linear solver-related fields in cvode_mem or   *
 *                 if the linear solver's init routine failed. In *
 *                 any case, the user should see the printed      *
 *                 error message for more details.                *
 *                                                                *
 * TOO_MUCH_WORK : The solver took mxstep internal steps but      *
 *                 could not reach tout. The default value for    *
 *                 mxstep is MXSTEP_DEFAULT = 500.                *
 *                                                                *
 * TOO_MUCH_ACC  : The solver could not satisfy the accuracy      *
 *                 demanded by the user for some internal step.   *
 *                                                                *
 * ERR_FAILURE   : Error test failures occurred too many times    *
 *                 (= MXNEF = 7) during one internal time step or *
 *                 occurred with |h| = hmin.                      *
 *                                                                *
 * CONV_FAILURE  : Convergence test failures occurred too many    *
 *                 times (= MXNCF = 10) during one internal time  *
 *                 step or occurred with |h| = hmin.              *
 *                                                                *
 * SETUP_FAILURE : The linear solver's setup routine failed in an *
 *                 unrecoverable manner.                          *
 *                                                                *
 * SOLVE_FAILURE : The linear solver's solve routine failed in an *
 *                 unrecoverable manner.                          *
 ******************************************************************/

int CVode(void *cvode_mem, realtype tout, N_Vector yout, 
          realtype *tret, int itask);

/* CVode return values */
enum { SUCCESS=0, TSTOP_RETURN=1, ROOT_RETURN=2,
       CVODE_NO_MEM=-1, CVODE_NO_MALLOC=-2, ILL_INPUT=-3,
       TOO_MUCH_WORK=-4, TOO_MUCH_ACC=-5, ERR_FAILURE=-6, 
       CONV_FAILURE=-7, SETUP_FAILURE=-8, SOLVE_FAILURE=-9 };

/*----------------------------------------------------------------*
 * Function : CVodeGetDky                                         *
 *----------------------------------------------------------------*
 * CVodeGetDky computes the kth derivative of the y function at   *
 * time t, where tn-hu <= t <= tn, tn denotes the current         *
 * internal time reached, and hu is the last internal step size   *
 * successfully used by the solver. The user may request          *
 * k=0, 1, ..., qu, where qu is the current order. The            *
 * derivative vector is returned in dky. This vector must be      *
 * allocated by the caller. It is only legal to call this         *
 * function after a successful return from CVode.                 *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *              CVodeCreate.                                      *
 *                                                                *
 * t   is the time at which the kth derivative of y is evaluated. *
 *        The legal range for t is [tn-hu,tn] as described above. *
 *                                                                *
 * k   is the order of the derivative of y to be computed. The    *
 *        legal range for k is [0,qu] as described above.         *
 *                                                                *
 * dky is the output derivative vector [(D_k)y](t).               *
 *                                                                *
 * The return values for CVodeGetDky are defined below.           *
 * Here is a brief description of each return value:              *
 *                                                                *
 * OKAY : CVodeGetDky succeeded.                                  *
 *                                                                *
 * BAD_K : k is not in the range 0, 1, ..., qu.                   *
 *                                                                *
 * BAD_T : t is not in the interval [tn-hu,tn].                   *
 *                                                                *
 * BAD_DKY : The dky argument was NULL.                           *
 *                                                                *
 * DKY_NO_MEM : The cvode_mem argument was NULL.                  *
 *----------------------------------------------------------------*/

int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky);

/*----------------------------------------------------------------*
 * Integrator optional output extraction functions                *
 *----------------------------------------------------------------*
 * The following functions can be called to get optional outputs  *
 * and statistics related to the main integrator.                 *
 * -------------------------------------------------------------- *
 *                                                                *
 * CVodeGetIntWorkSpace returns the CVODE integer workspace size  *
 * CVodeGetRealWorkSpace returns the CVODE real workspace size    *
 * CVodeGetNumSteps returns the cumulative number of internal     *
 *       steps taken by the solver                                *
 * CVodeGetNumRhsEvals returns the number of calls to the user's  *
 *       f function                                               *
 * CVodeGetNumLinSolvSetups returns the number of calls made to   *
 *       the linear solver's setup routine                        *
 * CVodeGetNumErrTestFails returns the number of local error test *
 *       failures that have occured                               *
 * CVodeGetLastOrder returns the order used during the last       *
 *       internal step                                            *
 * CVodeGetCurrentOrder returns the order to be used on the next  *
 *       internal step                                            *
 * CVodeGetNumStabLimOrderReds returns the number of order        *
 *       reductions due to stability limit detection              *
 * CVodeGetActualInitStep returns the actual initial step size    *
 *       used by CVODE                                            *
 * CVodeGetLastStep returns the step size for the last internal   *
 *       step                                                     *
 * CVodeGetCurrentStep returns the step size to be attempted on   *
 *       the next internal step                                   *
 * CVodeGetCurrentTime returns the current internal time reached  *
 *       by the solver                                            *
 * CVodeGetTolScaleFactor returns a suggested factor by which the *
 *       user's tolerances should be scaled when too much         *
 *       accuracy has been requested for some internal step       *
 * CVodeGetErrWeights returns the state error weight vector.      *
 *       The user need not allocate space for ewt.                *
 * CVodeGetEstLocalErrors returns the vector of estimated local   *
 *       errors. The user need not allocate space for ele.        *
 * CVodeGetNumGEvals returns the number of calls to the user's    *
 *       g function (for rootfinding)                             *
 * CVodeGetRootInfo returns an array of int's showing the indices *
 *       for which g_i was found to have a root.                  *
 *       For i = 0 ... nrtfn-1, rootsfound[i] = 1 if g_i has a    *
 *       root, and = 0 if not.                                    *
 *----------------------------------------------------------------*/

int CVodeGetIntWorkSpace(void *cvode_mem, long int *leniw);
int CVodeGetRealWorkSpace(void *cvode_mem, long int *lenrw);
int CVodeGetNumSteps(void *cvode_mem, long int *nsteps);
int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevals);
int CVodeGetNumLinSolvSetups(void *cvode_mem, long int *nlinsetups);
int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails);
int CVodeGetLastOrder(void *cvode_mem, int *qlast);
int CVodeGetCurrentOrder(void *cvode_mem, int *qcur);
int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred);
int CVodeGetActualInitStep(void *cvode_mem, realtype *hinused);
int CVodeGetLastStep(void *cvode_mem, realtype *hlast);
int CVodeGetCurrentStep(void *cvode_mem, realtype *hcur);
int CVodeGetCurrentTime(void *cvode_mem, realtype *tcur);
int CVodeGetTolScaleFactor(void *cvode_mem, realtype *tolsfact);
int CVodeGetErrWeights(void *cvode_mem, N_Vector *eweight);
int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector *ele);
int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals);
int CVodeGetRootInfo(void *cvode_mem, int **rootsfound);
 
/*----------------------------------------------------------------*
 * As a convenience, the following two functions provide the      *
 * optional outputs in groups.                                    *
 *----------------------------------------------------------------*/

int CVodeGetWorkSpace(void *cvode_mem, long int *leniw, long int *lenrw);
int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps, 
                            long int *nfevals, long int *nlinsetups, 
                            long int *netfails, int *qlast, int *qcur, 
                            realtype *hinused, realtype *hlast, 
                            realtype *hcur, realtype *tcur);

/*----------------------------------------------------------------*
 * Nonlinear solver optional output extraction functions          *
 *----------------------------------------------------------------*
 *                                                                *
 * The following functions can be called to get optional outputs  *
 * and statistics related to the nonlinear solver.                *
 * -------------------------------------------------------------- *
 *                                                                *
 * CVodeGetNumNonlinSolvIters returns the number of nonlinear     *
 *       solver iterations performed.                             *
 * CVodeGetNumNonlinSolvConvFails returns the number of nonlinear *
 *       convergence failures.                                    *
 *----------------------------------------------------------------*/

int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters);
int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails);

/*----------------------------------------------------------------*
 * As a convenience, the following function provides the          *
 * optional outputs in a group.                                   *
 *----------------------------------------------------------------*/

int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, 
                            long int *nncfails);

/* CVodeGet* return values */
enum { OKAY=0, CVG_NO_MEM=-1, CVG_NO_SLDET=-2, 
       BAD_K=-3, BAD_T=-4, BAD_DKY=-5 };
 
/******************************************************************
 * Function : CVodeFree                                           *
 *----------------------------------------------------------------*
 * CVodeFree frees the problem memory cvode_mem allocated by      *
 * CVodeMalloc.  Its only argument is the pointer cvode_mem       *
 * returned by CVodeMalloc.                                       *
 ******************************************************************/

void CVodeFree(void *cvode_mem);
 
/*================================================================*
 *                                                                *
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K   *
 *                                                                *
 *================================================================*/

/* Basic CVODE constants */

#define ADAMS_Q_MAX 12     /* max value of q for lmm == ADAMS     */
#define BDF_Q_MAX    5     /* max value of q for lmm == BDF       */
#define Q_MAX  ADAMS_Q_MAX /* max value of q for either lmm       */
#define L_MAX  (Q_MAX+1)   /* max value of L for either lmm       */
#define NUM_TESTS    5     /* number of error test quantities     */

/******************************************************************
 * Types : struct CVodeMemRec, CVodeMem                           *
 *----------------------------------------------------------------*
 * The type CVodeMem is type pointer to struct CVodeMemRec. This  *
 * structure contains fields to keep track of problem state.      *
 ******************************************************************/

typedef struct CVodeMemRec {

  realtype cv_uround;    /* machine unit roundoff */

  /*---------------------------- 
    Problem Specification Data 
  ----------------------------*/

  RhsFn cv_f;          /* y' = f(t,y(t))              */
  void *cv_f_data;     /* user pointer passed to f    */
  int cv_lmm;          /* lmm = ADAMS or BDF          */
  int cv_iter;         /* iter = FUNCTIONAL or NEWTON */
  int cv_itol;         /* itol = SS or SV             */
  realtype *cv_reltol; /* ptr to relative tolerance   */
  void *cv_abstol;     /* ptr to absolute tolerance   */

  /*-------------------------
    Nordsieck History Array 
  -------------------------*/

  N_Vector cv_zn[L_MAX];  /* Nordsieck array, of size N x (q+1).         */
                          /* zn[j] is a vector of length N (j=0,...,q)   */
                          /* zn[j] = [1/factorial(j)] * h^j * (jth       */ 
                          /* derivative of the interpolating polynomial  */

  /*---------------------
    Vectors of length N 
  ---------------------*/

  N_Vector cv_ewt;     /* error weight vector                          */
  N_Vector cv_y;       /* y is used as temporary storage by the solver */
                       /* The memory is provided by the user to CVode  */
                       /* where the vector is named yout.              */
  N_Vector cv_acor;    /* In the context of the solution of the        */
                       /* nonlinear equation, acor = y_n(m) - y_n(0).  */
                       /* On return, this vector is scaled to give     */
                       /* the estimated local error in y.              */
  N_Vector cv_tempv;   /* temporary storage vector                     */
  N_Vector cv_ftemp;   /* temporary storage vector                     */

  /*-----------------
    Tstop information
  -------------------*/
  booleantype cv_tstopset;
  realtype cv_tstop;

  /*-----------
    Step Data 
  -----------*/  

  int cv_q;         /* current order                           */
  int cv_qprime;    /* order to be used on the next step       */ 
                    /* = q-1, q, or q+1                        */
  int cv_next_q;    /* order to be used on the next step       */
  int cv_qwait;     /* number of internal steps to wait before */
                    /* considering a change in q               */
  int cv_L;         /* L = q + 1                               */

  realtype cv_hin;
  realtype cv_h;      /* current step size                     */
  realtype cv_hprime; /* step size to be used on the next step */ 
  realtype cv_next_h; /* step size to be used on the next step */ 
  realtype cv_eta;    /* eta = hprime / h                      */
  realtype cv_hscale; /* value of h used in zn                 */
  realtype cv_tn;     /* current internal value of t           */

  realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step     */
                               /* sizes indexed from 1 to q+1               */
  realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from     */
                               /* 1 to NUM_TESTS(=5)                        */
  realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)      */

  realtype cv_rl1;     /* 1 / l[1]                     */
  realtype cv_gamma;   /* gamma = h * rl1              */
  realtype cv_gammap;  /* gamma at the last setup call */
  realtype cv_gamrat;  /* gamma / gammap               */

  realtype cv_crate;   /* estimated corrector convergence rate     */
  realtype cv_acnrm;   /* | acor | wrms                            */
  realtype cv_nlscoef; /* coeficient in nonlinear convergence test */
  int  cv_mnewt;       /* Newton iteration counter                 */

  /*--------
    Limits 
  --------*/

  int cv_qmax;        /* q <= qmax                                          */
  long int cv_mxstep; /* maximum number of internal steps for one user call */
  int cv_maxcor;      /* maximum number of corrector iterations for the     */
                      /* solution of the nonlinear equation                 */
  int cv_mxhnil;      /* maximum number of warning messages issued to the   */
                      /* user that t + h == t for the next internal step    */
  int cv_maxnef;      /* maximum number of error test failures              */
  int cv_maxncf;      /* maximum number of nonlinear convergence failures   */

  realtype cv_hmin;     /* |h| >= hmin       */
  realtype cv_hmax_inv; /* |h| <= 1/hmax_inv */
  realtype cv_etamax;   /* eta <= etamax     */

  /*----------
    Counters 
  ----------*/

  long int cv_nst;              /* number of internal steps taken             */
  long int cv_nfe;              /* number of f calls                          */
  long int cv_ncfn;             /* number of corrector convergence failures   */
  long int cv_netf;             /* number of error test failures              */
  long int cv_nni;              /* number of Newton iterations performed      */
  long int cv_nsetups;          /* number of setup calls                      */
  int cv_nhnil;                 /* number of messages issued to the user that */
                                /* t + h == t for the next iternal step       */

  realtype cv_etaqm1;      /* ratio of new to old h for order q-1        */
  realtype cv_etaq;        /* ratio of new to old h for order q          */
  realtype cv_etaqp1;      /* ratio of new to old h for order q+1        */

  /*------------------------------- 
    Space requirements for CVODE 
  -------------------------------*/

  long int cv_lrw1;        /* no. of realtype words in 1 N_Vector         */ 
  long int cv_liw1;        /* no. of integer words in 1 N_Vector          */ 
  long int cv_lrw;         /* no. of realtype words in CVODE work vectors */
  long int cv_liw;         /* no. of integer words in CVODE work vectors  */

  /*--------------------
    Linear Solver Data 
  --------------------*/

  /* Linear Solver functions to be called */

  int (*cv_linit)(struct CVodeMemRec *cv_mem);

  int (*cv_lsetup)(struct CVodeMemRec *cv_mem, int convfail, N_Vector ypred,
                   N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                   N_Vector vtemp2, N_Vector vtemp3); 

  int (*cv_lsolve)(struct CVodeMemRec *cv_mem, N_Vector b, N_Vector weight,
                   N_Vector ycur, N_Vector fcur);

  void (*cv_lfree)(struct CVodeMemRec *cv_mem);

  /* Linear Solver specific memory */

  void *cv_lmem;           

  /*--------------
    Saved Values
  --------------*/

  int cv_qu;             /* last successful q value used   */
  long int cv_nstlp;          /* step number of last setup call */
  realtype cv_h0u;       /* actual initial stepsize        */
  realtype cv_hu;        /* last successful h value used   */
  realtype cv_saved_tq5; /* saved value of tq[5]           */
  booleantype cv_jcur;   /* Is the Jacobian info used by   */
                         /* linear solver current?         */
  realtype cv_tolsf;     /* tolerance scale factor         */
  booleantype cv_setupNonNull; /* Does setup do something? */

  booleantype cv_MallocDone;  

  /*------------
    Error File 
  ------------*/

  FILE *cv_errfp;       /* CVODE error messages are sent to errfp */

  /*-----------------------------------------------------
    Pointer to the vector specification structure for
    state N_Vectors
  -----------------------------------------------------*/

  NV_Spec cv_nvspec;

  /*-------------------------
    Stability Limit Detection
  ---------------------------*/

  booleantype cv_sldeton;     /* Is Stability Limit Detection on?          */
  realtype cv_ssdat[6][4];    /* scaled data array for STALD               */
  int cv_nscon;               /* counter for STALD method                  */
  long int cv_nor;            /* counter for number of order reductions    */

  /*----------------
    Rootfinding Data
  ------------------*/

  RootFn cv_gfun;       /* Function g for roots sought                     */
  int cv_nrtfn;         /* number of components of g                       */
  void *cv_g_data;      /* pointer to user data for g                      */
  int *cv_iroots;       /* int array for root information                  */
  realtype cv_tlo;      /* nearest endpoint of interval in root search     */
  realtype cv_thi;      /* farthest endpoint of interval in root search    */
  realtype cv_troot;    /* approximate root location                       */
  realtype *cv_glo;     /* saved array of g values at t = tlo              */
  realtype *cv_ghi;     /* saved array of g values at t = thi              */
  realtype *cv_groot;   /* array of g values at t = troot                  */
  realtype cv_tretlast; /* last value of t returned                        */
  realtype cv_toutc;    /* copy of tout (if NORMAL mode)                   */
  realtype cv_ttol;     /* tolerance on root location troot                */
  int cv_taskc;         /* copy of parameter task                          */
  int cv_irfnd;         /* flag showing whether last step had a root       */
  int cv_nge;           /* counter for g evaluations                       */


} *CVodeMem;


/******************************************************************
 * Communication between user and a CVODE Linear Solver           *
 *----------------------------------------------------------------*
 * Return values of the linear solver specification routine.      *
 * The values of these are given in the enum statement below.     *
 *                                                                *
 *    SUCCESS      : The routine was successful.                  *
 *                                                                *
 *    LIN_NO_MEM   : CVODE memory = NULL.                         *
 *                                                                *
 *    LMEM_FAIL    : A memory allocation failed.                  *
 *                                                                *
 *    LIN_ILL_INPUT: Some input was illegal (see message).        *
 *                                                                *
 *    LIN_NO_LMEM  : The linear solver's memory = NULL.           *
 *----------------------------------------------------------------*/

/* SUCCESS = 0  */
enum {LMEM_FAIL=-1, LIN_ILL_INPUT=-2, LIN_NO_MEM=-3, LIN_NO_LMEM=-4};

/******************************************************************
 * Communication between cvode.c and a CVODE Linear Solver        *
 *----------------------------------------------------------------*
 * (1) cv_linit return values                                     *
 *                                                                *
 * LINIT_OK    : The cv_linit routine succeeded.                  *
 *                                                                *
 * LINIT_ERR   : The cv_linit routine failed. Each linear solver  *
 *               init routine should print an appropriate error   *
 *               message to (cv_mem->errfp).                      *
 *                                                                *
 * (2) convfail (input to cv_lsetup)                              *
 *                                                                *
 * NO_FAILURES : Either this is the first cv_setup call for this  *
 *               step, or the local error test failed on the      *
 *               previous attempt at this step (but the Newton    *
 *               iteration converged).                            * 
 *                                                                *
 * FAIL_BAD_J  : This value is passed to cv_lsetup if             *
 *                                                                *
 *               (1) The previous Newton corrector iteration      *
 *                   did not converge and the linear solver's     *
 *                   setup routine indicated that its Jacobian-   *
 *                   related data is not current.                 *
 *                                   or                           *
 *               (2) During the previous Newton corrector         *
 *                   iteration, the linear solver's solve routine *
 *                   failed in a recoverable manner and the       *
 *                   linear solver's setup routine indicated that *
 *                   its Jacobian-related data is not current.    *
 *                                                                *
 * FAIL_OTHER  : During the current internal step try, the        *
 *               previous Newton iteration failed to converge     *
 *               even though the linear solver was using current  *
 *               Jacobian-related data.                           *
 *                                                                *
 * (3) Parameter documentation, as well as a brief description    *
 *     of purpose, for each CVODE linear solver routine to be     *
 *     called in cvode.c is given below the constant declarations *
 *     that follow.                                               *
 ******************************************************************/

/* cv_linit return values */

#define LINIT_OK        0
#define LINIT_ERR      -1

/* Constants for convfail (input to cv_lsetup) */

#define NO_FAILURES 0   
#define FAIL_BAD_J  1  
#define FAIL_OTHER  2  


/*******************************************************************
 * int (*cv_linit)(CVodeMem cv_mem);                               *
 *-----------------------------------------------------------------*
 * The purpose of cv_linit is to complete initializations for      *
 * specific linear solver, such as counters and statistics.        *
 * An LInitFn should return LINIT_OK (= 0) if it has successfully  *
 * initialized the CVODE linear solver and LINIT_ERR (= -1)        *
 * otherwise. These constants are defined above.  If an error does *
 * occur, an appropriate message should be sent to (cv_mem->errfp).*
 *******************************************************************/

/*******************************************************************
 * int (*cv_lsetup)(CVodeMem cv_mem, int convfail, N_Vector ypred, *
 *             N_Vector fpred, booleantype *jcurPtr,               *
 *             N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); *
 *-----------------------------------------------------------------*
 * The job of cv_lsetup is to prepare the linear solver for        *
 * subsequent calls to cv_lsolve. It may re-compute Jacobian-      *
 * related data is it deems necessary. Its parameters are as       *
 * follows:                                                        *
 *                                                                 *
 * cv_mem - problem memory pointer of type CVodeMem. See the big   *
 *          typedef earlier in this file.                          *
 *                                                                 *
 * convfail - a flag to indicate any problem that occurred during  *
 *            the solution of the nonlinear equation on the        *
 *            current time step for which the linear solver is     *
 *            being used. This flag can be used to help decide     *
 *            whether the Jacobian data kept by a CVODE linear     *
 *            solver needs to be updated or not.                   *
 *            Its possible values have been documented above.      *
 *                                                                 *
 * ypred - the predicted y vector for the current CVODE internal   *
 *         step.                                                   *
 *                                                                 *
 * fpred - f(tn, ypred).                                           *
 *                                                                 *
 * jcurPtr - a pointer to a boolean to be filled in by cv_lsetup.  *
 *           The function should set *jcurPtr=TRUE if its Jacobian *
 *           data is current after the call and should set         *
 *           *jcurPtr=FALSE if its Jacobian data is not current.   *
 *           Note: If cv_lsetup calls for re-evaluation of         *
 *           Jacobian data (based on convfail and CVODE state      *
 *           data), it should return *jcurPtr=TRUE unconditionally;*
 *           otherwise an infinite loop can result.                *
 *                                                                 *
 * vtemp1 - temporary N_Vector provided for use by cv_lsetup.      *
 *                                                                 *
 * vtemp3 - temporary N_Vector provided for use by cv_lsetup.      *
 *                                                                 *
 * vtemp3 - temporary N_Vector provided for use by cv_lsetup.      *
 *                                                                 *
 * The cv_lsetup routine should return 0 if successful,            *
 * a positive value for a recoverable error, and a negative value  *
 * for an unrecoverable error.                                     *
 *******************************************************************/

/*******************************************************************
 * int (*cv_lsolve)(CVodeMem cv_mem, N_Vector b, N_Vector ycur,    *
 *                  N_Vector fcur);                                *
 *-----------------------------------------------------------------*
 * cv_lsolve must solve the linear equation P x = b, where         *
 * P is some approximation to (I - gamma J), J = (df/dy)(tn,ycur)  *
 * and the RHS vector b is input. The N-vector ycur contains       *
 * the solver's current approximation to y(tn) and the vector      *
 * fcur contains the N-vector f(tn,ycur). The solution is to be    *
 * returned in the vector b. cv_lsolve returns a positive value    *
 * for a recoverable error and a negative value for an             *
 * unrecoverable error. Success is indicated by a 0 return value.  *
 *******************************************************************/

/*******************************************************************
 * void (*cv_lfree)(CVodeMem cv_mem);                              *
 *-----------------------------------------------------------------*
 * cv_lfree should free up any memory allocated by the linear      *
 * solver. This routine is called once a problem has been          *
 * completed and the linear solver is no longer needed.            *
 *******************************************************************/

#endif

#ifdef __cplusplus
}
#endif
