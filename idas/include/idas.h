/*******************************************************************
 * File          : idas.h                                          *
 * Programmers   : Alan C. Hindmarsh, Radu Serban and              *
 *                 Allan G. Taylor @ LLNL                          *
 * Version of    : 18 September 2003                               *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/idas/LICENSE                          *
 *-----------------------------------------------------------------*
 * This is the header (include) file for the main IDAS solver.     *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idas_h
#define _idas_h

#include <stdio.h>
#include "sundialstypes.h"
#include "nvector.h"

/******************************************************************
 * IDAS is used to solve numerically the initial value problem    *
 * for the differential algebraic equation (DAE) system           *
 *   F(t,y,y') = 0,                                               *
 * given initial conditions                                       *
 *   y(t0) = y0,   y'(t0) = yp0.                                  *
 * Here y and F are vectors of length N.                          *
 *                                                                *
 * Optionally, IDAS can perform forward sensitivity analysis      *
 * to find sensitivities of the solution y with respect to        *
 * parameters in the model function F and/or in the initial       *
 * conditions y0, yp0.                                            * 
 ******************************************************************/

/*----------------------------------------------------------------*
 * Enumerations for inputs to IDAMalloc, IDAReInit, IDACalcIC,    *
 * IDASolve, IDASensMalloc, IDASensReInit, IDAQuadMalloc,         *
 * IDAQuadReInit, and IDASet*.                                    *
 *----------------------------------------------------------------*
 *                                                                *
 * itol:  This parameter specifies the relative and absolute      *
 *        tolerance types to be used. The SS tolerance type means *
 *        a scalar relative and absolute tolerance, while the SV  *
 *        tolerance type means a scalar relative tolerance and a  *
 *        vector absolute tolerance (a potentially different      *
 *        absolute tolerance for each vector component).          *
 *                                                                *
 * itolQ: Same as itol for quadrature variables.                  *
 *                                                                *
 * ism:   This parameter specifies the sensitivity corrector type *
 *        to be used. In the SIMULTANEOUS case, the nonlinear     *
 *        systems for states and all sensitivities are solved     *
 *        simultaneously. In the STAGGERED case, the nonlinear    *
 *        system for states is solved first and then, the         *
 *        nonlinear systems for all sensitivities are solved      *
 *        at the same time. Finally, in the STAGGERED1 approach   *     
 *        all nonlinear systems are solved in a sequence.         *
 *                                                                *
 * ifS:   Type of the function returning the sensitivity right    *
 *        hand side. ifS can be either ALLSENS if the function    *
 *        (of type SensResFn) returns residuals for all           *
 *        sensitivity systems at once, or ONESENS if the function *
 *        (of type SensRes1Fn) returns the residual of one        *
 *        sensitivity system at a time.                           *
 *                                                                *
 * itask: The itask input parameter to IDASolve indicates the job *
 *        of the solver for the next user step. The NORMAL        *
 *        itask is to have the solver take internal steps until   *
 *        it has reached or just passed the user specified tout   *
 *        parameter. The solver then interpolates in order to     *
 *        return an approximate value of y(tout) and yp(tout).    *
 *        The ONE_STEP option tells the solver to just take one   *
 *        internal step and return the solution at the point      *
 *        reached by that step.                                   *
 *        The NORMAL_TSTOP and ONE_STEP_TSTOP modes are similar   *
 *        to NORMAL and ONE_STEP, respectively, except that the   *
 *        integration never proceeds past the value tstop         *
 *        (specified through the routine IDASetStopTime).         *
 *----------------------------------------------------------------*/

enum { SS, SV };                                  /* itol, itolQ */

enum { SIMULTANEOUS, STAGGERED, STAGGERED1 };             /* ism */

enum { NORMAL, ONE_STEP, NORMAL_TSTOP, ONE_STEP_TSTOP}; /* itask */

enum { ALLSENS, ONESENS };                                /* ifS */

enum { CALC_YA_YDP_INIT = 1 , CALC_Y_INIT = 2 };        /* icopt */


/******************************************************************
 * Type : ResFn                                                   *
 *----------------------------------------------------------------*        
 * The F function which defines the DAE system   F(t,y,y')=0      *
 * must have type ResFn.                                          *
 * Symbols are as follows: t  <-> tres     y <-> yy               *
 *                         y' <-> yp       F <-> res (type ResFn) *
 * A ResFn takes as input the independent variable value tres,    *
 * the dependent variable vector yy, and the derivative (with     *
 * respect to t) of the yy vector, yp.  It stores the result of   *
 * F(t,y,y') in the vector resval. The yy, yp, and resval         *
 * arguments are of type N_Vector. The rdata parameter is to be   *
 * of the same type as the rdata parameter passed by the user to  *
 * the IDASetRdata routine. This user-supplied pointer is passed  *
 * to the user's res function every time it is called, to provide *
 * access in res to user data.                                    *
 *                                                                *
 * A ResFn res will return the value ires, which has possible     *
 * values RES_ERROR_RECVR = 1, RES_ERROR_NONRECVR = -1,           *
 * and SUCCESS = 0. The file idas.h may be used to obtain these   *
 * values but is not required; returning 0, +1, or -1 suffices.   *
 * RES_ERROR_NONRECVR will ensure that the program halts.         *
 * RES_ERROR_RECVR should be returned if, say, a yy or other input*
 * value is illegal. IDAS will attempt to correct and retry.      *
 *                                                                *
 ******************************************************************/

typedef int (*ResFn)(realtype tres, 
                     N_Vector yy, N_Vector yp, N_Vector resval, 
                     void *rdata);

/******************************************************************
 * Type : SensResFn                                               *
 *----------------------------------------------------------------* 
 *                                                                *
 * A SensResFn resS will return the value ires, which has         *
 * possible values RES_ERROR_RECVR = 1, RES_ERROR_NONRECVR = -1,  *
 * and SUCCESS = 0 (defined above).                               *
 *                                                                *
 ******************************************************************/

typedef int (*SensResFn)(int Ns, realtype tres, 
                         N_Vector yy, N_Vector yp, N_Vector resval,
                         N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                         void *rdataS,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/******************************************************************
 * Type : SensRes1Fn                                              *
 *----------------------------------------------------------------* 
 *                                                                *
 * A SensRes1Fn resS1 will return the value ires, which has       *
 * possible values RES_ERROR_RECVR = 1, RES_ERROR_NONRECVR = -1,  *
 * and SUCCESS = 0 (defined above).                               *
 *                                                                *
 ******************************************************************/

typedef int (*SensRes1Fn)(int Ns, realtype tres, 
                          N_Vector yy, N_Vector yp, N_Vector resval,
                          int iS,
                          N_Vector yyS, N_Vector ypS, N_Vector resvalS,
                          void *rdataS,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/******************************************************************
 * Type : QuadRhsFn                                               *
 *----------------------------------------------------------------* 
 *                                                                *
 ******************************************************************/

typedef void (*QuadRhsFn)(realtype tres, 
                          N_Vector yy, N_Vector yp,
                          N_Vector ypQ,
                          void *rdataQ);

/*----------------------------------------------------------------*
 * Enumerations for res, lsetup and lsolve return values          *
 *----------------------------------------------------------------*/

enum {SUCCESS               = 0,
      RES_ERROR_RECVR       = 1, RES_ERROR_NONRECVR      = -1,
      LSETUP_ERROR_RECVR    = 2, LSETUP_ERROR_NONRECVR   = -2,
      LSOLVE_ERROR_RECVR    = 3, LSOLVE_ERROR_NONRECVR   = -3};

/*================================================================*
 *                                                                *
 *          U S E R - C A L L A B L E   R O U T I N E S           *
 *                                                                *
 *================================================================*/

/*----------------------------------------------------------------*
 * Function : IDACreate                                           *
 *----------------------------------------------------------------*
 * IDACreate creates an internal memory block for a problem to    *
 * be solved by IDAS.                                             *
 *                                                                *
 * If successful, IDACreate returns a pointer to initialized      *
 * problem memory. This pointer should be passed to IDAMalloc.    *
 * If an initialization error occurs, IDACreate prints an error   *
 * message to standard err and returns NULL.                      *
 *                                                                *
 *----------------------------------------------------------------*/

void *IDACreate(void);

/*----------------------------------------------------------------*
 * Integrator optional input specification functions              *
 *----------------------------------------------------------------*
 * The following functions can be called to set optional inputs   *
 * to values other than the defaults given below:                 *
 *                                                                *
 * Function             |  Optional input / [ default value ]     *
 *                      |                                         * 
 * -------------------------------------------------------------- *
 *                      |                                         * 
 * IDASetRdata          | a pointer to user data that will be     *
 *                      | passed to the user's res function every *
 *                      | time res is called.                     *
 *                      | [NULL]                                  *
 *                      |                                         * 
 * IDASetErrFile        | the file pointer for an error file      *
 *                      | where all IDAS warning and error        *
 *                      | messages will be written. This parameter*
 *                      | can be stdout (standard output), stderr *
 *                      | (standard error), a file pointer        *
 *                      | (corresponding to a user error file     *
 *                      | opened for writing) returned by fopen.  *
 *                      | If not called, then all messages will   *
 *                      | be written to standard output.          *
 *                      | [NULL]                                  *
 *                      |                                         * 
 * IDASetMaxOrd         | maximum lmm order to be used by the     *
 *                      | solver.                                 *
 *                      | [5]                                     * 
 *                      |                                         * 
 * IDASetMaxNumSteps    | maximum number of internal steps to be  *
 *                      | taken by the solver in its attempt to   *
 *                      | reach tout.                             *
 *                      | [500]                                   *
 *                      |                                         * 
 * IDASetInitStep       | initial step size.                      *
 *                      | [estimated by IDAS]                     * 
 *                      |                                         * 
 * IDASetMaxStep        | maximum absolute value of step size     *
 *                      | allowed.                                *
 *                      | [infinity]                              *
 *                      |                                         * 
 * IDASetStopTime       | the independent variable value past     *
 *                      | which the solution is not to proceed.   *
 *                      | [infinity]                              *
 *                      |                                         * 
 * IDASetNlinConvCoef   | Newton convergence test constant        *
 *                      | for use during integration.             *
 *                      | [0.33]                                  *
 *                      |                                         * 
 * IDASetMaxErrTestFails| Maximum number of error test failures   *
 *                      | in attempting one step.                 *
 *                      | [10]                                    *
 *                      |                                         *
 * IDASetMaxNonlinIters | Maximum number of nonlinear solver      *
 *                      | iterations at one solution.             *
 *                      | [4]                                     *
 *                      |                                         *
 * IDASetMaxConvFails   | Maximum number of allowable conv.       *
 *                      | failures in attempting one step.        *
 *                      | [10]                                    *
 *                      |                                         *
 * IDASetSuppressAlg    | flag to indicate whether or not to      *
 *                      | suppress algebraic variables in the     *
 *                      | local error tests:                      *
 *                      | FALSE = do not suppress;                * 
 *                      | TRUE = do suppress;                     *
 *                      | [FALSE]                                 *
 *                      | NOTE: if suppressed algebraic variables *
 *                      | is selected, the nvector 'id' must be   *
 *                      | supplied for identification of those    *
 *                      | algebraic components (see IDASetID).    *
 *                      |                                         * 
 * IDASetID             | an N_Vector, which states a given       *
 *                      | element to be either algebraic or       *
 *                      | differential.                           *
 *                      | A value of 1.0 indicates a differential *
 *                      | variable while a 0.0 indicates an       *
 *                      | algebraic variable. 'id' is required    *
 *                      | if optional input SUPPRESSALG is set,   *
 *                      | or if IDACalcIC is to be called with    *
 *                      | icopt = CALC_YA_YDP_INIT.               *
 *                      |                                         *
 * IDASetConstraints    | an N_Vector defining inequality         *
 *                      | constraints for each component of the   *
 *                      | solution vector y. If a given element   *
 *                      | of this vector has values +2 or -2,     *
 *                      | then the corresponding component of y   *
 *                      | will be constrained to be > 0.0 or      *
 *                      | <0.0, respectively, while if it is +1   *
 *                      | or -1, the y component is constrained   *
 *                      | to be >= 0.0 or <= 0.0, respectively.   *
 *                      | If a component of constraints is 0.0,   *
 *                      | then no constraint is imposed on the    *
 *                      | corresponding component of y.           *
 *                      | The presence of a non-NULL constraints  *
 *                      | vector that is not 0.0 (ZERO) in all    *
 *                      | components will cause constraint        *
 *                      | checking to be performed.               *
 *                      |                                         *
 * -------------------------------------------------------------- *
 * If successful, these functions return SUCCESS. If an argument  *
 * has an illegal value, they print an error message to the       *
 * file specified by errfp and return one of the error flags      *  
 * defined below.                                                 *
 *----------------------------------------------------------------*/

int IDASetRdata(void *ida_mem, void *rdata);
int IDASetErrFile(void *ida_mem, FILE *errfp);
int IDASetMaxOrd(void *ida_mem, int maxord);
int IDASetMaxNumSteps(void *ida_mem, int mxsteps);
int IDASetInitStep(void *ida_mem, realtype hin);
int IDASetMaxStep(void *ida_mem, realtype hmax);
int IDASetStopTime(void *ida_mem, realtype tstop);
int IDASetNlinConvCoef(void *ida_mem, realtype epcon);
int IDASetMaxErrTestFails(void *ida_mem, int maxnef);
int IDASetMaxNonlinIters(void *ida_mem, int maxcor);
int IDASetMaxConvFails(void *ida_mem, int maxncf);
int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg);
int IDASetID(void *ida_mem, N_Vector id);
int IDASetConstraints(void *ida_mem, N_Vector constraints);

/* Error return values for IDASet* functions */
/* SUCCESS = 0*/
enum {IDAS_NO_MEM = -1, IDAS_ILL_INPUT = -2};

/******************************************************************
 * Function : IDAMalloc                                           *
 *----------------------------------------------------------------*
 * IDAMalloc allocates and initializes memory for a problem to    *
 * to be solved by IDAS.                                          *
 *                                                                *
 * res     is the residual function F in F(t,y,y') = 0.           *          
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
 * reltol    is a pointer to the relative tolerance scalar.       *
 *                                                                *
 * abstol    is a pointer (void) to the absolute tolerance scalar *
 *            or an N_Vector tolerance.                           *
 * (ewt)                                                          *
 *         Both reltol and abstol are used to compute the error   *
 *         vector, ewt. The error test required of a correction   *
 *         delta is that the weighted-RMS norm of delta be less   *
 *         than or equal to 1.0. Other convergence tests use the  *
 *         same norm. The weighting vector used in this norm is   *
 *         ewt. The components of ewt are defined by              *
 *         ewt[i] = 1.0/(reltol*yy[i] + abstol[i]). Here, yy is   *
 *         the current approximate solution.  See the routine     *
 *         N_VWrmsNorm for the norm used in this error test.      *
 *                                                                *
 * nvspec  is a pointer to a vector specification structure       * 
 *                                                                *
 * Note: The tolerance values may be changed in between calls to  *
 *       IDASolve for the same problem. These values refer to     *
 *       (*reltol) and either (*abstol), for a scalar absolute    *
 *       tolerance, or the components of abstol, for a vector     *
 *       absolute tolerance.                                      *
 *                                                                * 
 * If successful, IDAMalloc returns SUCCESS.                      *
 * If an initialization error occurs, IDAMalloc prints an error   *
 * message to the file specified by errfp and returns one of the  *
 * error flags defined below.                                     *
 ******************************************************************/

int IDAMalloc(void *ida_mem, ResFn res,
              realtype t0, N_Vector y0, N_Vector yp0, 
              int itol, realtype *reltol, void *abstol,
              NV_Spec nvspec);

/* Error return values for IDAMalloc */
/* SUCCESS = 0 */
enum {IDAM_NO_MEM = -1, IDAM_MEM_FAIL=-2, IDAM_ILL_INPUT = -2};

/******************************************************************
 * Function : IDAReInit                                           *
 *----------------------------------------------------------------*
 * IDAReInit re-initializes IDAS for the solution of a problem,   *
 * where a prior call to IDAMalloc has been made.                 *
 * IDAReInit performs the same input checking and initializations *
 * that IDAMalloc does.                                           *
 * But it does no memory allocation, assuming that the existing   *
 * internal memory is sufficient for the new problem.             *
 *                                                                *
 * The use of IDAReInit requires that the maximum method order,   *
 * maxord, is no larger for the new problem than for the problem  *
 * specified in the last call to IDAMalloc.  This condition is    *
 * automatically fulfilled if the default value for maxord is     *
 * specified.                                                     *
 *                                                                *
 * Following the call to IDAReInit, a call to the linear solver   *
 * specification routine is necessary if a different linear solver*
 * is chosen, but may not be otherwise.  If the same linear solver*
 * is chosen, and there are no changes in its input parameters,   *
 * then no call to that routine is needed.                        *
 *                                                                *
 * The first argument to IDAReInit is:                            *
 *                                                                *
 * ida_mem = pointer to IDA memory returned by IDACreate.         *
 *                                                                *
 * All the remaining arguments to IDAReInit have names and        *
 * meanings identical to those of IDAMalloc.                      *
 *                                                                *
 * The return value of IDAReInit is equal to SUCCESS = 0 if there *
 * were no errors; otherwise it is a negative int equal to:       *
 *   IDAREI_NO_MEM     indicating ida_mem was NULL, or            *
 *   IDAREI_NO_MALLOC  indicating that ida_mem was not allocated. *
 *   IDAREI_ILL_INPUT  indicating an input argument was illegal   *
 *                     (including an attempt to increase maxord). *
 * In case of an error return, an error message is also printed.  *
 ******************************************************************/

int IDAReInit(void *ida_mem, ResFn res,
              realtype t0, N_Vector y0, N_Vector yp0,
              int itol, realtype *reltol, void *abstol);
 
/* Error return values for IDAReInit */
/* SUCCESS = 0 */ 
enum {IDAREI_NO_MEM = -1, IDAREI_NO_MALLOC = -2, IDAREI_ILL_INPUT = -3};

/*----------------------------------------------------------------*
 * Quadrature optional input specification functions              *
 *----------------------------------------------------------------*
 * The following functions can be called to set optional inputs   *
 * to values other than the defaults given below:                 *
 *                                                                *
 * Function             |  Optional input / [ default value ]     *
 *                      |                                         * 
 * -------------------------------------------------------------- *
 *                      |                                         *
 * IDASetQuadErrCon     | are quadratures considered in the error *
 *                      | control?                                *
 *                      | [FALSE]                                 *
 *                      |                                         *
 * IDASetQuadRdata      | a pointer to user data that will be     *
 *                      | passed to the user's rhsQ function      *
 *                      | every time rhsQ is called.              *
 *                      | [NULL]                                  *
 *                      |                                         *
 * IDASetQuadTolerances | set tolerances for quadrature           *
 *                      | integration. Only needed if errconQ=TRUE*
 *                      | [no default]
 *                      |                                         *
 * -------------------------------------------------------------- *
 * If successful, these functions return SUCCESS. If an argument  *
 * has an illegal value, they print an error message to the       *
 * file specified by errfp and return one of the error flags      *  
 * defined for the IDASet* routines.                              *
 *----------------------------------------------------------------*/

int IDASetQuadErrCon(void *ida_mem, booleantype errconQ);
int IDASetQuadRdata(void *ida_mem, void *rdataQ);
int IDASetQuadTolerances(void *ida_mem, int itolQ, 
                         realtype *reltolQ, void *abstolQ);


/*----------------------------------------------------------------*
 * Function : IDAQuadMalloc                                       *
 *----------------------------------------------------------------*
 * IDAQuadMalloc allocates and initializes memory related to      *
 * quadrature integration.                                        *
 *                                                                *
 * ida_mem is a pointer to IDAS memory returned by IDACreate      *
 *                                                                *
 * rhsQ    is the user-provided integrand routine.                *
 *                                                                *
 * nvspecQ   is a pointer to a vector specification structure     *
 *           for N_Vectors containing quadrature variables.       *
 *                                                                *
 *                                                                *
 * If successful, IDAQuadMalloc returns SUCCESS.                  *
 * If an initialization error occurs, CVodeQuadMalloc prints an   *
 * error message to the file specified by errfp and returns one   *
 * of the error flags defined below.                              *
 *----------------------------------------------------------------*/

int IDAQuadMalloc(void *ida_mem, QuadRhsFn rhsQ, NV_Spec nvspecQ);

enum {QIDAM_NO_MEM = -1, QIDAM_ILL_INPUT = -2, QIDAM_MEM_FAIL = -3};

/*----------------------------------------------------------------*
 * Function : IDAQuadReInit                                       *
 *----------------------------------------------------------------*
 * IDAQuadReInit re-initializes IDAS's quadrature related         *
 * memory for a problem, assuming it has already been allocated   *
 * in prior calls to IDAMalloc and IDAQuadMalloc.                 *
 *                                                                *
 * All problem specification inputs are checked for errors.       *
 * The number of quadratures Nq is assumed to be unchanged        *
 * since the previous call to IDAQuadMalloc.                      *
 * If any error occurs during initialization, it is reported to   *
 * the file whose file pointer is errfp, and one of the error     *
 * flags defined below is returned.                               *
 *                                                                *
 *----------------------------------------------------------------*/

int IDAQuadReInit(void *ida_mem, QuadRhsFn rhsQ);

enum {QIDAREI_NO_MEM = -1, QIDAREI_NO_QUAD = -2, QIDAREI_ILL_INPUT = -3};

/*----------------------------------------------------------------*
 * Forward sensitivity optional input specification functions     *
 *----------------------------------------------------------------*
 * The following functions can be called to set optional inputs   *
 * to other values than the defaults given below:                 *
 *                                                                *
 * Function             |  Optional input / [ default value ]     *
 *                      |                                         * 
 * -------------------------------------------------------------- *
 *                      |                                         *
 * IDASetSensResFn      | sensitivity residual function.          *
 *                      | This function must compute residuals    *
 *                      | for all sensitivity equations.          *
 *                      | [IDAS difference quotient approx.]      *
 *                      |                                         *
 * IDASetSensRes1Fn     | sensitivity residual function.          *
 *                      | This function must compute the residual *
 *                      | of one sensitivity equation at a time   *
 *                      | [IDAS difference quotient approx.]      *
 *                      |                                         *
 * IDASetSensErrCon     | are sensitivity variables considered in *
 *                      | the error control?                      *
 *                      | [TRUE]                                  *
 *                      |                                         *
 * IDASetSensRho        | controls the selection of finite        *
 *                      | difference schemes used in evaluating   *
 *                      | the sensitivity right hand sides.       *
 *                      | [0.0]                                   *
 *                      |                                         *
 * IDASetSensPbar       | a pointer to scaling factors used in    *
 *                      | computing sensitivity absolute          *
 *                      | tolerances as well as by the IDAS       *
 *                      | difference quotient routines for        *
 *                      | sensitivty right hand sides. pbar[i]    *
 *                      | must give the order of magnitude of     *
 *                      | parameter p[i]. Typically, if p[i] is   *
 *                      | nonzero, pbar[i]=p[i].                  *
 *                      | [NULL]                                  *
 *                      |                                         *
 * IDASetSensTolerances | a pointer to the sensitivity relative   *
 *                      | tolerance scalar and one for the        *
 *                      | absolute tolerance                      *
 *                      | [estimated by IDAS]                     *
 *                      |                                         *
 * IDASetSensRdata      | a pointer to user data that will be     *
 *                      | passed to the user's resS function      *
 *                      | every time resS is called.              *
 *                      | [NULL]                                  *
 *                          |                                     *
 * IDASetSensMaxNonlinIters | Maximum number of nonlinear solver  *
 *                          | iterations for sensitivity systems  *
 *                          | (staggered)                         *
 *                          | [4]                                 *
 *                          |                                     *
 * -------------------------------------------------------------- *
 * If successful, these functions return SUCCESS. If an argument  *
 * has an illegal value, they print an error message to the       *
 * file specified by errfp and return one of the error flags      *  
 * defined for the IDASet* routines.                              *
 *----------------------------------------------------------------*/

int IDASetSensResFn(void *ida_mem, SensResFn resS);
int IDASetSensRes1Fn(void *ida_mem, SensRes1Fn resS);
int IDASetSensErrCon(void *ida_mem, booleantype errconS);
int IDASetSensRho(void *ida_mem, realtype rho);
int IDASetSensPbar(void *ida_mem, realtype *pbar);
int IDASetSensTolerances(void *ida_mem, int itolS, 
                         realtype *reltolS, void *abstolS);
int IDASetSensRdata(void *ida_mem, void *rdataS);
int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS);

/*----------------------------------------------------------------*
 * Function : IDASensMalloc                                       *
 *----------------------------------------------------------------*
 * IDASensMalloc allocates and initializes memory related to      *
 * sensitivity computations.                                      *
 *                                                                *
 * ida_mem is a pointer to IDAS memory returned by IDAMalloc      *
 *                                                                *
 * Ns        is the number of sensitivities to be computed.       *
 *                                                                *
 * ism       is the type of corrector used in sensitivity         *
 *           analysis. The legal values are: SIMULTANEOUS,        *
 *           STAGGERED, and STAGGERED1 (see previous description) *
 *                                                                *
 * p         is a pointer to problem parameters with respect to   *
 *           which sensitivities may be computed (see description *
 *           of plist below). If the right hand sides of the      *
 *           sensitivity equations are to be evaluated by the     *
 *           difference quotient routines provided with IDAS,     *
 *           then p must also be a field in the user data         *
 *           structure pointed to by f_data.                      *
 *                                                                *
 * plist     is a pointer to a list of parameters with respect to *
 *           which sensitivities are to be computed.              *
 *           If plist[j]=i, then sensitivities with respect to    *
 *           the i-th parameter (i.e. p[i-1]) will be computed.   *
 *           A negative plist entry also indicates that the       *
 *           corresponding parameter affects only the initial     *
 *           conditions of the ODE and not its right hand side.   *
 *                                                                *
 * yS0       is the array of initial condition vectors for        *
 *           sensitivity variables.                               * 
 *                                                                *
 * ypS0      is the array of initial condition vectors for        *
 *           sensitivity derivatives.                             * 
 *                                                                *
 * If successful, IDASensMalloc returns SUCCESS. If an            *
 * initialization error occurs, IDASensMalloc prints an error     *
 * message to the file specified by errfp and returns one of      *
 * the error flags defined below.                                 *
 *                                                                *
 *----------------------------------------------------------------*/

int IDASensMalloc(void *ida_mem, int Ns, int ism, 
                  realtype *p, int *plist, 
                  N_Vector *yS0, N_Vector *ypS0);
    
enum {SIDAM_NO_MEM = -1, SIDAM_ILL_INPUT = -2, SIDAM_MEM_FAIL = -3};

/*----------------------------------------------------------------*
 * Function : IDASensReInit                                       *
 *----------------------------------------------------------------*
 * IDASensReInit re-initializes CVODES's sensitivity related      *
 * memory for a problem, assuming it has already been allocated   *
 * in prior calls to IDAMalloc and CvodeSensMalloc.               *
 *                                                                *
 * All problem specification inputs are checked for errors.       *
 * The number of sensitivities Ns is assumed to be unchanged      *
 * since the previous call to IDASensMalloc.                      *
 * If any error occurs during initialization, it is reported to   *
 * the file whose file pointer is errfp.                          *
 *                                                                *
 * IDASensReInit potentially does some minimal memory allocation  *
 * (for the sensitivity absolute tolerance and for arrays of      *
 * counters used by the STAGGERED1 method).                       *
 *                                                                *
 * The return value is equal to SUCCESS = 0 if there were no      *
 * errors; otherwise it is a negative int equal to:               *
 *   SIDAREI_NO_MEM    indicating cvode_mem was NULL, or          *
 *   SIDAREI_NO_SENSI  indicating there was not a prior call to   *
 *                     IDASensMalloc.                             *
 *   SIDAREI_ILL_INPUT indicating an input argument was illegal   *
 *   SIDAREI_MEM_FAIL  indicating a memory request failed.        *
 * In case of an error return, an error message is also printed.  *
 *----------------------------------------------------------------*/

int IDASensReInit(void *ida_mem, int ism,
                  realtype *p, int *plist, 
                  N_Vector *yS0, N_Vector *ypS0);
  
enum {SIDAREI_NO_MEM    = -1, SIDAREI_NO_SENSI = -2, 
      SIDAREI_ILL_INPUT = -3, SIDAREI_MEM_FAIL = -4};

/*----------------------------------------------------------------*
 * Initial Conditions optional input specification functions      *
 *----------------------------------------------------------------*
 * The following functions can be called to set optional inputs   *
 * to control the initial conditions calculations.                *
 *                                                                *
 * Function               |  Optional input / [ default value ]   *
 *                        |                                       * 
 * -------------------------------------------------------------- *
 *                        |                                       * 
 * IDASetNlinConvCoefIC   | positive scalar factor in the Newton  *
 *                        | convergence test.  This test uses a   *
 *                        | weighted RMS norm (with weights       *
 *                        | defined by the tolerances, as in      *
 *                        | IDASolve).  For new initial value     *
 *                        | vectors y and y' to be accepted, the  *
 *                        | norm of J-inverse F(t0,y,y') is       *
 *                        | required to be less than epiccon,     *
 *                        | where J is the system Jacobian.       *
 *                        | [0.01 * 0.33]                         * 
 *                        |                                       * 
 * IDASetMaxNumStepsIC    | maximum number of values of h allowed *
 *                        | when icopt = CALC_YA_YDP_INIT, where  *
 *                        | h appears in the system Jacobian,     *
 *                        | J = dF/dy + (1/h)dF/dy'.              *
 *                        | [5]                                   *
 *                        |                                       * 
 * IDASetMaxNumJacsIC     | maximum number of values of the       *
 *                        | approximate Jacobian or preconditioner*
 *                        | allowed, when the Newton iterations   *
 *                        | appear to be slowly converging.       *
 *                        | [4]                                   * 
 *                        |                                       * 
 * IDASetMaxNumItersIC    | maximum number of Newton iterations   *
 *                        | allowed in any one attempt to solve   *
 *                        | the IC problem.                       *
 *                        | [10]                                  *
 *                        |                                       * 
 * IDASetLineSearchOffIC  | a boolean flag to turn off the        *
 *                        | linesearch algorithm.                 *
 *                        | [FALSE]                               *
 *                        |                                       * 
 * IDASetStepToleranceIC  | positive lower bound on the norm of   *
 *                        | a Newton step.                        *
 *                        | [(unit roundoff)^(2/3)                *
 ******************************************************************/

int IDASetNlinConvCoefIC(void *ida_mem, realtype epiccon);
int IDASetMaxNumStepsIC(void *ida_mem, int maxnh);
int IDASetMaxNumJacsIC(void *ida_mem, int maxnj);
int IDASetMaxNumItersIC(void *ida_mem, int maxnit);
int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff);
int IDASetStepToleranceIC(void *ida_mem, realtype steptol);

/******************************************************************
 * Function : IDACalcIC                                           *
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
 * IDAMalloc or IDAReInit for the given DAE problem, and by a     *
 * successful call to the linear system solver specification      *
 * routine.                                                       *
 * In addition, IDACalcIC assumes that the vectors y0, yp0, and   *
 * (if relevant) id and constraints that were passed to IDAMalloc *
 * (or IDAReInit) remain unaltered since that call.               *
 *                                                                *
 * The call to IDACalcIC should precede the call(s) to IDASolve   *
 * for the given problem.                                         *  
 *                                                                *
 * The arguments to IDACalcIC are as follows.                     *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by IDACreate.    *
 *                                                                *
 * icopt  is the option of IDACalcIC to be used.                  *
 *        icopt = CALC_YA_YDP_INIT   directs IDACalcIC to compute *
 *                the algebraic components of y and differential  *
 *                components of y', given the differential        *
 *                components of y.  This option requires that the *
 *                N_Vector id was input to IDAMalloc or IDAReInit,*
 *                specifying the differential and algebraic       *
 *                components.                                     *
 *        icopt = CALC_Y_INIT   directs IDACalcIC to compute all  *
 *                components of y, given y'.  id is not required. *
 *                                                                *
 * tout1  is the first value of t at which a soluton will be      *
 *        requested (from IDASolve).  (This is needed here to     *
 *        determine the direction of integration and rough scale  *
 *        in the independent variable t.                          *
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
 ******************************************************************/

int IDACalcIC (void *ida_mem, int icopt, realtype tout1); 

/* IDACalcIC return values */

/* The following three values are IDASolve return values, repeated
   here for convenience. 
       SETUP_FAILURE=-7,  SOLVE_FAILURE=-8,  RES_NONRECOV_ERR=-11  */

enum { IC_IDA_NO_MEM =     -20,   IC_ILL_INPUT =       -21,
       IC_LINIT_FAIL =     -22,   IC_BAD_EWT =         -23,
       IC_FIRST_RES_FAIL = -24,   IC_NO_RECOVERY =     -25,
       IC_FAILED_CONSTR =  -26,   IC_FAILED_LINESRCH = -27,
       IC_CONV_FAILURE =   -28,   IC_NO_MALLOC =       -28  };

/******************************************************************
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
 * responsible for allocating the memory for this value.          *
 *                                                                *
 * IDA_mem is the pointer (void) to IDAS memory returned by       *
 *         IDAMalloc.                                             *
 *                                                                *
 * tout    is the next independent variable value at which a      *
 *         computed solution is desired.                          *
 *                                                                *
 * *tret   is the actual independent variable value corresponding *
 *         to the solution vector yret.                           *
 *                                                                *
 * yret    is the computed solution vector.  With no errors,      *
 *         yret = y(tret).                                        *
 *                                                                *
 * ypret   is the derivative of the computed solution at t = tret.*
 *                                                                *
 * Note: yret and ypret may be the same N_Vectors as y0 and yp0   *
 * in the call to IDAMalloc or IDAReInit.                         *
 *                                                                *
 * itask   is NORMAL, NORMAL_TSTOP, ONE_STEP, or ONE_STEP_TSTOP.  *
 *         These modes are described above.                       *
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

int IDASolve(void *ida_mem, realtype tout, realtype *tret,
             N_Vector yret, N_Vector ypret, int itask);

/* IDASolve return values */

enum { NORMAL_RETURN=0,     INTERMEDIATE_RETURN=1, TSTOP_RETURN=2,
       IDA_NO_MEM=-1,       NO_MALLOC=-2,          ILL_INPUT=-3,
       TOO_MUCH_WORK=-4,    TOO_MUCH_ACC=-5,       ERR_FAILURE=-6,
       CONV_FAILURE=-7,     SETUP_FAILURE=-8,      SOLVE_FAILURE=-9,
       CONSTR_FAILURE=-10,  REP_RES_REC_ERR=-11,   RES_NONRECOV_ERR=-12 };

/*----------------------------------------------------------------*
 * Function: IDAGetSolution                                       *
 *----------------------------------------------------------------*
 *                                                                *
 * This routine evaluates y(t) and y'(t) as the value and         *
 * derivative of the interpolating polynomial at the independent  *
 * variable t, and stores the results in the vectors yret and     *
 * ypret.  It uses the current independent variable value, tn,    *
 * and the method order last used, kused. This function is        *
 * called by IDASolve with t = tout, t = tn, or t = tstop.        *
 *                                                                *
 * If kused = 0 (no step has been taken), or if t = tn, then the  *
 * order used here is taken to be 1, giving yret = phi[0],        *
 * ypret = phi[1]/psi[0].                                         *
 *                                                                *
 * The return values are:                                         *
 *  OKAY  if t is legal, or                                       *
 *  BAD_T if t is not within the interval of the last step taken. *
 *                                                                *
 ******************************************************************/

int IDAGetSolution(void *ida_mem, realtype t, 
                   N_Vector yret, N_Vector ypret);

/*----------------------------------------------------------------*
 * Integrator optional output extraction functions                *
 *----------------------------------------------------------------*
 * The following functions can be called to get optional outputs  *
 * and statistics related to the main integrator.                 *
 * -------------------------------------------------------------- *
 *                                                                *
 * IDAGetIntWorkSpace returns the IDAS integer workspace size     *
 * IDAGetRealWorkSpace returns the IDAS real workspace size       *
 * IDAGetNumSteps returns the cumulative number of internal       *
 *       steps taken by the solver                                *
 * IDAGetNumRhsEvals returns the number of calls to the user's    *
 *       res function                                             *
 * IDAGetNumLinSolvSetups returns the number of calls made to     *
 *       the linear solver's setup routine                        *
 * IDAGetNumErrTestFails returns the number of local error test   *
 *       failures that have occured                               *
 * IDAGetNumBacktrackOps returns the number of backtrack          *
 *       operations done in the linesearch algorithm in IDACalcIC *
 * IDAGetLastOrder returns the order used during the last         *
 *       internal step                                            *
 * IDAGetCurrentOrder returns the order to be used on the next    *
 *       internal step                                            *
 * IDAGetActualInitStep returns the actual initial step size      *
 *       used by IDAS                                             *
 * IDAGetLastStep returns the step size for the last internal     *
 *       step (if from IDASolve), or the last value of the        *
 *       artificial step size h (if from IDACalcIC)               *
 * IDAGetCurrentStep returns the step size to be attempted on the *
 *       next internal step                                       *
 * IDAGetCurrentTime returns the current internal time reached    *
 *       by the solver                                            *
 * IDAGetTolScaleFactor returns a suggested factor by which the   *
 *       user's tolerances should be scaled when too much         *
 *       accuracy has been requested for some internal step       *
 * IDAGetErrWeights returns the state error weight vector.        *
 *       The user need not allocate space for ewt.                *
 * IDAGetEstLocalErrors returns the vector of estimated local     *
 *       errors. The user need not allocate space for ele.        *
 *----------------------------------------------------------------*/

int IDAGetIntWorkSpace(void *ida_mem, long int *leniw);
int IDAGetRealWorkSpace(void *ida_mem, long int *lenrw);
int IDAGetNumSteps(void *ida_mem, int *nsteps);
int IDAGetNumResEvals(void *ida_mem, int *nrevals);
int IDAGetNumLinSolvSetups(void *ida_mem, int *nlinsetups);
int IDAGetNumErrTestFails(void *ida_mem, int *netfails);
int IDAGetNumBacktrackOps(void *ida_mem, int *nbacktr);
int IDAGetLastOrder(void *ida_mem, int *klast);
int IDAGetCurrentOrder(void *ida_mem, int *kcur);
int IDAGetActualInitStep(void *ida_mem, realtype *hinused);
int IDAGetLastStep(void *ida_mem, realtype *hlast);
int IDAGetCurrentStep(void *ida_mem, realtype *hcur);
int IDAGetCurrentTime(void *ida_mem, realtype *tcur);
int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact);
int IDAGetErrWeights(void *ida_mem, N_Vector *eweight);

/*----------------------------------------------------------------*
 * As a convenience, the following two functions provide the      *
 * optional outputs in groups.                                    *
 *----------------------------------------------------------------*/

int IDAGetWorkSpace(void *ida_mem, long int *leniw, long int *lenrw);
int IDAGetIntegratorStats(void *ida_mem, int *nsteps, int *nrevals, 
                          int *nlinsetups, int *netfails,
                          int *qlast, int *qcur, realtype *hlast, 
                          realtype *hcur, realtype *tcur);

/*----------------------------------------------------------------*
 * Nonlinear solver optional output extraction functions          *
 *----------------------------------------------------------------*
 * The following functions can be called to get optional outputs  *
 * and statistics related to the nonlinear solver.                *
 * -------------------------------------------------------------- *
 *                                                                *
 * IDAGetNumNonlinSolvIters returns the number of nonlinear       *
 *       solver iterations performed.                             *
 * IDAGetNumNonlinSolvConvFails returns the number of nonlinear   *
 *       convergence failures.                                    *
 *----------------------------------------------------------------*/

int IDAGetNumNonlinSolvIters(void *ida_mem, int *nniters);
int IDAGetNumNonlinSolvConvFails(void *ida_mem, int *nncfails);

/*----------------------------------------------------------------*
 * As a convenience, the following function provides the          *
 * optional outputs in a group.                                   *
 *----------------------------------------------------------------*/

int IDAGetNonlinSolvStats(void *ida_mem, int *nniters, int *nncfails);

/*----------------------------------------------------------------*
 * Quadrature integration solution extraction routines            *
 *----------------------------------------------------------------*
 * The following functions can be called to obtain the quadrature *
 * variables after a successful integration step.                 *
 *                                                                *
 * Return values are similar to those of IDAGetSolution.          *
 * Additionally, IDAGetQuad can return IDAG_NO_QUAD if quadratures*
 * were not computed.                                             *
 *----------------------------------------------------------------*/

int IDAGetQuad(void *ida_mem, realtype t, N_Vector yretQ);

/*----------------------------------------------------------------*
 * Quadrature integration optional output extraction routines     *
 *----------------------------------------------------------------*
 * The following functions can be called to get optional outputs  *
 * and statistics related to the integration of quadratures.      *
 * -------------------------------------------------------------- *
 *                                                                * 
 * IDAGetNumQuadRhsEvals returns the number of calls to the       *
 *      user function rhsQ defining the right hand side of the    *
 *      quadrature variables.                                     *
 * IDAGetNumQuadErrTestFails returns the number of local error    *
 *      test failures for quadrature variables.                   *
 * IDAGetQuadErrWeights returns the vector of error weights for   *
 *      the quadrature variables. The user need not allocate      *
 *      space for ewtQ.                                           *
 *----------------------------------------------------------------*/

int IDAGetNumQuadRhsEvals(void *ida_mem, int *nrhsQevals);
int IDAGetNumQuadErrTestFails(void *ida_mem, int *nQetfails);
int IDAGetQuadErrWeights(void *ida_mem, N_Vector *eQweight);

/*----------------------------------------------------------------*
 * As a convenience, the following function provides the          *
 * optional outputs in a group.                                   *
 *----------------------------------------------------------------*/

int IDAGetQuadStats(void *ida_mem, int *nrhsQevals, int *nQetfails);

/*----------------------------------------------------------------*
 * Forward sensitivity solution extraction routines               *
 *----------------------------------------------------------------*
 * IDAGetSens returns sensitivities of the y function at          *
 * the time t. The argument yS must be a pointer to N_Vector      *
 * and must be allocated by the user to hold at least Ns vectors. *
 *                                                                *
 * Return values are similar to those of IDAGetSolution.          *
 * Additionally, IDAGetSens can return IDAG_NO_SENSI if           *
 * sensitivities were not computed and BAD_IS if                  *
 * is < 0 or is >= Ns.                                            *
 *----------------------------------------------------------------*/

int IDAGetSens(void *ida_mem, realtype t, 
               N_Vector *yretS, N_Vector *ypretS);

/*----------------------------------------------------------------*
 * Forward sensitivity optional output extraction routines        *
 *----------------------------------------------------------------*
 * The following functions can be called to get optional outputs  *
 * and statistics related to the integration of sensitivities.    *
 * -------------------------------------------------------------- *
 *                                                                * 
 * IDAGetNumSensResEvals returns the number of calls to the       *
 *     sensitivity residual routine.                              *
 * IDAGetNumResEvalsSens returns the number of calls to the       *
 *     user res routine due to finite difference evaluations of   *
 *     the sensitivity equations.                                 *
 * IDAGetNumSensErrTestFails returns the number of local error    *
 *     test failures for sensitivity variables.                   *
 * IDAGetNumSensLinSolvSetups returns the number of calls made    *
 *     to the linear solver's setup routine due to sensitivity    *
 *     computations.                                              *
 * IDAGetSensErrWeights returns the sensitivity error weight      *
 *     vectors. The user need not allocate space for ewtS.        *
 *----------------------------------------------------------------*/

int IDAGetNumSensRhsEvals(void *ida_mem, int *nresSevals);
int IDAGetNumRhsEvalsSens(void *ida_mem, int *nresevalsS);
int IDAGetNumSensErrTestFails(void *ida_mem, int *nSetfails);
int IDAGetNumSensLinSolvSetups(void *ida_mem, int *nlinsetupsS);
int IDAGetSensErrWeights(void *ida_mem, N_Vector_S *eSweight);

/*----------------------------------------------------------------*
 * As a convenience, the following function provides the          *
 * optional outputs in a group.                                   *
 *----------------------------------------------------------------*/

int IDAGetSensStats(void *ida_mem, int *nresSevals, int *nresevalsS, 
                      int *nSetfails, int *nlinsetupsS);

/*----------------------------------------------------------------*
 *                                                                *
 * Sensitivity nonlinear solver optional output extraction        *
 *----------------------------------------------------------------*
 *                                                                *
 * The following functions can be called to get optional outputs  *
 * and statistics related to the sensitivity nonlinear solver.    *
 * -------------------------------------------------------------- *
 *                                                                * 
 * IDAGetNumSensNonlinSolvIters returns the total number of       *
 *         nonlinear iterations for sensitivity variables.        *
 * IDAGetNumSensNonlinSolvConvFails returns the total number of   *
 *         nonlinear convergence failures for sensitivity         *
 *         variables                                              *
 * IDAGetNumStgrSensNonlinSolvIters returns a vector of Ns        *
 *         nonlinear iteration counters for sensitivity variables *
 *         in the STAGGERED1 method.                              *
 * IDAGetNumStgrSensNonlinSolvConvFails returns a vector of Ns    *
 *         nonlinear solver convergence failure counters for      *
 *         sensitivity variables in the STAGGERED1 method.        *
 *----------------------------------------------------------------*/

int IDAGetNumSensNonlinSolvIters(void *ida_mem, int *nSniters);
int IDAGetNumSensNonlinSolvConvFails(void *ida_mem, int *nSncfails);
int IDAGetNumStgrSensNonlinSolvIters(void *ida_mem, int *nSTGR1niters);
int IDAGetNumStgrSensNonlinSolvConvFails(void *ida_mem, int *nSTGR1ncfails);

/*----------------------------------------------------------------*
 * As a convenience, the following two functions provide the      *
 * optional outputs in groups.                                    *
 *----------------------------------------------------------------*/

int IDAGetSensNonlinSolvStats(void *ida_mem, int *nSniters, int *nSncfails);
int IDAGetStgrSensNonlinSolvStats(void *ida_mem, int *nSTGR1niters, 
                                    int *nSTGR1ncfails);

/* IDAGet* return values */
enum { OKAY = 0, IDAG_NO_MEM = -1, BAD_T = -2, 
       IDAG_NO_QUAD = -3, IDAG_NO_SENS = -4 };

/*----------------------------------------------------------------*
 * Function : IDAFree                                             *
 *----------------------------------------------------------------*
 * IDAFree frees the problem memory IDA_mem allocated by          *
 * IDAMalloc.  Its only argument is the pointer idamem            *
 * returned by IDAMalloc.                                         *
 *                                                                *
 *----------------------------------------------------------------*/

void IDAFree(void *ida_mem);

/*----------------------------------------------------------------*
 *                                                                *
 * Function : IDAQuadFree                                         *
 *----------------------------------------------------------------*
 * IDAQuadFree frees the problem memory in ida_mem allocated      *
 * for quadrature integration. Its only argument is the pointer   *
 * ida_mem returned by IDACreate.                                 *
 *                                                                *
 *----------------------------------------------------------------*/

void IDAQuadFree(void *ida_mem);

/*----------------------------------------------------------------*
 *                                                                *
 * Function : IDASensFree                                         *
 *----------------------------------------------------------------*
 * IDASensFree frees the problem memory in ida_mem allocated      *
 * for sensitivity analysis. Its only argument is the pointer     *
 * ida_mem returned by IDACreate.                                 *
 *                                                                *
 *----------------------------------------------------------------*/

void IDASensFree(void *ida_mem);


/*================================================================*
 *                                                                *
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K   *
 *                                                                *
 *================================================================*/

/* Basic IDAS constants */

#define MXORDP1 6 /* max. number of vectors kept in the phi array */

/******************************************************************
 * Types : struct IDAMemRec, IDAMem                               *
 *----------------------------------------------------------------*
 * The type IDAMem is type pointer to struct IDAMemRec. This      *
 * structure contains fields to keep track of problem state.      *
 *                                                                *
 ******************************************************************/

typedef struct IDAMemRec {

  realtype ida_uround;    /* machine unit roundoff */

  /*--------------------------
    Problem Specification Data 
    --------------------------*/

  ResFn          ida_res;            /* F(t,y(t),y'(t))=0; the function F  */
  void          *ida_rdata;          /* user pointer passed to res         */
  int            ida_itol;           /* itol = SS or SV                    */
  realtype      *ida_reltol;         /* ptr to relative tolerance          */
  void          *ida_abstol;         /* ptr to absolute tolerance          */  
  booleantype    ida_setupNonNull;   /* Does setup do something?           */
  booleantype    ida_constraintsSet; /* constraints vector present: 
                                        do constraints calc                */
  booleantype    ida_suppressalg;    /* true means suppress algebraic vars
                                        in local error tests               */

  /*-----------------------
    Quadrature Related Data 
    -----------------------*/

  booleantype    ida_quad;
  QuadRhsFn      ida_rhsQ;
  int            ida_itolQ;
  realtype      *ida_reltolQ;
  void          *ida_abstolQ;
  booleantype    ida_errconQ;
  void          *ida_rdataQ;

  /*------------------------
    Sensitivity Related Data
    ------------------------*/

  booleantype    ida_sensi;
  int            ida_Ns;
  SensResFn      ida_resS;
  SensRes1Fn     ida_resS1;
  booleantype    ida_resSDQ;
  int            ida_iresS;
  int            ida_ism;
  realtype      *ida_p;
  realtype      *ida_pbar;
  int           *ida_plist;
  int            ida_itolS;
  realtype      *ida_reltolS;
  void          *ida_abstolS;
  realtype       ida_rhomax;
  booleantype    ida_errconS;
  void          *ida_rdataS;

  /*-----------------------------------------------
    Divided differences array and associated arrays
    -----------------------------------------------*/

  int ida_maxcol;              /* Actual number of phi arrays allocated          */
  N_Vector ida_phi[MXORDP1];   /* phi = (maxord+1) arrays of divided differences */

  realtype ida_psi[MXORDP1];   /* differences in t (sums of recent step sizes)   */
  realtype ida_alpha[MXORDP1]; /* ratios of current stepsize to psi values       */
  realtype ida_beta[MXORDP1];  /* ratios of current to previous product of psi's */
  realtype ida_sigma[MXORDP1]; /* product successive alpha values and factorial  */
  realtype ida_gamma[MXORDP1]; /* sum of reciprocals of psi values               */

  /*-------------------------
    N_Vectors for integration
    -------------------------*/

  N_Vector ida_ewt;         /* error weight vector                           */
  N_Vector ida_y0;          /* initial y vector (user-supplied)              */
  N_Vector ida_yp0;         /* initial y' vector (user-supplied)             */
  N_Vector ida_yy;          /* work space for y vector (= user's yret)       */
  N_Vector ida_yp;          /* work space for y' vector (= user's ypret)     */
  N_Vector ida_delta;       /* residual vector                               */
  N_Vector ida_id;          /* bit vector for diff./algebraic components     */
  N_Vector ida_constraints; /* vector of inequality constraint options       */
  N_Vector ida_savres;      /* saved residual vector (= tempv1)              */
  N_Vector ida_ee;          /* accumulated corrections to y                  */
  N_Vector ida_mm;          /* mask vector in constraints tests (= tempv2)   */
  N_Vector ida_tempv1;      /* work space vector                             */
  N_Vector ida_tempv2;      /* work space vector                             */
  N_Vector ida_ynew;        /* work vector for y in IDACalcIC (= tempv2)     */
  N_Vector ida_ypnew;       /* work vector for yp in IDACalcIC (= ee)        */
  N_Vector ida_delnew;      /* work vector for delta in IDACalcIC (= phi[2]) */
  N_Vector ida_dtemp;       /* work vector in IDACalcIC (= phi[3])           */

  /*----------------------------
    Quadrature Related N_Vectors 
    ----------------------------*/

  N_Vector ida_phiQ[MXORDP1];
  N_Vector ida_yyQ;
  N_Vector ida_ypQ;
  N_Vector ida_ewtQ;
  N_Vector ida_eeQ;

  /*---------------------------
    Sensitivity Related Vectors 
    ---------------------------*/

  N_Vector *ida_phiS[MXORDP1];
  N_Vector *ida_ewtS;

  N_Vector *ida_yS0;        /* initial yS vector (user-supplied)            */
  N_Vector *ida_ypS0;       /* initial yS' vector (user-supplied)           */

  N_Vector *ida_eeS;        /* cumulative sensitivity corrections           */

  N_Vector *ida_yyS;        /* allocated and used for:                      */
  N_Vector *ida_ypS;        /*                 ism = SIMULTANEOUS           */
  N_Vector *ida_deltaS;     /*                 ism = STAGGERED              */

  N_Vector ida_yyS1;        /* allocated and used for:                      */
  N_Vector ida_ypS1;        /*                                              */
  N_Vector ida_deltaS1;     /*                 ism = STAGGERED1             */

  N_Vector ida_tmpS1;       /* work space vectors: tmpS1 = tempv1           */
  N_Vector ida_tmpS2;       /*                     tmpS2 = tempv2           */
  N_Vector ida_tmpS3;       /*                     tmpS3 = allocated        */    
  

  /*----------------------------
    Scalars for use by IDACalcIC
    ----------------------------*/

  int ida_icopt;            /* IC calculation user option                    */
  booleantype ida_lsoff;    /* IC calculation linesearch turnoff option      */
  int ida_maxnh;            /* max. number of h tries in IC calculation      */
  int ida_maxnj;            /* max. number of J tries in IC calculation      */
  int ida_maxnit;           /* max. number of Netwon iterations in IC calc.  */
  int ida_nbacktr;          /* number of IC linesearch backtrack operations  */
  int ida_sysindex;         /* computed system index (0 or 1)                */
  realtype ida_epiccon;     /* IC nonlinear convergence test constant        */
  realtype ida_steptol;     /* minimum Newton step size in IC calculation    */
  realtype ida_tscale;      /* time scale factor = abs(tout1 - t0)           */

  /*-------------------------------------------------
    Does IDASensMalloc allocate additional space?
  -------------------------------------------------*/  

  booleantype ida_tolSset;      /* tolerances set by IDAS?                   */
  booleantype ida_abstolSalloc; /* abstolS allocated by IDAS?                */
  booleantype ida_stgr1alloc;   /* ncfS1,ncfnS1,and nniS1 allocated by IDAS? */


  /*-----------------
    Tstop information
  -------------------*/
  booleantype ida_tstopset;
  realtype ida_tstop;

  /*---------
    Step Data
    ---------*/

  int ida_kk;        /* current BDF method order                          */
  int ida_kused;     /* method order used on last successful step         */
  int ida_knew;      /* order for next step from order decrease decision  */
  int ida_phase;     /* flag to trigger step doubling in first few steps  */
  int ida_ns;        /* counts steps at fixed stepsize and order          */

  realtype ida_hin;      /* initial step                                      */
  realtype ida_h0u;      /* actual initial stepsize                           */
  realtype ida_hh;       /* current step size h                               */
  realtype ida_hused;    /* step size used on last successful step            */
  realtype ida_rr;       /* rr = hnext / hused                                */
  realtype ida_tn;       /* current internal value of t                       */
  realtype ida_tretp;    /* value of tret previously returned by IDASolve     */
  realtype ida_cj;       /* current value of scalar (-alphas/hh) in Jacobian  */
  realtype ida_cjlast;   /* cj value saved from last successful step          */
  realtype ida_cjold;    /* cj value saved from last call to lsetup           */
  realtype ida_cjratio;  /* ratio of cj values: cj/cjold                      */
  realtype ida_ss;       /* scalar used in Newton iteration convergence test  */
  realtype ida_epsNewt;  /* test constant in Newton convergence test          */
  realtype ida_epcon;    /* Newton convergence test constant                  */
  realtype ida_toldel;   /* tolerance in direct test on Newton corrections    */
  
  realtype ida_ssS;      /* scalar ss for sensitivity variables (STAGGERED)   */
  realtype *ida_ssS1;    /* scalars ss for sensitivity variables (STAGGERED1) */

  /*------
    Limits
    ------*/

  int ida_maxncf;        /* max numer of convergence failures                 */
  int ida_maxcor;        /* max number of Newton corrections                  */
  int ida_maxnef;        /* max number of error test failures                 */

  int ida_maxord;        /* max value of method order k:                      */
  int ida_mxstep;        /* max number of internal steps for one user call    */
  realtype ida_hmax_inv; /* inverse of max. step size hmax (default = 0.0)    */

  int ida_maxcorS;       /* max number of Newton corrections for sensitivity
                            systems (staggered method)                        */

  /*--------
    Counters
    --------*/

  int ida_nst;           /* number of internal steps taken                    */

  int ida_nre;           /* number of function (res) calls                    */
  int ida_nrQe;
  int ida_nrSe;
  int ida_nreS;

  int ida_ncfn;          /* number of corrector convergence failures          */
  int ida_ncfnS;
  int *ida_ncfnS1;

  int ida_nni;           /* number of Newton iterations performed             */
  int ida_nniS;
  int *ida_nniS1;

  int ida_netf;          /* number of error test failures                     */
  int ida_netfQ;
  int ida_netfS;

  int ida_nsetups;       /* number of lsetup calls                            */
  int ida_nsetupsS;
  
  /*---------------------------
    Space requirements for IDAS
    ---------------------------*/

  long int ida_lrw1;     /* no. of realtype words in 1 N_Vector               */
  long int ida_liw1;     /* no. of integertype words in 1 N_Vector            */
  long int ida_lrw1Q;
  long int ida_liw1Q;
  long int ida_lrw;      /* number of realtype words in IDAS work vectors     */
  long int ida_liw;      /* no. of integertype words in IDAS work vectors     */

  realtype ida_tolsf;    /* tolerance scale factor (saved value)              */

  /*-------------------------------------------------------
    Flags to verify correct calling sequence
    Turned ON by IDAMalloc, IDASensMalloc, and IDAQuadMalloc 
    and read by IDAMalloc, IDASensReInit, and IDAQuadReInit
    --------------------------------------------------------*/

  booleantype ida_SetupDone;     /* set to FALSE by IDAMalloc and IDAReInit */
                                 /* set to TRUE by IDACalcIC or IDASolve    */

  booleantype ida_MallocDone;    /* set to FALSE by IDACreate               */
                                 /* set to TRUE by IDAMAlloc                */
                                 /* tested by IDAReInit and IDASolve        */

  booleantype ida_sensMallocDone;

  booleantype ida_quadMallocDone;
  
  /*-------------------------------------
    IDAS error messages are sent to errfp
    -------------------------------------*/

  FILE *ida_errfp;

  /*------------------
    Linear Solver Data
    ------------------*/

  /* Linear Solver functions to be called */

  int (*ida_linit)(struct IDAMemRec *idamem);

  int (*ida_lsetup)(struct IDAMemRec *idamem, N_Vector yyp, 
                    N_Vector ypp, N_Vector resp, 
                    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3); 

  int (*ida_lsolve)(struct IDAMemRec *idamem, N_Vector b, N_Vector ycur,
                    N_Vector ypcur, N_Vector rescur);

  int (*ida_lsolveS)(struct IDAMemRec *idamem, N_Vector b, N_Vector ycur,
                     N_Vector ypcur, N_Vector rescur, int is);

  int (*ida_lperf)(struct IDAMemRec *idamem, int perftask);

  int (*ida_lfree)(struct IDAMemRec *idamem);

  /* Linear Solver specific memory */

  void *ida_lmem;           

  /* Flag to indicate successful ida_linit call */

  booleantype ida_linitOK;

  /*-------------------------------------------------
    Pointer to the vector specification structure for
    state N_Vectors
    -------------------------------------------------*/

  NV_Spec ida_nvspec;

  /*-------------------------------------------------
    Pointer to the vector specification structure for 
    quadrature N_Vectors
    -------------------------------------------------*/

  NV_Spec ida_nvspecQ;

} *IDAMem;


/******************************************************************
 *                                                                *
 * Communication between user and an IDAS Linear Solver           *
 *----------------------------------------------------------------*
 * Return values of the linear solver specification routine.      *
 * The values of these are given in the enum statement below.     *
 * SUCCESS      : The routine was successful.                     *
 *                                                                *
 * LIN_NO_MEM   : IDAS memory = NULL.                             *
 *                                                                *
 * LMEM_FAIL    : A memory allocation failed.                     *
 *                                                                *
 * LIN_ILL_INPUT: Some input was illegal (see message).           *
 *                                                                *
 * LIN_NO_LMEM  : The linear solver's memory = NULL.              *
 *                                                                *
 ******************************************************************/

/* SUCCESS = 0  (defined above but listed here for completeness)  */
enum {LMEM_FAIL = -1, LIN_ILL_INPUT = -2, LIN_NO_MEM=-3, LIN_NO_LMEM=-4};


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
 *     of purpose, for each IDAS linear solver routine to be      *
 *     called in IDAS is given below the constant declarations    *
 *     that follow.                                               *
 *                                                                *
 ******************************************************************/

/* ida_linit return values */

#define LINIT_OK        0
#define LINIT_ERR      -1

/*******************************************************************
 *                                                                 *
 * int (*ida_linit)(IDAMem IDA_mem);                               *
 *-----------------------------------------------------------------*
 * The purpose of ida_linit is to allocate memory for the          *
 * solver-specific fields in the structure *(idamem->ida_lmem) and *
 * perform any needed initializations of solver-specific memory,   *
 * such as counters/statistics. An (*ida_linit) should return      *
 * LINIT_OK (== 0) if it has successfully initialized the IDAS     *
 * linear solver and LINIT_ERR (== -1) otherwise.                  *
 * These constants are defined above. If an error does occur, an   *
 * appropriate message should be sent to (idamem->errfp).          *
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * int (*ida_lsetup)(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,   *
 *                   N_Vector resp,                                *
 *                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); *
 *-----------------------------------------------------------------*
 * The job of ida_lsetup is to prepare the linear solver for       *
 * subsequent calls to ida_lsolve. Its parameters are as follows:  *
 *                                                                 *
 * idamem - problem memory pointer of type IDAMem. See the big     *
 *          typedef earlier in this file.                          *
 *                                                                 *
 *                                                                 *
 * yyp   - the predicted y vector for the current IDAS internal    *
 *         step.                                                   *
 *                                                                 *
 * ypp   - the predicted y' vector for the current IDAS internal   *
 *         step.                                                   *
 *                                                                 *
 * resp  - F(tn, yyp, ypp).                                        *
 *                                                                 *
 * tmp1, tmp2, tmp3 - temporary N_Vectors provided for use by      *
 *         ida_lsetup.                                             *
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
 *                   N_Vector ypcur, N_Vector rescur);             *
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
 * int (*ida_lsolveS)(IDAMem IDA_mem, N_Vector b, N_Vector ycur,   *
 *                    N_Vector ypcur, N_Vector rescur, int is);    *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * int (*ida_lperf)(IDAMem IDA_mem, int perftask);                 *
 *                                                                 *
 *-----------------------------------------------------------------*
 * ida_lperf is called two places in IDAS where linear solver      *
 * performance data is required by IDAS. For perftask = 0, an      *
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
