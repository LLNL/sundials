/*******************************************************************
 *                                                                 *
 * File          : idaspgmr.h                                      *
 * Programmers   : Alan C. Hindmarsh and Allan G. Taylor           *
 * Version of    : 17 July 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/ida/LICENSE                           *
 *-----------------------------------------------------------------*
 * This is the header file for the IDA Scaled Preconditioned       *
 * GMRES linear solver module, IDASPGMR.                           *
 *                                                                 *
 *******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idaspgmr_h
#define _idaspgmr_h


#include <stdio.h>
#include "ida.h"
#include "sundialstypes.h"
#include "spgmr.h"
#include "nvector.h"

/******************************************************************
 *                                                                *           
 * Type : IDASpgmrPrecSetupFn                                     *
 *----------------------------------------------------------------*
 * The optional user-supplied functions PrecSetup and PrecSolve   *
 * together must define the left preconditoner matrix P           *
 * approximating the system Jacobian matrix                       *
 *    J = dF/dy + c_j*dF/dy'                                      *
 * (where the DAE system is F(t,y,y') = 0), and solve the linear  *
 * systems P z = r.   PrecSetup is to do any necessary setup      *
 * operations, and PrecSolve is to compute the solution of        *
 * P z = r.                                                       *
 *                                                                *
 * The preconditioner setup function PrecSetup is to evaluate and *
 * preprocess any Jacobian-related data needed by the             *
 * preconditioner solve function PrecSolve.  This might include   *
 * forming a crude approximate Jacobian, and performing an LU     *
 * factorization on it.  This function will not be called in      *
 * advance of every call to PrecSolve, but instead will be called *
 * only as often as necessary to achieve convergence within the   *
 * Newton iteration in IDA.  If the PrecSolve function needs no   *
 * preparation, the PrecSetup function can be NULL.               *
 *                                                                *
 * Each call to the PrecSetup function is preceded by a call to   *
 * the system function res with the same (t,y,y') arguments.      *
 * Thus the PrecSetup function can use any auxiliary data that is *
 * computed and saved by the res function and made accessible     *
 * to PrecSetup.                                                  *
 *                                                                *
 * A preconditioner setup function PrecSetup must have the        *
 * prototype given below.  Its parameters are as follows:         *
 *                                                                *
 * tt  is the current value of the independent variable t.        *
 *                                                                *
 * yy  is the current value of the dependent variable vector,     *
 *        namely the predicted value of y(t).                     *
 *                                                                *
 * yp  is the current value of the derivative vector y',          *
 *        namely the predicted value of y'(t).                    *
 *                                                                *
 * rr  is the current value of the residual vector F(t,y,y').     *
 *                                                                *
 * cj  is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * pdata  is a pointer to user preconditioner data - the same as  *
 *        the pdata parameter passed to IDASpgmr.                 *
 *                                                                *
 * tmp1, tmp2, tmp3 are pointers to vectors of type N_Vector      *
 * which can be used by an IDASpgmrPrecSetupFn routine            *
 * as temporary storage or work space.                            *
 *                                                                *
 * NOTE: If the user's preconditioner needs other quantities,     *
 *     they are accessible as follows: hcur (the current stepsize)*
 *     and ewt (the error weight vector) are accessible through   *
 *     IDAGetCurrentStep and IDAGetErrWeights, respectively (see  *
 *     ida.h). The unit roundoff is available through a call to   *
 *     UnitRoundoff.                                              *
 *                                                                *
 * The IDASpgmrPrecSetupFn should return                          *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * In the case of a recoverable error return, IDA will attempt to *
 * recover by reducing the stepsize (which changes cj).           *
 ******************************************************************/
  
typedef int (*IDASpgmrPrecSetupFn)(realtype tt, 
                                   N_Vector yy, N_Vector yp, N_Vector rr, 
                                   realtype cj, void *pdata,
                                   N_Vector tmp1, N_Vector tmp2, 
                                   N_Vector tmp3);

/******************************************************************
 *                                                                *           
 * Type : IDASpgmrPrecSolveFn                                     *
 *----------------------------------------------------------------*
 * The optional user-supplied function PrecSolve must compute a   *
 * solution to the linear system P z = r, where P is the left     *
 * preconditioner defined by the user.  If no preconditioning     *
 * is desired, pass NULL for PrecSolve to IDASpgmr.               *
 *                                                                *
 * A preconditioner solve function PrecSolve must have the        *
 * prototype given below.  Its parameters are as follows:         *
 *                                                                *
 * tt  is the current value of the independent variable t.        *
 *                                                                *
 * yy  is the current value of the dependent variable vector y.   *
 *                                                                *
 * yp  is the current value of the derivative vector y'.          *
 *                                                                *
 * rr  is the current value of the residual vector F(t,y,y').     *
 *                                                                *
 * cj  is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * pdata  is a pointer to user preconditioner data - the same as  *
 *        the pdata parameter passed to IDASpgmr.                 *
 *                                                                *
 * delta  is an input tolerance for use by PrecSolve if it uses an*
 *        iterative method in its solution.   In that case, the   *
 *        the residual vector r - P z of the system should be     *
 *        made less than delta in weighted L2 norm, i.e.,         *
 *            sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta .            *
 *        Note: the error weight vector ewt can be obtained       *
 *        through a call to the routine IDAGetErrWeights.         *
 *                                                                *
 * rvec   is the input right-hand side vector r.                  *
 *                                                                *
 * zvec   is the computed solution vector z.                      *
 *                                                                *
 * tmp  is an N_Vector which can be used by the PrecSolve         *
 * routine as temporary storage or work space.                    *
 *                                                                *
 *                                                                *
 * The IDASpgmrPrecSolveFn should return                          *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * Following a recoverable error, IDA will attempt to recover by  *
 * updating the preconditioner and/or reducing the stepsize.      *
 *                                                                *
 ******************************************************************/
  
typedef int (*IDASpgmrPrecSolveFn)(realtype tt, 
                                   N_Vector yy, N_Vector yp, N_Vector rr, 
                                   N_Vector rvec, N_Vector zvec,
                                   realtype cj, realtype delta,
                                   void *pdata, N_Vector tmp);

/******************************************************************
 *                                                                *           
 * Type : IDASpgmrJacTimesVecFn                                   *
 *----------------------------------------------------------------*
 * The user-supplied function jtimes is to generate the product   *
 * J*v for given v, where J is the Jacobian matrix                *
 *    J = dF/dy + c_j*dF/dy'                                      *
 *  or an approximation to it, and v is a given vector.           *
 * It should return 0 if successful and a nonzero int otherwise.  *
 *                                                                *
 * A function jtimes must have the prototype given below. Its     *
 * parameters are as follows:                                     *
 *                                                                *
 *   v    is the N_Vector to be multiplied by J.                  *
 *                                                                *
 *   Jv   is the output N_Vector containing J*v.                  *
 *                                                                *
 *   t    is the current value of the independent variable.       *
 *                                                                *
 *   yy   is the current value of the dependent variable vector,  *
 *        namely the predicted value of y(t).                     *
 *                                                                *
 *   yp   is the current value of the derivative vector y',       *
 *        namely the predicted value of y'(t).                    *
 *                                                                *
 *   rr   is the current value of the residual vector F(t,y,y').  *
 *                                                                *
 *   cj   is the scalar in the system Jacobian, proportional      *
 *        to 1/hh.                                                *
 *                                                                *
 *   jac_data is a pointer to user Jacobian data, the same as the *
 *        pointer passed to CVSpgmr.                              *
 *                                                                *
 *   tmp1, tmp2 are two N_Vectors which can be used by Jtimes for *
 *         work space.                                            *
 *                                                                *
 ******************************************************************/

typedef int (*IDASpgmrJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
                                     N_Vector yy, N_Vector yp, N_Vector rr,
                                     realtype c_j, void *jac_data, 
                                     N_Vector tmp1, N_Vector tmp2);

/******************************************************************
 *                                                                *
 * Function : IDASpgmr                                            *
 *----------------------------------------------------------------*
 * A call to the IDASpgmr function links the main IDA integrator  *
 * with the IDASPGMR linear solver module.  Its parameters are    *
 * as follows:                                                    *
 *                                                                *
 * IDA_mem   is the pointer to IDA memory returned by IDACreate.  *
 *                                                                *
 * maxl      is the maximum Krylov subspace dimension, an         *
 *           optional input.  Pass 0 to use the default value,    *
 *           MIN(Neq, 5).  Otherwise pass a positive integer.     *
 *                                                                *
 * The return values of IDASpgmr are:                             *
 *    SUCCESS       = 0  if successful                            *
 *    LMEM_FAIL     = -1 if there was a memory allocation failure *
 *    LIN_ILL_INPUT = -2 if there was illegal input.              *
 *                                                                *
 ******************************************************************/

int IDASpgmr(void *ida_mem, int maxl);

/******************************************************************
 * Optional inputs to the IDASPGMR linear solver                  *
 *----------------------------------------------------------------*
 *                                                                *
 * IDASpgmrSetGSType specifies the type of Gram-Schmidt           *
 *           orthogonalization to be used. This must be one of    *
 *           the two enumeration constants MODIFIED_GS or         *
 *           CLASSICAL_GS defined in iterativ.h. These correspond *
 *           to using modified Gram-Schmidt and classical         *
 *           Gram-Schmidt, respectively.                          *
 *           Default value is MODIFIED_GS.                        *
 * IDASpgmrSetMaxRestarts specifies the maximum number of restarts*
 *           to be used in the GMRES algorithm.  maxrs must be a  *
 *           non-negative integer.  Pass 0 to specify no restarts.*
 *           Default is 5.                                        *
 * IDASpgmrSetEpsLin specifies the factor in the linear iteration *
 *           convergence test constant.                           *
 *           Default is 0.05                                      *
 * IDASpgmrSetIncrementFactor specifies a factor in the increments*
 *           to yy used in the difference quotient approximations *
 *           to matrix-vector products Jv.                        *
 *           Default is 1.0                                       *
 * IDASpgmrSetPrecSetupFn specifies the PrecSetup function.       *
 *           Default is NULL.                                     *
 * IDASpgmrSetPrecSolveFn specifies the PrecSolve function.       *
 *           Default is NULL.                                     *
 * IDASpgmrSetPrecData specifies a pointer to user preconditioner *
 *           data. This pointer is passed to PrecSetup and        *
 *           PrecSolve every time these routines are called.      *
 *           Default is NULL.                                     *
 * IDASpgmrSetJacTimesVecFn specifies the jtimes function.        *
 *           Default is to use an internal finite difference      *
 *           approximation routine.                               *
 * IDASpgmrSetJacData specifies a pointer to user Jacobian data.  *
 *           This pointer is passed to jtimes every time this     *
 *           routine is called.                                   *
 *           Default is NULL.                                     *
 *                                                                *
 ******************************************************************/

int IDASpgmrSetGSType(void *ida_mem, int gstype);
int IDASpgmrSetMaxRestarts(void *ida_mem, int maxrs);
int IDASpgmrSetEpsLin(void *ida_mem, realtype eplifac);
int IDASpgmrSetIncrementFactor(void *ida_mem, realtype dqincfac);
int IDASpgmrSetPrecSetupFn(void *ida_mem, IDASpgmrPrecSetupFn pset);
int IDASpgmrSetPrecSolveFn(void *ida_mem, IDASpgmrPrecSolveFn psolve);
int IDASpgmrSetPrecData(void *ida_mem, void *pdata);
int IDASpgmrSetJacTimesVecFn(void *ida_mem, IDASpgmrJacTimesVecFn jtimes);
int IDASpgmrSetJacData(void *ida_mem, void *jdata);

/******************************************************************
 * Optional outputs from the IDASPGMR linear solver               *
 *----------------------------------------------------------------*
 *                                                                *
 * IDASpgmrGetIntWorkSpace returns the integer workspace used by  *
 *     IDASPGMR.                                                  *
 * IDASpgmrGetRealWorkSpace returns the real workspace used by    *
 *     IDASPGMR.                                                  *
 * IDASpgmrGetNumPrecEvals returns the number of preconditioner   *
 *     evaluations, i.e. the number of calls made to PrecSetup    *
 *     with jok==FALSE.                                           *
 * IDASpgmrGetNumPrecSolves returns the number of calls made to   *
 *     PrecSolve.                                                 *
 * IDASpgmrGetNumLinIters returns the number of linear iterations.*
 * IDASpgmrGetNumConvFails returns the number of linear           *
 *     convergence failures.                                      *
 * IDASpgmrGetNumJtimesEvals returns the number of calls to jtimes*
 * IDASpgmrGetNumResEvals returns the number of calls to the user *
 *     res routine due to finite difference Jacobian times vector *
 *     evaluation.                                                *
 *                                                                *
 ******************************************************************/

int IDASpgmrGetIntWorkSpace(void *ida_mem, long int *leniwSG);
int IDASpgmrGetRealWorkSpace(void *ida_mem, long int *lenrwSG);
int IDASpgmrGetNumPrecEvals(void *ida_mem, int *npevals);
int IDASpgmrGetNumPrecSolves(void *ida_mem, int *npsolves);
int IDASpgmrGetNumLinIters(void *ida_mem, int *nliters);
int IDASpgmrGetNumConvFails(void *ida_mem, int *nlcfails);
int IDASpgmrGetNumJtimesEvals(void *ida_mem, int *njvevals);
int IDASpgmrGetNumResEvals(void *ida_mem, int *nrevalsSG); 

/******************************************************************
 *                                                                *           
 * Types : IDASpgmrMemRec, IDASpgmrMem                            *
 *----------------------------------------------------------------*
 * The type IDASpgmrMem is pointer to an IDASpgmrMemRec. This     *
 * structure contains IDASpgmr solver-specific data.              *
 *                                                                *
 ******************************************************************/


typedef struct {

  int  g_gstype;       /* type of Gram-Schmidt orthogonalization       */
  realtype g_sqrtN;    /* sqrt(N)                                      */
  int  g_maxl;         /* maxl = maximum dimension of the Krylov space */
  int  g_maxrs;        /* maxrs = max. number of GMRES restarts        */
  realtype g_eplifac;  /* eplifac = linear convergence factor          */
  realtype g_dqincfac; /* dqincfac = optional increment factor in Jv   */
  realtype g_epslin;   /* SpgrmSolve tolerance parameter               */

  int g_resflag;       /* flag from last res call                      */
  int g_npe;           /* npe = total number of precond calls          */   
  int g_nli;           /* nli = total number of linear iterations      */
  int g_nps;           /* nps = total number of psolve calls           */
  int g_ncfl;          /* ncfl = total number of convergence failures  */
  int g_nreSG;         /* nreSG = total number of calls to res         */    
  int g_njtimes;       /* njtimes = total number of calls to jtimes    */

  int g_nst0;          /* nst0 = saved nst (for performance monitor)   */   
  int g_nni0;          /* nni0 = saved nni (for performance monitor)   */   
  int g_nli0;          /* nli0 = saved nli (for performance monitor)   */   
  int g_ncfn0;         /* ncfn0 = saved ncfn (for performance monitor) */   
  int g_ncfl0;         /* ncfl0 = saved ncfl (for performance monitor) */   
  int g_nwarn;         /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector g_ytemp;    /* temp vector used by IDAAtimesDQ              */ 
  N_Vector g_yptemp;   /* temp vector used by IDAAtimesDQ              */ 
  N_Vector g_xx;       /* temp vector used by IDASpgmrSolve            */
  N_Vector g_ycur;     /* IDA current y vector in Newton iteration     */
  N_Vector g_ypcur;    /* IDA current yp vector in Newton iteration    */
  N_Vector g_rcur;     /* rcur = F(tn, ycur, ypcur)                    */

  IDASpgmrPrecSetupFn g_pset;     /* pset = user-supplied routine      */
                                  /* to compute a preconditioner       */

  IDASpgmrPrecSolveFn g_psolve;   /* psolve = user-supplied routine to */
                                  /* solve preconditioner linear system*/

  void *g_pdata;                  /* pdata passed to psolve and precond*/
  SpgmrMem g_spgmr_mem;           /* spgmr_mem is memory used by the   */
                                  /* generic Spgmr solver              */

  IDASpgmrJacTimesVecFn g_jtimes; /* Jacobian*vector routine           */ 
  void *g_jdata;                  /* data passed to Jtimes             */

} IDASpgmrMemRec, *IDASpgmrMem;

#endif

#ifdef __cplusplus
}
#endif
