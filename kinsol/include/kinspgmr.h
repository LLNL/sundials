/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2004-05-03 21:24:47 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                 Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/kinsol/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for the KINSol scaled, preconditioned
 * GMRES linear solver, KINSpgmr.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _kinspgmr_h
#define _kinspgmr_h

#include <stdio.h>
#include "kinsol.h"
#include "spgmr.h"
#include "sundialstypes.h"
#include "nvector.h"

/*
 * -----------------------------------------------------------------
 * KINSpgmr solver constants
 * -----------------------------------------------------------------
 * KINSPGMR_MAXL: default value for the maximum Krylov dimension
                  is KINSPGMR_MAXL.
 * -----------------------------------------------------------------
 */

#define KINSPGMR_MAXL 10

/*
 * -----------------------------------------------------------------
 * constants for error returns from KINSpgmr
 * -----------------------------------------------------------------
 */

#define KIN_MEM_NULL      -1
#define KINSPGMR_MEM_FAIL -2
#define SPGMR_MEM_FAIL    -3

/*
 * -----------------------------------------------------------------
 * Type : KINSpgmrPrecSetupFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner setup function
 * KINSpgmrPrecSetupFn and the user-supplied preconditioner solve
 * function KINSpgmrPrecSolveFn together must define the right
 * preconditoner matrix P chosen so as to provide an easier system
 * for the Krylov solver to solve. *PrecSetupFn is called to provide
 * any matrix data required by the subsequent call(s) to
 * *PrecSolveFn. The data is stored in the memory allocated to
 * P_data and the structuring of that memory is up to the user. More
 * specifically, the user-supplied preconditioner setup function
 * *PrecSetupFn is supposed to evaluate and preprocess any Jacobian-
 * related data needed by the preconditioner solve function
 * *PrecSolveFn. This might include forming a crude approximate
 * Jacobian, and performing an LU factorization on the resulting
 * approximation to J. This function will not be called in advance of
 * every call to *PrecSolveFn, but instead will be called only as
 * often as necessary to achieve convergence within the Newton
 * iteration in KINSol. If the *PrecSolveFn function needs no
 * preparation, then the *PrecSetupFn function can be NULL.
 *
 * *PrecSetupFn should not modify the contents of the arrays uu or
 * fval as those arrays are used elsewhere in the iteration process.
 *
 * Each call to the *PrecSetupFn function is preceded by a call to
 * the system function func. Thus the *PrecSetupFn function can use
 * any auxiliary data that is computed by the func function and
 * saved in a way accessible to *PrecSetupFn.
 *
 * The two scaling arrays, fscale and uscale, are provided to the
 * *PrecSetupFn function for possible use in approximating Jacobian
 * data, e.g. by difference quotients. Note: These arrays should also
 * not be altered.
 *
 * A *PrecSetupFn function must have the prototype given below and
 * must also have the following parameters:
 *
 * uu      an N_Vector giving the current iterate for the system
 *
 * uscale  an N_Vector giving the diagonal entries of the uu
 *         scaling matrix
 *
 * fval    an N_Vector giving the current function value
 *
 * fscale  an N_Vector giving the diagonal entries of the func
 *         scaling matrix
 *
 * P_data  is a pointer to user data (P_data)
 *
 * vtemp1, vtemp2  two N_Vector temporaries for use by pset
 *
 * Note: The value to be returned by the *PrecSetupFn function is
 * a flag indicating whether it was successful. This value should
 * be:
 *   0 if successful, or
 *   1 if failure - KINSol stops
 * -----------------------------------------------------------------
 */

typedef int (*KINSpgmrPrecSetupFn)(N_Vector uu, N_Vector uscale,
                                   N_Vector fval, N_Vector fscale,
                                   void *P_data,
                                   N_Vector vtemp1, N_Vector vtemp2);

/*
 * -----------------------------------------------------------------
 * Type : KINSpgmrPrecSolveFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve function
 * KINSpgmrPrecSolveFn is used to solve a linear system Px = r in
 * which the matrix P is the (right) preconditioner matrix P.
 *
 * *PrecSolveFn should not modify the contents of the iterate array
 * uu or the current function value array fval as those are used
 * elsewhere in the iteration process.
 *
 * A *PrecSolveFn function must have the prototype given below and
 * the following parameters:
 *
 * uu      an N_Vector giving the current iterate for the system
 *
 * uscale  an N_Vector giving the diagonal entries of the uu
 *         scaling matrix
 *
 * fval    an N_Vector giving the current function value
 *
 * fscale  an N_Vector giving the diagonal entries of the func
 *         scaling matrix
 *
 * vv      an N_Vector to hold the RHS vector r on input and
 *         the result vector x on return
 *
 * P_data  is a pointer to user data (P_data)
 *
 * vtemp   an N_Vector used for work space
 *
 * Note: The value to be returned by the *PrecSolveFn function is
 * a flag indicating whether it was successful. This value should
 * be:
 *   0 if successful, or
 *   1 if failure - KINSol stops
 * -----------------------------------------------------------------
 */
  
typedef int (*KINSpgmrPrecSolveFn)(N_Vector uu, N_Vector uscale, 
                                   N_Vector fval, N_Vector fscale, 
                                   N_Vector vv, void *P_data,
                                   N_Vector vtemp);

/*
 * -----------------------------------------------------------------
 * Type : KINSpgmrJacTimesVecFn
 * -----------------------------------------------------------------
 * The user-supplied J times v routine (optional) where J is the
 * Jacobian matrix dF/du, or an approximation to it, and v is a
 * given vector. This routine computes the product Jv = J*v. It
 * should return 0 if successful and a nonzero integer otherwise.
 *
 * v   is the N_Vector to be multiplied by J
 *     (preconditioned and unscaled as received)
 *
 * Jv  is the N_Vector resulting from the application of J to v
 *
 * uu  is the N_Vector equal to the current iterate u
 *
 * new_uu  is an input flag indicating whether or not the uu
 *         vector has been changed since the last call to this
 *         function. If this function computes and saves Jacobian
 *         data, then this computation can be skipped if
 *         new_uu = FALSE.
 *
 * J_data  is a pointer to the structure used to handle data for
 *         the evaluation of the jacobian of func
 * -----------------------------------------------------------------
 */

typedef int (*KINSpgmrJacTimesVecFn)(N_Vector v, N_Vector Jv,
                                     N_Vector uu, booleantype *new_uu, 
                                     void *J_data);
  
 
/*
 * -----------------------------------------------------------------
 * Function : KINSpgmr
 * -----------------------------------------------------------------
 * A call to the KINSpgmr function links the main KINSol solver
 * with the KINSpgmr linear solver. Among other things, it sets
 * the generic names linit, lsetup, lsolve, and lfree to the
 * specific names for this package:
 *   KINSpgmrInit
 *   KINSpgmrSetup
 *   KINSpgmrSolve
 *   KINSpgmrFree
 *
 * kinmem  is the pointer to KINSol memory returned by KINCreate
 *
 * maxl  is the maximum Krylov dimension
 *       Note: This is an optional input to the KINSpgmr solver and
 *       0 may be passed to use the default value KINSPGMR_MAXL.
 *
 * KINSpgmr returns an flag equal to one of the following:
 *   SUCCESS, KIN_MEM_NULL, KINSPGMR_MEM_FAIL, SPGMR_MEM_FAIL
 * -----------------------------------------------------------------
 */

int KINSpgmr(void *kinmem,  int maxl);

/*
 * -----------------------------------------------------------------
 * optional inputs to the KINSPGMR linear solver
 * -----------------------------------------------------------------
 * KINSpgmrSetMaxRestarts specifies the maximum number of restarts
 *                        to be used in the GMRES algorithm. maxrs
 *                        must be a non-negative integer. Pass 0 to
 *                        specify no restarts.
 *                        [Default: 0]
 *
 * KINSpgmrSetPrecSetupFn specifies the *PrecSetupFn function.
 *                        [Default: NULL]
 *
 * KINSpgmrSetPrecSolveFn specifies the *PrecSolveFn function
 *                        [Default: NULL]
 *
 * KINSpgmrSetPrecData specifies a pointer to user preconditioner
 *                     data. This pointer is passed to *PrecSetupFn
 *                     and *PrecSolveFn every time these routines
 *                     are called.
 *                     [Default: NULL]
 *
 * KINSpgmrSetJacTimesVecFn specifies the jtimes function
 *                          [Default: use an internal finite
 *                          difference approximation routine]
 *
 * KINSpgmrSetJAcData specifies a pointer to user Jac times vec
 *                    data. This pointer is passed to jtimes every
 *                    time this routine is called.
 *                    [Default: NULL]
 * -----------------------------------------------------------------
 */

int KINSpgmrSetMaxRestarts(void *kinmem, int maxrs);
int KINSpgmrSetPrecSetupFn(void *kinmem, KINSpgmrPrecSetupFn pset);
int KINSpgmrSetPrecSolveFn(void *kinmem, KINSpgmrPrecSolveFn psolve);
int KINSpgmrSetPrecData(void *kinmem, void *P_data);
int KINSpgmrSetJacTimesVecFn(void *kinmem, KINSpgmrJacTimesVecFn jtimes);
int KINSpgmrSetJacData(void *kinmem, void *J_data);

/*
 * -----------------------------------------------------------------
 * optional outputs from the KINSPGMR linear solver
 * -----------------------------------------------------------------
 * KINSpgmrGetIntWorkSpace returns the integer workspace used by
 *                         KINSPGMR
 *
 * KINSpgmrGetRealWorkSpace returns the real workspace used by
 *                          KINSPGMR
 *
 * KINSpgmrGetNumPrecEvals returns the number of preconditioner
 *                         evaluations, i.e. the number of calls
 *                         made to *PrecSetupFn
 *
 * KINSpgmrGetNumPrecSolves returns the number of calls made to
 *                          *PrecSolveFn
 *
 * KINSpgmrGetNumLinIters returns the number of linear iterations
 *
 * KINSpgmrGetNumConvFails returns the number of linear convergence
 *                         failures
 *
 * KINSpgmrGetNumJtimesEvals returns the number of calls to jtimes
 *
 * KINSpgmrGetNumFuncEvals returns the number of calls to the user
 *                         res routine due to finite difference
 *                         Jacobian times vector evaluation
 * -----------------------------------------------------------------
 */

int KINSpgmrGetIntWorkSpace(void *kinmem, long int *leniwSG);
int KINSpgmrGetRealWorkSpace(void *kinmem, long int *lenrwSG);
int KINSpgmrGetNumPrecEvals(void *kinmem, long int *npevals);
int KINSpgmrGetNumPrecSolves(void *kinmem, long int *npsolves);
int KINSpgmrGetNumLinIters(void *kinmem, long int *nliters);
int KINSpgmrGetNumConvFails(void *kinmem, long int *nlcfails);
int KINSpgmrGetNumJtimesEvals(void *kinmem, long int *njvevals);
int KINSpgmrGetNumFuncEvals(void *kinmem, long int *nfevalsSG); 

/*
 * -----------------------------------------------------------------
 * Types : KINSpgmrMemRec and KINSpgmrMem
 * -----------------------------------------------------------------
 * The type KINSpgmrMem is a pointer to a KINSpgmrMemRec data type.
 * This structure contains KINSpgmr solver-specific data.
 * -----------------------------------------------------------------
 */

typedef struct {

  int  g_maxl;           /* maxl = maximum dimension of the Krylov subspace */     
  int  g_pretype;        /* preconditioning type --for Spgmr call           */
  int  g_gstype;         /* gram-schmidt type --  for Spgmr call            */
  booleantype g_new_uu;  /* flag indicating that a new uu has been created --
                            indicating that a call to generate a new
                            user-supplied Jacobian is required              */
  int g_maxlrst;       /* max number of linear solver restarts allowed      */
  long int g_nli;      /* nli = total number of linear iterations           */
  long int g_npe;      /* npe = total number of pset calls                  */
  long int g_nps;      /* nps = total number of psolve calls                */
  long int g_ncfl;     /* ncfl = total number of convergence failures       */
  long int g_nfeSG;    /* nfeSG = total number of calls to func             */    
  long int g_njtimes;  /* njtimes = total number of calls to jtimes         */
  KINSpgmrPrecSetupFn g_pset;      /* *PrecSetupFn = user-supplied routine
                                      to compute a preconditioner           */
  KINSpgmrPrecSolveFn g_psolve;    /* *PrecSolveFn = user-supplied routine
                                      to solve preconditioned linear system */ 
  KINSpgmrJacTimesVecFn g_jtimes;  /* jtimes = user-supplied routine to
                                      optionally compute the product J*v
                                      as required                           */ 
  void *g_P_data;  /* P_data is a memory block passed to psolve and pset    */
  void *g_J_data;  /* J_data is a memory block passed to jtimes             */
  SpgmrMem g_spgmr_mem;  /* spgmr_mem is memory used by the generic Spgmr
                            solver                                          */
} KINSpgmrMemRec, *KINSpgmrMem;

#endif
#ifdef __cplusplus
}
#endif
