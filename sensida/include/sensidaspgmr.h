/******************************************************************
 *                                                                *
 * File          : sensidaspgmr.h                                 *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh            *
 * Version of    : 3 July 2001                                    *
 *----------------------------------------------------------------*
 *                                                                *
 * This is the header file for the IDA scaled, preconditioned     *
 * GMRES linear solver, IDASPGMR, for sensitivity analysis.       *
 * The linear system size is Ny, where Ny is the number of        *
 * equations contained in F(t,y,y',p) = 0.                        *
 *                                                                *
 * This file is simply the idaspgmr.h header file, along with a   *
 * function prototype for SensIDAspgmr().                         *
 *                                                                *
 * Note: The type integer must be large enough to store the value *
 * of the linear system size Ny.                                  *
 *                                                                *
 ******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _sensidaspgmr_h
#define _sensidaspgmr_h


#include <stdio.h>
#include "ida.h"
#include "llnltyps.h"
#include "spgmr.h"
#include "nvector.h"
#include "idaspgmr.h"

/*********************************************************************
 *                                                                   *
 * Function : SensIDASpgmr                                           *
 *-------------------------------------------------------------------*
 * A call to the SensIDASpgmr function links the main IDA integrator
 * with the SensIDASPGMR linear solver module.
 *
 * Its parameters are as follows:
 *
 * IDA_mem   is the pointer to IDA memory returned by SensIDAMalloc.
 *
 * precond   is the user's preconditioner setup routine. 
 *           It is used to evaluate and preprocess any Jacobian-related 
 *           data needed by the psolve routine.
 *           See the description of the type IDASpgmrPrecondFn.
 *           Pass NULL if no such data setup is required.
 *
 * psolve    is the user's preconditioner solve routine. 
 *           It is used to solve linear systems P z = r, where P is
 *           the preconditioner matrix.
 *           See the description of the type IDASpgmrPSolveFn.
 *           Pass NULL for psolve if no preconditioning is to be done.
 *           However, a preconditioner of some form is strongly encouraged.
 *
 * gstype    is the type of Gram-Schmidt orthogonalization to be used.
 *           This must be one of the two enumeration constants
 *           MODIFIED_GS or CLASSICAL_GS defined in iterativ.h.
 *           These correspond to using modified or classical
 *           Gram-Schmidt algorithms, respectively.
 *
 * maxl      is the maximum Krylov subspace dimension, an optional input.
 *           Pass 0 to use the default value, MIN(Ny, 5).
 *           Otherwise pass a positive integer.
 *
 * maxrs     is the maximum number of restarts to be used in the
 *           GMRES algorithm, an optional input.  
 *           maxrs must be a non-negative integer, or -1.
 *           Pass 0 to use the default value, which is 5.
 *           Pass -1 to use the value 0, meaning no restarts.
 *           In any case, maxrs will be restricted to the range 0 to
 *           Ny/maxl.
 *
 * eplifac   is a factor in the linear iteration convergence test 
 *           constant, an optional input.  
 *           Pass 0.0 to use the default, which is 1.0.  
 *           Otherwise eplifac must be a positive real number.
 *
 * dqincfac  is a factor in the increments to yy used in the difference
 *           quotient approximations to matrix-vector products Jv, an 
 *           optional input.
 *           Pass 0.0 to use the default, which is 1.0.  
 *           Otherwise dqincfac must be a positive real number.
 *
 * pdata     is a pointer to user preconditioner data.  
 *           This pointer is passed to precond and psolve every time 
 *           these routines are called.
 *
 * SensIDASpgmr returns either                                  
 *    SUCCESS = 0                 if successful, or
 *    IDA_SPGMR_FAIL = -1         if either IDA_mem was null or a
 *                                malloc failure occurred, or
 *    IDA_SPGMR_BAD_ARG = -2      if gstype was found illegal.
 *********************************************************************/

int SensIDASpgmr(void *IDA_mem, IDASpgmrPrecondFn precond, 
             IDASpgmrPSolveFn psolve, int gstype, int maxl, int maxrs,
             real eplifac, real dqincfac, void *pdata);


#endif

#ifdef __cplusplus
}
#endif
