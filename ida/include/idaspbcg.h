/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2005-04-27 21:38:01 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the public header file for the scaled preconditioned
 * Bi-CGSTAB linear solver module, IDASPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPBCG_H
#define _IDASPBCG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idaspils.h"
#include "spbcg.h"
#include "sundialstypes.h"
#include "nvector.h"

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcg
 * -----------------------------------------------------------------
 * A call to the IDASpbcg function links the main integrator with
 * the IDASPBCG linear solver module. Its parameters are as
 * follows:
 *
 * IDA_mem   is the pointer to memory block returned by IDACreate.
 *
 * maxl      is the maximum Krylov subspace dimension, an
 *           optional input. Pass 0 to use the default value.
 *           Otherwise pass a positive integer.
 *
 * The return values of IDASpbcg are:
 *    IDASPBCG_SUCCESS    if successful
 *    IDASPBCG_MEM_NULL   if the ida memory was NULL
 *    IDASPBCG_MEM_FAIL   if there was a memory allocation failure
 *    IDASPBCG_ILL_INPUT  if there was illegal input.
 * -----------------------------------------------------------------
 */

int IDASpbcg(void *ida_mem, int maxl);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDASPBCG linear solver
 * -----------------------------------------------------------------
 * IDASpbcgSetPreconditionr specifies the PrecSetup and PrecSolve
 *                   functions and the pointer to user data that
 *                   is pased to these functinos whenever they
 *                   are called.
 *                   Default is NULL for all three.
 * IDASpbcgSetJacTimesVecFn specifies the jtimes function and the 
 *                   pointer to user data that is passed to jtimes
 *                   whenever it is called.
 *                   Default is to use an internal finite
 *                   difference approximation routine.
 * IDASpbcgSetEpsLin specifies the factor in the linear iteration
 *                   convergence test constant.
 *                   Default is 0.05.
 * IDASpbcgSetIncrementFactor specifies a factor in the increments
 *                   to yy used in the difference quotient
 *                   approximations to matrix-vector products Jv.
 *                   Default is 1.0.
 *
 * The return value of IDASpbcgSet* is one of:
 *    IDASPBCG_SUCCESS   if successful
 *    IDASPBCG_MEM_NULL  if the ida memory was NULL
 *    IDASPBCG_LMEM_NULL if the idaspbcg memory was NULL
 * -----------------------------------------------------------------
 */

int IDASpbcgSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset,
                              IDASpilsPrecSolveFn psolve, void *prec_data);
int IDASpbcgSetJacTimesVecFn(void *ida_mem, 
                             IDASpilsJacTimesVecFn jtimes, void *jac_data);
int IDASpbcgSetEpsLin(void *ida_mem, realtype eplifac);
int IDASpbcgSetIncrementFactor(void *ida_mem, realtype dqincfac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDASPBCG linear solver
 * -----------------------------------------------------------------
 * IDASpbcgGetWorkSpace returns the real and integer workspace used
 *                      by IDASPBCG.
 * IDASpbcgGetNumPrecEvals returns the number of preconditioner
 *                         evaluations, i.e. the number of calls made
 *                         to PrecSetup with jok==FALSE.
 * IDASpbcgGetNumPrecSolves returns the number of calls made to
 *                          PrecSolve.
 * IDASpbcgGetNumLinIters returns the number of linear iterations.
 * IDASpbcgGetNumConvFails returns the number of linear
 *                         convergence failures.
 * IDASpbcgGetNumJtimesEvals returns the number of calls to jtimes
 * IDASpbcgGetNumResEvals returns the number of calls to the user
 *                        res routine due to finite difference Jacobian
 *                        times vector evaluation.
 * IDASpbcgGetLastFlag returns the last error flag set by any of
 *                     the IDASPBCG interface functions.
 *
 * The return value of IDASpbcgGet* is one of:
 *    IDASPBCG_SUCCESS   if successful
 *    IDASPBCG_MEM_NULL  if the ida memory was NULL
 *    IDASPBCG_LMEM_NULL if the idaspbcg memory was NULL
 * -----------------------------------------------------------------
 */

int IDASpbcgGetWorkSpace(void *ida_mem, long int *lenrwSG, long int *leniwSG);
int IDASpbcgGetNumPrecEvals(void *ida_mem, long int *npevals);
int IDASpbcgGetNumPrecSolves(void *ida_mem, long int *npsolves);
int IDASpbcgGetNumLinIters(void *ida_mem, long int *nliters);
int IDASpbcgGetNumConvFails(void *ida_mem, long int *nlcfails);
int IDASpbcgGetNumJtimesEvals(void *ida_mem, long int *njvevals);
int IDASpbcgGetNumResEvals(void *ida_mem, long int *nrevalsSG); 
int IDASpbcgGetLastFlag(void *ida_mem, int *flag);

/* IDASPBCG return values */

#define IDASPBCG_SUCCESS     0
#define IDASPBCG_MEM_NULL   -1 
#define IDASPBCG_LMEM_NULL  -2 
#define IDASPBCG_ILL_INPUT  -3
#define IDASPBCG_MEM_FAIL   -4

#ifdef __cplusplus
}
#endif

#endif
