/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2005-01-24 23:55:28 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSPBCG linear solver module header file (public version)
 * -----------------------------------------------------------------
 */

#ifndef _KINSPBCG_H
#define _KINSPBCG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type : KINSpbcgPrecSetupFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner setup subroutine should
 * compute the right-preconditioner matrix P (stored in memory
 * block referenced by P_data pointer) used to form the
 * scaled preconditioned linear system:
 *
 *  (Df*J(uu)*(P^-1)*(Du^-1)) * (Du*P*x) = Df*(-F(uu))
 *
 * where Du and Df denote the diagonal scaling matrices whose
 * diagonal elements are stored in the vectors uscale and
 * fscale, repsectively.
 *
 * The preconditioner setup routine (referenced by KINSPBCG module
 * via pset (type KINSpbcgPrecSetupFn)) will not be called prior
 * to every call made to the psolve function, but will instead be
 * called only as often as necessary to achieve convergence of the
 * Newton iteration.
 *
 * Note: If the psolve routine requires no preparation, then a
 * preconditioner setup function need not be given.
 *
 *  uu  current iterate (unscaled) [input]
 *
 *  uscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for vector uu [input]
 *
 *  fval  vector (type N_Vector) containing result of nonliear
 *        system function evaluated at current iterate:
 *        fval = F(uu) [input]
 *
 *  fscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for fval [input]
 *
 *  P_data  pointer to user-allocated system memory block used
 *          for storage of preconditioner matrix-related data
 *          [output]
 *
 *  vtemp1/vtemp2  available scratch vectors (temporary storage)
 *
 * If successful, the function should return 0 (zero). If an error
 * occurs, then the routine should return a non-zero integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINSpbcgPrecSetupFn)(N_Vector uu, N_Vector uscale,
                                   N_Vector fval, N_Vector fscale,
                                   void *P_data, N_Vector vtemp1,
				   N_Vector vtemp2);

/*
 * -----------------------------------------------------------------
 * Type : KINSpbcgPrecSolveFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve subroutine (referenced
 * by KINSPBCG module via psolve (type KINSpbcgPrecSolveFn)) should
 * solve a (scaled) preconditioned linear system of the generic form
 * Pz = r, where P denotes the right-preconditioner matrix computed
 * by the pset routine.
 *
 *  uu  current iterate (unscaled) [input]
 *
 *  uscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for vector uu [input]
 *
 *  fval  vector (type N_Vector) containing result of nonliear
 *        system function evaluated at current iterate:
 *        fval = F(uu) [input]
 *
 *  fscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for fval [input]
 *
 *  vv  vector initially set to the right-hand side vector r, but
 *      which upon return contains a solution of the linear system
 *      P*z = r [input/output]
 *
 *  P_data  pointer to user-allocated system memory block used
 *          for storage of preconditioner matrix-related data
 *          [output]
 *
 *  vtemp  available scratch vector (volatile storage)
 *
 * If successful, the function should return 0 (zero). If a
 * recoverable error occurs, then the subroutine should return
 * a positive integer value. However, if an unrecoverable error
 * occurs, then the function should return a negative integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINSpbcgPrecSolveFn)(N_Vector uu, N_Vector uscale, 
                                   N_Vector fval, N_Vector fscale, 
                                   N_Vector vv, void *P_data,
                                   N_Vector vtemp);

/*
 * -----------------------------------------------------------------
 * Type : KINSpbcgJacTimesVecFn
 * -----------------------------------------------------------------
 * The (optional) user-supplied matrix-vector product subroutine
 * (referenced internally via jtimes (type KINSpbcgJacTimesVecFn))
 * is used to compute Jv = J(uu)*v (system Jacobian applied to a
 * given vector). If a user-defined routine is not given, then the
 * private KINSPBCG module routine named KINSpbcgDQJtimes is used
 * (see kinspbcg.c).
 *
 *  v  unscaled variant of vector to be multiplied by J(uu) [input]
 *
 *  Jv  vector containing result of matrix-vector product J(uu)*v
 *      [output]
 *
 *  uu  current iterate (unscaled) [input]
 *
 *  new_uu  flag (reset by user) indicating if the iterate uu
 *          has been updated in the interim - Jacobian needs
 *          to be updated/reevaluated, if appropriate, unless
 *          new_uu = FALSE [input/output]
 *
 *  J_data  pointer to user-allocated memory block where J(uu) data
 *          is to be stored [input]
 *
 * If successful, the function should return 0 (zero). If an error
 * occurs, then the routine should return a non-zero integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINSpbcgJacTimesVecFn)(N_Vector v, N_Vector Jv,
                                     N_Vector uu, booleantype *new_uu, 
                                     void *J_data);

/*
 * -----------------------------------------------------------------
 * Function : KINSpbcg
 * -----------------------------------------------------------------
 * KINSpbcg links the main KINSOL solver module with the SPBCG
 * linear solver module. The routine establishes the inter-module
 * interface by setting the generic KINSOL pointers linit, lsetup,
 * lsolve, and lfree to KINSpbcgInit, KINSpbcgSetup, KINSpbcgSolve,
 * and KINSpbcgFree, respectively.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  maxl  maximum allowable dimension of Krylov subspace (passing
 *        a value of 0 (zero) will cause the default value
 *        KINSPBCG_MAXL (predefined constant) to be used)
 *
 * If successful, KINSpbcg returns KINSPBCG_SUCCESS. If an error
 * occurs, then KINSpbcg returns an error code (negative integer
 * value).
 * -----------------------------------------------------------------
 */

int KINSpbcg(void *kinmem, int maxl);

/*
 * -----------------------------------------------------------------
 * KINSpbcg Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpbcg subroutine are the
 * following:
 *
 * KINSPBCG_SUCCESS : means the KINSPBCG linear solver module
 *                    (implementation of the Bi-CGSTAB method) was
 *                    successfully initialized - allocated system
 *                    memory and set shared variables to default
 *                    values [0]
 *
 * KINSPBCG_MEM_NULL : means a NULL KINSOL memory block pointer
 *                     was given (must call the KINCreate and
 *                     KINMalloc memory allocation subroutines
 *                     prior to calling KINSpbcg) [-1]
 *
 * KINSPBCG_MEM_FAIL : means either insufficient system resources
 *                     were available to allocate memory for the
 *                     main KINSPBCG data structure (type
 *                     KINSpbcgMemRec), or the SpbcgMalloc subroutine
 *                     failed (unable to allocate enough system
 *                     memory for vector storate and/or the main
 *                     SPBCG data structure (type SpbcgMemRec)) [-4]
 *
 * KINSPBCG_ILL_INPUT : means either a supplied parameter was invalid,
 *                      or the NVECTOR implementation is NOT
 *                      compatible [-3]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Input Specification Functions (KINSPBCG)
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs:
 *
 *       Function Name       |   Optional Input  [Default Value]
 *                           |
 * -----------------------------------------------------------------
 *                           |
 * KINSpbcgSetPrecSetupFn    | name of user-supplied routine used to
 *                           | compute a preconditioner matrix for
 *                           | the given linear system (pset)
 *                           | [NULL]
 *                           |
 * KINSpbcgSetPrecSolveFn    | name of user-supplied routine used to
 *                           | apply preconditioner to linear system
 *                           | (psolve)
 *                           | [NULL]
 *                           |
 * KINSpbcgSetPrecData       | pointer to user-allocated system
 *                           | memory that is passed to the pset
 *                           | and psolve routines
 *                           | [NULL]
 *                           |
 * KINSpbcgSetJacTimesVecFn  | name of user-supplied subroutine used
 *                           | to compute the matrix-vector product
 *                           | J(u)*v, where J denotes the system
 *                           | Jacobian (jtimes)
 *                           | [KINSpbcgDQJtimes] (see kinspbcg.c)
 *                           |
 * KINSpbcgSetJacData        | pointer to a user-allocated memory
 *                           | block that is passed to the jtimes
 *                           | routine
 *                           | [NULL]
 * -----------------------------------------------------------------
 */

int KINSpbcgSetPrecSetupFn(void *kinmem, KINSpbcgPrecSetupFn pset);
int KINSpbcgSetPrecSolveFn(void *kinmem, KINSpbcgPrecSolveFn psolve);
int KINSpbcgSetPrecData(void *kinmem, void *P_data);
int KINSpbcgSetJacTimesVecFn(void *kinmem, KINSpbcgJacTimesVecFn jtimes);
int KINSpbcgSetJacData(void *kinmem, void *J_data);

/*
 * -----------------------------------------------------------------
 * KINSpbcgSet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpbcgSet* subroutines
 * are the following:
 *
 * KINSPBCG_SUCCESS : means the associated variable was successfully
 *                    set [0]
 *
 * KINSPBCG_ILL_INPUT : means the supplied parameter was invalid
 *                      (check error message) [-3]
 *
 * KINSPBCG_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                     given [-1]
 *
 * KINSPBCG_LMEM_NULL : means system memory has not yet been allocated
 *                      for SPBCG (lmem == NULL) [-2]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Output Extraction Functions (KINSPBCG)
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistical information related to the KINSPBCG linear
 * solver:
 *
 *        Function Name       |      Returned Value
 *                            |
 * -----------------------------------------------------------------
 *                            |
 * KINSpbcgGetIntWorkSpace    | returns both integer workspace size
 *                            | (total number of long int-sized blocks
 *                            | of memory allocated by KINSPBCG for
 *                            | vector storage), and real workspace
 *                            | size (total number of realtype-sized
 *                            | blocks of memory allocated by KINSPBCG
 *                            | for vector storage)
 *                            |
 * KINSpbcgGetNumPrecEvals    | total number of preconditioner
 *                            | evaluations (number of calls made
 *                            | to the user-defined pset routine)
 *                            |
 * KINSpbcgGetNumPrecSolves   | total number of times preconditioner
 *                            | was applied to linear system (number
 *                            | of calls made to the user-supplied
 *                            | psolve function)
 *                            |
 * KINSpbcgGetNumLinIters     | total number of linear iterations
 *                            | performed
 *                            |
 * KINSpbcgGetNumConvFails    | total number of linear convergence
 *                            | failures
 *                            |
 * KINSpbcgGetNumJtimesEvals  | total number of times the matrix-
 *                            | vector product J(u)*v was computed
 *                            | (number of calls made to the jtimes
 *                            | subroutine)
 *                            |
 * KINSpbcgGetNumFuncEvals    | total number of evaluations of the
 *                            | system function F(u) (number of
 *                            | calls made to the user-supplied
 *                            | func routine by the KINSPBCG module
 *                            | member subroutines)
 *                            |
 * KINSpbcgGetLastFlag        | returns last flag returned by the
 *                            | linear solver
 * -----------------------------------------------------------------
 */

int KINSpbcgGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG);
int KINSpbcgGetNumPrecEvals(void *kinmem, long int *npevals);
int KINSpbcgGetNumPrecSolves(void *kinmem, long int *npsolves);
int KINSpbcgGetNumLinIters(void *kinmem, long int *nliters);
int KINSpbcgGetNumConvFails(void *kinmem, long int *nlcfails);
int KINSpbcgGetNumJtimesEvals(void *kinmem, long int *njvevals);
int KINSpbcgGetNumFuncEvals(void *kinmem, long int *nfevalsSG); 
int KINSpbcgGetLastFlag(void *kinmem, int *flag);

/*
 * -----------------------------------------------------------------
 * KINSpbcgGet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpbcgGet* subroutines
 * are the following:
 *
 * KINSPBCG_SUCCESS : means the routine exited normally [0]
 *
 * KINSPBCG_ILL_INPUT : means at least one input parameter was
 *                      invalid (check error message(s)) [-3]
 *
 * KINSPBCG_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                     given [-1]
 *
 * KINSPBCG_LMEM : means a NULL KINSPBCG memory block pointer was
 *                 given [-2]
 * -----------------------------------------------------------------
 */

#define KINSPBCG_SUCCESS 0

#define KINSPBCG_MEM_NULL  -1
#define KINSPBCG_LMEM_NULL -2
#define KINSPBCG_ILL_INPUT -3
#define KINSPBCG_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
