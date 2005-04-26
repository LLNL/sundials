/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-04-26 20:31:27 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the common header file for the Scaled Preconditioned
 * Iterative Linear Solvers in KINSOL.
 * -----------------------------------------------------------------
 */

#ifndef _KINSPILS_H
#define _KINSPILS_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type : KINSpilsPrecSetupFn
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
 * The preconditioner setup routine (referenced by iterative linear
 * solver modules via pset (type KINSpilsPrecSetupFn)) will not be
 * called prior to every call made to the psolve function, but will
 * instead be called only as often as necessary to achieve convergence
 * of the Newton iteration.
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

typedef int (*KINSpilsPrecSetupFn)(N_Vector uu, N_Vector uscale,
                                   N_Vector fval, N_Vector fscale,
                                   void *P_data, N_Vector vtemp1,
				   N_Vector vtemp2);

/*
 * -----------------------------------------------------------------
 * Type : KINSpilsPrecSolveFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve subroutine (referenced
 * by iterative linear solver modules via psolve (type
 * KINSpilsPrecSolveFn)) should solve a (scaled) preconditioned
 * linear system of the generic form P*z = r, where P denotes the
 * right-preconditioner matrix computed by the pset routine.
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

typedef int (*KINSpilsPrecSolveFn)(N_Vector uu, N_Vector uscale, 
                                   N_Vector fval, N_Vector fscale, 
                                   N_Vector vv, void *P_data,
                                   N_Vector vtemp);

/*
 * -----------------------------------------------------------------
 * Type : KINSpilsJacTimesVecFn
 * -----------------------------------------------------------------
 * The (optional) user-supplied matrix-vector product subroutine
 * (referenced internally via jtimes (type KINSpilsJacTimesVecFn))
 * is used to compute Jv = J(uu)*v (system Jacobian applied to a
 * given vector). If a user-defined routine is not given, then the
 * private routine is used.
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

typedef int (*KINSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv,
                                     N_Vector uu, booleantype *new_uu, 
                                     void *J_data);

#ifdef __cplusplus
}
#endif

#endif
