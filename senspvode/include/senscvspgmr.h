/******************************************************************
 *                                                                *
 * File          : senscvspgmr.h                                  *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Last Modified : 25 August 2000                                 *
 *----------------------------------------------------------------*
 * This is the header file for the CVODE scaled, preconditioned   *
 * GMRES linear solver, CVSPGMR, for sensitivity analysis.        *
 * The linear system size is Ny, where Ny is the number of ODEs   *
 * contained in y' = f(t,y,p).                                    *
 *                                                                *
 * This file is simply the cvspgmr.h header file, along with a    *
 * function prototype for SensCVSpgmr().                          *
 *                                                                *
 * Note: The type integer must be large enough to store the value *
 * of the linear system size Ny.                                  *
 *                                                                *
 ******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _senscvspgmr_h
#define _senscvspgmr_h


#include <stdio.h>
#include "cvode.h"
#include "spgmr.h"
#include "llnltyps.h"
#include "nvector.h"
#include "cvspgmr.h"

/******************************************************************
 *                                                                *
 * Function : SensCVSpgmr                                         *
 *----------------------------------------------------------------*
 * A call to the SensCVSpgmr function links the main CVODE        *
 * integrator with the CVSPGMR linear solver for sensitivity      *
 * analysis.                                                      *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *             SensCVodeMalloc.                                   *
 *                                                                *
 * pretype   is the type of user preconditioning to be done.      *
 *             This must be one of the four enumeration constants *
 *             NONE, LEFT, RIGHT, or BOTH defined in iterativ.h.  *
 *             These correspond to no preconditioning,            *
 *             left preconditioning only, right preconditioning   *
 *             only, and both left and right preconditioning,     *
 *             respectively.                                      *
 *                                                                *
 * gstype    is the type of Gram-Schmidt orthogonalization to be  *
 *           used. This must be one of the two enumeration        *
 *           constants MODIFIED_GS or CLASSICAL_GS defined in     *
 *           iterativ.h. These correspond to using modified       *
 *           Gram-Schmidt and classical Gram-Schmidt,             *
 *           respectively.                                        *
 *                                                                *
 * maxl      is the maximum Krylov dimension. This is an          *
 *             optional input to the CVSPGMR solver. Pass 0 to    *
 *             use the default value MIN(Ny, CVSPGMR_MAXL=5).     *
 *                                                                *
 * delt      is the factor by which the tolerance on the          *
 *             nonlinear iteration is multiplied to get a         *
 *             tolerance on the linear iteration. This is an      *
 *             optional input to the CVSPGMR solver. Pass 0 to    *
 *             use the default value CVSPGMR_DELT = 0.05.         *
 *                                                                *
 * precond   is the user's preconditioner routine. It is used to  *
 *             evaluate and preprocess any Jacobian-related data  *
 *             needed by the psolve routine.  See the             *
 *             documentation for the type CVSpgmrPrecondFn for    *
 *             full details.  Pass NULL if no such setup of       *
 *             Jacobian data is required.  A precond routine is   *
 *             NOT required for any of the four possible values   *
 *             of pretype.                                        *
 *                                                                *
 * psolve    is the user's preconditioner solve routine. It is    *
 *             used to solve Pz=r, where P is a preconditioner    *
 *             matrix.  See the documentation for the type        *
 *             CVSpgmrPSolveFn for full details.  The only case   *
 *             in which psolve is allowed to be NULL is when      *
 *             pretype is NONE.  A valid psolve function must be  *
 *             supplied when any preconditioning is to be done.   *
 *                                                                *
 * P_data    is a pointer to user preconditioner data. This       *
 *             pointer is passed to precond and psolve every time *
 *             these routines are called.                         *
 *                                                                *
 ******************************************************************/
  
void SensCVSpgmr(void *cvode_mem, int pretype, int gstype, 
		 int maxl, real delt, CVSpgmrPrecondFn precond,
		 CVSpgmrPSolveFn psolve, void *P_data);
 
#endif

#ifdef __cplusplus
}
#endif
