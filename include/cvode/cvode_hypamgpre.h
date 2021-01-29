/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Michael Wittman, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CVHYPRE_BOOMERAMG module, for a
 * interface to hypre's BoomerAMG preconditioner, for use with CVODE, a CVSPILS linear
 * solver, and the ParHyp implementation of NVECTOR.
 *
 * Summary:
 *
 *
 * The user's calling program should have the following form:
 *
 *   #include <nvector_parhyp.h>
 *   #include <cvode/cvode_hypamgpre.h>
 *   ...
 *   void *cvode_mem;
 *   ...
 *   Set y0
 *   ...
 *   cvode_mem = CVodeCreate(...);
 *   ier = CVodeInit(...);
 *   ...
 *   flag = CVSpgmr(cvode_mem, pretype, maxl);
 *      -or-
 *   flag = CVSpbcg(cvode_mem, pretype, maxl);
 *      -or-
 *   flag = CVSptfqmr(cvode_mem, pretype, maxl);
 *   ...
 *   flag = CVBoomerAMGInit(cvode_mem, );
 *   ...
 *   ier = CVode(...);
 *   ...
 *   CVodeFree(&cvode_mem);
 * 
 *   Free y0
 *
 * The user-supplied routines required are:
 *
 *   f    = function defining the ODE right-hand side f(t,y).
 *
 *   gloc = function defining the approximation g(t,y).
 *
 *   cfn  = function to perform communication need for gloc.
 *
 * Notes:
 *
 * 1) This header file is included by the user for the definition
 *    of the CVBBDData type and for needed function prototypes.
 *
 * 2) The CVBoomerAMGInit call includes half-bandwiths mudq and mldq
 *    to be used in the difference quotient calculation of the
 *    approximate Jacobian. They need not be the true
 *    half-bandwidths of the Jacobian of the local block of g,
 *    when smaller values may provide a greater efficiency.
 *    Also, the half-bandwidths mukeep and mlkeep of the retained
 *    banded approximate Jacobian block may be even smaller,
 *    to reduce storage and computation costs further.
 *    For all four half-bandwidths, the values need not be the
 *    same on every processor.
 *
 * 3) The actual name of the user's f function is passed to
 *    CVodeInit, and the names of the user's gloc and cfn
 *    functions are passed to CVBoomerAMGInit.
 *
 * 4) The pointer to the user-defined data block user_data, which is
 *    set through CVodeSetUserData is also available to the user in
 *    gloc and cfn.
 *
 * 5) Optional outputs specific to this module are available by
 *    way of routines listed below. These include work space sizes
 *    and the cumulative number of gloc calls. The costs
 *    associated with this module also include nsetups banded LU
 *    factorizations, nlinsetups cfn calls, and npsolves banded
 *    backsolve calls, where nlinsetups and npsolves are
 *    integrator/CVSPGMR/CVSPBCG/CVSPTFQMR optional outputs.
 * -----------------------------------------------------------------
 */

#ifndef _CVHYPRE_BOOMERAMG_H
#define _CVHYPRE_BOOMERAMG_H

#include <sundials/sundials_nvector.h>
#include <mpi.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"
#include <nvector/nvector_parhyp.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Type: CVParCsrJacFn
 * -----------------------------------------------------------------
 *
 * A dense Jacobian approximation function Jac must be of type 
 * CVDlsDenseJacFn. Its parameters are:
 *
 * N   is the problem size.
 *
 * Jac is the dense matrix (of type DlsMat) that will be loaded
 *     by a CVDlsDenseJacFn with an approximation to the Jacobian 
 *     matrix J = (df_i/dy_j) at the point (t,y). 
 *
 * t   is the current value of the independent variable.
 *
 * y   is the current value of the dependent variable vector,
 *     namely the predicted value of y(t).
 *
 * fy  is the vector f(t,y).
 *
 * user_data is a pointer to user data - the same as the user_data
 *     parameter passed to CVodeSetFdata.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 * vectors of length N which can be used by a CVDlsDenseJacFn
 * as temporary storage or work space.
 *
 * A CVParCsrJacFn should return 0 if successful, a positive 
 * value if a recoverable error occurred, and a negative value if 
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
 *     (see cvode.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*CVParCsrJacFn)(sunindextype N, sunindextype ilower, sunindextype iupper, sunindextype jlower, sunindextype jupper, realtype gamma, realtype t,
			       N_Vector y, N_Vector fy, 
			       HYPRE_IJMatrix* Jac, void *user_data,
			       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
typedef int (*CVJacobianIJUpdateFn)(HYPRE_IJMatrix* pointer_A, int ilower, int iupper, int jlower, int jupper, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Function : CVBoomerAMGInit
 * -----------------------------------------------------------------
 * CVBoomerAMGInit allocates and initializes the BBD preconditioner.
 *
 * The parameters of CVBoomerAMGInit are as follows:
 *
 * cvode_mem is the pointer to the integrator memory.
 *
 * Nlocal is the length of the local block of the vectors y etc.
 *        on the current processor.
 *
 * mudq, mldq are the upper and lower half-bandwidths to be used
 *            in the difference quotient computation of the local
 *            Jacobian block.
 *
 * mukeep, mlkeep are the upper and lower half-bandwidths of the
 *                retained banded approximation to the local Jacobian
 *                block.
 *
 * dqrely is an optional input. It is the relative increment
 *        in components of y used in the difference quotient
 *        approximations. To specify the default, pass 0.
 *        The default is dqrely = sqrt(unit roundoff).
 *
 * gloc is the name of the user-supplied function g(t,y) that
 *      approximates f and whose local Jacobian blocks are
 *      to form the preconditioner.
 *
 * cfn is the name of the user-defined function that performs
 *     necessary interprocess communication for the
 *     execution of gloc.
 *
 * The return value of CVBoomerAMGInit is one of:
 *   CVSPILS_SUCCESS if no errors occurred
 *   CVSPILS_MEM_NULL if the integrator memory is NULL
 *   CVSPILS_LMEM_NULL if the linear solver memory is NULL
 *   CVSPILS_ILL_INPUT if an input has an illegal value
 *   CVSPILS_MEM_FAIL if a memory allocation request failed
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBoomerAMGInit(void *cvode_mem, int ilower, int iupper, int jlower, int jupper, int N);

/*
SUNDIALS_EXPORT int CVBoomerAMGSetup(realtype t, N_Vector y, N_Vector fy, 
                          booleantype jok, booleantype *jcurPtr, 
                          realtype gamma, void *bbd_data, 
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

SUNDIALS_EXPORT int CVBoomerAMGSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *bbd_data, N_Vector tmp);*/

/*
 * -----------------------------------------------------------------
 * Function : CVBoomerAMGReInit
 * -----------------------------------------------------------------
 * CVBoomerAMGReInit re-initializes the HYPRE_BOOMERAMG module when solving a
 * sequence of problems of the same size with CVSPGMR/CVHYPRE_BOOMERAMG or
 * CVSPBCG/CVHYPRE_BOOMERAMG or CVSPTFQMR/CVHYPRE_BOOMERAMG provided there is no change 
 * in Nlocal, mukeep, or mlkeep. After solving one problem, and after 
 * calling CVodeReInit to re-initialize the integrator for a subsequent 
 * problem, call CVBoomerAMGReInit.
 *
 * All arguments have the same names and meanings as those
 * of CVBoomerAMGInit.
 *
 * The return value of CVBoomerAMGReInit is one of:
 *   CVSPILS_SUCCESS if no errors occurred
 *   CVSPILS_MEM_NULL if the integrator memory is NULL
 *   CVSPILS_LMEM_NULL if the linear solver memory is NULL
 *   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBoomerAMGReInit(void *cvode_mem, sunindextype mudq, sunindextype mldq,
				    realtype dqrely);

/*
 * -----------------------------------------------------------------
 * HYPRE_BOOMERAMG optional output extraction routines
 * -----------------------------------------------------------------
 * CVBoomerAMGGetWorkSpace returns the HYPRE_BOOMERAMG real and integer work space
 *                       sizes.
 * CVBoomerAMGGetNumGfnEvals returns the number of calls to gfn.
 *
 * The return value of CVBoomerAMGGet* is one of:
 *   CVSPILS_SUCCESS if no errors occurred
 *   CVSPILS_MEM_NULL if the integrator memory is NULL
 *   CVSPILS_LMEM_NULL if the linear solver memory is NULL
 *   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBoomerAMGGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CVBoomerAMGGetNumGfnEvals(void *cvode_mem, long int *ngevalsBBDP);

#ifdef __cplusplus
}
#endif

#endif
