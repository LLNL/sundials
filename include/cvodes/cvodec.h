/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmers: Radu Serban @ LLNL                              
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
 * This is the interface file for CVODES derivative calculations
 * using the complex step method.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvodec_h
#define _cvodec_h

/*----------------------------------------------------------------*
 *                                                                *
 * Type : RhsCSFn                                                 *
 *----------------------------------------------------------------*
 *                                                                *
 * A RhsCSFn f_cs does not have a return value.                   *
 *                                                                *
 *----------------------------------------------------------------*/

typedef void (*RhsCSFn1)(realtype t_re, realtype t_im,
                         N_Vector y_re, N_Vector y_im,
                         N_Vector ydot_re, N_Vector ydot_im,
                         void *f_data);

typedef void (*RhsCSFn2)(realtype t_re, realtype t_im,
                         N_Vector y_re, N_Vector y_im,
                         N_Vector ydot_re, N_Vector ydot_im,
                         void *f_data_re, void *f_data_im);


/*----------------------------------------------------------------*
 *                                                                *
 * User-callable setup routines                                   *
 *----------------------------------------------------------------*
 *                                                                *
 * CVodeSetCSDerivs - enables complex step derivative calculation *
 *                                                                *
 *          f_cs is a user-defined function, of type RhCSFn which *
 *          should compute the ODE right hand side in complex     *
 *          arithmetic.
 *                                                                *
 *          f_data_im is a pointer to a structure of the same     *
 *          type as the f_data pointer passed to CVode through    *
 *          CVodeSetFdata. It is the user's responsibility to     *
 *          allocate and initialize to 0.0 all fields in f_data_im*
 *          that are flagged.                                     *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * CVodeSetCSStep - sets the perturbation del_cs used in the      *
 *                  complex step method (default = unit roundoff) *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * CVodeSetSensCSRhs - enables complex step calculation of the    *
 *                     sensitivity right hand sides.              *
 *                                                                *
 *           p_im is a pointer in f_data_im data structure which  *
 *           corresponds to the field p in f_data. p_im must be   *
 *           allocated and initialized to 0.0 by the user.        *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * CVDenseSetCSJac - enables complex step calculation of the      *
 *                   dense Jacobian for use with CVDENSE.         *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * CVBandSetCSJac - enables complex step calculation of the       *
 *                  band Jacobian approximation for use with      *
 *                  CVBAND.                                       *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * CVSpgmrSetCSJacTimesVec - enables complex step calculation of  *
 *                  the Jacobian x Vector for use with CVSPGMR.   *
 *                                                                *
 *----------------------------------------------------------------*/

SUNDIALS_EXPORT int CVodeSetCSDerivs(void *cvode_mem, void *f_cs, void *f_data_im);
SUNDIALS_EXPORT int CVodeSetCSStep(void *cvode_mem, realtype del_cs);
SUNDIALS_EXPORT int CVodeSetSensCSRhs(void *cvode_mem, realtype *p_im);
SUNDIALS_EXPORT int CVDenseSetCSJac(void *cvode_mem);
SUNDIALS_EXPORT int CVBandSetCSJac(void *cvode_mem);
SUNDIALS_EXPORT int CVSpgmrSetCSJacTimesVec(void *cvode_mem);

/* Return values */
/* SUCCESS */
enum {CVCS_NO_MEM=-1,   CVCS_MEM_FAIL=-2, CVCS_ILL_INPUT=-3,
      CVCS_NO_CSMEM=-4, CVCS_NO_LMEM=-5                     };

/*----------------------------------------------------------------*
 *                                                                *
 * Complex step approximation routines                            *
 *----------------------------------------------------------------*
 *                                                                *
 * CVSensRhsCS     - computes sensitivity right hand side         *
 * CVDenseCSJac    - computes the dense Jacobian                  *
 * CVBandCSJac     - computes the band Jacobian approximation     *
 * CVSpgmrCSJtimes - computes the Jacobian x Vector product       *
 *                                                                *
 *----------------------------------------------------------------*/

SUNDIALS_EXPORT void CVSensRhsCS(int Ns, realtype t,
                                 N_Vector y, N_Vector ydot,
                                 int iS, N_Vector yS, N_Vector ySdot,
                                 void *fS_data,
                                 N_Vector tmp1, N_Vector tmp2);

SUNDIALS_EXPORT void CVDenseCSJac(sunindextype N, DenseMat J, realtype t,
                                  N_Vector y, N_Vector fy, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

SUNDIALS_EXPORT void CVBandCSJac(sunindextype N, sunindextype mupper,
                                 sunindextype mlower, BandMat J, realtype t,
                                 N_Vector y, N_Vector fy, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

SUNDIALS_EXPORT int CVSpgmrCSJtimes(N_Vector v, N_Vector Jv, realtype t,
                                    N_Vector y, N_Vector fy,
                                    void *jac_data, N_Vector work);

/*----------------------------------------------------------------*
 *                                                                *
 * Types : struct CVCSMemRec, CVCSMem                             *
 *----------------------------------------------------------------*
 *                                                                *
 *                                                                *
 *----------------------------------------------------------------*/

typedef struct {

  booleantype cvcs_type1;

  RhsCSFn1 cvcs_f_cs1;
  RhsCSFn2 cvcs_f_cs2;

  void *cvcs_f_data_im;

  realtype cvcs_del_cs;

  realtype *cvcs_p_im;

} CVCSMemRec, *CVCSMem;

#endif

#ifdef __cplusplus
}
#endif
