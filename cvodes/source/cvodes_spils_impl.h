/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-12 22:53:38 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the iterative linear solvers CVSP*
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPILS_IMPL_H
#define _CVSSPILS_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_spils.h"

  /*
   * ------------------------------------------------
   * Wrapper functions for using the iterative linear 
   * solvers on adjoint (backward) problems
   * ------------------------------------------------
   */

  /* 
   * CVAspilsPrecSetup has type CVSpilsPrecSetupFn
   * It wraps around the user-provided function of type CVSpilsPrecSetupFnB
   */

  int CVAspilsPrecSetup(realtype t, N_Vector yB, 
                        N_Vector fyB, booleantype jokB, 
                        booleantype *jcurPtrB, realtype gammaB,
                        void *cvadj_mem,
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

  /* 
   * CVAspilsPrecSolve has type CVSpilsPrecSolveFn 
   * It wraps around the user-provided function of type CVSpilsPrecSolveFnB
   */

  int CVAspilsPrecSolve(realtype t, N_Vector yB, N_Vector fyB,
                        N_Vector rB, N_Vector zB,
                        realtype gammaB, realtype deltaB,
                        int lrB, void *cvadj_mem, N_Vector tmpB);
  
  /* 
   * CVAspilsJacTimesVec has type CVSpilsJacTimesVecFn 
   * It wraps around the user-provided function of type CVSpilsJacTimesVecFnB
   */

  int CVAspilsJacTimesVec(N_Vector vB, N_Vector JvB, realtype t, 
                          N_Vector yB, N_Vector fyB, 
                          void *cvadj_mem, N_Vector tmpB);

  /*
   * -----------------------------------------------------------------
   * Types : CVSpilsMemRecB, CVSpilsMemB       
   * -----------------------------------------------------------------
   * CVSpgmrB, CVSpbcgB, and CVSptfqmr attach such a structure to the 
   * lmemB filed of CVadjMem
   * -----------------------------------------------------------------
   */

  typedef struct {

    CVSpilsJacTimesVecFnB s_jtimesB;
    CVSpilsPrecSetupFnB s_psetB;
    CVSpilsPrecSolveFnB s_psolveB;
    void *s_P_dataB;
    void *s_jac_dataB;

  } CVSpilsMemRecB, *CVSpilsMemB;


#ifdef __cplusplus
}
#endif

#endif
