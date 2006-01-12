/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-12 22:53:38 $
 * ----------------------------------------------------------------- 
 * Programmer(s):Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Wrappers for using the CVODES iterative linear solvers on adjoint
 * (backward) problems.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_spils_impl.h"
#include "cvodea_impl.h"

/* Readability replacements */

#define ytmp        (ca_mem->ca_ytmp)
#define getY        (ca_mem->ca_getY)
#define lmemB       (ca_mem->ca_lmemB)

#define pset_B      (cvspilsB_mem->s_psetB)
#define psolve_B    (cvspilsB_mem->s_psolveB)
#define jtimes_B    (cvspilsB_mem->s_jtimesB)
#define P_data_B    (cvspilsB_mem->s_P_dataB)
#define jac_data_B  (cvspilsB_mem->s_jac_dataB)

/*
 * CVAspilsPrecSetup
 *
 * This routine interfaces to the CVSpilsPrecSetupFnB routine 
 * provided by the user.
 * NOTE: p_data actually contains cvadj_mem
 */

int CVAspilsPrecSetup(realtype t, N_Vector yB, 
                      N_Vector fyB, booleantype jokB, 
                      booleantype *jcurPtrB, realtype gammaB,
                      void *cvadj_mem,
                      N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVadjMem ca_mem;
  CVSpilsMemB cvspilsB_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvspilsB_mem = (CVSpilsMemB) lmemB;

  /* Forward solution from interpolation */
  flag = getY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint precondB routine */
  flag = pset_B(t, ytmp, yB, fyB, jokB, jcurPtrB, gammaB,
                P_data_B, tmp1B, tmp2B, tmp3B);

  return(flag);
}


/*
 * CVAspilsPrecSolve
 *
 * This routine interfaces to the CVSpilsPrecSolveFnB routine 
 * provided by the user.
 * NOTE: p_data actually contains cvadj_mem
 */

int CVAspilsPrecSolve(realtype t, N_Vector yB, N_Vector fyB,
                      N_Vector rB, N_Vector zB,
                      realtype gammaB, realtype deltaB,
                      int lrB, void *cvadj_mem, N_Vector tmpB)
{
  CVadjMem ca_mem;
  CVSpilsMemB cvspilsB_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvspilsB_mem = (CVSpilsMemB) lmemB;

  /* Forward solution from interpolation */
  flag = getY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint psolveB routine */
  flag = psolve_B(t, ytmp, yB, fyB, rB, zB, gammaB, deltaB, 
                  lrB, P_data_B, tmpB);

  return(flag);
}


/*
 * CVAspilsJacTimesVec
 *
 * This routine interfaces to the CVSpilsJacTimesVecFnB routine 
 * provided by the user.
 * NOTE: jac_data actually contains cvadj_mem
 */

int CVAspilsJacTimesVec(N_Vector vB, N_Vector JvB, realtype t, 
                        N_Vector yB, N_Vector fyB, 
                        void *cvadj_mem, N_Vector tmpB)
{
  CVadjMem ca_mem;
  CVSpilsMemB cvspilsB_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvspilsB_mem = (CVSpilsMemB) lmemB;

  /* Forward solution from interpolation */
  flag = getY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint jtimesB routine */
  flag = jtimes_B(vB, JvB, t, ytmp, yB, fyB, jac_data_B, tmpB);

  return(flag);
}

