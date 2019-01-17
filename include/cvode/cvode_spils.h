/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * Header file for the deprecated Scaled, Preconditioned Iterative
 * Linear Solver interface in CVODE; these routines now just wrap 
 * the updated CVODE generic linear solver interface in cvode_ls.h.
 * -----------------------------------------------------------------*/

#ifndef _CVSPILS_H
#define _CVSPILS_H

#include <cvode/cvode_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  CVSPILS user-supplied function prototypes (typedefs for 
  equivalent types in cvode_ls.h)
  ===============================================================*/

typedef CVLsPrecSetupFn CVSpilsPrecSetupFn;
typedef CVLsPrecSolveFn CVSpilsPrecSolveFn;
typedef CVLsJacTimesSetupFn CVSpilsJacTimesSetupFn;
typedef CVLsJacTimesVecFn CVSpilsJacTimesVecFn;

/*=================================================================
  CVSPILS Exported functions (wrappers for equivalent routines in 
  cvode_ls.h)
  =================================================================*/

int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS);

int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac);
  
int CVSpilsSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset, CVSpilsPrecSolveFn psolve);

int CVSpilsSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup, CVSpilsJacTimesVecFn jtimes);
  
int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
  
int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals);
 
int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves);
  
int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters);
  
int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails);
  
int CVSpilsGetNumJTSetupEvals(void *cvode_mem, long int *njtsetups);
  
int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
  
int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);
  
int CVSpilsGetLastFlag(void *cvode_mem, long int *flag);
  
char *CVSpilsGetReturnFlagName(long int flag);
  

#ifdef __cplusplus
}
#endif

#endif
