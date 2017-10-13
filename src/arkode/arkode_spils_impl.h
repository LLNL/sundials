/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 *---------------------------------------------------------------
 * Implementation header file for the scaled, preconditioned
 * linear solver interface.
 *--------------------------------------------------------------*/

#ifndef _ARKSPILS_IMPL_H
#define _ARKSPILS_IMPL_H

#include <arkode/arkode_spils.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 Types: ARKSpilsMemRec, ARKSpilsMem

 The type ARKSpilsMem is pointer to a ARKSpilsMemRec.
---------------------------------------------------------------*/
typedef struct ARKSpilsMemRec {

  realtype sqrtN;     /* sqrt(N)                                      */
  realtype eplifac;   /* eplifac = user specified or EPLIN_DEFAULT    */
  realtype deltar;    /* deltar = delt * LTE                          */
  realtype delta;     /* delta = deltar * sqrtN                       */

  booleantype jbad;   /* heuristic suggestion for pset/JTimes         */
  long int nstlpre;   /* value of nst at the last pset call           */
  long int npe;       /* npe = total number of pset calls             */
  long int nli;       /* nli = total number of linear iterations      */
  long int nps;       /* nps = total number of psolve calls           */
  long int ncfl;      /* ncfl = total number of convergence failures  */
  long int njtsetup;  /* njtsetup = total number of calls to jtsetup  */
  long int njtimes;   /* njtimes = total number of calls to jtimes    */
  long int nfes;      /* nfeSG = total number of calls to f for     
                         difference quotient Jacobian-vector products */

  SUNLinearSolver LS; /* generic iterative linear solver object       */
  
  N_Vector ytemp;     /* temp vector passed to jtimes and psolve      */
  N_Vector x;         /* solution vector used by SUNLinearSolver      */
  N_Vector ycur;      /* ARKODE current y vector in Newton Iteration  */
  N_Vector fcur;      /* fcur = f(tn, ycur)                           */

  /* Preconditioner computation
    (a) user-provided:
        - P_data == user_data
        - pfree == NULL (the user dealocates memory for user_data)
    (b) internal preconditioner module
        - P_data == arkode_mem
        - pfree == set by the prec. module and called in ARKodeFree  */
  ARKSpilsPrecSetupFn pset;
  ARKSpilsPrecSolveFn psolve;
  int (*pfree)(ARKodeMem ark_mem);
  void *P_data;

  /* Jacobian times vector compuation
    (a) jtimes function provided by the user:
        - j_data == user_data
        - jtimesDQ == SUNFALSE
    (b) internal jtimes
        - j_data == arkode_mem
        - jtimesDQ == SUNTRUE   */
  booleantype jtimesDQ;
  ARKSpilsJacTimesSetupFn jtsetup;
  ARKSpilsJacTimesVecFn jtimes;
  void *j_data;

  long int last_flag; /* last error flag returned by any function */

} *ARKSpilsMem;


/*---------------------------------------------------------------
 Types: ARKSpilsMassMemRec, ARKSpilsMassMem

 The type ARKSpilsMassMem is pointer to a ARKSpilsMassMemRec.
---------------------------------------------------------------*/
typedef struct ARKSpilsMassMemRec {

  realtype sqrtN;     /* sqrt(N)                                      */
  realtype eplifac;   /* eplifac = user specified or EPLIN_DEFAULT    */
  realtype deltar;    /* deltar = delt * LTE                          */
  realtype delta;     /* delta = deltar * sqrtN                       */

  long int npe;       /* npe = total number of pset calls             */
  long int nli;       /* nli = total number of linear iterations      */
  long int nps;       /* nps = total number of psolve calls           */
  long int ncfl;      /* ncfl = total number of convergence failures  */
  long int nmtsetup;  /* nmtsetup = total number of calls to mtsetup  */
  long int nmtimes;   /* nmtimes = total number of calls to mtimes    */

  SUNLinearSolver LS; /* generic iterative linear solver object       */

  booleantype time_dependent;  /* flag stating whether M depends on t */
  
  N_Vector x;         /* solution vector used by SUNLinearSolver      */
  N_Vector ycur;      /* ARKODE current y vector                      */

  /* Preconditioner computation
    (a) user-provided:
        - P_data == user_data
        - pfree == NULL (the user dealocates memory for user_data)
    (b) internal preconditioner module
        - P_data == arkode_mem
        - pfree == set by the prec. module and called in ARKodeFree */
  ARKSpilsMassPrecSetupFn pset;
  ARKSpilsMassPrecSolveFn psolve;
  int (*pfree)(ARKodeMem ark_mem);
  void *P_data;

  /* Mass matrix times vector setup and product routines, data */
  ARKSpilsMassTimesSetupFn mtsetup;
  ARKSpilsMassTimesVecFn mtimes;
  void *mt_data;
  
  long int last_flag; /* last error flag returned by any function     */

} *ARKSpilsMassMem;


/*---------------------------------------------------------------
 Prototypes of internal functions
---------------------------------------------------------------*/

/* Interface routines called by system SUNLinearSolver */
int ARKSpilsATimes(void *ark_mem, N_Vector v, N_Vector z);
int ARKSpilsPSetup(void *ark_mem);
int ARKSpilsPSolve(void *ark_mem, N_Vector r, N_Vector z,
                   realtype tol, int lr);

/* Interface routines called by mass SUNLinearSolver */
int ARKSpilsMTimes(void *ark_mem, N_Vector v, N_Vector z);
int ARKSpilsMPSetup(void *ark_mem);
int ARKSpilsMPSolve(void *ark_mem, N_Vector r, N_Vector z,
                    realtype tol, int lr);

/* Difference quotient approximation for Jac times vector */
int ARKSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                     N_Vector y, N_Vector fy, void *data,
                     N_Vector work);

/* Generic linit/lsetup/lsolve/lfree interface routines for ARKode to call */
int arkSpilsInitialize(ARKodeMem ark_mem);

int arkSpilsSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                  N_Vector fpred, booleantype *jcurPtr, 
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 

int arkSpilsSolve(ARKodeMem ark_mem, N_Vector b, N_Vector ycur, N_Vector fcur);

int arkSpilsFree(ARKodeMem ark_mem);

/* Generic minit/msetup/mmult/msolve/mfree routines for ARKode to call */  
int arkSpilsMassInitialize(ARKodeMem ark_mem);
  
int arkSpilsMassSetup(ARKodeMem ark_mem, N_Vector vtemp1,
                      N_Vector vtemp2, N_Vector vtemp3); 

int arkSpilsMassMult(ARKodeMem ark_mem, N_Vector v, N_Vector Mv);

int arkSpilsMassSolve(ARKodeMem ark_mem, N_Vector b);

int arkSpilsMassFree(ARKodeMem ark_mem);

/* Auxilliary functions */
int arkSpilsInitializeCounters(ARKSpilsMem arkspils_mem);

int arkSpilsInitializeMassCounters(ARKSpilsMassMem arkspils_mem);


/*---------------------------------------------------------------
 Error Messages -- REMOVE SOME???
---------------------------------------------------------------*/
#define MSGS_ARKMEM_NULL   "Integrator memory is NULL."
#define MSGS_MEM_FAIL      "A memory request failed."
#define MSGS_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE    "Incompatible linear solver type."
#define MSGS_LMEM_NULL     "Linear solver memory is NULL."
#define MSGS_MASSMEM_NULL  "Mass matrix solver memory is NULL."
#define MSGS_BAD_EPLIN     "eplifac < 0 illegal."

#define MSGS_PSET_FAILED   "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTSETUP_FAILED "The Jacobian x vector setup routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."
#define MSGS_MTSETUP_FAILED "The mass matrix x vector setup routine failed in an unrecoverable manner."
#define MSGS_MTIMES_FAILED "The mass matrix x vector routine failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
