/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * Common implementation header file for the scaled, preconditioned
 * linear solver modules.
 *--------------------------------------------------------------*/

#ifndef _ARKSPILS_IMPL_H
#define _ARKSPILS_IMPL_H

#include <arkode/arkode_spils.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Types of iterative linear solvers */
#define SPILS_SPGMR   1
#define SPILS_SPBCG   2
#define SPILS_SPTFQMR 3
#define SPILS_PCG     4
#define SPILS_SPFGMR  5


/*---------------------------------------------------------------
 Types: ARKSpilsMemRec, ARKSpilsMem

 The type ARKSpilsMem is pointer to a ARKSpilsMemRec.
---------------------------------------------------------------*/
typedef struct ARKSpilsMemRec {

  int s_type;           /* type of scaled preconditioned iterative LS   */

  int  s_pretype;       /* type of preconditioning                      */
  int  s_gstype;        /* type of Gram-Schmidt orthogonalization       */
  realtype s_sqrtN;     /* sqrt(N)                                      */
  realtype s_eplifac;   /* eplifac = user specified or EPLIN_DEFAULT    */
  realtype s_deltar;    /* deltar = delt * LTE                          */
  realtype s_delta;     /* delta = deltar * sqrtN                       */
  int  s_maxl;          /* maxl = maximum dimension of the Krylov space */

  long int s_nstlpre;   /* value of nst at the last pset call           */
  long int s_npe;       /* npe = total number of pset calls             */
  long int s_nli;       /* nli = total number of linear iterations      */
  long int s_nps;       /* nps = total number of psolve calls           */
  long int s_ncfl;      /* ncfl = total number of convergence failures  */
  long int s_njtimes;   /* njtimes = total number of calls to jtimes    */
  long int s_nfes;      /* nfeSG = total number of calls to f for     
                           difference quotient Jacobian-vector products */

  N_Vector s_ytemp;     /* temp vector passed to jtimes and psolve      */
  N_Vector s_x;         /* temp vector used by ARKSpilsSolve            */
  N_Vector s_ycur;      /* ARKODE current y vector in Newton Iteration  */
  N_Vector s_fcur;      /* fcur = f(tn, ycur)                           */

  void* s_spils_mem;    /* memory used by the generic solver            */

  /* Preconditioner computation
    (a) user-provided:
        - P_data == user_data
        - pfree == NULL (the user dealocates memory for user_data)
    (b) internal preconditioner module
        - P_data == arkode_mem
        - pfree == set by the prec. module and called in ARKodeFree  */
  ARKSpilsPrecSetupFn s_pset;
  ARKSpilsPrecSolveFn s_psolve;
  void (*s_pfree)(ARKodeMem ark_mem);
  void *s_P_data;

  /* Jacobian times vector compuation
    (a) jtimes function provided by the user:
        - j_data == user_data
        - jtimesDQ == FALSE
    (b) internal jtimes
        - j_data == arkode_mem
        - jtimesDQ == TRUE   */
  booleantype s_jtimesDQ;
  ARKSpilsJacTimesVecFn s_jtimes;
  void *s_j_data;

  long int s_last_flag; /* last error flag returned by any function     */

} *ARKSpilsMem;


/*---------------------------------------------------------------
 Types: ARKSpilsMassMemRec, ARKSpilsMassMem

 The type ARKSpilsMassMem is pointer to a ARKSpilsMassMemRec.
---------------------------------------------------------------*/
typedef struct ARKSpilsMassMemRec {

  int s_type;           /* type of scaled preconditioned iterative LS   */

  int  s_pretype;       /* type of preconditioning                      */
  int  s_gstype;        /* type of Gram-Schmidt orthogonalization       */
  realtype s_sqrtN;     /* sqrt(N)                                      */
  realtype s_eplifac;   /* eplifac = user specified or EPLIN_DEFAULT    */
  realtype s_deltar;    /* deltar = delt * LTE                          */
  realtype s_delta;     /* delta = deltar * sqrtN                       */
  int  s_maxl;          /* maxl = maximum dimension of the Krylov space */

  long int s_npe;       /* npe = total number of pset calls             */
  long int s_nli;       /* nli = total number of linear iterations      */
  long int s_nps;       /* nps = total number of psolve calls           */
  long int s_ncfl;      /* ncfl = total number of convergence failures  */

  N_Vector s_ytemp;     /* temp vector passed to mtimes and psolve      */
  N_Vector s_x;         /* temp vector used by ARKSpilsSolve            */
  N_Vector s_ycur;      /* ARKODE current y vector                      */

  void* s_spils_mem;    /* memory used by the generic solver            */

  /* Preconditioner computation
    (a) user-provided:
        - P_data == user_data
        - pfree == NULL (the user dealocates memory for user_data)
    (b) internal preconditioner module
        - P_data == arkode_mem
        - pfree == set by the prec. module and called in ARKodeFree */
  ARKSpilsMassPrecSetupFn s_pset;
  ARKSpilsMassPrecSolveFn s_psolve;
  void (*s_pfree)(ARKodeMem ark_mem);
  void *s_P_data;

  long int s_last_flag; /* last error flag returned by any function     */

} *ARKSpilsMassMem;


/*---------------------------------------------------------------
 Prototypes of internal functions
---------------------------------------------------------------*/

/* Atimes and PSolve routines called by generic solver */
int ARKSpilsAtimes(void *ark_mem, N_Vector v, N_Vector z);
int ARKSpilsPSolve(void *ark_mem, N_Vector r, N_Vector z, int lr);

/* Mtimes and MPSolve routines called by mass matrix solver */
int ARKSpilsMtimes(void *ark_mem, N_Vector v, N_Vector z);
int ARKSpilsMPSolve(void *ark_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */
int ARKSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
		     N_Vector y, N_Vector fy, void *data,
		     N_Vector work);


/*---------------------------------------------------------------
 Error Messages
---------------------------------------------------------------*/
#define MSGS_ARKMEM_NULL   "Integrator memory is NULL."
#define MSGS_MEM_FAIL      "A memory request failed."
#define MSGS_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE    "Incompatible linear solver type."
#define MSGS_BAD_PRETYPE   "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."
#define MSGS_PSOLVE_REQ    "pretype != PREC_NONE, but PSOLVE = NULL is illegal."
#define MSGS_LMEM_NULL     "Linear solver memory is NULL."
#define MSGS_MASSMEM_NULL  "Mass matrix solver memory is NULL."
#define MSGS_BAD_GSTYPE    "Illegal value for gstype. Legal values are MODIFIED_GS and CLASSICAL_GS."
#define MSGS_BAD_EPLIN     "eplifac < 0 illegal."

#define MSGS_PSET_FAILED   "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."
#define MSGS_MTIMES_FAILED "The mass matrix * vector routine failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
