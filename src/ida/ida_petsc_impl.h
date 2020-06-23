/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmers: Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is implementation of the interface to PETSc libraries.
 * So far, only linear solvers are interfaced.
 * -----------------------------------------------------------------
 */

#ifndef _IDA_PETSC_IMPL_H
#define _IDA_PETSC_IMPL_H

#include <petscksp.h>
#include <ida/ida_petsc.h>
#include "ida_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Types : IDAPETScMemRec, IDAPETScMem
 *
 * TODO: Some of these variables are not used as this 
 * code is still under heavy development. Once we settle down on 
 * the features we want to support here, the redundant declarations 
 * will be removed. For now, they are left here as a reminder.                              
 * -----------------------------------------------------------------
 */

typedef struct IDAPETScMemRec {

  int s_type;          /* type of scaled preconditioned iterative LS   */

  int  s_gstype;       /* type of Gram-Schmidt orthogonalization       */
  realtype s_sqrtN;    /* sqrt(N)                                      */
  int  s_maxl;         /* maxl = maximum dimension of the Krylov space */
  int  s_maxrs;        /* maxrs = max. number of GMRES restarts        */
  realtype s_eplifac;  /* eplifac = linear convergence factor          */
  realtype s_dqincfac; /* dqincfac = optional increment factor in Jv   */
  realtype s_epslin;   /* SpgrmSolve tolerance parameter               */

  long int s_npe;      /* npe = total number of precond calls          */   
  long int s_nli;      /* nli = total number of linear iterations      */
  long int s_nps;      /* nps = total number of psolve calls           */
  long int s_ncfl;     /* ncfl = total number of convergence failures  */
  long int s_nres;     /* nres = total number of calls to res          */
  long int s_njtimes;  /* njtimes = total number of calls to jtimes    */
  long int s_nje;      /* nje = total number of Jacobian calls         */

  long int s_nst0;     /* nst0 = saved nst (for performance monitor)   */   
  long int s_nni0;     /* nni0 = saved nni (for performance monitor)   */   
  long int s_nli0;     /* nli0 = saved nli (for performance monitor)   */   
  long int s_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor) */   
  long int s_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor) */   
  long int s_nwarn;    /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector s_ytemp;    /* temp vector used only to compute sqrt(N)     */ 

  KSP *s_ksp_mem;      /* memory used by the KSP solver                */

  long int s_last_flag; /* last error return flag                      */
  
  void *s_udata;       /* pointer to user defined data                 */

  /* 
   * Jacobian evaluation function that stores Jacobian in PETSc matrix 
   * 
   */
  IDAPETScJacFn s_jaceval;

  /* 
   * Pointer to PETSc matrix storing Jacobian
   * 
   */
  Mat *JacMat;

} *IDAPETScMem;


/*
 * -----------------------------------------------------------------
 * Error and Warning Messages
 * -----------------------------------------------------------------
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGS_TIME "at t = %Lg, "
#define MSGS_FRMT "%Le."

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGS_TIME "at t = %lg, "
#define MSGS_FRMT "%le."

#else

#define MSGS_TIME "at t = %g, "
#define MSGS_FRMT "%e."

#endif


/* Error Messages */

#define MSGS_IDAMEM_NULL   "Integrator memory is NULL."
#define MSGS_MEM_FAIL      "A memory request failed."
#define MSGS_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE    "Incompatible linear solver type."
#define MSGS_LMEM_NULL     "Linear solver memory is NULL."
#define MSGS_BAD_GSTYPE    "gstype has an illegal value."
#define MSGS_NEG_MAXRS     "maxrs < 0 illegal."
#define MSGS_NEG_EPLIFAC   "eplifac < 0.0 illegal."
#define MSGS_NEG_DQINCFAC  "dqincfac < 0.0 illegal."

#define MSGS_PSET_FAILED  "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."
#define MSGS_JAC_FAILED    "The Jacobian setup routine failed in an unrecoverable manner."
#define MSGS_JAC_NULL      "The pointer to Jacobian setup routine is NULL."
#define MSGS_JAC_MAT_NULL  "The pointer to Jacobian matrix is NULL."

/* Warning Messages */

#define MSGS_WARN  "Warning: " MSGS_TIME "poor iterative algorithm performance. "

#define MSGS_AVD_WARN  MSGS_WARN "Average number of linear iterations is " MSGS_FRMT
#define MSGS_CFN_WARN  MSGS_WARN "Nonlinear convergence failure rate is " MSGS_FRMT
#define MSGS_CFL_WARN  MSGS_WARN "Linear convergence failure rate is " MSGS_FRMT



#ifdef __cplusplus
}
#endif

#endif /* _IDA_PETSC_IMPL_H */
