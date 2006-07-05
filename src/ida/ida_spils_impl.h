/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:35 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the common header file (private version) for the Scaled
 * Preconditioned Iterative Linear Solver modules.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPILS_IMPL_H
#define _IDASPILS_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <ida/ida_spils.h>

/* Types of iterative linear solvers */

#define SPILS_SPGMR   1
#define SPILS_SPBCG   2
#define SPILS_SPTFQMR 3

/* Constants */

#define IDA_SPILS_MAXL    5
#define IDA_SPILS_MAXRS   5

/*
 * -----------------------------------------------------------------
 * Types : IDASpilsMemRec, IDASpilsMem                             
 * -----------------------------------------------------------------
 */

typedef struct {

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

  long int s_nst0;     /* nst0 = saved nst (for performance monitor)   */   
  long int s_nni0;     /* nni0 = saved nni (for performance monitor)   */   
  long int s_nli0;     /* nli0 = saved nli (for performance monitor)   */   
  long int s_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor) */   
  long int s_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor) */   
  long int s_nwarn;    /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector s_ytemp;    /* temp vector used by IDAAtimesDQ              */ 
  N_Vector s_yptemp;   /* temp vector used by IDAAtimesDQ              */ 
  N_Vector s_xx;       /* temp vector used by the solve function       */
  N_Vector s_ycur;     /* current y vector in Newton iteration         */
  N_Vector s_ypcur;    /* current yp vector in Newton iteration        */
  N_Vector s_rcur;     /* rcur = F(tn, ycur, ypcur)                    */

  IDASpilsPrecSetupFn s_pset;     /* pset = user-supplied routine      */
                                  /* to compute a preconditioner       */

  IDASpilsPrecSolveFn s_psolve;   /* psolve = user-supplied routine to */
                                  /* solve preconditioner linear system*/

  void *s_pdata;                  /* pdata passed to psolve and precond*/

  void *s_spils_mem;              /* memory used by the generic solver */

  IDASpilsJacTimesVecFn s_jtimes; /* Jacobian*vector routine           */ 

  void *s_jdata;                  /* data passed to Jtimes             */

  int s_last_flag;                /* last error return flag            */

} IDASpilsMemRec, *IDASpilsMem;


/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

/* Atimes and PSolve routines called by generic solver */

int IDASpilsAtimes(void *ida_mem, N_Vector v, N_Vector z);

int IDASpilsPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

int IDASpilsDQJtimes(realtype tt,
                     N_Vector yy, N_Vector yp, N_Vector rr,
                     N_Vector v, N_Vector Jv, 
                     realtype c_j, void *jac_data, 
                     N_Vector work1, N_Vector work2);



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

#define MSGS_PSET_FAILED "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."

/* Warning Messages */

#define MSGS_WARN  "Warning: " MSGS_TIME "poor iterative algorithm performance. "

#define MSGS_AVD_WARN  MSGS_WARN "Average number of linear iterations is " MSGS_FRMT
#define MSGS_CFN_WARN  MSGS_WARN "Nonlinear convergence failure rate is " MSGS_FRMT
#define MSGS_CFL_WARN  MSGS_WARN "Linear convergence failure rate is " MSGS_FRMT



#ifdef __cplusplus
}
#endif

#endif
