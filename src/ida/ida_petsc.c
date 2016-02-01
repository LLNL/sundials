/*
 * -----------------------------------------------------------------
 * $Revision:  $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmers: Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the IDA interface to linear
 * solvers from PETSc library.            
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <petscksp.h>

#include <ida/ida_petsc_ksp.h>
#include "ida_petsc_impl.h"
#include "ida_impl.h"

#include <sundials/sundials_petsc_ksp.h>
#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define PT9          RCONST(0.9)
#define PT05         RCONST(0.05)

/* IDAKSP linit, lsetup, lsolve, lperf, and lfree routines */

static int IDAKSPInit(IDAMem IDA_mem);

static int IDAKSPSetup(IDAMem IDA_mem, 
                         N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int IDAKSPSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector yy_now, N_Vector yp_now, N_Vector rr_now);

static int IDAKSPPerf(IDAMem IDA_mem, int perftask);

static int IDAKSPFree(IDAMem IDA_mem);


/* Readability Replacements */

#define nst          (IDA_mem->ida_nst)
#define tn           (IDA_mem->ida_tn)
#define cj           (IDA_mem->ida_cj)
#define epsNewt      (IDA_mem->ida_epsNewt)
#define res          (IDA_mem->ida_res)
#define user_data    (IDA_mem->ida_user_data)
#define ewt          (IDA_mem->ida_ewt)
#define errfp        (IDA_mem->ida_errfp)
// #define linit        (IDA_mem->ida_linit)
// #define lsetup       (IDA_mem->ida_lsetup)
// #define lsolve       (IDA_mem->ida_lsolve)
// #define lperf        (IDA_mem->ida_lperf)
// #define lfree        (IDA_mem->ida_lfree)
// #define lmem         (IDA_mem->ida_lmem)
#define nni          (IDA_mem->ida_nni)
#define ncfn         (IDA_mem->ida_ncfn)
//#define setupNonNull (IDA_mem->ida_setupNonNull)
//#define vec_tmpl     (IDA_mem->ida_tempv1)

//#define sqrtN     (idapetsc_mem->s_sqrtN)
#define epslin    (idapetsc_mem->s_epslin)
#define ytemp     (idapetsc_mem->s_ytemp)
#define yptemp    (idapetsc_mem->s_yptemp)
#define xx        (idapetsc_mem->s_xx)
#define ycur      (idapetsc_mem->s_ycur)
#define ypcur     (idapetsc_mem->s_ypcur)
#define rcur      (idapetsc_mem->s_rcur)
// #define npe       (idapetsc_mem->s_npe)
// #define nli       (idapetsc_mem->s_nli)
// #define nps       (idapetsc_mem->s_nps)
// #define ncfl      (idapetsc_mem->s_ncfl)
#define nst0      (idapetsc_mem->s_nst0)
#define nni0      (idapetsc_mem->s_nni0)
#define nli0      (idapetsc_mem->s_nli0)
#define ncfn0     (idapetsc_mem->s_ncfn0)
#define ncfl0     (idapetsc_mem->s_ncfl0)
#define nwarn     (idapetsc_mem->s_nwarn)
// #define njtimes   (idapetsc_mem->s_njtimes)
// #define nres      (idapetsc_mem->s_nres)
//#define spils_mem (idapetsc_mem->s_spils_mem)

#define jtimesDQ  (idapetsc_mem->s_jtimesDQ)
#define jtimes    (idapetsc_mem->s_jtimes)
#define jdata     (idapetsc_mem->s_jdata)

// #define last_flag (idapetsc_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * IDAKSP
 * -----------------------------------------------------------------
 *
 * This routine initializes the memory record and sets various function
 * fields specific to the IDAKSP linear solver module.  
 *
 * IDAKSP first calls the existing lfree routine if this is not NULL.
 * It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDAKSPInit, IDAKSPSetup,
 * IDAKSPSolve, IDAKSPPerf, and IDAKSPFree, respectively.
 * It allocates memory for a structure of type IDAPETScMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem).  It then various fields in the
 * IDAPETScMemRec structure. Finally, IDAKSP allocates memory for 
 * ytemp, yptemp, and xx, and calls KSPMalloc to allocate memory
 * for the KSP solver.
 *
 * The return value of IDAKSP is:
 *   IDASPILS_SUCCESS   =  0  if successful
 *   IDASPILS_MEM_FAIL  = -1 if IDA_mem is NULL or a memory allocation failed
 *   IDASPILS_ILL_INPUT = -2 if the gstype argument is illegal.
 *
 * -----------------------------------------------------------------
 */

int IDAKSP(void *ida_mem, int maxl, MPI_Comm comm)
{
  IDAMem IDA_mem;
  IDAPETScMem idapetsc_mem;
  // KSPMem spgmr_mem;
  int flag, maxl1;
  KSP *solver;
  PetscErrorCode ierr;
  
  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDAKSP", "IDAKSP", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if(IDA_mem->ida_tempv1->ops->nvdotprod == NULL) {
    IDAProcessError(NULL, IDASPILS_ILL_INPUT, "IDAKSP", "IDAKSP", MSGS_BAD_NVECTOR);
    return(IDASPILS_ILL_INPUT);
  }

  if (IDA_mem->ida_lfree != NULL) 
    flag = IDA_mem->ida_lfree((IDAMem) ida_mem);

  /* Set five main function fields in ida_mem */
  IDA_mem->ida_linit  = IDAKSPInit;
  IDA_mem->ida_lsetup = IDAKSPSetup;
  IDA_mem->ida_lsolve = IDAKSPSolve;
  IDA_mem->ida_lperf  = IDAKSPPerf;
  IDA_mem->ida_lfree  = IDAKSPFree;

  /* Allocate memory for KSP solver */
  solver = NULL;
  solver = (KSP*) malloc(sizeof(KSP));
  if (solver == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDAKSP", "IDAKSP", MSGS_MEM_FAIL);
    return(IDASPILS_MEM_FAIL);
  }

  /* Create KSP solver */
  ierr = KSPCreate(comm, solver);
  CHKERRQ(ierr);

  /* Set ILS type */
  //idapetsc_mem->s_type = SPILS_KSP;

  /* Set KSP parameters that were passed in call sequence */
  //maxl1 = (maxl <= 0) ? IDA_SPILS_MAXL : maxl;
  //idapetsc_mem->s_maxl     = maxl1;

  /* Set defaults for Jacobian-related fileds */
  //jtimesDQ = TRUE;
  //jtimes   = NULL;
  //jdata    = NULL;

  /* Set defaults for preconditioner-related fields */
  //idapetsc_mem->s_pset   = NULL;
  //idapetsc_mem->s_psolve = NULL;
  //idapetsc_mem->s_pfree  = NULL;
  //idapetsc_mem->s_pdata  = IDA_mem->ida_user_data;

  /* Set default values for the rest of the KSP parameters */
  //idapetsc_mem->s_gstype   = MODIFIED_GS;
  //idapetsc_mem->s_maxrs    = IDA_SPILS_MAXRS;
  //idapetsc_mem->s_eplifac  = PT05;
  //idapetsc_mem->s_dqincfac = ONE;

  idapetsc_mem->s_last_flag  = IDASPILS_SUCCESS;

  /* Set setupNonNull to FALSE */
  IDA_mem->ida_setupNonNull = FALSE;
  //setupNonNull = FALSE;

  /* Allocate memory for ytemp, yptemp, and xx */

  ytemp = N_VClone(IDA_mem->ida_tempv1);
  if (ytemp == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDAKSP", "IDAKSP", MSGS_MEM_FAIL);
    free(idapetsc_mem); idapetsc_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  yptemp = N_VClone(IDA_mem->ida_tempv1);
  if (yptemp == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDAKSP", "IDAKSP", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    free(idapetsc_mem); idapetsc_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  xx = N_VClone(IDA_mem->ida_tempv1);
  if (xx == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDAKSP", "IDAKSP", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    free(idapetsc_mem); idapetsc_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  idapetsc_mem->s_sqrtN = SUNRsqrt( N_VDotProd(ytemp, ytemp) );

  /* Call KSPMalloc to allocate workspace for KSP */
//   spgmr_mem = NULL;
//   spgmr_mem = KSPMalloc(maxl1, vec_tmpl);
//   if (spgmr_mem == NULL) {
//     IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDAKSP", "IDAKSP", MSGS_MEM_FAIL);
//     N_VDestroy(ytemp);
//     N_VDestroy(yptemp);
//     N_VDestroy(xx);
//     free(idapetsc_mem); idapetsc_mem = NULL;
//     return(IDASPILS_MEM_FAIL);
//   }

  /* Attach KSP solver to its SUNDIALS memory structure */
  idapetsc_mem->s_ksp_mem = (void *) solver;

  /* Attach linear solver memory to the integrator memory */
  IDA_mem->ida_lmem = idapetsc_mem;

  return(IDASPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDAKDP interface routines
 * -----------------------------------------------------------------
 */

/* Additional readability Replacements */

#define gstype   (idapetsc_mem->s_gstype)
#define maxl     (idapetsc_mem->s_maxl)
#define maxrs    (idapetsc_mem->s_maxrs)
#define eplifac  (idapetsc_mem->s_eplifac)
#define psolve   (idapetsc_mem->s_psolve)
#define pset     (idapetsc_mem->s_pset)
#define pdata    (idapetsc_mem->s_pdata)

static int IDAKSPInit(IDAMem IDA_mem)
{
  IDAPETScMem idapetsc_mem;

  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

  /* Initialize counters */
  idapetsc_mem->s_npe  = 0;
  idapetsc_mem->s_nli  = 0;
  idapetsc_mem->s_nps  = 0;
  idapetsc_mem->s_ncfl = 0;
  idapetsc_mem->s_njtimes = 0;
  idapetsc_mem->s_nres = 0;

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  IDA_mem->ida_setupNonNull = (psolve != NULL) && (pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (jtimesDQ) {
    //jtimes = IDAPETScDQJtimes;
    //jdata = IDA_mem;
  } else {
    jdata = user_data;
  }

  idapetsc_mem->s_last_flag = IDASPILS_SUCCESS;
  return(0);
}

static int IDAKSPSetup(IDAMem IDA_mem, 
                       N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  IDAPETScMem idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

  /* Call user setup routine pset and update counter npe. */
  retval = idapetsc_mem->s_pset(IDA_mem->ida_tn, yy_p, yp_p, rr_p, IDA_mem->ida_cj, idapetsc_mem->s_pdata,
                                tmp1, tmp2, tmp3);
  (idapetsc_mem->s_npe)++;

  /* Return flag showing success or failure of pset. */
  if (retval < 0) {
    IDAProcessError(IDA_mem, KSP_PSET_FAIL_UNREC, "IDAKSP", "IDAKSPSetup", MSGS_PSET_FAILED);
    idapetsc_mem->s_last_flag = KSP_PSET_FAIL_UNREC;
    return(-1);
  }
  if (retval > 0) {
    idapetsc_mem->s_last_flag = KSP_PSET_FAIL_REC;
    return(+1);
  }

  idapetsc_mem->s_last_flag = KSP_SUCCESS;
  return(0);
}


/*
 * The x-scaling and b-scaling arrays are both equal to weight.
 *  
 * We set the initial guess, x = 0, then call KSPSolve.  
 * We copy the solution x into b, and update the counters nli, nps, ncfl.
 * If KSPSolve returned nli_inc = 0 (hence x = 0), we take the KSP
 * vtemp vector (= P_inverse F) as the correction vector instead.
 *  Finally, we set the return value according to the success of KSPSolve.
 */

static int IDAKSPSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                       N_Vector yy_now, N_Vector yp_now, N_Vector rr_now)
{
  IDAPETScMem idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;
  KSP *solver = (KSP*) idapetsc_mem->s_ksp_mem;
  int pretype;
//   int nli_inc, nps_inc, retval;
  realtype res_norm;
  PetscErrorCode ierr;

  /* Set KSPSolve convergence test constant epslin, in terms of the
    Newton convergence test constant epsNewt and safety factors.  The factor 
    sqrt(Neq) assures that the GMRES convergence test is applied to the
    WRMS norm of the residual vector, rather than the weighted L2 norm. */
  epslin = (idapetsc_mem->s_sqrtN)*(idapetsc_mem->s_eplifac)*(IDA_mem->ida_epsNewt);

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  idapetsc_mem->s_ycur  = yy_now;
  idapetsc_mem->s_ypcur = yp_now;
  idapetsc_mem->s_rcur  = rr_now;

  /* Set KSPSolve inputs pretype and initial guess xx = 0. */  
  pretype = (psolve == NULL) ? PREC_NONE : PREC_LEFT;
//   N_VConst(ZERO, xx);
  
  /* Call KSPSolve and copy xx to bb. */
  ierr = KSPSolve(*solver, *(NV_PVEC_PTC(bb)), *(NV_PVEC_PTC(bb)));

//   if (nli_inc == 0) N_VScale(ONE, KSP_VTEMP(spgmr_mem), bb);
//   else N_VScale(ONE, xx, bb);
//   
//   /* Increment counters nli, nps, and return if successful. */
//   idapetsc_mem->s_nli += nli_inc;
//   idapetsc_mem->s_nps += nps_inc;
//   if (retval != KSP_SUCCESS) 
//     (idapetsc_mem->s_ncfl)++;

  /* Interpret return value from KSPSolve */

  idapetsc_mem->s_last_flag = ierr;
  CHKERRQ(ierr);
  
//   switch(retval) {
// 
//   case KSP_SUCCESS:
//     return(0);
//     break;
//   case KSP_RES_REDUCED:
//     return(1);
//     break;
//   case KSP_CONV_FAIL:
//     return(1);
//     break;
//   case KSP_QRFACT_FAIL:
//     return(1);
//     break;
//   case KSP_PSOLVE_FAIL_REC:
//     return(1);
//     break;
//   case KSP_ATIMES_FAIL_REC:
//     return(1);
//     break;
//   case KSP_MEM_NULL:
//     return(-1);
//     break;
//   case KSP_ATIMES_FAIL_UNREC:
//     IDAProcessError(IDA_mem, KSP_ATIMES_FAIL_UNREC, "IDAKSP", "IDAKSPSolve", MSGS_JTIMES_FAILED);    
//     return(-1);
//     break;
//   case KSP_PSOLVE_FAIL_UNREC:
//     IDAProcessError(IDA_mem, KSP_PSOLVE_FAIL_UNREC, "IDAKSP", "IDAKSPSolve", MSGS_PSOLVE_FAILED);
//     return(-1);
//     break;
//   case KSP_GS_FAIL:
//     return(-1);
//     break;
//   case KSP_QRSOL_FAIL:
//     return(-1);
//     break;
//   }

  return(0);
}

/*
 * This routine handles performance monitoring specific to the IDAKSP
 * linear solver.  When perftask = 0, it saves values of various counters.
 * When perftask = 1, it examines difference quotients in these counters,
 * and depending on their values, it prints up to three warning messages.
 * Messages are printed up to a maximum of 10 times.
 * 
 * TODO: Need to figure out how to use PETSc built-in stats. Disable for now!
 */

static int IDAKSPPerf(IDAMem IDA_mem, int perftask)
{
  IDAPETScMem idapetsc_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

//   if (perftask == 0) {
//     nst0 = nst;  nni0 = nni;  nli0 = nli;
//     ncfn0 = ncfn;  ncfl0 = ncfl;  
//     nwarn = 0;
//     return(0);
//   }
// 
//   nstd = nst - nst0;  nnid = nni - nni0;
//   if (nstd == 0 || nnid == 0) return(0);
//   avdim = (realtype) ((nli - nli0)/((realtype) nnid));
//   rcfn = (realtype) ((ncfn - ncfn0)/((realtype) nstd));
//   rcfl = (realtype) ((ncfl - ncfl0)/((realtype) nnid));
//   lavd = (avdim > ((realtype) maxl ));
//   lcfn = (rcfn > PT9);
//   lcfl = (rcfl > PT9);
//   if (!(lavd || lcfn || lcfl)) return(0);
//   nwarn++;
//   if (nwarn > 10) return(1);
//   if (lavd) 
//     IDAProcessError(IDA_mem, IDA_WARNING, "IDAKSP", "IDAKSPPerf", MSGS_AVD_WARN, tn, avdim);
//   if (lcfn) 
//     IDAProcessError(IDA_mem, IDA_WARNING, "IDAKSP", "IDAKSPPerf", MSGS_CFN_WARN, tn, rcfn);
//   if (lcfl) 
//     IDAProcessError(IDA_mem, IDA_WARNING, "IDAKSP", "IDAKSPPerf", MSGS_CFL_WARN, tn, rcfl);

  return(0);
}

static int IDAKSPFree(IDAMem IDA_mem)
{
  IDAPETScMem idapetsc_mem;
  KSP *solver;
  PetscErrorCode ierr;

  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;
  
  N_VDestroy(ytemp);
  N_VDestroy(yptemp);
  N_VDestroy(xx);

  solver = (KSP*) idapetsc_mem->s_ksp_mem;
  ierr = KSPDestroy(solver);
  CHKERRQ(ierr);

  if (idapetsc_mem->s_pfree != NULL) 
    (idapetsc_mem->s_pfree)(IDA_mem);

  free(idapetsc_mem); idapetsc_mem = NULL;

  return(0);
}

