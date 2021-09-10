/*
 * -----------------------------------------------------------------
 * $Revision:  $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmers: Slaven Peles @ LLNL
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
 * This is the implementation file for the IDA interface to linear
 * solvers from PETSc library.            
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

//#include <petscksp.h>

#include <ida/ida_petsc.h>
#include <nvector/nvector_petsc.h>
#include "ida_petsc_impl.h"
#include "ida_impl.h"

#include <sundials/sundials_math.h>

#define NV_CONTENT_PTC(v)    ( (N_VectorContent_Petsc)(v->content) )
#define NV_PVEC_PTC(v)       ( NV_CONTENT_PTC(v)->pvec )


/* Return values for KSPSolve */

#define KSP_SUCCESS            0  /* Converged                     */
#define KSP_RES_REDUCED        1  /* Did not converge, but reduced
                                     norm of residual              */
#define KSP_CONV_FAIL          2  /* Failed to converge            */
#define KSP_QRFACT_FAIL        3  /* QRfact found singular matrix  */
#define KSP_PSOLVE_FAIL_REC    4  /* psolve failed recoverably     */
#define KSP_ATIMES_FAIL_REC    5  /* atimes failed recoverably     */
#define KSP_JAC_FAIL_REC       6  /* Jacobian failedrecoverably    */

#define KSP_MEM_NULL          -1  /* mem argument is NULL          */
#define KSP_ATIMES_FAIL_UNREC -2  /* atimes returned failure flag  */
#define KSP_PSOLVE_FAIL_UNREC -3  /* psolve failed unrecoverably   */
#define KSP_GS_FAIL           -4  /* Gram-Schmidt routine faiuled  */        
#define KSP_QRSOL_FAIL        -5  /* QRsol found singular R        */
#define KSP_JAC_FAIL_UNREC    -6  /* Jacobian failed unrecoverably */


/* Error message defaults */

#define SUNMODULE "IDA PETSc KSP" /* SUNDIALS module ID            */
#define SUNFUNC   "IDAPETScKSP"   /* SUNDIALS function prefix      */

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)
#define PT9          RCONST(0.9)
#define PT05         RCONST(0.05)

/* IDAPETScKSP linit, lsetup, lsolve, lperf, and lfree routines */

static int IDAPETScKSPInit(IDAMem IDA_mem);

static int IDAPETScKSPSetup(IDAMem IDA_mem, 
                         N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int IDAPETScKSPSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector yy_now, N_Vector yp_now, N_Vector rr_now);

static int IDAPETScKSPPerf(IDAMem IDA_mem, int perftask);

static int IDAPETScKSPFree(IDAMem IDA_mem);


/*
 * --------------------------------------------------------------------------
 * IDAPETScKSP
 * --------------------------------------------------------------------------
 *
 * This routine initializes the memory record and sets various function
 * fields specific to the IDAPETScKSP linear solver module.  
 *
 * IDAPETScKSP first calls the existing lfree routine if this is not NULL.
 * It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDAPETScKSPInit, IDAPETScKSPSetup,
 * IDAPETScKSPSolve, IDAPETScKSPPerf, and IDAPETScKSPFree, respectively.
 * It allocates memory for a structure of type IDAPETScMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem).  It then various fields in the
 * IDAPETScMemRec structure. Finally, IDAPETScKSP allocates memory for 
 * ytemp, yptemp, and xx, and calls KSPMalloc to allocate memory
 * for the KSP solver.
 *
 * The return value of IDAPETScKSP is:
 *   PETSC_KSP_SUCCESS   =  0 if successful
 *   PETSC_KSP_MEM_NULL  = -1 if IDA_mem is NULL or a memory allocation failed
 *   PETSC_KSP_LMEM_NULL = -2 if linear solver memory is NULL
 *   PETSC_KSP_ILL_INPUT = -3 if the argument is illegal.
 *
 * --------------------------------------------------------------------------
 */
int IDAPETScKSP(void *ida_mem, MPI_Comm comm, Mat *JacMat)
{
  IDAMem IDA_mem;
  IDAPETScMem idapetsc_mem;
  KSP *solver;
  int flag;
  PetscErrorCode ierr;
  
  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_NULL, SUNMODULE, SUNFUNC, MSGS_IDAMEM_NULL);
    return(PETSC_KSP_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if(IDA_mem->ida_tempv1->ops->nvdotprod == NULL) {
    IDAProcessError(NULL, PETSC_KSP_ILL_INPUT, SUNMODULE, SUNFUNC, MSGS_BAD_NVECTOR);
    return(PETSC_KSP_ILL_INPUT);
  }

  /* If there is an instance of linear solver associated with IDA, delete it. */
  if (IDA_mem->ida_lfree != NULL) 
    flag = IDA_mem->ida_lfree((IDAMem) ida_mem);

  /* Set five main function fields in ida_mem */
  IDA_mem->ida_linit  = IDAPETScKSPInit;
  IDA_mem->ida_lsetup = IDAPETScKSPSetup;
  IDA_mem->ida_lsolve = IDAPETScKSPSolve;
  IDA_mem->ida_lperf  = IDAPETScKSPPerf;
  IDA_mem->ida_lfree  = IDAPETScKSPFree;

  /* Set setupNonNull to SUNTRUE */
  /* IDA_mem->ida_setupNonNull = SUNTRUE; */

  /* Get memory for IDAPETScMemRec. */
  idapetsc_mem = (IDAPETScMem) malloc(sizeof(struct IDAPETScMemRec));
  if (idapetsc_mem == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_MEM_FAIL, SUNMODULE, SUNFUNC, MSGS_MEM_FAIL);
    return(PETSC_KSP_MEM_FAIL);
  }

  /* Set pointer to user data attached to the solver */
  idapetsc_mem->s_udata  = IDA_mem->ida_user_data;
  idapetsc_mem->s_last_flag  = PETSC_KSP_SUCCESS;

  /* Allocate memory for ytemp, yptemp, and xx */

  idapetsc_mem->s_ytemp = N_VClone(IDA_mem->ida_tempv1);
  if (idapetsc_mem->s_ytemp == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_FAIL, SUNMODULE, SUNFUNC, MSGS_MEM_FAIL);
    free(idapetsc_mem); idapetsc_mem = NULL;
    return(PETSC_KSP_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, idapetsc_mem->s_ytemp);
  idapetsc_mem->s_sqrtN = SUNRsqrt( N_VDotProd(idapetsc_mem->s_ytemp, idapetsc_mem->s_ytemp) );
  
  /* Set linear solver tolerance factor */
  idapetsc_mem->s_eplifac = PT05;

  /* Initialize counters */
  idapetsc_mem->s_npe  = 0;
  idapetsc_mem->s_nli  = 0;
  idapetsc_mem->s_nps  = 0;
  idapetsc_mem->s_ncfl = 0;
  idapetsc_mem->s_njtimes = 0;
  idapetsc_mem->s_nres = 0;
  idapetsc_mem->s_nje  = 0;

  /* Allocate memory for KSP solver */
  solver = NULL;
  solver = (KSP*) malloc(sizeof(KSP));
  if (solver == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_FAIL, SUNMODULE, SUNFUNC, MSGS_MEM_FAIL);
    N_VDestroy(idapetsc_mem->s_ytemp);
    free(idapetsc_mem);
    idapetsc_mem = NULL;
    return(PETSC_KSP_MEM_FAIL);
  }

  /* Create KSP solver */
  ierr = KSPCreate(comm, solver);
  CHKERRQ(ierr);

  /* Attach Jacobian matrix to the KSP solver */
  ierr = KSPSetOperators(*solver, *JacMat, *JacMat);
  CHKERRQ(ierr);

  /* store pointer to Jacobian matrix */
  idapetsc_mem->JacMat = JacMat;  
    
  /* Attach KSP solver to its SUNDIALS memory structure */
  idapetsc_mem->s_ksp_mem = solver;

  /* Attach linear solver memory to the integrator memory */
  IDA_mem->ida_lmem = idapetsc_mem;

  return PETSC_KSP_SUCCESS;
}


/*
 * -----------------------------------------------------------------
 * IDAPETScKSP interface routines
 * -----------------------------------------------------------------
 */

static int IDAPETScKSPInit(IDAMem IDA_mem)
{
  IDAPETScMem idapetsc_mem;
  KSP *solver;
  PetscErrorCode ierr;

  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;
  solver = idapetsc_mem->s_ksp_mem;
  
  /* Initialize counters */
  idapetsc_mem->s_npe  = 0;
  idapetsc_mem->s_nli  = 0;
  idapetsc_mem->s_nps  = 0;
  idapetsc_mem->s_ncfl = 0;
  idapetsc_mem->s_njtimes = 0;
  idapetsc_mem->s_nres = 0;
  idapetsc_mem->s_nje  = 0;

  /* Set options for the PETSc linear solver */
  ierr = KSPSetFromOptions(*solver);
  CHKERRQ(ierr);

  idapetsc_mem->s_last_flag = PETSC_KSP_SUCCESS;
  return(0);
}

/*
 * --------------------------------------------------------------------------
 * IDAPETScKSPSetup
 * --------------------------------------------------------------------------
 *
 * Calls user supplied function for Jacobian evaluation.
 *
 * The return value of IDAPETScKSP is:
 *   PETSC_KSP_SUCCESS   =  0 if successful
 *   PETSC_KSP_MEM_NULL  = -1 if IDA_mem is NULL or a memory allocation failed
 *   PETSC_KSP_LMEM_NULL = -2 if linear solver memory is NULL
 *   PETSC_KSP_ILL_INPUT = -3 if the argument is illegal.
 * 
 * TODO: Consolidate error messages!
 *
 * --------------------------------------------------------------------------
 */
static int IDAPETScKSPSetup(IDAMem IDA_mem, 
                            N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  IDAPETScMem idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;
  
  if (idapetsc_mem->s_jaceval == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_ILL_INPUT, SUNMODULE, SUNFUNC "Setup", MSGS_JAC_NULL);
    return(PETSC_KSP_ILL_INPUT);
  }
  if (idapetsc_mem->JacMat == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_ILL_INPUT, SUNMODULE, SUNFUNC "Setup", MSGS_JAC_MAT_NULL);
    return(PETSC_KSP_ILL_INPUT);
  }
  
  /* Call user setup routine jaceval and update counter nje. */
  retval = idapetsc_mem->s_jaceval(IDA_mem->ida_tn, IDA_mem->ida_cj,
                                   yy_p, yp_p, rr_p, 
                                   *(idapetsc_mem->JacMat), idapetsc_mem->s_udata,
                                   tmp1, tmp2, tmp3);
  (idapetsc_mem->s_nje)++;

  /* Return flag showing success or failure of jaceval. */
  if (retval < 0) {
    IDAProcessError(IDA_mem, KSP_JAC_FAIL_UNREC, SUNMODULE, SUNFUNC "Setup", MSGS_JAC_FAILED);
    idapetsc_mem->s_last_flag = KSP_JAC_FAIL_UNREC;
    return(-1);
  }
  if (retval > 0) {
    idapetsc_mem->s_last_flag = KSP_JAC_FAIL_REC;
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
 * Finally, we set the return value according to the success of KSPSolve.
 * 
 * Use KSPGetConvergedReason() to find the outcome of the solve. 
 * typedef enum {// converged 
              KSP_CONVERGED_RTOL_NORMAL        =  1,
              KSP_CONVERGED_ATOL_NORMAL        =  9,
              KSP_CONVERGED_RTOL               =  2,
              KSP_CONVERGED_ATOL               =  3,
              KSP_CONVERGED_ITS                =  4,
              KSP_CONVERGED_CG_NEG_CURVE       =  5,
              KSP_CONVERGED_CG_CONSTRAINED     =  6,
              KSP_CONVERGED_STEP_LENGTH        =  7,
              KSP_CONVERGED_HAPPY_BREAKDOWN    =  8,
              // diverged 
              KSP_DIVERGED_NULL                = -2,
              KSP_DIVERGED_ITS                 = -3,
              KSP_DIVERGED_DTOL                = -4,
              KSP_DIVERGED_BREAKDOWN           = -5,
              KSP_DIVERGED_BREAKDOWN_BICG      = -6,
              KSP_DIVERGED_NONSYMMETRIC        = -7,
              KSP_DIVERGED_INDEFINITE_PC       = -8,
              KSP_DIVERGED_NANORINF            = -9,
              KSP_DIVERGED_INDEFINITE_MAT      = -10,
              KSP_DIVERGED_PCSETUP_FAILED      = -11,

              KSP_CONVERGED_ITERATING          =  0} KSPConvergedReason;
 */

static int IDAPETScKSPSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                            N_Vector yy_now, N_Vector yp_now, N_Vector rr_now)
{
  IDAPETScMem idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;
  KSP *solver = idapetsc_mem->s_ksp_mem;
  int pretype;
  int nli_inc;
  realtype res_norm;
  PetscErrorCode ierr;
  
  realtype cjratio = IDA_mem->ida_cjratio;
  /* printf("cjratio should be 1, is: %g.\n", cjratio); */

  /* Set KSPSolve convergence test constant epslin, in terms of the
    Newton convergence test constant epsNewt and safety factors.  The factor 
    sqrt(Neq) assures that the GMRES convergence test is applied to the
    WRMS norm of the residual vector, rather than the weighted L2 norm. */
  idapetsc_mem->s_epslin = (idapetsc_mem->s_sqrtN)*(idapetsc_mem->s_eplifac)*(IDA_mem->ida_epsNewt);
//   printf("Estimated linear solver tolerance: %g\n",idapetsc_mem->s_epslin );
//   printf("sqrtN: %g\n", idapetsc_mem->s_sqrtN);
//   printf("eplifac: %g\n", idapetsc_mem->s_eplifac);
//   printf("epsNewt: %g\n", IDA_mem->ida_epsNewt);
  ierr = KSPSetTolerances(*solver, idapetsc_mem->s_epslin, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  
  /* Call KSPSolve and copy xx to bb. */
  ierr = KSPSolve(*solver, *(NV_PVEC_PTC(bb)), *(NV_PVEC_PTC(bb)));

  KSPGetIterationNumber(*solver, &nli_inc);
  
/*   if (nli_inc == 0) N_VScale(ONE, KSP_VTEMP(spgmr_mem), bb);
//   else N_VScale(ONE, xx, bb);
*/   
  
  /* Increment counters nli, nps, and return if successful. */
  idapetsc_mem->s_nli += nli_inc;
  if (ierr != KSP_SUCCESS) 
    (idapetsc_mem->s_ncfl)++;
  //printf("nli = %d\n", idapetsc_mem->s_nli);

  /* Interpret return value from KSPSolve */

  idapetsc_mem->s_last_flag = ierr;
  //CHKERRQ(ierr);
  
  
  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), bb, bb);

  idapetsc_mem->s_last_flag = PETSC_KSP_SUCCESS;

  
/*   switch(retval) {
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
//     IDAProcessError(IDA_mem, KSP_ATIMES_FAIL_UNREC, SUNMODULE, "IDAKSPSolve", MSGS_JTIMES_FAILED);    
//     return(-1);
//     break;
//   case KSP_PSOLVE_FAIL_UNREC:
//     IDAProcessError(IDA_mem, KSP_PSOLVE_FAIL_UNREC, SUNMODULE, "IDAKSPSolve", MSGS_PSOLVE_FAILED);
//     return(-1);
//     break;
//   case KSP_GS_FAIL:
//     return(-1);
//     break;
//   case KSP_QRSOL_FAIL:
//     return(-1);
//     break;
//   }
*/

  return(0);
}

/*
 * This routine handles performance monitoring specific to the IDAPETScKSP
 * linear solver.  When perftask = 0, it saves values of various counters.
 * When perftask = 1, it examines difference quotients in these counters,
 * and depending on their values, it prints up to three warning messages.
 * Messages are printed up to a maximum of 10 times.
 * 
 * TODO: Need to figure out how to use PETSc built-in stats. Disable for now!
 */

static int IDAPETScKSPPerf(IDAMem IDA_mem, int perftask)
{
  IDAPETScMem idapetsc_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

/*
//   if (perftask == 0) {
//     nst0 = nst;  nni0 = nni;  nli0 = idapetsc_mem->s_nli;
//     ncfn0 = ncfn;  ncfl0 = idapetsc_mem->s_ncfl;  
//     nwarn = 0;
//     return(0);
//   }
// 
//   nstd = nst - nst0;  nnid = nni - nni0;
//   if (nstd == 0 || nnid == 0) return(0);
//   avdim = (realtype) ((idapetsc_mem->s_nli - nli0)/((realtype) nnid));
//   rcfn = (realtype) ((ncfn - ncfn0)/((realtype) nstd));
//   rcfl = (realtype) ((idapetsc_mem->s_ncfl - ncfl0)/((realtype) nnid));
//   lavd = (avdim > ((realtype) maxl ));
//   lcfn = (rcfn > PT9);
//   lcfl = (rcfl > PT9);
//   if (!(lavd || lcfn || lcfl)) return(0);
//   nwarn++;
//   if (nwarn > 10) return(1);
//   if (lavd) 
//     IDAProcessError(IDA_mem, IDA_WARNING, SUNMODULE, "IDAKSPPerf", MSGS_AVD_WARN, tn, avdim);
//   if (lcfn) 
//     IDAProcessError(IDA_mem, IDA_WARNING, SUNMODULE, "IDAKSPPerf", MSGS_CFN_WARN, tn, rcfn);
//   if (lcfl) 
//     IDAProcessError(IDA_mem, IDA_WARNING, SUNMODULE, "IDAKSPPerf", MSGS_CFL_WARN, tn, rcfl);
*/

  return(0);
}

/*
 * Delete KSP context and idapetsc_mem.
 */
static int IDAPETScKSPFree(IDAMem IDA_mem)
{
  IDAPETScMem idapetsc_mem;
  PetscErrorCode ierr;

  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;
  
  N_VDestroy(idapetsc_mem->s_ytemp);

  ierr = KSPDestroy(idapetsc_mem->s_ksp_mem);
  idapetsc_mem->s_ksp_mem = NULL;
  
  free(idapetsc_mem); 
  idapetsc_mem = NULL;

  CHKERRQ(ierr);
  return(0);
}


/*
 * IDAPETScSetJacFn specifies the PETSc Jacobian function.
 */
int IDAPETScSetJacFn(void* ida_mem, IDAPETScJacFn jac)
{
  IDAMem IDA_mem;
  IDAPETScMem idapetsc_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_NULL, SUNMODULE, "IDAPETScSetJacFn", MSGS_IDAMEM_NULL);
    return(PETSC_KSP_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_LMEM_NULL, SUNMODULE, 
            "IDAPETScSetJacFn", MSGS_LMEM_NULL);
    return(PETSC_KSP_LMEM_NULL);
  }
  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

  idapetsc_mem->s_jaceval = jac;

  return(PETSC_KSP_SUCCESS);
}

/*
 * Get number of PETSc linear solver iterations.
 */
int IDAPETScGetNumLinIters(void* ida_mem, long int* nliters)
{
  IDAMem IDA_mem;
  IDAPETScMem idapetsc_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_NULL, SUNMODULE, "IDAPETScGetNumLinIters", MSGS_IDAMEM_NULL);
    return(PETSC_KSP_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_LMEM_NULL, SUNMODULE, "IDAPETScGetNumLinIters", MSGS_LMEM_NULL);
    return(PETSC_KSP_LMEM_NULL);
  }
  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

  *nliters = idapetsc_mem->s_nli;

  return(PETSC_KSP_SUCCESS);
}

/*
 * Get number of PETSc linear solver convergence fails.
 */
int IDAPETScGetNumConvFails(void *ida_mem, long int *nlcfails)
{
  IDAMem IDA_mem;
  IDAPETScMem idapetsc_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_NULL, SUNMODULE, "IDAPETScGetNumConvFails", MSGS_IDAMEM_NULL);
    return(PETSC_KSP_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_LMEM_NULL, SUNMODULE, "IDAPETScGetNumConvFails", MSGS_LMEM_NULL);
    return(PETSC_KSP_LMEM_NULL);
  }
  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

  *nlcfails = idapetsc_mem->s_ncfl;

  return(PETSC_KSP_SUCCESS);
}

/*
 * Get number of Jacobian evaluations.
 */
int IDAPETScGetNumJacEvals(void* ida_mem, long int* njacevals)
{
  IDAMem IDA_mem;
  IDAPETScMem idapetsc_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, PETSC_KSP_MEM_NULL, SUNMODULE, "IDAPETScGetNumJtimesEvals", MSGS_IDAMEM_NULL);
    return(PETSC_KSP_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, PETSC_KSP_LMEM_NULL, SUNMODULE, "IDAPETScGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(PETSC_KSP_LMEM_NULL);
  }
  idapetsc_mem = (IDAPETScMem) IDA_mem->ida_lmem;

  *njacevals = idapetsc_mem->s_nje;

  return(PETSC_KSP_SUCCESS);
}

