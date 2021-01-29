/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with CPODES, a CPSPILS linear
 * solver, and the parallel implementation of NVECTOR.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpodes/cpodes_sptfqmr.h>
#include <cpodes/cpodes_spbcgs.h>
#include <cpodes/cpodes_spgmr.h>

#include "cpodes_bbdpre_impl.h"
#include "cpodes_private.h"
#include "cpodes_spils_impl.h"

#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define MIN_INC_MULT RCONST(1000.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* cpBBDPrecSetupExpl and cpBBDPrecSetupImpl */

static int cpBBDPrecSetupExpl(realtype t, N_Vector y, N_Vector fy, 
                              booleantype jok, booleantype *jcurPtr, 
                              realtype gamma, void *bbd_data, 
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cpBBDPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              realtype gamma, void *bbd_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* cpBBDPrecSolveExpl and cpBBDPrecSolveImpl */

static int cpBBDPrecSolveExpl(realtype t, N_Vector y, N_Vector fy, 
                              N_Vector b, N_Vector x, 
                              realtype gamma, realtype delta,
                              int lr, void *bbd_data, N_Vector tmp);

static int cpBBDPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              N_Vector b, N_Vector x,
                              realtype gamma, realtype delta, 
                              void *bbd_data, N_Vector tmp);

/* Prototype for cpBBDPrecFree */
static int cpBBDPrecFree(CPodeMem cp_mem);

/* Difference quotient Jacobian calculation routines */

static int cpBBDDQJacExpl(CPBBDPrecData pdata, realtype t, 
                          N_Vector y, N_Vector gy, 
                          N_Vector ytemp, N_Vector gtemp);

static int cpBBDDQJacImpl(CPBBDPrecData pdata, realtype tt, realtype gamma,
                          N_Vector yy, N_Vector yp, N_Vector gref, 
                          N_Vector ytemp, N_Vector yptemp, N_Vector gtemp);

/* Redability replacements */

#define ode_type (cp_mem->cp_ode_type)
#define uround   (cp_mem->cp_uround)
#define vec_tmpl (cp_mem->cp_tempv)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

int CPBBDPrecInit(void *cpode_mem, int Nlocal, 
                  int mudq, int mldq,
                  int mukeep, int mlkeep, 
                  realtype dqrely, 
                  void *gloc, CPBBDCommFn cfn)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBBDPrecData pdata;
  N_Vector tmp4;
  int muk, mlk, storage_mu;
  int flag=0;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  /* Test if the NVECTOR package is compatible with the BLOCK BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_BAD_NVECTOR);
    return(CPSPILS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (CPBBDPrecData) malloc(sizeof *pdata);  
  if (pdata == NULL) {
    cpProcessError(cp_mem, CPSPILS_ILL_INPUT, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CPSPILS_ILL_INPUT);
  }

  /* Set pointers to gloc and cfn; load half-bandwidths */

  pdata->cpode_mem = cpode_mem;

  switch (ode_type) {
  case CP_EXPL:
    pdata->glocE = (CPBBDLocalRhsFn) gloc;
    pdata->glocI = NULL;
    break;
  case CP_IMPL:
    pdata->glocI = (CPBBDLocalResFn) gloc;
    pdata->glocE = NULL;
    break;
  }

  pdata->cfn = cfn;

  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));

  muk = SUNMIN(Nlocal-1, SUNMAX(0,mukeep));
  mlk = SUNMIN(Nlocal-1, SUNMAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = NewBandMat(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) { 
    free(pdata); pdata = NULL; 
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CPSPILS_MEM_FAIL); 
  }

  /* Allocate memory for preconditioner matrix */
  storage_mu = SUNMIN(Nlocal-1, muk + mlk);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Allocate memory for pivots */
  pdata->pivots = NULL;
  pdata->pivots = NewIndexArray(Nlocal);
  if (pdata->savedJ == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }

  /* Allocate tmp4 for use by cpBBDDQJacImpl */
  tmp4 = NULL;
  tmp4 = N_VClone(vec_tmpl); 
  if (tmp4 == NULL){
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    DestroyArray(pdata->pivots);
    free(pdata); pdata = NULL;
    cpProcessError(cp_mem, CPSPILS_MEM_FAIL, "CPBBDPRE", "CPBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CPSPILS_MEM_FAIL);
  }
  pdata->tmp4 = tmp4;

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : SUNRsqrt(uround);

  /* Store Nlocal to be used in cpBBDPrecSetupExpl and cpBBDPrecSetupImpl */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(muk + 2*mlk + storage_mu + 2);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  /* Overwrite the P_data field in the SPILS memory */
  cpspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  cpspils_mem->s_pfree = cpBBDPrecFree;

  /* Attach preconditioner solve and setup functions */
  switch (ode_type) {
  case CP_EXPL:
    flag = CPSpilsSetPrecFnExpl(cpode_mem, cpBBDPrecSetupExpl, cpBBDPrecSolveExpl);
    break;
  case CP_IMPL:
    flag = CPSpilsSetPrecFnImpl(cpode_mem, cpBBDPrecSetupImpl, cpBBDPrecSolveImpl);
    break;
  }

  return(flag);
}


int CPBBDPrecReInit(void *cpode_mem, int mudq, int mldq, realtype dqrely)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBBDPrecData pdata;
  int Nlocal;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBBDPRE", "CPBBDPrecReInit", MSGBBD_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBBDPRE", "CPBBDPrecReInit", MSGBBD_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  /* Test if the preconditioner data is non-NULL */
  if (cpspils_mem->s_P_data == NULL) {
    cpProcessError(cp_mem, CPSPILS_PMEM_NULL, "CPBBDPRE", "CPBBDPrecReInit", MSGBBD_PMEM_NULL);
    return(CPSPILS_PMEM_NULL);
  } 
  pdata = (CPBBDPrecData) cpspils_mem->s_P_data;

  /* Load half-bandwidths */
  Nlocal = pdata->n_local;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : SUNRsqrt(uround);

  /* Re-initialize nge */
  pdata->nge = 0;

  return(CPSPILS_SUCCESS);
}


int CPBBDPrecGetWorkSpace(void *cpode_mem, long int *lenrwBBDP, long int *leniwBBDP)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBBDPrecData pdata;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBBDPRE", "CPBBDPrecGetWorkSpace", MSGBBD_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBBDPRE", "CPBBDPrecGetWorkSpace", MSGBBD_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  if (cpspils_mem->s_P_data == NULL) {
    cpProcessError(cp_mem, CPSPILS_PMEM_NULL, "CPBBDPRE", "CPBBDPrecGetWorkSpace", MSGBBD_PMEM_NULL);
    return(CPSPILS_PMEM_NULL);
  } 
  pdata = (CPBBDPrecData) cpspils_mem->s_P_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(CPSPILS_SUCCESS);
}


int CPBBDPrecGetNumGfnEvals(void *cpode_mem, long int *ngevalsBBDP)
{
  CPodeMem cp_mem;
  CPSpilsMem cpspils_mem;
  CPBBDPrecData pdata;

  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPSPILS_MEM_NULL, "CPBBDPRE", "CPBBDPrecGetNumGfnEvals", MSGBBD_CPMEM_NULL);
    return(CPSPILS_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (cp_mem->cp_lmem == NULL) {
    cpProcessError(cp_mem, CPSPILS_LMEM_NULL, "CPBBDPRE", "CPBBDPrecGetNumGfnEvals", MSGBBD_LMEM_NULL);
    return(CPSPILS_LMEM_NULL);
  }
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;

  if (cpspils_mem->s_P_data == NULL) {
    cpProcessError(cp_mem, CPSPILS_PMEM_NULL, "CPBBDPRE", "CPBBDPrecGetNumGfnEvals", MSGBBD_PMEM_NULL);
    return(CPSPILS_PMEM_NULL);
  } 
  pdata = (CPBBDPrecData) cpspils_mem->s_P_data;

  *ngevalsBBDP = pdata->nge;

  return(CPSPILS_SUCCESS);
}


/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/* Readability Replacements */

#define Nlocal (pdata->n_local)
#define mudq   (pdata->mudq)
#define mldq   (pdata->mldq)
#define mukeep (pdata->mukeep)
#define mlkeep (pdata->mlkeep)
#define dqrely (pdata->dqrely)
#define glocE  (pdata->glocE)
#define glocI  (pdata->glocI)
#define cfn    (pdata->cfn)
#define savedJ (pdata->savedJ)
#define savedP (pdata->savedP)
#define pivots (pdata->pivots)
#define nge    (pdata->nge)

/*
 * -----------------------------------------------------------------
 * Function: cpBBDPrecSetupExpl                                      
 * -----------------------------------------------------------------
 * cpBBDPrecSetupExpl generates and factors a banded block of the
 * preconditioner matrix on each processor, via calls to the
 * user-supplied glocE and cfn functions. It uses difference
 * quotient approximations to the Jacobian elements.
 *
 * cpBBDPrecSetupExpl calculates a new J,if necessary, then
 * calculates P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of cpBBDPrecSetupExpl used here are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == SUNFALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == SUNTRUE  means that Jacobian data from the
 *                  previous cpBBDPrecSetupExpl call can be reused
 *                  (with the current value of gamma).
 *         A cpBBDPrecSetupExpl call with jok == SUNTRUE should only
 *         occur after a call with jok == SUNFALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         set by cpBBDPrecSetupExpl as follows:
 *           *jcurPtr = SUNTRUE if Jacobian data was recomputed.
 *           *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
 *                      but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bbd_data  is a pointer to user data returned by CPBBDPrecInit.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 *           NVectors which are be used by cpBBDPrecSetupExpl
 *           as temporary storage or work space.
 *
 * Return value:
 * The value returned by this cpBBDPrecSetupExpl function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried).
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSetupExpl(realtype t, N_Vector y, N_Vector fy, 
                              booleantype jok, booleantype *jcurPtr, 
                              realtype gamma, void *bbd_data, 
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int ier;
  CPBBDPrecData pdata;
  CPodeMem cp_mem;
  int retval;

  pdata = (CPBBDPrecData) bbd_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  if (jok) {

    /* If jok = SUNTRUE, use saved copy of J */
    *jcurPtr = SUNFALSE;
    BandCopy(savedJ, savedP, mukeep, mlkeep);

  } else {

    /* Otherwise call cpBBDDQJacExpl for new J value */
    *jcurPtr = SUNTRUE;
    retval = cpBBDDQJacExpl(pdata, t, y, tmp1, tmp2, tmp3);
    if (retval < 0) {
      cpProcessError(cp_mem, -1, "CPBBDPRE", "cpBBDPrecSetup", MSGBBD_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mukeep, mlkeep);

  }
  
  /* Scale and add I to get P = I - gamma*J */
  BandScale(-gamma, savedP);
  AddIdentity(savedP);
 
  /* Do LU factorization of P in place */
  ier = BandGBTRF(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function: cpBBDPrecSolveExpl
 * -----------------------------------------------------------------
 * cpBBDPrecSolveExpl solves a linear system P z = r, with the
 * band-block-diagonal preconditioner matrix P generated and
 * factored by cpBBDPrecSetupExpl.
 *
 * The parameters of cpBBDPrecSolveExpl used here are as follows:
 *
 *   r - right-hand side vector of the linear system.
 *
 *   bbd_data - pointer to the preconditioner data returned by
 *              CPBBDPrecInit.
 *
 *   z - output vector computed by cpBBDPrecSolveExpl.
 *
 * The value returned by the cpBBDPrecSolveExpl function is always
 * 0, indicating success.
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSolveExpl(realtype t, N_Vector y, N_Vector fy, 
                              N_Vector b, N_Vector x, 
                              realtype gamma, realtype delta,
                              int lr, void *bbd_data, N_Vector tmp)
{
  CPBBDPrecData pdata;
  realtype *xd;

  pdata = (CPBBDPrecData) bbd_data;

  /* Copy b to x, then do backsolve and return */
  N_VScale(ONE, b, x);
  
  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : cpBBDPrecSetupImpl
 * -----------------------------------------------------------------
 * cpBBDPrecSetupImpl generates a band-block-diagonal preconditioner
 * matrix, where the local block (on this processor) is a band
 * matrix. Each local block is computed by a difference quotient
 * scheme via calls to the user-supplied routines glocal, gcomm.
 * After generating the block in the band matrix savedP, this routine
 * does an LU factorization in place in savedP.
 *
 * The cpBBDPrecSetupImpl parameters used here are as follows:
 *
 * t is the current value of the independent variable t.
 *
 * yy is the current value of the dependent variable vector,
 *    namely the predicted value of y(t).
 *
 * yp is the current value of the derivative vector y',
 *    namely the predicted value of y'(t).
 *
 * gamma is the scalar in the system Jacobian, proportional to h.
 *
 * bbd_data is a pointer to user preconditioner data returned by
 *    CPBBDPrecInit.
 *
 * tmp1, tmp2, tmp3 are pointers to vectors of type N_Vector, 
 *    used for temporary storage or work space.
 *
 * Return value:
 * The value returned by cpBBDPrecSetupImpl function is a int
 * flag indicating whether it was successful. This value is
 *    0    if successful,
 *  > 0    for a recoverable error (step will be retried), or
 *  < 0    for a nonrecoverable error (step fails).
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSetupImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              realtype gamma, void *bbd_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int ier, retval;
  CPBBDPrecData pdata;
  CPodeMem cp_mem;

  pdata =(CPBBDPrecData) bbd_data;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Call cpBBDDQJacImpl for a new Jacobian calculation and store in savedP. */
  retval = cpBBDDQJacImpl(pdata, t, gamma, y, yp,
                          tmp1, tmp2, tmp3, pdata->tmp4);
  if (retval < 0) {
    cpProcessError(cp_mem, -1, "CPBBDPRE", "cpBBDPrecSetupImpl", MSGBBD_FUNC_FAILED);
    return(-1);
  }
  if (retval > 0) {
    return(+1);
  } 

  /* Do LU factorization of preconditioner block in place (in savedP). */
  ier = BandGBTRF(savedP, pivots);

  /* Return 0 if the LU was complete, or +1 otherwise. */
  if (ier > 0) return(+1);
  return(0);
}


/*
 * -----------------------------------------------------------------
 * Function: cpBBDPrecSolveImpl
 * -----------------------------------------------------------------
 * The function cpBBDPrecSolve computes a solution to the linear
 * system P x = b, where P is the left preconditioner defined by
 * the routine cpBBDPrecSetupImpl.
 *
 * The cpBBDPrecSolveImpl parameters used here are as follows:
 *
 * b is the input right-hand side vector r.
 *
 * x is the computed solution vector.
 *
 * bbd_data is a pointer to user preconditioner data returned
 *     by CPBBDPrecInit.
 *
 * The arguments t, y, yp, r, gamma, delta, and tmp are NOT used.
 *
 * cpBBDPrecSolveImpl always returns 0, indicating success.
 * -----------------------------------------------------------------
 */

static int cpBBDPrecSolveImpl(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                              N_Vector b, N_Vector x,
                              realtype gamma, realtype delta, 
                              void *bbd_data, N_Vector tmp)
{
  CPBBDPrecData pdata;
  realtype *xd;

  pdata = (CPBBDPrecData) bbd_data;

  /* Copy b to x, do the backsolve, and return. */
  N_VScale(ONE, b, x);

  xd = N_VGetArrayPointer(x);

  BandGBTRS(savedP, pivots, xd);

  return(0);
}

/* 
 * -----------------------------------------------------------------
 * Function: cpBBDPrecFree
 * -----------------------------------------------------------------
 * This is the pfree function (in the SPILS memory) defined by the
 * BBD preconditioner.
 */

static int cpBBDPrecFree(CPodeMem cp_mem)
{
  CPSpilsMem cpspils_mem;
  CPBBDPrecData pdata;
  
  if (cp_mem->cp_lmem == NULL) return(0);
  cpspils_mem = (CPSpilsMem) cp_mem->cp_lmem;
  
  if (cpspils_mem->s_P_data == NULL) return(0);
  pdata = (CPBBDPrecData) cpspils_mem->s_P_data;

  DestroyMat(savedJ);
  DestroyMat(savedP);
  DestroyArray(pivots);

  free(pdata);
  pdata = NULL;

  return(0);
}


/* 
 * =================================================================
 * DQ LOCAL JACOBIAN APROXIMATIONS
 * =================================================================
 */

#define ewt       (cp_mem->cp_ewt)
#define h         (cp_mem->cp_h)
#define user_data (cp_mem->cp_user_data)

/*
 * -----------------------------------------------------------------
 * Function: cpBBDDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the local block of the Jacobian of g(t,y). It assumes that a
 * band matrix of type BandMat is stored columnwise, and that elements
 * within each column are contiguous. All matrix elements are generated
 * as difference quotients, by way of calls to the user routine glocE.
 * By virtue of the band structure, the number of these calls is
 * bandwidth + 1, where bandwidth = mldq + mudq + 1.
 * But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 * This routine also assumes that the local elements of a vector are
 * stored contiguously.
 * -----------------------------------------------------------------
 */

static int cpBBDDQJacExpl(CPBBDPrecData pdata, realtype t, 
                          N_Vector y, N_Vector gy, 
                          N_Vector ytemp, N_Vector gtemp)
{
  CPodeMem cp_mem;
  realtype gnorm, minInc, inc, inc_inv;
  int group, i, j, width, ngroups, i1, i2;
  realtype *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;
  int retval;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and glocE to get base value of g(t,y) */
  if (cfn != NULL) {
    retval = cfn(Nlocal, t, y, NULL, user_data);
    if (retval != 0) return(retval);
  }

  retval = glocE(Nlocal, t, ytemp, gy, user_data);
  nge++;
  if (retval != 0) return(retval);

  /* Obtain pointers to the data for various vectors */
  y_data     =  N_VGetArrayPointer(y);
  gy_data    =  N_VGetArrayPointer(gy);
  ewt_data   =  N_VGetArrayPointer(ewt);
  ytemp_data =  N_VGetArrayPointer(ytemp);
  gtemp_data =  N_VGetArrayPointer(gtemp);

  /* Set minimum increment based on uround and norm of g */
  gnorm = N_VWrmsNorm(gy, ewt);
  minInc = (gnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(h) * uround * Nlocal * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mldq + mudq + 1;
  ngroups = SUNMIN(width, Nlocal);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = SUNMAX(dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    retval = glocE(Nlocal, t, ytemp, gtemp, user_data);
    nge++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < Nlocal; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = SUNMAX(dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mukeep);
      i2 = SUNMIN(j+mlkeep, Nlocal-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }

  return(0);
}


/*
 * -----------------------------------------------------------------
 * cpBBDDQJacImpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the local block of the Jacobian of G(t,y,y'). It assumes that
 * a band matrix of type BandMat is stored column-wise, and that
 * elements within each column are contiguous.
 *
 * All matrix elements are generated as difference quotients, by way
 * of calls to the user routine glocI. By virtue of the band
 * structure, the number of these calls is bandwidth + 1, where
 * bandwidth = mldq + mudq + 1. But the band matrix kept has
 * bandwidth = mlkeep + mukeep + 1. This routine also assumes that
 * the local elements of a vector are stored contiguously.
 *
 * Return values are: 0 (success), > 0 (recoverable error),
 * or < 0 (nonrecoverable error).
 * -----------------------------------------------------------------
 */

static int cpBBDDQJacImpl(CPBBDPrecData pdata, realtype t, realtype gamma,
                          N_Vector y, N_Vector yp, N_Vector gref, 
                          N_Vector ytemp, N_Vector yptemp, N_Vector gtemp)
{
  CPodeMem cp_mem;
  realtype inc, inc_inv;
  int  retval;
  int group, i, j, width, ngroups, i1, i2;
  realtype *ydata, *ypdata, *ytempdata, *yptempdata, *grefdata, *gtempdata;
  realtype *ewtdata;
  realtype *col_j, yj, ypj, ewtj;

  cp_mem = (CPodeMem) pdata->cpode_mem;

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, y, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Obtain pointers as required to the data array of vectors. */
  ydata     = N_VGetArrayPointer(y);
  ypdata    = N_VGetArrayPointer(yp);
  gtempdata = N_VGetArrayPointer(gtemp);
  ewtdata   = N_VGetArrayPointer(ewt);
  ytempdata = N_VGetArrayPointer(ytemp);
  yptempdata= N_VGetArrayPointer(yptemp);
  grefdata = N_VGetArrayPointer(gref);

  /* Call cfn and glocI to get base value of G(t,y,y'). */
  if (cfn != NULL) {
    retval = cfn(Nlocal, t, y, yp, user_data);
    if (retval != 0) return(retval);
  }

  retval = glocI(Nlocal, t, y, yp, gref, user_data); 
  nge++;
  if (retval != 0) return(retval);

  /* Set bandwidth and number of column groups for band differencing. */
  width = mldq + mudq + 1;
  ngroups = SUNMIN(width, Nlocal);

  /* Loop over groups. */
  for(group = 1; group <= ngroups; group++) {
    
    /* Loop over the components in this group. */
    for(j = group-1; j < Nlocal; j += width) {
      yj = ydata[j];
      ypj = ypdata[j];
      ewtj = ewtdata[j];
      
      /* Set increment inc to yj based on rel_yy*abs(yj), with
         adjustments using ypj and ewtj if this is small, and a further
         adjustment to give it the same sign as hh*ypj. */
      inc = dqrely*SUNMAX(SUNRabs(yj), SUNMAX( SUNRabs(h*ypj), ONE/ewtj));
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      
      /* Increment yj and ypj. */
      ytempdata[j]  += gamma*inc;
      yptempdata[j] += inc;
      
    }

    /* Evaluate G with incremented y and yp arguments. */
    retval = glocI(Nlocal, t, ytemp, yptemp, gtemp, user_data); 
    nge++;
    if (retval != 0) return(retval);

    /* Loop over components of the group again; restore ytemp and yptemp. */
    for(j = group-1; j < Nlocal; j += width) {
      yj  = ytempdata[j]  = ydata[j];
      ypj = yptempdata[j] = ypdata[j];
      ewtj = ewtdata[j];

      /* Set increment inc as before .*/
      inc = dqrely*SUNMAX(SUNRabs(yj), SUNMAX( SUNRabs(h*ypj), ONE/ewtj));
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;

      /* Form difference quotients and load into savedP. */
      inc_inv = ONE/inc;
      col_j = BAND_COL(savedP,j);
      i1 = SUNMAX(0, j-mukeep);
      i2 = SUNMIN(j+mlkeep, Nlocal-1);
      for(i = i1; i <= i2; i++) 
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (gtempdata[i] - grefdata[i]);
    }
  }
  
  return(0);
}
