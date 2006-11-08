/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:05 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CPBAND linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_band_impl.h"
#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

/* CPBAND linit, lsetup, lsolve, and lfree routines */
static int cpBandInit(CPodeMem cp_mem);
static int cpBandSetup(CPodeMem cp_mem, int convfail, 
                       N_Vector yP, N_Vector ypP, N_Vector fctP,
                        booleantype *jcurPtr,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                       N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpBandFree(CPodeMem cp_mem);

/* CPBAND DQ integration Jacobian functions */
static int cpBandDQJacExpl(long int N, long int mupper, long int mlower,
                           realtype t, N_Vector y, N_Vector fy, 
                           BandMat J, void *jac_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpBandDQJacImpl(long int N, long int mupper, long int mlower,
                           realtype t, realtype gm,
                           N_Vector y, N_Vector yp, N_Vector r,
                           BandMat J, void *jac_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define ode_type       (cp_mem->cp_ode_type)
#define lmm_type       (cp_mem->cp_lmm_type)
#define fe             (cp_mem->cp_fe)
#define fi             (cp_mem->cp_fi)
#define f_data         (cp_mem->cp_f_data)
#define uround         (cp_mem->cp_uround)
#define nst            (cp_mem->cp_nst)
#define tn             (cp_mem->cp_tn)
#define h              (cp_mem->cp_h)
#define gamma          (cp_mem->cp_gamma)
#define gammap         (cp_mem->cp_gammap)
#define gamrat         (cp_mem->cp_gamrat)
#define ewt            (cp_mem->cp_ewt)
#define tempv          (cp_mem->cp_tempv)

#define linit          (cp_mem->cp_linit)
#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)
#define lfree          (cp_mem->cp_lfree)
#define lmem           (cp_mem->cp_lmem)
#define lsetup_exists  (cp_mem->cp_lsetup_exists)

#define n              (cpband_mem->b_n)
#define jacE           (cpband_mem->b_jacE)
#define jacI           (cpband_mem->b_jacI)
#define M              (cpband_mem->b_M)
#define mu             (cpband_mem->b_mu)
#define ml             (cpband_mem->b_ml)
#define storage_mu     (cpband_mem->b_storage_mu)
#define pivots         (cpband_mem->b_pivots)
#define savedJ         (cpband_mem->b_savedJ)
#define nstlj          (cpband_mem->b_nstlj)
#define nje            (cpband_mem->b_nje)
#define nfeB           (cpband_mem->b_nfeB)
#define J_data         (cpband_mem->b_J_data)
#define last_flag      (cpband_mem->b_last_flag)

/*
 * -----------------------------------------------------------------
 * CPBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  CPBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cp_linit, cp_lsetup, cp_lsolve, and cp_lfree fields in (*cpode_mem)
 * to be CPBandInit, CPBandSetup, CPBandSolve, and CPBandFree,
 * respectively.  It allocates memory for a structure of type
 * CPBandMemRec and sets the cp_lmem field in (*cpode_mem) to the
 * address of this structure.  It sets lsetup_exists in (*cpode_mem) to be
 * TRUE, b_mu to be mupper, b_ml to be mlower, and the b_jac field to be 
 * CPBandDQJac.
 * Finally, it allocates memory for M, savedJ, and pivot.  The CPBand
 * return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int CPBand(void *cpode_mem, long int N,
           long int mupper, long int mlower)
{
  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPBAND_MEM_NULL, "CPBAND", "CPBand", MSGB_CPMEM_NULL);
    return(CPBAND_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, CPBAND_ILL_INPUT, "CPBAND", "CPBand", MSGB_BAD_NVECTOR);
    return(CPBAND_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */  
  linit  = cpBandInit;
  lsetup = cpBandSetup;
  lsolve = cpBandSolve;
  lfree  = cpBandFree;
  
  /* Get memory for CPBandMemRec */
  cpband_mem = NULL;
  cpband_mem = (CPBandMem) malloc(sizeof(CPBandMemRec));
  if (cpband_mem == NULL) {
    cpProcessError(cp_mem, CPBAND_MEM_FAIL, "CPBAND", "CPBand", MSGB_MEM_FAIL);
    return(CPBAND_MEM_FAIL);
  }
  
  /* Set default Jacobian routine and Jacobian data */
  jacE = NULL;
  jacI = NULL;
  J_data = NULL;

  last_flag = CPBAND_SUCCESS;
  lsetup_exists = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in cpband_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    cpProcessError(cp_mem, CPBAND_ILL_INPUT, "CPBAND", "CPBand", MSGB_BAD_SIZES);
    return(CPBAND_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  storage_mu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = BandAllocMat(N, mu, ml, storage_mu);
  if (M == NULL) {
    cpProcessError(cp_mem, CPBAND_MEM_FAIL, "CPBAND", "CPBand", MSGB_MEM_FAIL);
    free(cpband_mem);
    return(CPBAND_MEM_FAIL);
  }  
  pivots = BandAllocPiv(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPBAND_MEM_FAIL, "CPBAND", "CPBand", MSGB_MEM_FAIL);
    BandFreeMat(M);
    free(cpband_mem);
    return(CPBAND_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = BandAllocMat(N, mu, ml, mu);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPBAND_MEM_FAIL, "CPBAND", "CPBand", MSGB_MEM_FAIL);
      BandFreeMat(M);
      BandFreePiv(pivots);
      free(cpband_mem);
      return(CPBAND_MEM_FAIL);
    }
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cpband_mem;

  return(CPBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandSetJacFn
 * -----------------------------------------------------------------
 */

int CPBandSetJacFn(void *cpode_mem, void *bjac, void *jac_data)
{
  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPBAND_MEM_NULL, "CPBAND", "CPBandSetJacFn", MSGB_CPMEM_NULL);
    return(CPBAND_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPBAND_LMEM_NULL, "CPBAND", "CPBandSetJacFn", MSGB_LMEM_NULL);
    return(CPBAND_LMEM_NULL);
  }
  cpband_mem = (CPBandMem) lmem;

  if (ode_type == CP_EXPL) jacE = (CPBandJacExplFn) bjac;
  else                     jacI = (CPBandJacImplFn) bjac;
  J_data = jac_data;

  return(CPBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandGetWorkSpace
 * -----------------------------------------------------------------
 */

int CPBandGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS)
{
  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPBAND_MEM_NULL, "CPBAND", "CPBandGetWorkSpace", MSGB_CPMEM_NULL);
    return(CPBAND_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPBAND_LMEM_NULL, "CPBAND", "CPBandGetWorkSpace", MSGB_LMEM_NULL);
    return(CPBAND_LMEM_NULL);
  }
  cpband_mem = (CPBandMem) lmem;

  if (ode_type == CP_EXPL) *lenrwLS = n*(storage_mu + mu + 2*ml + 2);
  else                     *lenrwLS = n*(storage_mu + ml + 1);                     
  *leniwLS = n;

  return(CPBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandGetNumJacEvals
 * -----------------------------------------------------------------
 */

int CPBandGetNumJacEvals(void *cpode_mem, long int *njevals)
{
  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPBAND_MEM_NULL, "CPBAND", "CPBandGetNumJacEvals", MSGB_CPMEM_NULL);
    return(CPBAND_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPBAND_LMEM_NULL, "CPBAND", "CPBandGetNumJacEvals", MSGB_LMEM_NULL);
    return(CPBAND_LMEM_NULL);
  }
  cpband_mem = (CPBandMem) lmem;

  *njevals = nje;

  return(CPBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandGetNumFctEvals
 * -----------------------------------------------------------------
 */

int CPBandGetNumFctEvals(void *cpode_mem, long int *nfevalsLS)
{
  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPBAND_MEM_NULL, "CPBAND", "CPBandGetNumRhsEvals", MSGB_CPMEM_NULL);
    return(CPBAND_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPBAND_LMEM_NULL, "CPBAND", "CPBandGetNumRhsEvals", MSGB_LMEM_NULL);
    return(CPBAND_LMEM_NULL);
  }
  cpband_mem = (CPBandMem) lmem;

  *nfevalsLS = nfeB;

  return(CPBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandGetLastFlag
 * -----------------------------------------------------------------
 */

int CPBandGetLastFlag(void *cpode_mem, int *flag)
{
  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPBAND_MEM_NULL, "CPBAND", "CPBandGetLastFlag", MSGB_CPMEM_NULL);
    return(CPBAND_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPBAND_LMEM_NULL, "CPBAND", "CPBandGetLastFlag", MSGB_LMEM_NULL);
    return(CPBAND_LMEM_NULL);
  }
  cpband_mem = (CPBandMem) lmem;

  *flag = last_flag;

  return(CPBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPBandGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CPBandGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPBAND_SUCCESS:
    sprintf(name,"CPBAND_SUCCESS");
    break;   
  case CPBAND_MEM_NULL:
    sprintf(name,"CPBAND_MEM_NULL");
    break;
  case CPBAND_LMEM_NULL:
    sprintf(name,"CPBAND_LMEM_NULL");
    break;
  case CPBAND_ILL_INPUT:
    sprintf(name,"CPBAND_ILL_INPUT");
    break;
  case CPBAND_MEM_FAIL:
    sprintf(name,"CPBAND_MEM_FAIL");
    break;
  case CPBAND_JACFUNC_UNRECVR:
    sprintf(name,"CPBAND_JACFUNC_UNRECVR");
    break;
  case CPBAND_JACFUNC_RECVR:
    sprintf(name,"CPBAND_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * cpBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpBandInit(CPodeMem cp_mem)
{
  CPBandMem cpband_mem;

  cpband_mem = (CPBandMem) lmem;

  nje   = 0;
  nfeB  = 0;
  nstlj = 0;

  if (ode_type == CP_EXPL && jacE == NULL) {
    jacE = cpBandDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && jacI == NULL) {
    jacI = cpBandDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPBAND_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the band LU factorization routine.
 * -----------------------------------------------------------------
 */

static int cpBandSetup(CPodeMem cp_mem, int convfail, 
                       N_Vector yP, N_Vector ypP, N_Vector fctP,
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  CPBandMem cpband_mem;
  int retval;

  cpband_mem = (CPBandMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CPB_MSBJ) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPB_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    jok = !jbad;
    
    if (jok) {
      
      /* If jok = TRUE, use saved copy of J */
      *jcurPtr = FALSE;
      BandCopy(savedJ, M, mu, ml);
      
    } else {
      
      /* If jok = FALSE, call jac routine for new J value */
      nje++;
      nstlj = nst;
      *jcurPtr = TRUE;

      BandZero(M);
      retval = jacE(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      if (retval == 0) {
        BandCopy(M, savedJ, mu, ml);
      }else if (retval < 0) {
        cpProcessError(cp_mem, CPBAND_JACFUNC_UNRECVR, "CPBAND", "CPBandSetup", MSGB_JACFUNC_FAILED);
        last_flag = CPBAND_JACFUNC_UNRECVR;
        return(-1);
      }else if (retval > 0) {
        last_flag = CPBAND_JACFUNC_RECVR;
        return(1);
      }

    }
  
    /* Scale and add I to get M = I - gamma*J */
    BandScale(-gamma, M);
    BandAddI(M);

    break;

  case CP_IMPL:

    BandZero(M);
    retval = jacI(n, mu, ml, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval < 0) {
      cpProcessError(cp_mem, CPBAND_JACFUNC_UNRECVR, "CPBAND", "CPBandSetup", MSGB_JACFUNC_FAILED);
      last_flag = CPBAND_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CPBAND_JACFUNC_RECVR;
      return(+1);
    }

    break;

  }

  /* Do LU factorization of M */
  ier = BandGBTRF(M, pivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier>0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int cpBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                       N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPBandMem cpband_mem;
  realtype *bd;

  cpband_mem = (CPBandMem) lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(M, pivots, bd);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }

  last_flag = CPBAND_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void cpBandFree(CPodeMem cp_mem)
{
  CPBandMem cpband_mem;

  cpband_mem = (CPBandMem) lmem;

  BandFreeMat(M);
  BandFreePiv(pivots);
  if (ode_type == CP_EXPL) BandFreeMat(savedJ);
  free(cpband_mem); 
  cpband_mem = NULL;
}

/* 
 * =================================================================
 *  DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpBandDQJacExpl
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y).  It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

static int cpBandDQJacExpl(long int N, long int mupper, long int mlower,
                           realtype t, N_Vector y, N_Vector fy, 
                           BandMat J, void *jac_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, srur;
  N_Vector ftemp, ytemp;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval = 0;

  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* jac_dat points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpband_mem = (CPBandMem) lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    retval = fe(tn, ytemp, ftemp, f_data);
    nfeB++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(J,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  return(retval);

}


/*
 */
static int cpBandDQJacImpl(long int N, long int mupper, long int mlower,
                           realtype t, realtype gm,
                           N_Vector y, N_Vector yp, N_Vector r,
                           BandMat J, void *jac_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, ewtj;
  realtype *y_data, *yp_data, *ewt_data;
  realtype *ytemp_data, *yptemp_data, *ftemp_data, *r_data, *col_j;
  int group;
  
  N_Vector ftemp, ytemp, yptemp;
  long int i, j, i1, i2, width, ngroups;
  int retval = 0;

  CPodeMem cp_mem;
  CPBandMem cpband_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpband_mem = (CPBandMem) lmem;

  ftemp = tmp1; /* Rename work vector for use as the perturbed residual. */
  ytemp = tmp2; /* Rename work vector for use as a temporary for yy. */
  yptemp= tmp3; /* Rename work vector for use as a temporary for yp. */

  /* Obtain pointers to the data for all vectors used.  */
  ewt_data    = N_VGetArrayPointer(ewt);
  r_data      = N_VGetArrayPointer(r);
  y_data      = N_VGetArrayPointer(y);
  yp_data     = N_VGetArrayPointer(yp);
  ftemp_data  = N_VGetArrayPointer(ftemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, y, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */
  srur = RSqrt(uround);
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y[j] and yp[j] for j in this group. */
    for (j=group-1; j<N; j+=width) {
        yj = y_data[j];
        ypj = yp_data[j];
        ewtj = ewt_data[j];

        /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
           adjustments using ypj and ewtj if this is small, and a further
           adjustment to give it the same sign as h*ypj. */
        inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );

        if (h*ypj < ZERO) inc = -inc;
        inc = (yj + inc) - yj;

        /* Increment yj and ypj. */
        ytemp_data[j]  += gamma*inc;
        yptemp_data[j] += inc;
    }

    /* Call ODE fct. with incremented arguments. */
    retval = fi(tn, ytemp, yptemp, ftemp, f_data);
    nfeB++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */
    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */
      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = BAND_COL(J, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */
      inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      
      /* Load the difference quotient Jacobian elements for column j. */
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower,N-1);
      for (i=i1; i<=i2; i++) 
        BAND_COL_ELEM(col_j,i,j) = inc_inv*(ftemp_data[i]-r_data[i]);
    }
    
  }
  
  return(retval);

}
