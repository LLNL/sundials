/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-10-11 16:34:19 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the KINBAND linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "kinsol_band_impl.h"
#include "kinsol_impl.h"

#include <sundials/sundials_math.h>

/* Other Constants */

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* KINBAND linit, lsetup, lsolve, and lfree routines */

static int KINBandInit(KINMem kin_mem);

static int KINBandSetup(KINMem kin_mem);

static int KINBandSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                        realtype *res_norm);

static void KINBandFree(KINMem kin_mem);

/* KINBAND DQJac routine */

static int KINBandDQJac(long int N, long int mupper, long int mlower,
                        BandMat J, N_Vector u, N_Vector fu, void *jac_data,
                        N_Vector tmp1, N_Vector tmp2);

/* Readability Replacements */

#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define uround         (kin_mem->kin_uround)
#define nfe            (kin_mem->kin_nfe)
#define nni            (kin_mem->kin_nni)
#define nnilset        (kin_mem->kin_nnilset)
#define func           (kin_mem->kin_func)
#define f_data         (kin_mem->kin_f_data)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)

#define n          (kinband_mem->b_n)
#define jac        (kinband_mem->b_jac)
#define J          (kinband_mem->b_J)
#define mu         (kinband_mem->b_mu)
#define ml         (kinband_mem->b_ml)
#define storage_mu (kinband_mem->b_storage_mu)
#define pivots     (kinband_mem->b_pivots)
#define nje        (kinband_mem->b_nje)
#define nfeB       (kinband_mem->b_nfeB)
#define J_data     (kinband_mem->b_J_data)
#define last_flag  (kinband_mem->b_last_flag)

/*
 * -----------------------------------------------------------------
 * KINBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  KINBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be KINBandInit, KINBandSetup, KINBandSolve, and KINBandFree,
 * respectively.  It allocates memory for a structure of type
 * KINBandMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to be
 * TRUE, b_mu to be mupper, b_ml to be mlower, and the b_jac field to be 
 * KINBandDQJac.
 * Finally, it allocates memory for M, savedJ, and pivot.  The KINBand
 * return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int KINBand(void *kinmem, long int N,
           long int mupper, long int mlower)
{
  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINBAND_MEM_NULL, "KINBAND", "KINBand", MSGB_KINMEM_NULL);
    return(KINBAND_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, KINBAND_ILL_INPUT, "KINBAND", "KINBand", MSGB_BAD_NVECTOR);
    return(KINBAND_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */  
  linit  = KINBandInit;
  lsetup = KINBandSetup;
  lsolve = KINBandSolve;
  lfree  = KINBandFree;
  
  /* Get memory for KINBandMemRec */
  kinband_mem = NULL;
  kinband_mem = (KINBandMem) malloc(sizeof(KINBandMemRec));
  if (kinband_mem == NULL) {
    KINProcessError(kin_mem, KINBAND_MEM_FAIL, "KINBAND", "KINBand", MSGB_MEM_FAIL);
    return(KINBAND_MEM_FAIL);
  }
  
  /* Set default Jacobian routine and Jacobian data */
  jac = KINBandDQJac;
  J_data = kinmem;
  last_flag = KINBAND_SUCCESS;

  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in kinband_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    KINProcessError(kin_mem, KINBAND_MEM_FAIL, "KINBAND", "KINBand", MSGB_MEM_FAIL);
    return(KINBAND_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  storage_mu = MIN(N-1, mu + ml);

  /* Allocate memory for J and pivot array */
  J = NULL;
  J = BandAllocMat(N, mu, ml, storage_mu);
  if (J == NULL) {
    KINProcessError(kin_mem, KINBAND_MEM_FAIL, "KINBAND", "KINBand", MSGB_MEM_FAIL);
    free(kinband_mem); kinband_mem = NULL;
    return(KINBAND_MEM_FAIL);
  }

  pivots = NULL;
  pivots = BandAllocPiv(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINBAND_MEM_FAIL, "KINBAND", "KINBand", MSGB_MEM_FAIL);
    BandFreeMat(J);
    free(kinband_mem); kinband_mem = NULL;
    return(KINBAND_MEM_FAIL);
  }

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kinband_mem;

  return(KINBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandSetJacFn
 * -----------------------------------------------------------------
 */

int KINBandSetJacFn(void *kinmem, KINBandJacFn bjac, void *jac_data)
{
  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINBAND_MEM_NULL, "KINBAND", "KINBandSetJacFn", MSGB_KINMEM_NULL);
    return(KINBAND_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINBAND_LMEM_NULL, "KINBAND", "KINBandSetJacFn", MSGB_LMEM_NULL);
    return(KINBAND_LMEM_NULL);
  }
  kinband_mem = (KINBandMem) lmem;

  jac = bjac;
  if (bjac != NULL) J_data = jac_data;

  return(KINBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINBandGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB)
{
  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINBAND_MEM_NULL, "KINBAND", "KINBandGetWorkSpace", MSGB_KINMEM_NULL);
    return(KINBAND_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINBAND_LMEM_NULL, "KINBAND", "KINBandGetWorkSpace", MSGB_LMEM_NULL);
    return(KINBAND_LMEM_NULL);
  }
  kinband_mem = (KINBandMem) lmem;

  *lenrwB = n*(storage_mu + mu + 2*ml + 2);
  *leniwB = n;

  return(KINBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandGetNumJacEvals
 * -----------------------------------------------------------------
 */

int KINBandGetNumJacEvals(void *kinmem, long int *njevalsB)
{
  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINBAND_MEM_NULL, "KINBAND", "KINBandGetNumJacEvals", MSGB_KINMEM_NULL);
    return(KINBAND_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINBAND_LMEM_NULL, "KINBAND", "KINBandGetNumJacEvals", MSGB_LMEM_NULL);
    return(KINBAND_LMEM_NULL);
  }
  kinband_mem = (KINBandMem) lmem;

  *njevalsB = nje;

  return(KINBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINBandGetNumFuncEvals(void *kinmem, long int *nfevalsB)
{
  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINBAND_MEM_NULL, "KINBAND", "KINBandGetNumFuncEvals", MSGB_KINMEM_NULL);
    return(KINBAND_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINBAND_LMEM_NULL, "KINBAND", "KINBandGetNumGuncEvals", MSGB_LMEM_NULL);
    return(KINBAND_LMEM_NULL);
  }
  kinband_mem = (KINBandMem) lmem;

  *nfevalsB = nfeB;

  return(KINBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandGetLastFlag
 * -----------------------------------------------------------------
 */

int KINBandGetLastFlag(void *kinmem, int *flag)
{
  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINBAND_MEM_NULL, "KINBAND", "KINBandGetLastFlag", MSGB_KINMEM_NULL);
    return(KINBAND_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINBAND_LMEM_NULL, "KINBAND", "KINBandGetLastFlag", MSGB_LMEM_NULL);
    return(KINBAND_LMEM_NULL);
  }
  kinband_mem = (KINBandMem) lmem;

  *flag = last_flag;

  return(KINBAND_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINBandGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *KINBandGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINBAND_SUCCESS:
    sprintf(name, "KINBAND_SUCCESS");
    break;
  case KINBAND_MEM_NULL:
    sprintf(name, "KINBAND_MEM_NULL");
    break;
  case KINBAND_LMEM_NULL:
    sprintf(name, "KINBAND_LMEM_NULL");
    break;
  case KINBAND_ILL_INPUT:
    sprintf(name, "KINBAND_ILL_INPUT");
    break;
  case KINBAND_MEM_FAIL:
    sprintf(name, "KINBAND_MEM_FAIL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * KINBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int KINBandInit(KINMem kin_mem)
{
  KINBandMem kinband_mem;

  kinband_mem = (KINBandMem) lmem;

  nje   = 0;
  nfeB  = 0;

  if (jac == NULL) {
    jac = KINBandDQJac;
    J_data = kin_mem;
  }

  last_flag = KINBAND_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int KINBandSetup(KINMem kin_mem)
{
  KINBandMem kinband_mem;
  long int ier;
  int retval;

  kinband_mem = (KINBandMem) lmem;

  nje++;
  BandZero(J); 
  retval = jac(n, mu, ml, J, uu, fval, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }
  
  /* Do LU factorization of J */
  ier = BandGBTRF(J, pivots);

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int KINBandSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINBandMem kinband_mem;
  realtype *xd;

  kinband_mem = (KINBandMem) lmem;

  /* Copy the right-hand side into x */

  N_VScale(ONE, b, x);
  
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */

  BandGBTRS(J, pivots, xd);

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.
     
     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */

  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);

  last_flag = KINBAND_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void KINBandFree(KINMem kin_mem)
{
  KINBandMem kinband_mem;

  kinband_mem = (KINBandMem) lmem;

  BandFreeMat(J);
  BandFreePiv(pivots);
  free(kinband_mem); kinband_mem = NULL;
}

/*
 * -----------------------------------------------------------------
 * KINBandDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of F(u).  It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 *
 * NOTE: Any type of failure of the system function her leads to an
 *       unrecoverable failure of the Jacobian function and thus
 *       of the linear solver setup function, stopping KINSOL.
 * -----------------------------------------------------------------
 */

#undef n
#undef J

static int KINBandDQJac(long int n, long int mupper, long int mlower,
                        BandMat J, N_Vector u, N_Vector fu, void *jac_data,
                        N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv;
  N_Vector futemp, utemp;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *fu_data, *futemp_data, *u_data, *utemp_data, *uscale_data;
  int retval;

  KINMem kin_mem;
  KINBandMem kinband_mem;

  /* jac_dat points to kinmem */
  kin_mem = (KINMem) jac_data;
  kinband_mem = (KINBandMem) lmem;

  /* Rename work vectors for use as temporary values of u and fu */
  futemp = tmp1;
  utemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, futemp, y, ytemp */
  fu_data    = N_VGetArrayPointer(fu);
  futemp_data = N_VGetArrayPointer(futemp);
  u_data     = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(uscale);
  utemp_data = N_VGetArrayPointer(utemp);

  /* Load utemp with u */
  N_VScale(ONE, u, utemp);

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, n);
  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all utemp components in group */
    for(j=group-1; j < n; j+=width) {
      inc = sqrt_relfunc*MAX(ABS(u_data[j]), ABS(uscale_data[j]));
      utemp_data[j] += inc;
    }

    /* Evaluate f with incremented u */
    retval = func(utemp, futemp, f_data);
    if (retval != 0) return(-1); 

    /* Restore utemp components, then form and load difference quotients */
    for (j=group-1; j < n; j+=width) {
      utemp_data[j] = u_data[j];
      col_j = BAND_COL(J,j);
      inc = sqrt_relfunc*MAX(ABS(u_data[j]), ABS(uscale_data[j]));
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, n-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (futemp_data[i] - fu_data[i]);
    }
  }
  
  /* Increment counter nfeB */
  nfeB += ngroups;

  return(0);
}
