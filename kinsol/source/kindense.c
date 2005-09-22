/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2005-09-22 23:12:33 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the KINDENSE linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "kindense_impl.h"
#include "kinsol_impl.h"

#include "sundialsmath.h"

/* Other Constants */

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* KINDENSE linit, lsetup, lsolve, and lfree routines */
 
static int KINDenseInit(KINMem kin_mem);

static int KINDenseSetup(KINMem kin_mem);

static int KINDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                         realtype *res_norm);

static int KINDenseFree(KINMem kin_mem);

/* KINDENSE DQJac routine */

static int KINDenseDQJac(long int n, DenseMat J,
                         N_Vector u, N_Vector fu, void *jac_data,
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

#define n         (kindense_mem->d_n)
#define jac       (kindense_mem->d_jac)
#define J         (kindense_mem->d_J)
#define pivots    (kindense_mem->d_pivots)
#define nje       (kindense_mem->d_nje)
#define nfeD      (kindense_mem->d_nfeD)
#define J_data    (kindense_mem->d_J_data)
#define last_flag (kindense_mem->d_last_flag)
                  
/*
 * -----------------------------------------------------------------
 * KINDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module. 
 * KINDense sets the kin_linit, kin_lsetup, kin_lsolve, kin_lfree fields 
 * in *kinmem to be KINDenseInit, KINDenseSetup, KINDenseSolve, and 
 * KINDenseFree, respectively.  
 * It allocates memory for a structure of type KINDenseMemRec and sets 
 * the kin_lmem field in *kinmem to the address of this structure.  
 * It sets setupNonNull in *kinmem to TRUE, and the d_jac field to the 
 * default KINDenseDQJac.
 * Finally, it allocates memory for J and pivots.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int KINDense(void *kinmem, long int N)
{
  KINMem kin_mem;
  KINDenseMem kindense_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stderr, MSGDS_KINMEM_NULL);
    return(KINDENSE_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_BAD_NVECTOR);
    return(KINDENSE_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */
  linit  = KINDenseInit;
  lsetup = KINDenseSetup;
  lsolve = KINDenseSolve;
  lfree  = KINDenseFree;

  /* Get memory for KINDenseMemRec */
  kindense_mem = (KINDenseMem) malloc(sizeof(KINDenseMemRec));
  if (kindense_mem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_MEM_FAIL);
    return(KINDENSE_MEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  jac = KINDenseDQJac;
  J_data = kin_mem;
  last_flag = KINDENSE_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for J and pivot array */
  
  J = DenseAllocMat(N);
  if (J == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_MEM_FAIL);
    return(KINDENSE_MEM_FAIL);
  }
  pivots = DenseAllocPiv(N);
  if (pivots == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_MEM_FAIL);
    DenseFreeMat(J);
    return(KINDENSE_MEM_FAIL);
  }

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kindense_mem;

  return(KINDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINDenseSetJacFn
 * -----------------------------------------------------------------
 */

int KINDenseSetJacFn(void *kinmem, KINDenseJacFn djac, void *jac_data)
{
  KINMem kin_mem;
  KINDenseMem kindense_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stderr, MSGDS_SETGET_KINMEM_NULL);
    return(KINDENSE_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_SETGET_LMEM_NULL);
    return(KINDENSE_LMEM_NULL);
  }
  kindense_mem = (KINDenseMem) lmem;

  jac = djac;
  if (djac != NULL) J_data = jac_data;

  return(KINDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINDenseGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINDenseGetWorkSpace(void *kinmem, long int *lenrwD, long int *leniwD)
{
  KINMem kin_mem;
  KINDenseMem kindense_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stderr, MSGDS_SETGET_KINMEM_NULL);
    return(KINDENSE_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_SETGET_LMEM_NULL);
    return(KINDENSE_LMEM_NULL);
  }
  kindense_mem = (KINDenseMem) lmem;

  *lenrwD = 2*n*n;
  *leniwD = n;

  return(KINDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINDenseGetNumJacEvals
 * -----------------------------------------------------------------
 */

int KINDenseGetNumJacEvals(void *kinmem, long int *njevalsD)
{
  KINMem kin_mem;
  KINDenseMem kindense_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stderr, MSGDS_SETGET_KINMEM_NULL);
    return(KINDENSE_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_SETGET_LMEM_NULL);
    return(KINDENSE_LMEM_NULL);
  }
  kindense_mem = (KINDenseMem) lmem;

  *njevalsD = nje;

  return(KINDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINDenseGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINDenseGetNumFuncEvals(void *kinmem, long int *nfevalsD)
{
  KINMem kin_mem;
  KINDenseMem kindense_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stderr, MSGDS_SETGET_KINMEM_NULL);
    return(KINDENSE_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_SETGET_LMEM_NULL);
    return(KINDENSE_LMEM_NULL);
  }
  kindense_mem = (KINDenseMem) lmem;

  *nfevalsD = nfeD;

  return(KINDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINDenseGetLastFlag
 * -----------------------------------------------------------------
 */

int KINDenseGetLastFlag(void *kinmem, int *flag)
{
   KINMem kin_mem;
  KINDenseMem kindense_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stderr, MSGDS_SETGET_KINMEM_NULL);
    return(KINDENSE_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSGDS_SETGET_LMEM_NULL);
    return(KINDENSE_LMEM_NULL);
  }
  kindense_mem = (KINDenseMem) lmem;

  *flag = last_flag;

  return(KINDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int KINDenseInit(KINMem kin_mem)
{
  KINDenseMem kindense_mem;

  kindense_mem = (KINDenseMem) lmem;
  
  nje   = 0;
  nfeD  = 0;
  
  if (jac == NULL) {
    jac = KINDenseDQJac;
    J_data = kin_mem;
  }

  last_flag = KINDENSE_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It calls the dense LU factorization routine.
 * -----------------------------------------------------------------
 */

static int KINDenseSetup(KINMem kin_mem)
{
  KINDenseMem kindense_mem;
  long int ier;
  
  kindense_mem = (KINDenseMem) lmem;
 
  nje++;
  DenseZero(J); 
  jac(n, J, uu, fval, J_data, vtemp1, vtemp2);

  /* Do LU factorization of J */
  ier = DenseFactor(J, pivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int KINDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINDenseMem kindense_mem;
  realtype *xd;

  kindense_mem = (KINDenseMem) lmem;

  /* Copy the right-hand side into x */

  N_VScale(ONE, b, x);
  
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */
  
  DenseBacksolve(J, pivots, xd);

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

  last_flag = KINDENSE_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static int KINDenseFree(KINMem kin_mem)
{
  KINDenseMem  kindense_mem;

  kindense_mem = (KINDenseMem) lmem;
  
  DenseFreeMat(J);
  DenseFreePiv(pivots);
  free(kindense_mem);
  
  return(0);
}

/*
 * -----------------------------------------------------------------
 * KINDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of F(u). It assumes that a dense matrix of type
 * DenseMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 *
 * The increment used in the finitie-difference approximation
 *   J_ij = ( F_i(u+sigma_j * e_j) - F_i(u)  ) / sigma_j
 * is
 *  sigma_j = max{|u_j|, |1/uscale_j|} * sqrt(uround)
 *
 * Note: uscale_j = 1/typ(u_j)
 * -----------------------------------------------------------------
 */

#undef n
#undef J
 
static int KINDenseDQJac(long int n, DenseMat J,
                         N_Vector u, N_Vector fu, void *jac_data,
                         N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv, ujsaved, ujscale, sign;
  realtype *tmp2_data, *u_data, *uscale_data;
  N_Vector ftemp, jthCol;
  long int j;

  KINMem kin_mem;
  KINDenseMem  kindense_mem;

  /* jac_data points to kin_mem */
  kin_mem = (KINMem) jac_data;
  kindense_mem = (KINDenseMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for u and uscale */
  u_data   = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(uscale);

  /* This is the only for loop for 0..N-1 in KINSOL */

  for (j = 0; j < n; j++) {

    /* Generate the jth col of J(u) */

    N_VSetArrayPointer(DENSE_COL(J,j), jthCol);

    ujsaved = u_data[j];
    ujscale = ONE/uscale_data[j];
    sign = (ujsaved >= ZERO) ? ONE : -ONE;
    inc = sqrt_relfunc*MAX(ABS(ujsaved), ujscale)*sign;
    u_data[j] += inc;
    func(u, ftemp, f_data);
    u_data[j] = ujsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fu, jthCol);

  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  /* Increment counter nfeD */
  nfeD += n;

  return(0);
}
