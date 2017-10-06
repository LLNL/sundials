/*
 * -----------------------------------------------------------------
 * $Revision: 4951 $
 * $Date: 2016-09-22 10:21:00 -0700 (Thu, 22 Sep 2016) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the implementation file for the IDADENSE linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <ida/ida_dense.h>
#include "ida_direct_impl.h"
#include "ida_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* IDADENSE linit, lsetup, lsolve, and lfree routines */
 
static int IDADenseInit(IDAMem IDA_mem);

static int IDADenseSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector rrp, N_Vector tmp1,
                         N_Vector tmp2, N_Vector tmp3);

static int IDADenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                         N_Vector ycur, N_Vector ypcur, N_Vector rrcur);

static int IDADenseFree(IDAMem IDA_mem);

/* Readability Replacements */

#define res          (IDA_mem->ida_res)
#define tn           (IDA_mem->ida_tn)
#define hh           (IDA_mem->ida_hh)
#define cj           (IDA_mem->ida_cj)
#define cjratio      (IDA_mem->ida_cjratio)
#define ewt          (IDA_mem->ida_ewt)
#define constraints  (IDA_mem->ida_constraints)
#define linit        (IDA_mem->ida_linit)
#define lsetup       (IDA_mem->ida_lsetup)
#define lsolve       (IDA_mem->ida_lsolve)
#define lperf        (IDA_mem->ida_lperf)
#define lfree        (IDA_mem->ida_lfree)
#define lmem         (IDA_mem->ida_lmem)
#define setupNonNull (IDA_mem->ida_setupNonNull)
#define vec_tmpl     (IDA_mem->ida_tempv1)

#define mtype        (idadls_mem->d_type)
#define neq          (idadls_mem->d_n)
#define jacDQ        (idadls_mem->d_jacDQ)
#define djac         (idadls_mem->d_djac)
#define JJ           (idadls_mem->d_J)
#define lpivots      (idadls_mem->d_lpivots)
#define nje          (idadls_mem->d_nje)
#define nreDQ        (idadls_mem->d_nreDQ)
#define jacdata      (idadls_mem->d_J_data)
#define last_flag    (idadls_mem->d_last_flag)
                  
/*
 * -----------------------------------------------------------------
 * IDADense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDADENSE linear solver module.  
 * IDADense first calls the existing lfree routine if this is not NULL.
 * Then it sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDADenseInit, IDADenseSetup,
 * IDADenseSolve, NULL, and IDADenseFree, respectively.
 * It allocates memory for a structure of type IDADlsMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem) to TRUE, sets the d_jdata field
 * in the IDADlsMemRec structure to be the input parameter jdata,
 * and sets the d_jac field to be:
 *   (1) the input parameter djac, if djac != NULL, or                
 *   (2) IDADenseDQJac, if djac == NULL.                             
 * Finally, it allocates memory for JJ and lpivots.
 * The return value is IDADLS_SUCCESS = 0, IDADLS_LMEM_FAIL = -1,
 * or IDADLS_ILL_INPUT = -2.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDADense will first 
 *       test for a compatible N_Vector internal
 *       representation by checking that the functions N_VGetArrayPointer
 *       and N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int IDADense(void *ida_mem, long int Neq)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADENSE", "IDADense", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if(vec_tmpl->ops->nvgetarraypointer == NULL ||
     vec_tmpl->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDADENSE", "IDADense", MSGD_BAD_NVECTOR);
    return(IDADLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  linit  = IDADenseInit;
  lsetup = IDADenseSetup;
  lsolve = IDADenseSolve;
  lperf  = NULL;
  lfree  = IDADenseFree;

  /* Get memory for IDADlsMemRec. */
  idadls_mem = NULL;
  idadls_mem = (IDADlsMem) malloc(sizeof(struct IDADlsMemRec));
  if (idadls_mem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    return(IDADLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  jacDQ   = TRUE;
  djac    = NULL;
  jacdata = NULL;

  last_flag = IDADLS_SUCCESS;

  idaDlsInitializeCounters(idadls_mem);

  setupNonNull = TRUE;

  /* Store problem size */
  neq = Neq;

  /* Allocate memory for JJ and pivot array. */
  JJ = NULL;
  JJ = NewDenseMat(Neq, Neq);
  if (JJ == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }

  lpivots = NULL;
  lpivots = NewLintArray(Neq);
  if (lpivots == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    DestroyMat(JJ);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }

  /* Attach linear solver memory to the integrator memory */
  lmem = idadls_mem;

  return(IDADLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDADENSE interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the IDADENSE
  linear solver module.  It returns 0.
*/

static int IDADenseInit(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;
  
  idadls_mem = (IDADlsMem) lmem;

  idaDlsInitializeCounters(idadls_mem);
  /*
   nje   = 0;
   nreDQ = 0;
  */

  if (jacDQ) {
    djac = idaDlsDenseDQJac;
    jacdata = IDA_mem;
  } else {
    jacdata = IDA_mem->ida_user_data;
  }
  
  last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the IDADENSE linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the dense LU factorization routine.
  The return value is either
     IDADLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int IDADenseSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector rrp, N_Vector tmp1, N_Vector tmp2,
                         N_Vector tmp3)
{
  int retval;
  long int retfac;
  IDADlsMem idadls_mem;
  
  idadls_mem = (IDADlsMem) lmem;

  /* Increment nje counter. */
  nje++;

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  SetToZero(JJ);
  retval = djac(neq, tn, cj, yyp, ypp, rrp, JJ, jacdata, 
                tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDADLS_JACFUNC_UNRECVR, "IDADENSE", "IDADenseSetup", MSGD_JACFUNC_FAILED);
    last_flag = IDADLS_JACFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    last_flag = IDADLS_JACFUNC_RECVR;
    return(+1);
  }

  /* Do LU factorization of JJ; return success or fail flag. */
  retfac = DenseGETRF(JJ, lpivots);

  if (retfac != 0) {
    last_flag = retfac;
    return(+1);
  }
  last_flag = IDADLS_SUCCESS;
  return(0);
}

/*
  This routine handles the solve operation for the IDADENSE linear
  solver module.  It calls the dense backsolve routine, scales the
  solution vector according to cjratio, then returns IDADLS_SUCCESS = 0.
*/

static int IDADenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                         N_Vector ycur, N_Vector ypcur, N_Vector rrcur)
{
  IDADlsMem idadls_mem;
  realtype *bd;
  
  idadls_mem = (IDADlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(JJ, lpivots, bd);

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), b, b);

  last_flag = 0;
  return(0);
}

/*
  This routine frees memory specific to the IDADENSE linear solver.
*/

static int IDADenseFree(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;

  idadls_mem = (IDADlsMem) lmem;
  
  DestroyMat(JJ);
  DestroyArray(lpivots);
  free(lmem); lmem = NULL;

  return(0);
}
