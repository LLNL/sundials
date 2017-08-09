/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
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
 * Finally, it allocates memory for J and lpivots.
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

int IDADense(void *ida_mem, sunindextype Neq)
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
  if(IDA_mem->ida_tempv1->ops->nvgetarraypointer == NULL ||
     IDA_mem->ida_tempv1->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDADENSE", "IDADense", MSGD_BAD_NVECTOR);
    return(IDADLS_ILL_INPUT);
  }

  if (IDA_mem->ida_lfree != NULL) IDA_mem->ida_lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  IDA_mem->ida_linit  = IDADenseInit;
  IDA_mem->ida_lsetup = IDADenseSetup;
  IDA_mem->ida_lsolve = IDADenseSolve;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = IDADenseFree;

  /* Get memory for IDADlsMemRec. */
  idadls_mem = NULL;
  idadls_mem = (IDADlsMem) malloc(sizeof(struct IDADlsMemRec));
  if (idadls_mem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    return(IDADLS_MEM_FAIL);
  }

  /* Set matrix type */
  idadls_mem->d_type = SUNDIALS_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  idadls_mem->d_jacDQ   = TRUE;
  idadls_mem->d_djac    = NULL;
  idadls_mem->d_J_data = NULL;

  idadls_mem->d_last_flag = IDADLS_SUCCESS;

  idaDlsInitializeCounters(idadls_mem);

  IDA_mem->ida_setupNonNull = TRUE;

  /* Store problem size */
  idadls_mem->d_n = Neq;

  /* Allocate memory for J and pivot array. */
  idadls_mem->d_J = NULL;
  idadls_mem->d_J = NewDenseMat(Neq, Neq);
  if (idadls_mem->d_J == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }

  idadls_mem->d_lpivots = NULL;
  idadls_mem->d_lpivots = NewIndexArray(Neq);
  if (idadls_mem->d_lpivots == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    DestroyMat(idadls_mem->d_J);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }

  /* Attach linear solver memory to the integrator memory */
  IDA_mem->ida_lmem = idadls_mem;

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
  
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  idaDlsInitializeCounters(idadls_mem);
  /*
   nje   = 0;
   nreDQ = 0;
  */

  if (idadls_mem->d_jacDQ) {
    idadls_mem->d_djac = idaDlsDenseDQJac;
    idadls_mem->d_J_data = IDA_mem;
  } else {
    idadls_mem->d_J_data = IDA_mem->ida_user_data;
  }
  
  idadls_mem->d_last_flag = 0;
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
  sunindextype retfac;
  IDADlsMem idadls_mem;
  
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* Increment nje counter. */
  idadls_mem->d_nje++;

  /* Zero out J; call Jacobian routine jac; return if it failed. */
  SetToZero(idadls_mem->d_J);
  retval = idadls_mem->d_djac(idadls_mem->d_n, IDA_mem->ida_tn, IDA_mem->ida_cj, yyp, ypp, rrp, idadls_mem->d_J, idadls_mem->d_J_data, 
                tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDADLS_JACFUNC_UNRECVR, "IDADENSE", "IDADenseSetup", MSGD_JACFUNC_FAILED);
    idadls_mem->d_last_flag = IDADLS_JACFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    idadls_mem->d_last_flag = IDADLS_JACFUNC_RECVR;
    return(+1);
  }

  /* Do LU factorization of J; return success or fail flag. */
  retfac = DenseGETRF(idadls_mem->d_J, idadls_mem->d_lpivots);

  if (retfac != 0) {
    idadls_mem->d_last_flag = (long int) retfac;
    return(+1);
  }
  idadls_mem->d_last_flag = IDADLS_SUCCESS;
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
  
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(idadls_mem->d_J, idadls_mem->d_lpivots, bd);

  /* Scale the correction to account for change in cj. */
  if (IDA_mem->ida_cjratio != ONE) N_VScale(TWO/(ONE + IDA_mem->ida_cjratio), b, b);

  idadls_mem->d_last_flag = 0;
  return(0);
}

/*
  This routine frees memory specific to the IDADENSE linear solver.
*/

static int IDADenseFree(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;

  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;
  
  DestroyMat(idadls_mem->d_J);
  DestroyArray(idadls_mem->d_lpivots);
  free(IDA_mem->ida_lmem); IDA_mem->ida_lmem = NULL;

  return(0);
}
