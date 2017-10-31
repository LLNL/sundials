/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Radu Serban @ LLNL
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
 * This is the implementation file for the KINDLS linear solver
 * interface
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "kinsol_impl.h"
#include "kinsol_direct_impl.h"
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>


/*==================================================================
  FUNCTION SPECIFIC CONSTANTS
  ==================================================================*/

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/*==================================================================
  KINDLS Exported fuctions -- Required
  ==================================================================*/

/*------------------------------------------------------------------
  KINDlsSetLinearSolver speficies the direct linear solver
  ------------------------------------------------------------------*/
int KINDlsSetLinearSolver(void *kinmem, SUNLinearSolver LS,
                          SUNMatrix A)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if any input is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS",
                    "KINDlsSetLinearSolver", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  if ( (LS == NULL) || (A == NULL) ) {
    KINProcessError(NULL, KINDLS_ILL_INPUT, "KINDLS", 
                    "KINDlsSetLinearSolver",
                    "Both LS and A must be non-NULL");
    return(KINDLS_ILL_INPUT);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if solver and vector are compatible with DLS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_DIRECT) {
    KINProcessError(kin_mem, KINDLS_ILL_INPUT, "KINDLS", 
                    "KINDlsSetLinearSolver", 
                    "Non-direct LS supplied to KINDls interface");
    return(KINDLS_ILL_INPUT);
  }
  if (kin_mem->kin_vtemp1->ops->nvgetarraypointer == NULL ||
      kin_mem->kin_vtemp1->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINDLS_ILL_INPUT, "KINDLS", 
                    "KINDlsSetLinearSolver", MSGD_BAD_NVECTOR);
    return(KINDLS_ILL_INPUT);
  }

  /* free any existing system solver attached to KINSOL */
  if (kin_mem->kin_lfree)  kin_mem->kin_lfree(kin_mem);

  /* Set four main system linear solver function fields in kin_mem */
  kin_mem->kin_linit  = kinDlsInitialize;
  kin_mem->kin_lsetup = kinDlsSetup;
  kin_mem->kin_lsolve = kinDlsSolve;
  kin_mem->kin_lfree  = kinDlsFree;

  /* Get memory for KINDlsMemRec */
  kindls_mem = NULL;
  kindls_mem = (KINDlsMem) malloc(sizeof(struct KINDlsMemRec));
  if (kindls_mem == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINDLS",
                    "KINDlsSetLinearSolver", MSGD_MEM_FAIL);
    return(KINDLS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  kindls_mem->LS = LS;

  /* set SUNMatrix pointer */
  kindls_mem->J = A;

  /* Initialize Jacobian-related data */
  kindls_mem->jacDQ     = SUNTRUE;
  kindls_mem->jac       = kinDlsDQJac;
  kindls_mem->J_data    = kin_mem;
  kindls_mem->last_flag = KINDLS_SUCCESS;

  /* Initialize counters */
  kinDlsInitializeCounters(kindls_mem);

  /* Attach linear solver memory to integrator memory */
  kin_mem->kin_lmem = kindls_mem;

  return(KINDLS_SUCCESS);
}


/*==================================================================
  KINDLS Exported fuctions -- Optional input/output
  ==================================================================*/

/*------------------------------------------------------------------
  KINDlsSetJacFn specifies the Jacobian function
  ------------------------------------------------------------------*/
int KINDlsSetJacFn(void *kinmem, KINDlsJacFn jac)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS",
                    "KINDlsSetJacFn", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS",
                    "KINDlsSetJacFn", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  if (jac != NULL) {
    kindls_mem->jacDQ  = SUNFALSE;
    kindls_mem->jac    = jac;
    kindls_mem->J_data = kin_mem->kin_user_data;
  } else {
    kindls_mem->jacDQ  = SUNTRUE;
    kindls_mem->jac    = kinDlsDQJac;
    kindls_mem->J_data = kin_mem;
  }

  return(KINDLS_SUCCESS);
}


/*------------------------------------------------------------------
  KINDlsGetWorkSpace returns the lenght of workspace allocated for
  the KINDls linear solver interface
  ------------------------------------------------------------------*/
int KINDlsGetWorkSpace(void *kinmem, long int *lenrwLS,
                       long int *leniwLS)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;
  long int lrw, liw;
  int flag;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS", 
                    "KINDlsGetWorkSpace", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS",
                    "KINGetWorkSpace", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  /* initialize outputs with requirements from KINDlsMem structure */
  *lenrwLS = 0;
  *leniwLS = 3;

  /* add LS sizes */
  if (kindls_mem->LS->ops->space) {
    flag = SUNLinSolSpace(kindls_mem->LS, &lrw, &liw);
    *lenrwLS += lrw;
    *leniwLS += liw;
  }

  return(KINDLS_SUCCESS);
}


/*------------------------------------------------------------------
  KINDlsGetNumJacEvals returns the number of Jacobian evaluations
  ------------------------------------------------------------------*/
int KINDlsGetNumJacEvals(void *kinmem, long int *njevals)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS",
                    "KINDlsGetNumJacEvals", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS",
                    "KINDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  *njevals = kindls_mem->nje;

  return(KINDLS_SUCCESS);
}


/*------------------------------------------------------------------
  KINDlsGetNumFuncEvals returns the number of calls to the user's
  F routine needed for the DQ Jacobian approximation
  ------------------------------------------------------------------*/
int KINDlsGetNumFuncEvals(void *kinmem, long int *nfevals)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS",
                    "KINDlsGetNumFuncEvals", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL,
                    "KINDLS", "KINDlsGetNumGuncEvals", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  *nfevals = kindls_mem->nfeDQ;

  return(KINDLS_SUCCESS);
}


/*------------------------------------------------------------------
  KINDlsGetLastFlag returns the last flag set in the KINDLS function
  ------------------------------------------------------------------*/
int KINDlsGetLastFlag(void *kinmem, long int *flag)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS",
                    "KINDlsGetLastFlag", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS",
                    "KINDlsGetLastFlag", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  *flag = kindls_mem->last_flag;

  return(KINDLS_SUCCESS);
}


/*------------------------------------------------------------------
  KINDlsGetReturnFlagName
  ------------------------------------------------------------------*/
char *KINDlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINDLS_SUCCESS:
    sprintf(name, "KINDLS_SUCCESS");
    break;
  case KINDLS_MEM_NULL:
    sprintf(name, "KINDLS_MEM_NULL");
    break;
  case KINDLS_LMEM_NULL:
    sprintf(name, "KINDLS_LMEM_NULL");
    break;
  case KINDLS_ILL_INPUT:
    sprintf(name, "KINDLS_ILL_INPUT");
    break;
  case KINDLS_MEM_FAIL:
    sprintf(name, "KINDLS_MEM_FAIL");
    break;
  case KINDLS_JACFUNC_ERR:
    sprintf(name,"KINDLS_JACFUNC_ERR");
    break;
  case KINDLS_SUNMAT_FAIL:
    sprintf(name,"KINDLS_SUNMAT_FAIL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}


/*==================================================================
  KINDLS Private function 
  ==================================================================*/

/*------------------------------------------------------------------
  kinDlsDenseDQJac 

  This routine is a wrapper for the Dense and Band implementations
  of the difference quotient Jacobian approximation routines. 
  ------------------------------------------------------------------*/
int kinDlsDQJac(N_Vector u, N_Vector fu, SUNMatrix Jac, 
                void *kinmem, N_Vector tmp1, N_Vector tmp2)
{
  int retval;
  KINMem kin_mem;
  kin_mem = (KINMem) kinmem;

  /* verify that Jac is non-NULL */
  if (Jac == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS", 
		    "kinDlsDQJac", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }

  if (SUNMatGetID(Jac) == SUNMATRIX_DENSE) {
    retval = kinDlsDenseDQJac(u, fu, Jac, kin_mem, tmp1, tmp2);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_BAND) {
    retval = kinDlsBandDQJac(u, fu, Jac, kin_mem, tmp1, tmp2);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_SPARSE) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINDLS",
		    "kinDlsDQJac",
                    "kinDlsDQJac not implemented for SUNMATRIX_SPARSE");
    retval = KIN_ILL_INPUT;
  } else {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINDLS",
		    "kinDlsDQJac",
                    "unrecognized matrix type for kinDlsDQJac");
    retval = KIN_ILL_INPUT;
  }
  return(retval);
}


/*------------------------------------------------------------------
  kinDlsDenseDQJac

  This routine generates a dense difference quotient approximation
  to the Jacobian of F(u). It assumes a dense SUNMatrix input
  stored column-wise, and that elements within each column are 
  contiguous. The address of the jth column of J is obtained via
  the function SUNDenseMatrix_Column() and this pointer is
  associated with an N_Vector using the N_VGetArrayPointer and
  N_VSetArrayPointer functions. Finally, the actual computation of
  the jth column of the Jacobian is done with a call to N_VLinearSum.
 
  The increment used in the finite-difference approximation
    J_ij = ( F_i(u+sigma_j * e_j) - F_i(u)  ) / sigma_j
  is
   sigma_j = max{|u_j|, |1/uscale_j|} * sqrt(uround)
 
  Note: uscale_j = 1/typ(u_j)
 
  NOTE: Any type of failure of the system function here leads to an
        unrecoverable failure of the Jacobian function and thus of 
        the linear solver setup function, stopping KINSOL.
  ------------------------------------------------------------------*/
int kinDlsDenseDQJac(N_Vector u, N_Vector fu, SUNMatrix Jac, 
                     KINMem kin_mem, N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv, ujsaved, ujscale, sign;
  realtype *tmp2_data, *u_data, *uscale_data;
  N_Vector ftemp, jthCol;
  sunindextype j, N;
  int retval = 0;
  KINDlsMem kindls_mem;

  /* access DlsMem interface structure */
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  /* access matrix dimension */
  N = SUNDenseMatrix_Rows(Jac);

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp  = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for u and uscale */
  u_data      = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(kin_mem->kin_uscale);

  /* This is the only for loop for 0..N-1 in KINSOL */

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(u) */

    /* Set data address of jthCol, and save u_j values and scaling */
    N_VSetArrayPointer(SUNDenseMatrix_Column(Jac,j), jthCol);
    ujsaved = u_data[j];
    ujscale = ONE/uscale_data[j];

    /* Compute increment */
    sign = (ujsaved >= ZERO) ? ONE : -ONE;
    inc = kin_mem->kin_sqrt_relfunc*SUNMAX(SUNRabs(ujsaved), ujscale)*sign;

    /* Increment u_j, call F(u), and return if error occurs */
    u_data[j] += inc;

    retval = kin_mem->kin_func(u, ftemp, kin_mem->kin_user_data);
    kindls_mem->nfeDQ++;
    if (retval != 0) break;

    /* reset u_j */
    u_data[j] = ujsaved;

    /* Construct difference quotient in jthCol */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fu, jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}


/*------------------------------------------------------------------
  kinDlsBandDQJac

  This routine generates a banded difference quotient approximation
  to the Jacobian of F(u).  It assumes a SUNBandMatrix input stored
  column-wise, and that elements within each column are contiguous.
  This makes it possible to get the address of a column of J via the
  function SUNBandMatrix_Column() and to write a simple for loop to
  set each of the elements of a column in succession.
 
  NOTE: Any type of failure of the system function her leads to an
        unrecoverable failure of the Jacobian function and thus of
        the linear solver setup function, stopping KINSOL.
  ------------------------------------------------------------------*/
int kinDlsBandDQJac(N_Vector u, N_Vector fu, SUNMatrix Jac,
                    KINMem kin_mem, N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv;
  N_Vector futemp, utemp;
  sunindextype group, i, j, width, ngroups, i1, i2;
  sunindextype N, mupper, mlower;
  realtype *col_j, *fu_data, *futemp_data, *u_data, *utemp_data, *uscale_data;
  int retval = 0;
  KINDlsMem kindls_mem;

  /* access DlsMem interface structure */
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  /* access matrix dimensions */
  N = SUNBandMatrix_Columns(Jac);
  mupper = SUNBandMatrix_UpperBandwidth(Jac);
  mlower = SUNBandMatrix_LowerBandwidth(Jac);

  /* Rename work vectors for use as temporary values of u and fu */
  futemp = tmp1;
  utemp  = tmp2;

  /* Obtain pointers to the data for ewt, fy, futemp, y, ytemp */
  fu_data     = N_VGetArrayPointer(fu);
  futemp_data = N_VGetArrayPointer(futemp);
  u_data      = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(kin_mem->kin_uscale);
  utemp_data  = N_VGetArrayPointer(utemp);

  /* Load utemp with u */
  N_VScale(ONE, u, utemp);

  /* Set bandwidth and number of column groups for band differencing */
  width   = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);
  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all utemp components in group */
    for(j=group-1; j < N; j+=width) {
      inc = kin_mem->kin_sqrt_relfunc*SUNMAX(SUNRabs(u_data[j]),
                                             ONE/SUNRabs(uscale_data[j]));
      utemp_data[j] += inc;
    }

    /* Evaluate f with incremented u */
    retval = kin_mem->kin_func(utemp, futemp, kin_mem->kin_user_data);
    if (retval != 0) return(retval); 

    /* Restore utemp components, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      utemp_data[j] = u_data[j];
      col_j = SUNBandMatrix_Column(Jac, j);
      inc = kin_mem->kin_sqrt_relfunc*SUNMAX(SUNRabs(u_data[j]),
                                             ONE/SUNRabs(uscale_data[j]));
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) = inc_inv * (futemp_data[i] - fu_data[i]);
    }
  }
  
  /* Increment counter nfeDQ */
  kindls_mem->nfeDQ += ngroups;

  return(0);
}


/*------------------------------------------------------------------
  kinDlsInitialize performs remaining initializations specific to 
  the direct linear solver interface (and solver itself)
  ------------------------------------------------------------------*/
int kinDlsInitialize(KINMem kin_mem)
{
  KINDlsMem kindls_mem;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS", 
                    "kinDlsInitialize", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS", 
                    "kinDlsInitialize", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;
  
  kinDlsInitializeCounters(kindls_mem);

  /* Set Jacobian function and data, depending on jacDQ (in case 
     it has changed based on user input) */
  if (kindls_mem->jacDQ) {
    kindls_mem->jac    = kinDlsDQJac;
    kindls_mem->J_data = kin_mem;
  } else {
    kindls_mem->J_data = kin_mem->kin_user_data;
  }

  if ( (kin_mem->kin_globalstrategy == KIN_PICARD) && kindls_mem->jacDQ ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL",
                    "kinDlsInitialize", MSG_NOL_FAIL);
    return(KIN_ILL_INPUT);
  }

  /* Call LS initialize routine */
  kindls_mem->last_flag = SUNLinSolInitialize(kindls_mem->LS);
  return(kindls_mem->last_flag);
}


/*------------------------------------------------------------------
  kinDlsSetup does the setup operations for the KINDLS linear solver
  interface. It calls the Jacobian evaluation routine, updates 
  counters, and calls the LS 'setup' routine to prepare for 
  subsequent calls to the LS 'solve' routine.

  The return value is either
     KINLS_SUCCESS = 0  if successful,
      1  if the LS 'setup' routine failed recoverably, or
     -1  if the jac or LS 'setup' routines failed unrecoverably.
  ------------------------------------------------------------------*/
int kinDlsSetup(KINMem kin_mem)
{
  KINDlsMem kindls_mem;
  int retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS", 
                    "kinDlsSetup", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS", 
                    "kinDlsSetup", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  /* Increment nje counter. */
  kindls_mem->nje++;

  /* Zero out J; call Jacobian routine jac; return if it failed. */
  retval = SUNMatZero(kindls_mem->J);
  if (retval != 0) {
    KINProcessError(kin_mem, KINDLS_SUNMAT_FAIL, "KINDLS",
                    "kinDlsSetup", MSGD_MATZERO_FAILED);
    kindls_mem->last_flag = KINDLS_SUNMAT_FAIL;
    return(-1);
  }

  retval = kindls_mem->jac(kin_mem->kin_uu, kin_mem->kin_fval, kindls_mem->J,
                           kindls_mem->J_data, kin_mem->kin_vtemp1, kin_mem->kin_vtemp2);
  if (retval != 0) {
    KINProcessError(kin_mem, KINDLS_JACFUNC_ERR, "KINDLS",
                    "kinDlsSetup", MSGD_JACFUNC_FAILED);
    kindls_mem->last_flag = KINDLS_JACFUNC_ERR;
    return(-1);
  }

  /* Call generic linear solver 'setup' with this system matrix, and
     return success/failure flag */
  kindls_mem->last_flag = SUNLinSolSetup(kindls_mem->LS, kindls_mem->J);
  return(kindls_mem->last_flag);
}


/*------------------------------------------------------------------
  kinDlsSolve interfaces between KINSOL and the generic 
  SUNLinearSolver object LS, by calling the LS 'solve' routine

  The return value is KINLS_SUCCESS = 0 if successful. The argument
  *sJpnorm is ignored.
  ------------------------------------------------------------------*/
int kinDlsSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                realtype *sJpnorm, realtype *sFdotJp)
{
  KINDlsMem kindls_mem;
  int retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDLS", 
		    "kinDlsSolve", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINDLS_LMEM_NULL, "KINDLS", 
		    "kinDlsSolve", MSGD_LMEM_NULL);
    return(KINDLS_LMEM_NULL);
  }
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  /* call the generic linear system solver */
  retval = SUNLinSolSolve(kindls_mem->LS, kindls_mem->J, x, b, ZERO);

  /*
    Compute the term sFdotJp for use in the linesearch routine.
    This term is subsequently corrected if the step is reduced by
    constraints or the linesearch.
    
    sFdotJp is the dot product of the scaled f vector and the scaled
    vector J*p, where the scaling uses fscale.
  */
  N_VProd(b, kin_mem->kin_fscale, b);
  N_VProd(b, kin_mem->kin_fscale, b);
  *sFdotJp = N_VDotProd(kin_mem->kin_fval, b);
  
  /* store solver return value and return */
  kindls_mem->last_flag = retval;
  return(retval);
}


/*------------------------------------------------------------------
  kinDlsFree frees memory associates with the KINDls system solver
  interface.
  ------------------------------------------------------------------*/
int kinDlsFree(KINMem kin_mem)
{
  KINDlsMem kindls_mem;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) return(KINDLS_SUCCESS);
  if (kin_mem->kin_lmem == NULL) return(KINDLS_SUCCESS);
  kindls_mem = (KINDlsMem) kin_mem->kin_lmem;

  /* Nullify SUNMatrix pointer */
  kindls_mem->J = NULL;

  /* free KINDls interface structure */
  free(kin_mem->kin_lmem);
  kindls_mem = NULL;
  
  return(KINDLS_SUCCESS);
}


/*------------------------------------------------------------------
  kinDlsInitializeCounters resets counters for the DLS interface
  ------------------------------------------------------------------*/
int kinDlsInitializeCounters(KINDlsMem kindls_mem)
{
  kindls_mem->nje   = 0;
  kindls_mem->nfeDQ = 0;
  
  return(0);
}
