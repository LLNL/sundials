/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for an IDADLS linear solver 
 * interface
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "ida_direct_impl.h"
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

/*=================================================================
  FUNCTION SPECIFIC CONSTANTS
  =================================================================*/

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/*===============================================================
  IDADLS Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  IDADlsSetLinearSolver specifies the direct linear solver.
  ---------------------------------------------------------------*/
int IDADlsSetLinearSolver(void *ida_mem, SUNLinearSolver LS,
                          SUNMatrix A)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if any input is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", 
                    "IDADlsSetLinearSolver", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  if ( (LS == NULL)  || (A == NULL) ) {
    IDAProcessError(NULL, IDADLS_ILL_INPUT, "IDADLS", 
                    "IDADlsSetLinearSolver",
                    "Both LS and A must be non-NULL");
    return(IDADLS_ILL_INPUT);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if solver and vector are compatible with DLS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_DIRECT) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDADLS", 
                    "IDADlsSetLinearSolver", 
                    "Non-direct LS supplied to IDADls interface");
    return(IDADLS_ILL_INPUT);
  }
  if (IDA_mem->ida_tempv1->ops->nvgetarraypointer == NULL ||
      IDA_mem->ida_tempv1->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDADLS", 
                    "IDADlsSetLinearSolver", MSGD_BAD_NVECTOR);
    return(IDADLS_ILL_INPUT);
  }

  /* free any existing system solver attached to IDA */
  if (IDA_mem->ida_lfree)  IDA_mem->ida_lfree(IDA_mem);

  /* Set four main system linear solver function fields in IDA_mem */
  IDA_mem->ida_linit  = idaDlsInitialize;
  IDA_mem->ida_lsetup = idaDlsSetup;
  IDA_mem->ida_lsolve = idaDlsSolve;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = idaDlsFree;
  
  /* Get memory for IDADlsMemRec */
  idadls_mem = NULL;
  idadls_mem = (IDADlsMem) malloc(sizeof(struct IDADlsMemRec));
  if (idadls_mem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADLS",
                    "IDADlsSetLinearSolver", MSGD_MEM_FAIL);
    return(IDADLS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  idadls_mem->LS = LS;
  
  /* set SUNMatrix pointer */
  idadls_mem->J = A;

  /* Initialize Jacobian-related data */
  idadls_mem->jacDQ     = SUNTRUE;
  idadls_mem->jac       = idaDlsDQJac;
  idadls_mem->J_data    = IDA_mem;
  idadls_mem->last_flag = IDADLS_SUCCESS;

  /* Initialize counters */
  idaDlsInitializeCounters(idadls_mem);

  /* Allocate memory for x */
  idadls_mem->x = N_VClone(IDA_mem->ida_tempv1);
  if (idadls_mem->x == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDADLS", 
                    "IDADlsSetLinearSolver", MSGD_MEM_FAIL);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }
  
  /* Attach linear solver memory to integrator memory */
  IDA_mem->ida_lmem = idadls_mem;

  return(IDADLS_SUCCESS);
}


/*===============================================================
  IDADLS Exported functions -- Optional input/output
  ===============================================================*/

/*---------------------------------------------------------------
  IDADlsSetJacFn specifies the Jacobian function.
  ---------------------------------------------------------------*/
int IDADlsSetJacFn(void *ida_mem, IDADlsJacFn jac)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS",
                    "IDADlsSetJacFn", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS",
                    "IDADlsSetJacFn", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  if (jac != NULL) {
    idadls_mem->jacDQ  = SUNFALSE;
    idadls_mem->jac    = jac;
    idadls_mem->J_data = IDA_mem->ida_user_data;
  } else {
    idadls_mem->jacDQ  = SUNTRUE;
    idadls_mem->jac    = idaDlsDQJac;
    idadls_mem->J_data = IDA_mem;
  }

  return(IDADLS_SUCCESS);
}

/*---------------------------------------------------------------
  IDADlsGetWorkSpace returns the length of workspace allocated 
  for the IDADLS linear solver interface
  ---------------------------------------------------------------*/
int IDADlsGetWorkSpace(void *ida_mem, long int *lenrwLS,
                       long int *leniwLS)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS",
                    "IDADlsGetWorkSpace", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS",
                    "IDADlsGetWorkSpace", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* initialize outputs with requirements from IDADlsMem structure */
  *lenrwLS = 0;
  *leniwLS = 3;

  /* add vector sizes */
  if (idadls_mem->x->ops->nvspace) {
    N_VSpace(idadls_mem->x, &lrw1, &liw1);
    *lenrwLS += lrw1;
    *leniwLS += liw1;
  }

  /* add LS sizes */
  if (idadls_mem->LS->ops->space) {
    flag = SUNLinSolSpace(idadls_mem->LS, &lrw, &liw);
    *lenrwLS += lrw;
    *leniwLS += liw;
  }

  return(IDADLS_SUCCESS);
}


/*---------------------------------------------------------------
  IDADlsGetNumJacEvals returns the number of Jacobian evaluations.
  ---------------------------------------------------------------*/
int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS",
                    "IDADlsGetNumJacEvals", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS",
                    "IDADlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  *njevals = idadls_mem->nje;

  return(IDADLS_SUCCESS);
}


/*---------------------------------------------------------------
  IDADlsGetNumResEvals returns the number of calls to the DAE 
  function needed for the DQ Jacobian approximation.
  ---------------------------------------------------------------*/
int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS",
                    "IDADlsGetNumResEvals", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS",
                    "IDADlsGetNumResEvals", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  *nrevalsLS = idadls_mem->nreDQ;

  return(IDADLS_SUCCESS);
}


/*---------------------------------------------------------------
  IDADlsGetReturnFlagName returns the name associated with a 
  IDADLS return value.
  ---------------------------------------------------------------*/
char *IDADlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDADLS_SUCCESS:
    sprintf(name,"IDADLS_SUCCESS");
    break;   
  case IDADLS_MEM_NULL:
    sprintf(name,"IDADLS_MEM_NULL");
    break;
  case IDADLS_LMEM_NULL:
    sprintf(name,"IDADLS_LMEM_NULL");
    break;
  case IDADLS_ILL_INPUT:
    sprintf(name,"IDADLS_ILL_INPUT");
    break;
  case IDADLS_MEM_FAIL:
    sprintf(name,"IDADLS_MEM_FAIL");
    break;
  case IDADLS_JACFUNC_UNRECVR:
    sprintf(name,"IDADLS_JACFUNC_UNRECVR");
    break;
  case IDADLS_JACFUNC_RECVR:
    sprintf(name,"IDADLS_JACFUNC_RECVR");
    break;
  case IDADLS_SUNMAT_FAIL:
    sprintf(name,"IDADLS_SUNMAT_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
  IDADlsGetLastFlag returns the last flag set in a IDADLS function.
  ---------------------------------------------------------------*/
int IDADlsGetLastFlag(void *ida_mem, long int *flag)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS",
                    "IDADlsGetLastFlag", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS",
                    "IDADlsGetLastFlag", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  *flag = idadls_mem->last_flag;

  return(IDADLS_SUCCESS);
}


/*===============================================================
  IDADLS Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  idaDlsDQJac:

  This routine is a wrapper for the Dense and Band 
  implementations of the difference quotient Jacobian 
  approximation routines.
---------------------------------------------------------------*/
int idaDlsDQJac(realtype t, realtype c_j, N_Vector y,  
                N_Vector yp, N_Vector r, SUNMatrix Jac, 
                void *ida_mem, N_Vector tmp1, N_Vector tmp2,
                N_Vector tmp3)
{
  int retval;
  IDAMem IDA_mem;
  IDA_mem = (IDAMem) ida_mem;

  /* verify that Jac is non-NULL */
  if (Jac == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", 
		    "idaDlsDQJac", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }

  if (SUNMatGetID(Jac) == SUNMATRIX_DENSE) {
    retval = idaDlsDenseDQJac(t, c_j, y, yp, r, Jac, IDA_mem, tmp1);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_BAND) {
    retval = idaDlsBandDQJac(t, c_j, y, yp, r, Jac, IDA_mem, tmp1, tmp2, tmp3);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_SPARSE) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDADLS", 
		    "idaDlsDQJac", 
                    "idaDlsDQJac not implemented for SUNMATRIX_SPARSE");
    retval = IDA_ILL_INPUT;
  } else {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDADLS", 
		    "idaDlsDQJac", 
                    "unrecognized matrix type for idaDlsDQJac");
    retval = IDA_ILL_INPUT;
  }
  return(retval);
}


/*---------------------------------------------------------------
  idaDlsDenseDQJac 

  This routine generates a dense difference quotient approximation
  to the Jacobian F_y + c_j*F_y'. It assumes a dense SUNmatrix 
  input (stored column-wise, and that elements within each column
  are contiguous). The address of the jth column of J is obtained 
  via the function SUNDenseMatrix_Column() and this pointer is 
  associated with an N_Vector using the 
  N_VGetArrayPointer/N_VSetArrayPointer functions.  Finally, the 
  actual computation of the jth column of the Jacobian is 
  done with a call to N_VLinearSum.
---------------------------------------------------------------*/
int idaDlsDenseDQJac(realtype tt, realtype c_j, N_Vector yy, 
                     N_Vector yp, N_Vector rr, SUNMatrix Jac,
                     IDAMem IDA_mem, N_Vector tmp1)
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  sunindextype j, N;
  int retval = 0;
  IDADlsMem idadls_mem;

  /* access DlsMem interface structure */
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* access matrix dimension */
  N = SUNDenseMatrix_Rows(Jac);

  /* Rename work vectors for readibility */
  rtemp = tmp1;

  /* Create an empty vector for matrix column calculations */
  jthCol = N_VCloneEmpty(tmp1);
  
  /* Obtain pointers to the data for ewt, yy, yp. */
  ewt_data = N_VGetArrayPointer(IDA_mem->ida_ewt);
  y_data   = N_VGetArrayPointer(yy);
  yp_data  = N_VGetArrayPointer(yp);
  if(IDA_mem->ida_constraints!=NULL)
    cns_data = N_VGetArrayPointer(IDA_mem->ida_constraints);

  srur = SUNRsqrt(IDA_mem->ida_uround);

  for (j=0; j < N; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetArrayPointer(SUNDenseMatrix_Column(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as hh*yp_j. */

    inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(IDA_mem->ida_hh*ypj) ),
                  ONE/ewt_data[j] );

    if (IDA_mem->ida_hh*ypj < ZERO) inc = -inc;
    inc = (yj + inc) - yj;

    /* Adjust sign(inc) again if y_j has an inequality constraint. */
    if (IDA_mem->ida_constraints != NULL) {
      conj = cns_data[j];
      if (SUNRabs(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
      else if (SUNRabs(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
    }

    /* Increment y_j and yp_j, call res, and break on error return. */
    y_data[j] += inc;
    yp_data[j] += c_j*inc;

    retval = IDA_mem->ida_res(tt, yy, yp, rtemp, IDA_mem->ida_user_data);
    idadls_mem->nreDQ++;
    if (retval != 0) break;

    /* Construct difference quotient in jthCol */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, rr, jthCol);

    /* DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);  /\* UNNECESSARY *\/ */

    /*  reset y_j, yp_j */     
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  /* Destroy jthCol vector */
  N_VSetArrayPointer(NULL, jthCol);  /* SHOULDN'T BE NEEDED */
  N_VDestroy(jthCol);

  return(retval);

}


/*---------------------------------------------------------------
  idaDlsBandDQJac 

  This routine generates a banded difference quotient approximation 
  JJ to the DAE system Jacobian J.  It assumes a band SUNMatrix 
  input (stored column-wise, and that elements within each column
  are contiguous).  This makes it possible to get the address 
  of a column of JJ via the function SUNBandMatrix_Column(). The 
  columns of the Jacobian are constructed using mupper + mlower + 1 
  calls to the res routine, and appropriate differencing.
  The return value is either IDABAND_SUCCESS = 0, or the nonzero 
  value returned by the res routine, if any.
  ---------------------------------------------------------------*/
int idaDlsBandDQJac(realtype tt, realtype c_j, N_Vector yy,
                    N_Vector yp, N_Vector rr, SUNMatrix Jac,
                    IDAMem IDA_mem, N_Vector tmp1, N_Vector tmp2,
                    N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj, ewtj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  realtype *ytemp_data, *yptemp_data, *rtemp_data, *r_data, *col_j;
  N_Vector rtemp, ytemp, yptemp;
  sunindextype i, j, i1, i2, width, ngroups, group;
  sunindextype N, mupper, mlower;
  int retval = 0;
  IDADlsMem idadls_mem;

  /* access DlsMem interface structure */
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* access matrix dimensions */
  N = SUNBandMatrix_Columns(Jac);
  mupper = SUNBandMatrix_UpperBandwidth(Jac);
  mlower = SUNBandMatrix_LowerBandwidth(Jac);

  /* Rename work vectors for use as temporary values of r, y and yp */
  rtemp = tmp1;
  ytemp = tmp2;
  yptemp= tmp3;

  /* Obtain pointers to the data for all eight vectors used.  */
  ewt_data    = N_VGetArrayPointer(IDA_mem->ida_ewt);
  r_data      = N_VGetArrayPointer(rr);
  y_data      = N_VGetArrayPointer(yy);
  yp_data     = N_VGetArrayPointer(yp);
  rtemp_data  = N_VGetArrayPointer(rtemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);
  if (IDA_mem->ida_constraints != NULL)
    cns_data = N_VGetArrayPointer(IDA_mem->ida_constraints);

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, yy, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */
  srur = SUNRsqrt(IDA_mem->ida_uround);
  width = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all yy[j] and yp[j] for j in this group. */
    for (j=group-1; j<N; j+=width) {
        yj = y_data[j];
        ypj = yp_data[j];
        ewtj = ewt_data[j];

        /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
        adjustments using ypj and ewtj if this is small, and a further
        adjustment to give it the same sign as hh*ypj. */
        inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(IDA_mem->ida_hh*ypj) ),
                      ONE/ewtj );
        if (IDA_mem->ida_hh*ypj < ZERO)  inc = -inc;
        inc = (yj + inc) - yj;

        /* Adjust sign(inc) again if yj has an inequality constraint. */
        if (IDA_mem->ida_constraints != NULL) {
          conj = cns_data[j];
          if (SUNRabs(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
          else if (SUNRabs(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
        }

        /* Increment yj and ypj. */
        ytemp_data[j] += inc;
        yptemp_data[j] += IDA_mem->ida_cj*inc;
    }

    /* Call res routine with incremented arguments. */
    retval = IDA_mem->ida_res(tt, ytemp, yptemp, rtemp, IDA_mem->ida_user_data);
    idadls_mem->nreDQ++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */
    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */
      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = SUNBandMatrix_Column(Jac, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */
      inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(IDA_mem->ida_hh*ypj) ),
                    ONE/ewtj );
      if (IDA_mem->ida_hh*ypj < ZERO)  inc = -inc;
      inc = (yj + inc) - yj;
      if (IDA_mem->ida_constraints != NULL) {
        conj = cns_data[j];
        if (SUNRabs(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
        else if (SUNRabs(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
      }
      
      /* Load the difference quotient Jacobian elements for column j */
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower,N-1);
      for (i=i1; i<=i2; i++) 
        SM_COLUMN_ELEMENT_B(col_j,i,j) = inc_inv * (rtemp_data[i]-r_data[i]);
    }
  }
  
  return(retval);
}


/*---------------------------------------------------------------
 idaDlsInitialize performs remaining initializations specific
 to the direct linear solver interface (and solver itself)
---------------------------------------------------------------*/
int idaDlsInitialize(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", 
                    "idaDlsInitialize", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", 
                    "idaDlsInitialize", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;
  
  idaDlsInitializeCounters(idadls_mem);

  /* Set Jacobian function and data, depending on jacDQ (in case 
     it has changed based on user input) */
  if (idadls_mem->jacDQ) {
    idadls_mem->jac    = idaDlsDQJac;
    idadls_mem->J_data = IDA_mem;
  } else {
    idadls_mem->J_data = IDA_mem->ida_user_data;
  }

  /* Call LS initialize routine */
  idadls_mem->last_flag = SUNLinSolInitialize(idadls_mem->LS);
  return(idadls_mem->last_flag);
}


/*---------------------------------------------------------------
  idaDlsSetup does the setup operations for the IDADLS linear 
  solver interface.  It calls the Jacobian evaluation routine, 
  updates counters, and calls the LS 'setup' routine to prepare
  for subsequent calls to the LS 'solve' routine.

  The return value is either
     IDADLS_SUCCESS = 0  if successful,
      1  if the jac or LS 'setup' routine failed recoverably, or
     -1  if the jac or LS 'setup' routines failed unrecoverably.
---------------------------------------------------------------*/
int idaDlsSetup(IDAMem IDA_mem, N_Vector y, N_Vector yp, N_Vector r, 
                N_Vector vt1, N_Vector vt2, N_Vector vt3)
{
  int retval;
  IDADlsMem idadls_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", 
                    "idaDlsSetup", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", 
                    "idaDlsSetup", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* Increment nje counter. */
  idadls_mem->nje++;

  /* Zero out J; call Jacobian routine jac; return if it failed. */
  retval = SUNMatZero(idadls_mem->J);
  if (retval != 0) {
    IDAProcessError(IDA_mem, IDADLS_SUNMAT_FAIL, "IDADLS",
                    "idaDlsSetup", MSGD_MATZERO_FAILED);
    idadls_mem->last_flag = IDADLS_SUNMAT_FAIL;
    return(-1);
  }
  
  retval = idadls_mem->jac(IDA_mem->ida_tn, IDA_mem->ida_cj, y,
                           yp, r, idadls_mem->J,
                           idadls_mem->J_data, vt1, vt2, vt3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDADLS_JACFUNC_UNRECVR, "IDADLS",
                    "idaDlsSetup", MSGD_JACFUNC_FAILED);
    idadls_mem->last_flag = IDADLS_JACFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    idadls_mem->last_flag = IDADLS_JACFUNC_RECVR;
    return(1);
  }

  /* Call generic linear solver 'setup' with this system matrix, and
     return success/failure flag */
  idadls_mem->last_flag = SUNLinSolSetup(idadls_mem->LS, idadls_mem->J);
  return(idadls_mem->last_flag);
}


/*---------------------------------------------------------------
 idaDlsSolve interfaces between IDA and the generic 
 SUNLinearSolver object LS, by calling the LS 'solve' routine 
 and scaling the result if cjratio does not equal one.
---------------------------------------------------------------*/
int idaDlsSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                N_Vector ycur, N_Vector ypcur, N_Vector rescur)
{
  int retval;
  IDADlsMem idadls_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", 
		    "idaDlsSolve", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", 
		    "idaDlsSolve", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* call the generic linear system solver, and copy x to b */
  retval = SUNLinSolSolve(idadls_mem->LS, idadls_mem->J, idadls_mem->x, b, ZERO);
  N_VScale(ONE, idadls_mem->x, b);
  
  /* scale the correction to account for change in cj */
  if (IDA_mem->ida_cjratio != ONE) 
    N_VScale(TWO/(ONE + IDA_mem->ida_cjratio), b, b);
  
  /* store solver return value and return */
  idadls_mem->last_flag = retval;
  return(retval);
}


/*---------------------------------------------------------------
 idaDlsFree frees memory associates with the IDADls system
 solver interface.
---------------------------------------------------------------*/
int idaDlsFree(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL)  return (IDADLS_SUCCESS);
  if (IDA_mem->ida_lmem == NULL)  return(IDADLS_SUCCESS);
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* Free x vector */
  if (idadls_mem->x) {
    N_VDestroy(idadls_mem->x);
    idadls_mem->x = NULL;
  }

  /* Nullify SUNMatrix pointer */
  idadls_mem->J = NULL;

  /* free IDADls interface structure */
  free(IDA_mem->ida_lmem);
  
  return(IDADLS_SUCCESS);
}


/*---------------------------------------------------------------
  idaDlsInitializeCounters resets counters for the DLS interface
  ---------------------------------------------------------------*/
int idaDlsInitializeCounters(IDADlsMem idadls_mem)
{
  idadls_mem->nje   = 0;
  idadls_mem->nreDQ = 0;
  return(0);
}
