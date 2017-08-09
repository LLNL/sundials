/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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
 * This is the implementation file for an IDADLS linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "ida_direct_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * IDADlsSetDenseJacFn specifies the dense Jacobian function.
 */
int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn jac)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", "IDADlsSetDenseJacFn", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", "IDADlsSetDenseJacFn", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  if (jac != NULL) {
    idadls_mem->d_jacDQ = FALSE;
    idadls_mem->d_djac = jac;
  } else {
    idadls_mem->d_jacDQ = TRUE;
  }

  return(IDADLS_SUCCESS);
}

/*
 * IDADlsSetBandJacFn specifies the band Jacobian function.
 */
int IDADlsSetBandJacFn(void *ida_mem, IDADlsBandJacFn jac)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", "IDADlsSetBandJacFn", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", "IDADlsSetBandJacFn", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  if (jac != NULL) {
    idadls_mem->d_jacDQ = FALSE;
    idadls_mem->d_bjac = jac;
  } else {
    idadls_mem->d_jacDQ = TRUE;
  }

  return(IDADLS_SUCCESS);
}

/*
 * IDADlsGetWorkSpace returns the length of workspace allocated for the
 * IDALAPACK linear solver.
 */
int IDADlsGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", "IDADlsGetWorkSpace", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", "IDADlsGetWorkSpace", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  if (idadls_mem->d_type == SUNDIALS_DENSE) {
    *lenrwLS = idadls_mem->d_n*idadls_mem->d_n;
    *leniwLS = idadls_mem->d_n;
  } else if (idadls_mem->d_type == SUNDIALS_BAND) {
    *lenrwLS = idadls_mem->d_n*(idadls_mem->d_smu + idadls_mem->d_ml + 1);
    *leniwLS = idadls_mem->d_n;
  }
    
  return(IDADLS_SUCCESS);
}

/*
 * IDADlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", "IDADlsGetNumJacEvals", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", "IDADlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  *njevals = idadls_mem->d_nje;

  return(IDADLS_SUCCESS);
}

/*
 * IDADlsGetNumResEvals returns the number of calls to the DAE function
 * needed for the DQ Jacobian approximation.
 */
int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", "IDADlsGetNumFctEvals", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", "IDADlsGetNumFctEvals", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  *nrevalsLS = idadls_mem->d_nreDQ;

  return(IDADLS_SUCCESS);
}

/*
 * IDADlsGetReturnFlagName returns the name associated with a IDALAPACK
 * return value.
 */
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
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * IDADlsGetLastFlag returns the last flag set in a IDALAPACK function.
 */
int IDADlsGetLastFlag(void *ida_mem, long int *flag)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDADLS", "IDADlsGetLastFlag", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_LMEM_NULL, "IDADLS", "IDADlsGetLastFlag", MSGD_LMEM_NULL);
    return(IDADLS_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  *flag = idadls_mem->d_last_flag;

  return(IDADLS_SUCCESS);
}

/* 
 * =================================================================
 * DQ JACOBIAN APPROXIMATIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * idaDlsDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian F_y + c_j*F_y'. It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro LAPACK_DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 
int idaDlsDenseDQJac(sunindextype N, realtype tt, realtype c_j,
                     N_Vector yy, N_Vector yp, N_Vector rr, 
                     DlsMat Jac, void *data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  sunindextype j;
  int retval = 0;

  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* data points to IDA_mem */
  IDA_mem = (IDAMem) data;
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  rtemp  = tmp1;
  jthCol = tmp2;

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
    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as hh*yp_j. */

    inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(IDA_mem->ida_hh*ypj) ) , ONE/ewt_data[j] );

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
    idadls_mem->d_nreDQ++;
    if (retval != 0) break;

    /* Construct difference quotient in jthCol */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, rr, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);

    /*  reset y_j, yp_j */     
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);

}

/*
 * -----------------------------------------------------------------
 * idaDlsBandDQJac 
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation JJ
 * to the DAE system Jacobian J.  It assumes that a band matrix of type
 * BandMat is stored column-wise, and that elements within each column
 * are contiguous.  The address of the jth column of JJ is obtained via
 * the macros BAND_COL and BAND_COL_ELEM. The columns of the Jacobian are 
 * constructed using mupper + mlower + 1 calls to the res routine, and
 * appropriate differencing.
 * The return value is either IDABAND_SUCCESS = 0, or the nonzero value returned
 * by the res routine, if any.
 */

int idaDlsBandDQJac(sunindextype N, sunindextype mupper, sunindextype mlower,
                    realtype tt, realtype c_j, 
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    DlsMat Jac, void *data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj, ewtj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  realtype *ytemp_data, *yptemp_data, *rtemp_data, *r_data, *col_j;
  N_Vector rtemp, ytemp, yptemp;
  sunindextype i, j, i1, i2, width, ngroups, group;
  int retval = 0;

  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* data points to IDA_mem */
  IDA_mem = (IDAMem) data;
  idadls_mem = (IDADlsMem) IDA_mem->ida_lmem;

  rtemp = tmp1; /* Rename work vector for use as the perturbed residual. */

  ytemp = tmp2; /* Rename work vector for use as a temporary for yy. */


  yptemp= tmp3; /* Rename work vector for use as a temporary for yp. */

  /* Obtain pointers to the data for all eight vectors used.  */

  ewt_data = N_VGetArrayPointer(IDA_mem->ida_ewt);
  r_data   = N_VGetArrayPointer(rr);
  y_data   = N_VGetArrayPointer(yy);
  yp_data  = N_VGetArrayPointer(yp);

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

        inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(IDA_mem->ida_hh*ypj) ) , ONE/ewtj );

        if (IDA_mem->ida_hh*ypj < ZERO) inc = -inc;
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
    idadls_mem->d_nreDQ++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */

    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */

      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = BAND_COL(Jac, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */

      inc = SUNMAX( srur * SUNMAX( SUNRabs(yj), SUNRabs(IDA_mem->ida_hh*ypj) ) , ONE/ewtj );
      if (IDA_mem->ida_hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      if (IDA_mem->ida_constraints != NULL) {
        conj = cns_data[j];
        if (SUNRabs(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
        else if (SUNRabs(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
      }
      
      /* Load the difference quotient Jacobian elements for column j. */

      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower,N-1);
      
      for (i=i1; i<=i2; i++) 
            BAND_COL_ELEM(col_j,i,j) = inc_inv*(rtemp_data[i]-r_data[i]);
    }
    
  }
  
  return(retval);
  
}

int idaDlsInitializeCounters(IDADlsMem idadls_mem)
{
  idadls_mem->d_nje   = 0;
  idadls_mem->d_nreDQ = 0;
  return(0);
}
