/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-12-01 22:48:59 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for an IDADIRECT linear solver.
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
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define res            (IDA_mem->ida_res)
#define rdata          (IDA_mem->ida_rdata)
#define uround         (IDA_mem->ida_uround)
#define nst            (IDA_mem->ida_nst)
#define tn             (IDA_mem->ida_tn)
#define hh             (IDA_mem->ida_hh)
#define cj             (IDA_mem->ida_cj)
#define cjratio        (IDA_mem->ida_cjratio)
#define ewt            (IDA_mem->ida_ewt)
#define constraints    (IDA_mem->ida_constraints)

#define linit          (IDA_mem->ida_linit)
#define lsetup         (IDA_mem->ida_lsetup)
#define lsolve         (IDA_mem->ida_lsolve)
#define lfree          (IDA_mem->ida_lfree)
#define lperf          (IDA_mem->ida_lperf)
#define lmem           (IDA_mem->ida_lmem)
#define tempv          (IDA_mem->ida_tempv1)
#define setupNonNull   (IDA_mem->ida_setupNonNull)

#define mtype          (idadls_mem->d_type)
#define n              (idadls_mem->d_n)
#define ml             (idadls_mem->d_ml)
#define mu             (idadls_mem->d_mu)
#define smu            (idadls_mem->d_smu)
#define djac           (idadls_mem->d_djac)
#define bjac           (idadls_mem->d_bjac)
#define M              (idadls_mem->d_J)
#define pivots         (idadls_mem->d_pivots)
#define nje            (idadls_mem->d_nje)
#define nreDQ          (idadls_mem->d_nreDQ)
#define J_data         (idadls_mem->d_J_data)
#define last_flag      (idadls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * IDADlsSetJacFn specifies the (dense or band) Jacobian function.
 */
int IDADlsSetJacFn(void *ida_mem, void *jac, void *jac_data)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADIRECT_MEM_NULL, "IDADIRECT", "IDADlsSetJacFn", MSGD_IDAMEM_NULL);
    return(IDADIRECT_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADIRECT_LMEM_NULL, "IDADIRECT", "IDADlsSetJacFn", MSGD_LMEM_NULL);
    return(IDADIRECT_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) lmem;

  if (mtype == SUNDIALS_DENSE) 
    djac = (IDADlsDenseJacFn) jac;
  else if (mtype == SUNDIALS_BAND)
    bjac = (IDADlsBandJacFn) jac;

  J_data = jac_data;

  return(IDADIRECT_SUCCESS);
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
    IDAProcessError(NULL, IDADIRECT_MEM_NULL, "IDADIRECT", "IDADlsGetWorkSpace", MSGD_IDAMEM_NULL);
    return(IDADIRECT_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADIRECT_LMEM_NULL, "IDADIRECT", "IDADlsGetWorkSpace", MSGD_LMEM_NULL);
    return(IDADIRECT_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) lmem;

  if (mtype == SUNDIALS_DENSE) {
    *lenrwLS = n*n;
    *leniwLS = n;
  } else if (mtype == SUNDIALS_BAND) {
    *lenrwLS = n*(smu + ml + 1);
    *leniwLS = n;
  }
    
  return(IDADIRECT_SUCCESS);
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
    IDAProcessError(NULL, IDADIRECT_MEM_NULL, "IDADIRECT", "IDADlsGetNumJacEvals", MSGD_IDAMEM_NULL);
    return(IDADIRECT_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADIRECT_LMEM_NULL, "IDADIRECT", "IDADlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(IDADIRECT_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) lmem;

  *njevals = nje;

  return(IDADIRECT_SUCCESS);
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
    IDAProcessError(NULL, IDADIRECT_MEM_NULL, "IDADIRECT", "IDADlsGetNumFctEvals", MSGD_IDAMEM_NULL);
    return(IDADIRECT_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADIRECT_LMEM_NULL, "IDADIRECT", "IDADlsGetNumFctEvals", MSGD_LMEM_NULL);
    return(IDADIRECT_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) lmem;

  *nrevalsLS = nreDQ;

  return(IDADIRECT_SUCCESS);
}

/*
 * IDADlsGetReturnFlagName returns the name associated with a IDALAPACK
 * return value.
 */
char *IDADlsGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDADIRECT_SUCCESS:
    sprintf(name,"IDADIRECT_SUCCESS");
    break;   
  case IDADIRECT_MEM_NULL:
    sprintf(name,"IDADIRECT_MEM_NULL");
    break;
  case IDADIRECT_LMEM_NULL:
    sprintf(name,"IDADIRECT_LMEM_NULL");
    break;
  case IDADIRECT_ILL_INPUT:
    sprintf(name,"IDADIRECT_ILL_INPUT");
    break;
  case IDADIRECT_MEM_FAIL:
    sprintf(name,"IDADIRECT_MEM_FAIL");
    break;
  case IDADIRECT_JACFUNC_UNRECVR:
    sprintf(name,"IDADIRECT_JACFUNC_UNRECVR");
    break;
  case IDADIRECT_JACFUNC_RECVR:
    sprintf(name,"IDADIRECT_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * IDADlsGetLastFlag returns the last flag set in a IDALAPACK function.
 */
int IDADlsGetLastFlag(void *ida_mem, int *flag)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADIRECT_MEM_NULL, "IDADIRECT", "IDADlsGetLastFlag", MSGD_IDAMEM_NULL);
    return(IDADIRECT_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADIRECT_LMEM_NULL, "IDADIRECT", "IDADlsGetLastFlag", MSGD_LMEM_NULL);
    return(IDADIRECT_LMEM_NULL);
  }
  idadls_mem = (IDADlsMem) lmem;

  *flag = last_flag;

  return(IDADIRECT_SUCCESS);
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
int idaDlsDenseDQJac(int N, realtype tt, realtype c_j,
                     N_Vector yy, N_Vector yp, N_Vector rr, 
                     DlsMat Jac, void *jac_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  int j;
  int retval = 0;

  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* jac_data points to IDA_mem */
  IDA_mem = (IDAMem) jac_data;
  idadls_mem = (IDADlsMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  rtemp  = tmp1;
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, yy, yp. */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(yy);
  yp_data  = N_VGetArrayPointer(yp);
  if(constraints!=NULL) cns_data = N_VGetArrayPointer(constraints);

  srur = RSqrt(uround);

  for (j=0; j < N; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as hh*yp_j. */

    inc = MAX( srur * MAX( ABS(yj), ABS(hh*ypj) ) , ONE/ewt_data[j] );

    if (hh*ypj < ZERO) inc = -inc;
    inc = (yj + inc) - yj;

    /* Adjust sign(inc) again if y_j has an inequality constraint. */
    if (constraints != NULL) {
      conj = cns_data[j];
      if (ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
      else if (ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
    }

    /* Increment y_j and yp_j, call res, and break on error return. */
    y_data[j] += inc;
    yp_data[j] += c_j*inc;

    retval = res(tt, yy, yp, rtemp, rdata);
    nreDQ++;
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

int idaDlsBandDQJac(int N, int mupper, int mlower,
                    realtype tt, realtype c_j, 
                    N_Vector yy, N_Vector yp, N_Vector rr,
                    DlsMat Jac, void *jac_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj, ewtj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  realtype *ytemp_data, *yptemp_data, *rtemp_data, *r_data, *col_j;
  int group;
  
  N_Vector rtemp, ytemp, yptemp;
  int i, j, i1, i2, width, ngroups;
  int retval = 0;

  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* jac_data points to IDA_mem */
  IDA_mem = (IDAMem) jac_data;
  idadls_mem = (IDADlsMem) lmem;

  rtemp = tmp1; /* Rename work vector for use as the perturbed residual. */

  ytemp = tmp2; /* Rename work vector for use as a temporary for yy. */


  yptemp= tmp3; /* Rename work vector for use as a temporary for yp. */

  /* Obtain pointers to the data for all eight vectors used.  */

  ewt_data = N_VGetArrayPointer(ewt);
  r_data   = N_VGetArrayPointer(rr);
  y_data   = N_VGetArrayPointer(yy);
  yp_data  = N_VGetArrayPointer(yp);

  rtemp_data  = N_VGetArrayPointer(rtemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);

  if (constraints != NULL) cns_data = N_VGetArrayPointer(constraints);

  /* Initialize ytemp and yptemp. */

  N_VScale(ONE, yy, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */

  srur = RSqrt(uround);
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

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

        inc = MAX( srur * MAX( ABS(yj), ABS(hh*ypj) ) , ONE/ewtj );

        if (hh*ypj < ZERO) inc = -inc;
        inc = (yj + inc) - yj;

        /* Adjust sign(inc) again if yj has an inequality constraint. */

        if (constraints != NULL) {
          conj = cns_data[j];
          if (ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
          else if (ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
        }

        /* Increment yj and ypj. */

        ytemp_data[j] += inc;
        yptemp_data[j] += cj*inc;
    }

    /* Call res routine with incremented arguments. */

    retval = res(tt, ytemp, yptemp, rtemp, rdata);
    nreDQ++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */

    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */

      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = BAND_COL(Jac, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */

      inc = MAX( srur * MAX( ABS(yj), ABS(hh*ypj) ) , ONE/ewtj );
      if (hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      if (constraints != NULL) {
        conj = cns_data[j];
        if (ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
        else if (ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
      }
      
      /* Load the difference quotient Jacobian elements for column j. */

      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower,N-1);
      
      for (i=i1; i<=i2; i++) 
            BAND_COL_ELEM(col_j,i,j) = inc_inv*(rtemp_data[i]-r_data[i]);
    }
    
  }
  
  return(retval);
  
}
