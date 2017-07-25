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
 * This is the implementation file for the CVDLS linear solvers
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_direct_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_band.h>
#include <sundials/sundials_dense.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PRIVATE FUNCTIONS
 * =================================================================
 */

/*---------------------------------------------------------------
 CVDlsSetupMatrix determines whether to create a new dense/band 
 Jacobian matrix (or use a stored version), based on heuristics
 regarding previous converence issues, the number of time steps 
 since it was last updated, etc.; it then creates the system
 matrix from this, the 'gamma' factor and the identity 
 matrix, A = I-gamma*J.
---------------------------------------------------------------*/
int CVDlsSetupMatrix(void *cvode_mem, N_Vector vtemp1,
                     N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  sunindextype ier;
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;
  int retval;
  booleantype DENSE;
  N_Vector ypred, fpred;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", 
		    "CVDlsCallSetup", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", 
		    "CVDlsCallSetup", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;
  ypred = cv_mem->cv_zn[0];
  fpred = cv_mem->cv_ftemp;

  /* set flag to indicate DENSE vs BAND */
  DENSE = (cvdls_mem->d_savedJ->type == SUNDIALS_DENSE) ? TRUE : FALSE;
  
  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  jbad = (cv_mem->cv_nst == 0) || (cv_mem->cv_nst > cvdls_mem->d_nstlj + CVD_MSBJ) ||
         ((cv_mem->cv_convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (cv_mem->cv_convfail == CV_FAIL_OTHER);
  jok = !jbad;
 
  /* If jok = TRUE, use saved copy of J */
  if (jok) {
    cv_mem->cv_jcur = FALSE;
    if (DENSE) {
      DenseCopy(cvdls_mem->d_savedJ, cvdls_mem->d_M);
    } else {
      BandCopy(cvdls_mem->d_savedJ, cvdls_mem->d_M, cvdls_mem->d_mu, cvdls_mem->d_ml);
    }

  /* If jok = FALSE, call jac routine for new J value */
  } else {
    cvdls_mem->d_nje++;
    cvdls_mem->d_nstlj = cv_mem->cv_nst;
    cv_mem->cv_jcur = TRUE;
    SetToZero(cvdls_mem->d_M);

    if (DENSE) {
      retval = cvdls_mem->d_djac(cvdls_mem->d_n, cv_mem->cv_tn, ypred,
                                 fpred, cvdls_mem->d_M, cvdls_mem->d_J_data,
                                 vtemp1, vtemp2, vtemp3);
    } else {
      retval = cvdls_mem->d_bjac(cvdls_mem->d_n, cvdls_mem->d_mu,
                               cvdls_mem->d_ml, cv_mem->cv_tn, ypred,
                               fpred, cvdls_mem->d_M, cvdls_mem->d_J_data,
                               vtemp1, vtemp2, vtemp3);
    }
    if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVDLS",
                     "CVDLSSetupMatrix", MSGD_JACFUNC_FAILED);
      cvdls_mem->d_last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      cvdls_mem->d_last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }


    if (DENSE) {
      DenseCopy(cvdls_mem->d_M, cvdls_mem->d_savedJ);
    } else {
      BandCopy(cvdls_mem->d_M, cvdls_mem->d_savedJ, cvdls_mem->d_mu, cvdls_mem->d_ml);
    }

  }
  
  /* Scale and add I to get M = I - gamma*J */
  if (DENSE) {
    DenseScale(-cv_mem->cv_gamma, cvdls_mem->d_M);
  } else {
    BandScale(-cv_mem->cv_gamma, cvdls_mem->d_M);
  }
  AddIdentity(cvdls_mem->d_M);
    
  return(CVDLS_SUCCESS);
}


/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */
              
/*
 * CVDlsSetDenseJacFn specifies the dense Jacobian function.
 */
int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn jac)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsSetDenseJacFn", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsSetDenseJacFn", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  if (jac != NULL) {
    cvdls_mem->d_jacDQ = FALSE;
    cvdls_mem->d_djac = jac;
  } else {
    cvdls_mem->d_jacDQ = TRUE;
  }

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsSetBandJacFn specifies the band Jacobian function.
 */
int CVDlsSetBandJacFn(void *cvode_mem, CVDlsBandJacFn jac)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsSetBandJacFn", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsSetBandJacFn", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  if (jac != NULL) {
    cvdls_mem->d_jacDQ = FALSE;
    cvdls_mem->d_bjac = jac;
  } else {
    cvdls_mem->d_jacDQ = TRUE;
  }

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetWorkSpace returns the length of workspace allocated for the
 * CVDLS linear solver.
 */
int CVDlsGetWorkSpace(void *cvode_mem, sunindextype *lenrwLS, sunindextype *leniwLS)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetWorkSpace", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetWorkSpace", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  if (cvdls_mem->d_type == SUNDIALS_DENSE) {
    *lenrwLS = 2*cvdls_mem->d_n*cvdls_mem->d_n;
    *leniwLS = cvdls_mem->d_n;
  } else if (cvdls_mem->d_type == SUNDIALS_BAND) {
    *lenrwLS = cvdls_mem->d_n*(cvdls_mem->d_smu + cvdls_mem->d_mu + 2*cvdls_mem->d_ml + 2);
    *leniwLS = cvdls_mem->d_n;
  }

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetNumJacEvals", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  *njevals = cvdls_mem->d_nje;

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetNumRhsEvals returns the number of calls to the ODE function
 * needed for the DQ Jacobian approximation.
 */
int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetNumRhsEvals", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetNumRhsEvals", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  *nfevalsLS = cvdls_mem->d_nfeDQ;

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetReturnFlagName returns the name associated with a CVDLS
 * return value.
 */
char *CVDlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVDLS_SUCCESS:
    sprintf(name,"CVDLS_SUCCESS");
    break;   
  case CVDLS_MEM_NULL:
    sprintf(name,"CVDLS_MEM_NULL");
    break;
  case CVDLS_LMEM_NULL:
    sprintf(name,"CVDLS_LMEM_NULL");
    break;
  case CVDLS_ILL_INPUT:
    sprintf(name,"CVDLS_ILL_INPUT");
    break;
  case CVDLS_MEM_FAIL:
    sprintf(name,"CVDLS_MEM_FAIL");
    break;
  case CVDLS_JACFUNC_UNRECVR:
    sprintf(name,"CVDLS_JACFUNC_UNRECVR");
    break;
  case CVDLS_JACFUNC_RECVR:
    sprintf(name,"CVDLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CVDlsGetLastFlag returns the last flag set in a CVDLS function.
 */
int CVDlsGetLastFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetLastFlag", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetLastFlag", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  *flag = cvdls_mem->d_last_flag;

  return(CVDLS_SUCCESS);
}

/* 
 * =================================================================
 * DQ JACOBIAN APPROXIMATIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cvDlsDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 

int cvDlsDenseDQJac(sunindextype N, realtype t,
                    N_Vector y, N_Vector fy, 
                    DlsMat Jac, void *data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  sunindextype j;
  int retval = 0;

  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* data points to cvode_mem */
  cv_mem = (CVodeMem) data;
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(cv_mem->cv_ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(cv_mem->cv_uround);
  fnorm = N_VWrmsNorm(fy, cv_mem->cv_ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */

    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = SUNMAX(srur*SUNRabs(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = cv_mem->cv_f(t, y, ftemp, cv_mem->cv_user_data);
    cvdls_mem->d_nfeDQ++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}

/*
 * -----------------------------------------------------------------
 * cvDlsBandDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation to
 * the Jacobian of f(t,y).  It assumes that a band matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. This makes it possible to get the address of a column
 * of J via the macro BAND_COL and to write a simple for loop to set
 * each of the elements of a column in succession.
 * -----------------------------------------------------------------
 */

int cvDlsBandDQJac(sunindextype N, sunindextype mupper, sunindextype mlower,
                   realtype t, N_Vector y, N_Vector fy, 
                   DlsMat Jac, void *data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  sunindextype group, i, j, width, ngroups, i1, i2;
  int retval = 0;

  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* data points to cvode_mem */
  cv_mem = (CVodeMem) data;
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(cv_mem->cv_ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(cv_mem->cv_uround);
  fnorm = N_VWrmsNorm(fy, cv_mem->cv_ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    retval = cv_mem->cv_f(cv_mem->cv_tn, ytemp, ftemp, cv_mem->cv_user_data);
    cvdls_mem->d_nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(Jac,j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  return(retval);
}

int cvDlsInitializeCounters(CVDlsMem cvdls_mem)
{
  cvdls_mem->d_nje   = 0;
  cvdls_mem->d_nfeDQ = 0;
  cvdls_mem->d_nstlj = 0;
  return(0);
}
