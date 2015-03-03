/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 *---------------------------------------------------------------
 * This is the implementation file for the ARKDLS linear solvers
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_direct_impl.h"
#include <sundials/sundials_math.h>

/*===============================================================
 FUNCTION SPECIFIC CONSTANTS
===============================================================*/

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/*===============================================================
 EXPORTED FUNCTIONS
===============================================================*/
              
/*---------------------------------------------------------------
 ARKDlsSetDenseJacFn specifies the dense Jacobian function.
---------------------------------------------------------------*/
int ARKDlsSetDenseJacFn(void *arkode_mem, ARKDlsDenseJacFn jac)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetDenseJacFn", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsSetDenseJacFn", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  if (jac != NULL) {
    arkdls_mem->d_jacDQ = FALSE;
    arkdls_mem->d_djac  = jac;
  } else {
    arkdls_mem->d_jacDQ = TRUE;
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsSetBandJacFn specifies the band Jacobian function.
---------------------------------------------------------------*/
int ARKDlsSetBandJacFn(void *arkode_mem, ARKDlsBandJacFn jac)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetBandJacFn", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsSetBandJacFn", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  if (jac != NULL) {
    arkdls_mem->d_jacDQ = FALSE;
    arkdls_mem->d_bjac  = jac;
  } else {
    arkdls_mem->d_jacDQ = TRUE;
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsSetDenseMassFn specifies the dense mass matrix function.
---------------------------------------------------------------*/
int ARKDlsSetDenseMassFn(void *arkode_mem, ARKDlsDenseMassFn mass)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetDenseMassFn", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsSetDenseMassFn", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  if (mass != NULL) {
    arkdls_mem->d_dmass = mass;
  } else {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
		    "ARKDlsSetDenseMassFn", "DenseMassFn must be non-NULL");
    return(ARKDLS_ILL_INPUT);
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsSetBandMassFn specifies the band mass matrix function.
---------------------------------------------------------------*/
int ARKDlsSetBandMassFn(void *arkode_mem, ARKDlsBandMassFn mass)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetBandMassFn", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsSetBandMassFn", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  if (mass != NULL) {
    arkdls_mem->d_bmass = mass;
  } else {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
		    "ARKDlsSetBandMassFn", "BandMassFn must be non-NULL");
    return(ARKDLS_MASSMEM_NULL);
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetWorkSpace returns the length of workspace allocated for 
 the ARKDLS linear solver.
---------------------------------------------------------------*/
int ARKDlsGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetWorkSpace", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetWorkSpace", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  if (arkdls_mem->d_type == SUNDIALS_DENSE) {
    *lenrwLS = 2*arkdls_mem->d_n*arkdls_mem->d_n;
    *leniwLS = arkdls_mem->d_n;
  } else if (arkdls_mem->d_type == SUNDIALS_BAND) {
    *lenrwLS = arkdls_mem->d_n*(arkdls_mem->d_smu + arkdls_mem->d_mu + 2*arkdls_mem->d_ml + 2);
    *leniwLS = arkdls_mem->d_n;
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetMassWorkSpace returns the length of workspace allocated 
 for the ARKDLS mass matrix linear solver.
---------------------------------------------------------------*/
int ARKDlsGetMassWorkSpace(void *arkode_mem, long int *lenrwMLS, 
			   long int *leniwMLS)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetMassWorkSpace", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetMassWorkSpace", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  if (arkdls_mem->d_type == SUNDIALS_DENSE) {
    *lenrwMLS = 2*arkdls_mem->d_n*arkdls_mem->d_n;
    *leniwMLS = arkdls_mem->d_n;
  } else if (arkdls_mem->d_type == SUNDIALS_BAND) {
    *lenrwMLS = arkdls_mem->d_n*(arkdls_mem->d_smu + arkdls_mem->d_mu + 2*arkdls_mem->d_ml + 2);
    *leniwMLS = arkdls_mem->d_n;
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumJacEvals returns the number of Jacobian evaluations.
---------------------------------------------------------------*/
int ARKDlsGetNumJacEvals(void *arkode_mem, long int *njevals)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumJacEvals", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  *njevals = arkdls_mem->d_nje;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumMassEvals returns the number of mass matrix evaluations.
---------------------------------------------------------------*/
int ARKDlsGetNumMassEvals(void *arkode_mem, long int *nmevals)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassEvals", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassEvals", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  *nmevals = arkdls_mem->d_nme;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumRhsEvals returns the number of calls to the ODE function
 needed for the DQ Jacobian approximation.
---------------------------------------------------------------*/
int ARKDlsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumRhsEvals", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumRhsEvals", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  *nfevalsLS = arkdls_mem->d_nfeDQ;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetReturnFlagName returns the name associated with a ARKDLS
 return value.
---------------------------------------------------------------*/
char *ARKDlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKDLS_SUCCESS:
    sprintf(name,"ARKDLS_SUCCESS");
    break;   
  case ARKDLS_MEM_NULL:
    sprintf(name,"ARKDLS_MEM_NULL");
    break;
  case ARKDLS_LMEM_NULL:
    sprintf(name,"ARKDLS_LMEM_NULL");
    break;
  case ARKDLS_MASSMEM_NULL:
    sprintf(name,"ARKDLS_MASSMEM_NULL");
    break;
  case ARKDLS_ILL_INPUT:
    sprintf(name,"ARKDLS_ILL_INPUT");
    break;
  case ARKDLS_MEM_FAIL:
    sprintf(name,"ARKDLS_MEM_FAIL");
    break;
  case ARKDLS_JACFUNC_UNRECVR:
    sprintf(name,"ARKDLS_JACFUNC_UNRECVR");
    break;
  case ARKDLS_JACFUNC_RECVR:
    sprintf(name,"ARKDLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
 ARKDlsGetLastFlag returns the last flag set in a ARKDLS function.
---------------------------------------------------------------*/
int ARKDlsGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastFlag", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastFlag", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  *flag = arkdls_mem->d_last_flag;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetLastMassFlag returns the last flag set in a ARKDLS mass
  matrix function.
---------------------------------------------------------------*/
int ARKDlsGetLastMassFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastMassFlag", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastMassFlag", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  *flag = arkdls_mem->d_last_flag;

  return(ARKDLS_SUCCESS);
}


/*===============================================================
 DQ JACOBIAN APPROXIMATIONS
===============================================================*/

/*---------------------------------------------------------------
 arkDlsDenseDQJac:

 This routine generates a dense difference quotient approximation 
 to the Jacobian of f(t,y). It assumes that a dense matrix of 
 type DlsMat is stored column-wise, and that elements within each 
 column are contiguous. The address of the jth column of J is 
 obtained via the macro DENSE_COL and this pointer is associated 
 with an N_Vector using the N_VGetArrayPointer/N_VSetArrayPointer 
 functions.  Finally, the actual computation of the jth column of
 the Jacobian is done with a call to N_VLinearSum.
---------------------------------------------------------------*/
int arkDlsDenseDQJac(long int N, realtype t, N_Vector y, 
		     N_Vector fy, DlsMat Jac, void *data,
		     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  long int j;
  int retval = 0;

  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* data points to arkode_mem */
  ark_mem = (ARKodeMem) data;
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ark_mem->ark_ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(ark_mem->ark_uround);
  /* fnorm = N_VWrmsNorm(fy, ark_mem->ark_ewt); */
  fnorm = N_VWrmsNorm(fy, ark_mem->ark_rwt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(ark_mem->ark_h) * ark_mem->ark_uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */
    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = SUNMAX(srur*SUNRabs(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = ark_mem->ark_fi(t, y, ftemp, ark_mem->ark_user_data);
    arkdls_mem->d_nfeDQ++;
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


/*---------------------------------------------------------------
 arkDlsBandDQJac:

 This routine generates a banded difference quotient approximation 
 to the Jacobian of f(t,y).  It assumes that a band matrix of type
 DlsMat is stored column-wise, and that elements within each 
 column are contiguous. This makes it possible to get the address 
 of a column of J via the macro BAND_COL and to write a simple for 
 loop to set each of the elements of a column in succession.
---------------------------------------------------------------*/
int arkDlsBandDQJac(long int N, long int mupper, long int mlower,
                   realtype t, N_Vector y, N_Vector fy, 
                   DlsMat Jac, void *data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  long int group, i, j, width, ngroups, i1, i2;
  int retval = 0;

  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* data points to arkode_mem */
  ark_mem = (ARKodeMem) data;
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ark_mem->ark_ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(ark_mem->ark_uround);
  /* fnorm = N_VWrmsNorm(fy, ark_mem->ark_ewt); */
  fnorm = N_VWrmsNorm(fy, ark_mem->ark_rwt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(ark_mem->ark_h) * ark_mem->ark_uround * N * fnorm) : ONE;

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
    retval = ark_mem->ark_fi(ark_mem->ark_tn, ytemp, ftemp, 
			     ark_mem->ark_user_data);
    arkdls_mem->d_nfeDQ++;
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


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
