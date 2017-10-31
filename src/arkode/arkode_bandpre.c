/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * Based off of cvode_bandpre.c by Scott D. Cohen, 
 *      Alan C. Hindmarsh, Radu Serban, and Aaron Collier @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the ARKSPILS linear solver interface.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_bandpre_impl.h"
#include "arkode_spils_impl.h"
#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/* Prototypes of ARKBandPrecSetup and ARKBandPrecSolve */
static int ARKBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
			    booleantype jok, booleantype *jcurPtr, 
			    realtype gamma, void *bp_data);
static int ARKBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
			    N_Vector r, N_Vector z, 
			    realtype gamma, realtype delta,
			    int lr, void *bp_data);

/* Prototype for ARKBandPrecFree */
static int ARKBandPrecFree(ARKodeMem ark_mem);

/* Prototype for difference quotient Jacobian calculation routine */
static int ARKBandPDQJac(ARKBandPrecData pdata,
			 realtype t, N_Vector y, N_Vector fy, 
			 N_Vector ftemp, N_Vector ytemp);


/*---------------------------------------------------------------
 Initialization, Free, and Get Functions
 NOTE: The band linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKBandPrecInit will
       first test for a compatible N_Vector internal 
       representation by checking that the function 
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKBandPrecInit(void *arkode_mem, sunindextype N, 
		    sunindextype mu, sunindextype ml)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBandPrecData pdata;
  sunindextype mup, mlp, storagemu;
  int flag;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the SPILS linear solver interface has been created */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Test compatibility of NVECTOR package with the BAND preconditioner */
  if(ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (ARKBandPrecData) malloc(sizeof *pdata);
  if (pdata == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->arkode_mem = arkode_mem;
  pdata->N = N;
  pdata->mu = mup = SUNMIN(N-1, SUNMAX(0,mu));
  pdata->ml = mlp = SUNMIN(N-1, SUNMAX(0,ml));

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = NULL;
  pdata->savedJ = SUNBandMatrix(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = SUNMIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = SUNBandMatrix(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Allocate memory for banded linear solver */
  pdata->LS = NULL;
  pdata->LS = SUNBandLinearSolver(ark_mem->ark_tempv, pdata->savedP);
  if (pdata->LS == NULL) {
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }
  
  /* allocate memory for temporary N_Vectors */
  pdata->tmp1 = NULL;
  pdata->tmp1 = N_VClone(ark_mem->ark_tempv);
  if (pdata->tmp1 == NULL) {
    SUNLinSolFree(pdata->LS);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }
  pdata->tmp2 = NULL;
  pdata->tmp2 = N_VClone(ark_mem->ark_tempv);
  if (pdata->tmp2 == NULL) {
    SUNLinSolFree(pdata->LS);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    N_VDestroy(pdata->tmp1);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* initialize band linear solver object */
  flag = SUNLinSolInitialize(pdata->LS);
  if (flag != SUNLS_SUCCESS) {
    SUNLinSolFree(pdata->LS);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    N_VDestroy(pdata->tmp1);
    N_VDestroy(pdata->tmp2);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSGBP_SUNLS_FAIL);
    return(ARKSPILS_SUNLS_FAIL);
  }
  
  /* make sure s_P_data is free from any previous allocations */
  if (arkspils_mem->pfree)
    arkspils_mem->pfree(ark_mem);

  /* Point to the new P_data field in the SPILS memory */
  arkspils_mem->P_data = pdata;

  /* Attach the pfree function */
  arkspils_mem->pfree = ARKBandPrecFree;

  /* Attach preconditioner solve and setup functions */
  flag = ARKSpilsSetPreconditioner(arkode_mem, 
                                   ARKBandPrecSetup, 
                                   ARKBandPrecSolve);
  return(flag);
}


int ARKBandPrecGetWorkSpace(void *arkode_mem, long int *lenrwBP, 
                            long int *leniwBP)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBandPrecData pdata;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;
  
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetWorkSpace", MSGBP_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetWorkSpace", MSGBP_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (arkspils_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetWorkSpace", MSGBP_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBandPrecData) arkspils_mem->P_data;

  /* sum space requirements for all objects in pdata */
  *leniwBP = 4;
  *lenrwBP = 0;
  if (ark_mem->ark_tempv->ops->nvspace) {
    N_VSpace(ark_mem->ark_tempv, &lrw1, &liw1);
    *leniwBP += 2*liw1;
    *lenrwBP += 2*lrw1;
  }
  if (pdata->savedJ->ops->space) {
    flag = SUNMatSpace(pdata->savedJ, &lrw, &liw);
    *leniwBP += liw;
    *lenrwBP += lrw;
  }
  if (pdata->savedP->ops->space) {
    flag = SUNMatSpace(pdata->savedP, &lrw, &liw);
    *leniwBP += liw;
    *lenrwBP += lrw;
  }
  if (pdata->LS->ops->space) {
    flag = SUNLinSolSpace(pdata->LS, &lrw, &liw);
    *leniwBP += liw;
    *lenrwBP += lrw;
  }

  return(ARKSPILS_SUCCESS);
}


int ARKBandPrecGetNumRhsEvals(void *arkode_mem, long int *nfevalsBP)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBandPrecData pdata;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetNumRhsEvals", MSGBP_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetNumRhsEvals", MSGBP_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (arkspils_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetNumRhsEvals", MSGBP_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBandPrecData) arkspils_mem->P_data;

  *nfevalsBP = pdata->nfeBP;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKBandPrecSetup:

 Together ARKBandPrecSetup and ARKBandPrecSolve use a banded
 difference quotient Jacobian to create a preconditioner.
 ARKBandPrecSetup calculates a new J, if necessary, then
 calculates P = I - gamma*J, and does an LU factorization of P.

 The parameters of ARKBandPrecSetup are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
         namely the predicted value of y(t).

 fy      is the vector f(t,y).

 jok     is an input flag indicating whether Jacobian-related
         data needs to be recomputed, as follows:
           jok == SUNFALSE means recompute Jacobian-related data
                  from scratch.
           jok == SUNTRUE means that Jacobian data from the
                  previous PrecSetup call will be reused
                  (with the current value of gamma).
         A ARKBandPrecSetup call with jok == SUNTRUE should only
         occur after a call with jok == SUNFALSE.

 *jcurPtr is a pointer to an output integer flag which is
          set by ARKBandPrecond as follows:
            *jcurPtr = SUNTRUE if Jacobian data was recomputed.
            *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
                       but saved data was reused.

 gamma   is the scalar appearing in the Newton matrix.

 bp_data is a pointer to preconditoner data (set by ARKBandPrecInit)

 The value to be returned by the ARKBandPrecSetup function is
   0  if successful, or
   1  if the band factorization failed.
---------------------------------------------------------------*/
static int ARKBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data)
{
  ARKBandPrecData pdata;
  ARKodeMem ark_mem;
  int retval;
  sunindextype ier;

  /* Assume matrix and lpivots have already been allocated. */
  pdata = (ARKBandPrecData) bp_data;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  if (jok) {

    /* If jok = SUNTRUE, use saved copy of J. */
    *jcurPtr = SUNFALSE;
    retval = SUNMatCopy(pdata->savedJ, pdata->savedP);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSGBP_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

  } else {

    /* If jok = SUNFALSE, call ARKBandPDQJac for new J value. */
    *jcurPtr = SUNTRUE;
    retval = SUNMatZero(pdata->savedJ);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSGBP_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    retval = ARKBandPDQJac(pdata, t, y, fy, 
                           pdata->tmp1, pdata->tmp2);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSGBP_RHSFUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    retval = SUNMatCopy(pdata->savedJ, pdata->savedP);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSGBP_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

  }
  
  /* Scale and add identity to get savedP = I - gamma*J. */
  retval = SUNMatScaleAddI(-gamma, pdata->savedP);
  if (retval) {
    arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                    "ARKBandPrecSetup", MSGBP_SUNMAT_FAIL);
    return(-1);
  }

  /* Do LU factorization of matrix and return error flag */
  ier = SUNLinSolSetup_Band(pdata->LS, pdata->savedP);
  return(ier);
}


/*---------------------------------------------------------------
 ARKBandPrecSolve:

 ARKBandPrecSolve solves a linear system P z = r, where P is the
 matrix computed by ARKBandPrecond.

 The parameters of ARKBandPrecSolve used here are as follows:

 r is the right-hand side vector of the linear system.

 bp_data is a pointer to preconditoner data (set by ARKBandPrecInit)

 z is the output vector computed by ARKBandPrecSolve.

 The value returned by the ARKBandPrecSolve function is always 0,
 indicating success.
---------------------------------------------------------------*/ 
static int ARKBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                            N_Vector r, N_Vector z, 
                            realtype gamma, realtype delta,
                            int lr, void *bp_data)
{
  ARKBandPrecData pdata;
  int retval;

  /* Assume matrix and linear solver have already been allocated. */
  pdata = (ARKBandPrecData) bp_data;

  /* Call banded solver object to do the work */
  retval = SUNLinSolSolve(pdata->LS, pdata->savedP, z, r, ZERO);
  return(retval);
}


/*---------------------------------------------------------------
 ARKBandPrecFree:

 Frees data associated with the ARKBand preconditioner.
---------------------------------------------------------------*/ 
static int ARKBandPrecFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  ARKBandPrecData pdata;

  if (ark_mem->ark_lmem == NULL) return(0);
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  
  if (arkspils_mem->P_data == NULL) return(0);
  pdata = (ARKBandPrecData) arkspils_mem->P_data;

  SUNLinSolFree(pdata->LS);
  SUNMatDestroy(pdata->savedP);
  SUNMatDestroy(pdata->savedJ);
  N_VDestroy(pdata->tmp1);
  N_VDestroy(pdata->tmp2);

  free(pdata);
  pdata = NULL;

  return(0);
}


/*---------------------------------------------------------------
 ARKBandPDQJac:

 This routine generates a banded difference quotient approximation to
 the Jacobian of f(t,y). It assumes that a band matrix of type
 DlsMat is stored column-wise, and that elements within each column
 are contiguous. This makes it possible to get the address of a column
 of J via the macro BAND_COL and to write a simple for loop to set
 each of the elements of a column in succession.
---------------------------------------------------------------*/
static int ARKBandPDQJac(ARKBandPrecData pdata, 
			 realtype t, N_Vector y, N_Vector fy, 
			 N_Vector ftemp, N_Vector ytemp)
{
  ARKodeMem ark_mem;
  realtype fnorm, minInc, inc, inc_inv, srur;
  sunindextype group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetArrayPointer(ark_mem->ark_ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = SUNRsqrt(ark_mem->ark_uround);
  fnorm = N_VWrmsNorm(fy, ark_mem->ark_rwt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(ark_mem->ark_h) * 
     ark_mem->ark_uround * pdata->N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = pdata->ml + pdata->mu + 1;
  ngroups = SUNMIN(width, pdata->N);
  
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < pdata->N; j += width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */
    retval = ark_mem->ark_fi(t, ytemp, ftemp, ark_mem->ark_user_data);
    pdata->nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < pdata->N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = SUNBandMatrix_Column(pdata->savedJ,j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-pdata->mu);
      i2 = SUNMIN(j+pdata->ml, pdata->N-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) =
          inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(0);
}


/*---------------------------------------------------------------
   EOF
---------------------------------------------------------------*/ 
