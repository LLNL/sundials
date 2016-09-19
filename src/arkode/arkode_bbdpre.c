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
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with ARKODE, a ARKSPILS 
 * linear solver, and the parallel implementation of NVECTOR.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_bbdpre_impl.h"
#include "arkode_spils_impl.h"
#include <sundials/sundials_math.h>


#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/* Prototypes of functions ARKBBDPrecSetup and ARKBBDPrecSolve */
static int ARKBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
			   booleantype jok, booleantype *jcurPtr, 
			   realtype gamma, void *bbd_data, 
			   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int ARKBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
			   N_Vector r, N_Vector z, 
			   realtype gamma, realtype delta,
			   int lr, void *bbd_data, N_Vector tmp);

/* Prototype for ARKBBDPrecFree */
static int ARKBBDPrecFree(ARKodeMem ark_mem);

/* Prototype for difference quotient Jacobian calculation routine */
static int ARKBBDDQJac(ARKBBDPrecData pdata, realtype t, 
		       N_Vector y, N_Vector gy, 
		       N_Vector ytemp, N_Vector gtemp);


/*---------------------------------------------------------------
 User-Callable Functions: initialization, reinit and free
---------------------------------------------------------------*/
int ARKBBDPrecInit(void *arkode_mem, long int Nlocal, 
                   long int mudq, long int mldq,
                   long int mukeep, long int mlkeep, 
                   realtype dqrely, 
                   ARKLocalFn gloc, ARKCommFn cfn)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;
  long int muk, mlk, storage_mu;
  int flag;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Test if the NVECTOR package is compatible with the BLOCK 
     BAND preconditioner */
  if(ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (ARKBBDPrecData) malloc(sizeof *pdata);  
  if (pdata == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set pointers to gloc and cfn; load half-bandwidths */
  pdata->arkode_mem = arkode_mem;
  pdata->gloc = gloc;
  pdata->cfn = cfn;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));
  muk = SUNMIN(Nlocal-1, SUNMAX(0,mukeep));
  mlk = SUNMIN(Nlocal-1, SUNMAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = NewBandMat(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) { 
    free(pdata); pdata = NULL; 
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL); 
  }

  /* Allocate memory for preconditioner matrix */
  storage_mu = SUNMIN(Nlocal-1, muk + mlk);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }
  /* Allocate memory for lpivots */
  pdata->lpivots = NULL;
  pdata->lpivots = NewLintArray(Nlocal);
  if (pdata->lpivots == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", 
		    "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : SUNRsqrt(ark_mem->ark_uround);

  /* Store Nlocal to be used in ARKBBDPrecSetup */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(muk + 2*mlk + storage_mu + 2);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  /* make sure s_P_data is free from any previous allocations */
  if (arkspils_mem->s_pfree != NULL) {
    arkspils_mem->s_pfree(ark_mem);
  }

  /* Point to the new P_data field in the SPILS memory */
  arkspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  arkspils_mem->s_pfree = ARKBBDPrecFree;

  /* Attach preconditioner solve and setup functions */
  flag = ARKSpilsSetPreconditioner(arkode_mem, ARKBBDPrecSetup, 
				   ARKBBDPrecSolve);

  return(flag);
}


/*-------------------------------------------------------------*/
int ARKBBDPrecReInit(void *arkode_mem, long int mudq, 
		     long int mldq, realtype dqrely)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;
  long int Nlocal;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecReInit", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecReInit", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Test if the preconditioner data is non-NULL */
  if (arkspils_mem->s_P_data == NULL) {
    arkProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecReInit", MSGBBD_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  /* Load half-bandwidths */
  Nlocal = pdata->n_local;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : SUNRsqrt(ark_mem->ark_uround);

  /* Re-initialize nge */
  pdata->nge = 0;

  return(ARKSPILS_SUCCESS);
}


/*-------------------------------------------------------------*/
int ARKBBDPrecGetWorkSpace(void *arkode_mem, long int *lenrwBBDP, 
			   long int *leniwBBDP)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecGetWorkSpace", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecGetWorkSpace", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (arkspils_mem->s_P_data == NULL) {
    arkProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecGetWorkSpace", MSGBBD_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(ARKSPILS_SUCCESS);
}


/*-------------------------------------------------------------*/
int ARKBBDPrecGetNumGfnEvals(void *arkode_mem, long int *ngevalsBBDP)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecGetNumGfnEvals", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecGetNumGfnEvals", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (arkspils_mem->s_P_data == NULL) {
    arkProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBBDPRE", 
		    "ARKBBDPrecGetNumGfnEvals", MSGBBD_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  *ngevalsBBDP = pdata->nge;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKBBDPrecSetup:

 ARKBBDPrecSetup generates and factors a banded block of the
 preconditioner matrix on each processor, via calls to the
 user-supplied gloc and cfn functions. It uses difference
 quotient approximations to the Jacobian elements.

 ARKBBDPrecSetup calculates a new J,if necessary, then calculates
 P = I - gamma*J, and does an LU factorization of P.

 The parameters of ARKBBDPrecSetup used here are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
         namely the predicted value of y(t).

 fy      is the vector f(t,y).

 jok     is an input flag indicating whether Jacobian-related
         data needs to be recomputed, as follows:
           jok == FALSE means recompute Jacobian-related data
                  from scratch.
           jok == TRUE  means that Jacobian data from the
                  previous ARKBBDPrecon call can be reused
                  (with the current value of gamma).
         A ARKBBDPrecon call with jok == TRUE should only occur
         after a call with jok == FALSE.

 jcurPtr is a pointer to an output integer flag which is
         set by ARKBBDPrecon as follows:
           *jcurPtr = TRUE if Jacobian data was recomputed.
           *jcurPtr = FALSE if Jacobian data was not recomputed,
                      but saved data was reused.

 gamma   is the scalar appearing in the Newton matrix.

 bbd_data is a pointer to the preconditioner data set by
          ARKBBDPrecInit

 tmp1, tmp2, and tmp3 are pointers to memory allocated
           for NVectors which are be used by ARKBBDPrecSetup
           as temporary storage or work space.

 Return value:
 The value returned by this ARKBBDPrecSetup function is the int
   0  if successful,
   1  for a recoverable error (step will be retried).
---------------------------------------------------------------*/
static int ARKBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
			   booleantype jok, booleantype *jcurPtr, 
			   realtype gamma, void *bbd_data, 
			   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int ier;
  ARKBBDPrecData pdata;
  ARKodeMem ark_mem;
  int retval;

  pdata = (ARKBBDPrecData) bbd_data;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* If jok = TRUE, use saved copy of J */
  if (jok) {
    *jcurPtr = FALSE;
    BandCopy(pdata->savedJ, pdata->savedP, pdata->mukeep, pdata->mlkeep);

  /* Otherwise call ARKBBDDQJac for new J value */
  } else {
    *jcurPtr = TRUE;
    SetToZero(pdata->savedJ);

    retval = ARKBBDDQJac(pdata, t, y, tmp1, tmp2, tmp3);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBBDPRE", "ARKBBDPrecSetup", 
		      MSGBBD_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }
    BandCopy(pdata->savedJ, pdata->savedP, pdata->mukeep, pdata->mlkeep);

  }
  
  /* Scale and add I to get P = I - gamma*J */
  BandScale(-gamma, pdata->savedP);
  AddIdentity(pdata->savedP);
 
  /* Do LU factorization of P in place */
  ier = BandGBTRF(pdata->savedP, pdata->lpivots);
 
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}


/*---------------------------------------------------------------
 ARKBBDPrecSolve:

 ARKBBDPrecSolve solves a linear system P z = r, with the
 band-block-diagonal preconditioner matrix P generated and
 factored by ARKBBDPrecSetup.

 The parameters of ARKBBDPrecSolve used here are as follows:

 r is the right-hand side vector of the linear system.

 bbd_data is a pointer to the preconditioner data set by
   ARKBBDPrecInit.

 z is the output vector computed by ARKBBDPrecSolve.

 The value returned by the ARKBBDPrecSolve function is always 0,
 indicating success.
---------------------------------------------------------------*/
static int ARKBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
			   N_Vector r, N_Vector z, 
			   realtype gamma, realtype delta,
			   int lr, void *bbd_data, N_Vector tmp)
{
  ARKBBDPrecData pdata;
  realtype *zd;

  pdata = (ARKBBDPrecData) bbd_data;

  /* Copy r to z, then do backsolve and return */
  N_VScale(ONE, r, z);
  zd = N_VGetArrayPointer(z);
  BandGBTRS(pdata->savedP, pdata->lpivots, zd);

  return(0);
}


/*-------------------------------------------------------------*/
static int ARKBBDPrecFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;
  
  if (ark_mem->ark_lmem == NULL) return(0);
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  
  if (arkspils_mem->s_P_data == NULL) return(0);
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  DestroyMat(pdata->savedJ);
  DestroyMat(pdata->savedP);
  DestroyArray(pdata->lpivots);

  free(pdata);
  pdata = NULL;

  return(0);
}


/*---------------------------------------------------------------
 ARKBBDDQJac:

 This routine generates a banded difference quotient approximation
 to the local block of the Jacobian of g(t,y). It assumes that a
 band matrix of type DlsMat is stored columnwise, and that elements
 within each column are contiguous. All matrix elements are generated
 as difference quotients, by way of calls to the user routine gloc.
 By virtue of the band structure, the number of these calls is
 bandwidth + 1, where bandwidth = mldq + mudq + 1.
 But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 This routine also assumes that the local elements of a vector are
 stored contiguously.
---------------------------------------------------------------*/
static int ARKBBDDQJac(ARKBBDPrecData pdata, realtype t, 
		       N_Vector y, N_Vector gy, 
		       N_Vector ytemp, N_Vector gtemp)
{
  ARKodeMem ark_mem;
  realtype gnorm, minInc, inc, inc_inv;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;
  int retval;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and gloc to get base value of g(t,y) */
  if (pdata->cfn != NULL) {
    retval = pdata->cfn(pdata->n_local, t, y, ark_mem->ark_user_data);
    if (retval != 0) return(retval);
  }

  retval = pdata->gloc(pdata->n_local, t, ytemp, gy, ark_mem->ark_user_data);
  pdata->nge++;
  if (retval != 0) return(retval);

  /* Obtain pointers to the data for various vectors */
  y_data     =  N_VGetArrayPointer(y);
  gy_data    =  N_VGetArrayPointer(gy);
  ewt_data   =  N_VGetArrayPointer(ark_mem->ark_ewt);
  ytemp_data =  N_VGetArrayPointer(ytemp);
  gtemp_data =  N_VGetArrayPointer(gtemp);

  /* Set minimum increment based on uround and norm of g */
  /* gnorm = N_VWrmsNorm(gy, ark_mem->ark_ewt); */
  gnorm = N_VWrmsNorm(gy, ark_mem->ark_rwt);
  minInc = (gnorm != ZERO) ? (MIN_INC_MULT * SUNRabs(ark_mem->ark_h) *
			      ark_mem->ark_uround * pdata->n_local * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = pdata->mldq + pdata->mudq + 1;
  ngroups = SUNMIN(width, pdata->n_local);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < pdata->n_local; j+=width) {
      inc = SUNMAX(pdata->dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    retval = pdata->gloc(pdata->n_local, t, ytemp, gtemp, 
			 ark_mem->ark_user_data);
    pdata->nge++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < pdata->n_local; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(pdata->savedJ,j);
      inc = SUNMAX(pdata->dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-pdata->mukeep);
      i2 = SUNMIN(j+pdata->mlkeep, pdata->n_local-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }

  return(0);
}



/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
