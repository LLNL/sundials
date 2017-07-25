/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
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
 * This is the implementation file for the CVBAND linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvode/cvode_band.h>
#include "cvode_direct_impl.h"
#include "cvode_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* CVBAND linit, lsetup, lsolve, and lfree routines */

static int cvBandInit(CVodeMem cv_mem);

static int cvBandSetup(CVodeMem cv_mem, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);

static int cvBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);

static int cvBandFree(CVodeMem cv_mem);

/*
 * -----------------------------------------------------------------
 * CVBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  CVBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be cvBandInit, cvBandSetup, cvBandSolve, and cvBandFree,
 * respectively.  It allocates memory for a structure of type
 * CVDlsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to be
 * TRUE, d_mu to be mupper, d_ml to be mlower, and the d_bjac field to be 
 * cvDlsBandDQJac.
 * Finally, it allocates memory for M, savedJ, and pivot.  The CVBand
 * return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int CVBand(void *cvode_mem, sunindextype N, sunindextype mupper, sunindextype mlower)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVBAND", "CVBand", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (cv_mem->cv_tempv->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVBAND", "CVBand", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }

  if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem);

  /* Set four main function fields in cv_mem */  
  cv_mem->cv_linit  = cvBandInit;
  cv_mem->cv_lsetup = cvBandSetup;
  cv_mem->cv_lsolve = cvBandSolve;
  cv_mem->cv_lfree  = cvBandFree;
  
  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVBAND", "CVBand", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  cvdls_mem->d_type = SUNDIALS_BAND;
  
  /* Initialize Jacobian-related data */
  cvdls_mem->d_jacDQ = TRUE;
  cvdls_mem->d_bjac = NULL;
  cvdls_mem->d_J_data = NULL;

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;

  cvDlsInitializeCounters(cvdls_mem);

  cv_mem->cv_setupNonNull = TRUE;
  
  /* Load problem dimension */
  cvdls_mem->d_n = N;

  /* Load half-bandwiths in cvdls_mem */
  cvdls_mem->d_ml = mlower;
  cvdls_mem->d_mu = mupper;

  /* Test ml and mu for legality */
  if ((cvdls_mem->d_ml < 0) || (cvdls_mem->d_mu < 0) ||
      (cvdls_mem->d_ml >= N) || (cvdls_mem->d_mu >= N)) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVBAND", "CVBand", MSGD_BAD_SIZES);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  cvdls_mem->d_smu = SUNMIN(N-1, cvdls_mem->d_mu + cvdls_mem->d_ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  cvdls_mem->d_M = NULL;
  cvdls_mem->d_M = NewBandMat(N, cvdls_mem->d_mu, cvdls_mem->d_ml, cvdls_mem->d_smu);
  if (cvdls_mem->d_M == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVBAND", "CVBand", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  cvdls_mem->d_savedJ = NULL;
  cvdls_mem->d_savedJ = NewBandMat(N, cvdls_mem->d_mu, cvdls_mem->d_ml, cvdls_mem->d_mu);
  if (cvdls_mem->d_savedJ == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVBAND", "CVBand", MSGD_MEM_FAIL);
    DestroyMat(cvdls_mem->d_M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  cvdls_mem->d_lpivots = NULL;
  cvdls_mem->d_lpivots = NewIndexArray(N);
  if (cvdls_mem->d_lpivots == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVBAND", "CVBand", MSGD_MEM_FAIL);
    DestroyMat(cvdls_mem->d_M);
    DestroyMat(cvdls_mem->d_savedJ);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * cvBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cvBandInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  cvDlsInitializeCounters(cvdls_mem);

  /* Set Jacobian function and data, depending on jacDQ */
  if (cvdls_mem->d_jacDQ) {
    cvdls_mem->d_bjac = cvDlsBandDQJac;
    cvdls_mem->d_J_data = cv_mem;
  } else {
    cvdls_mem->d_J_data = cv_mem->cv_user_data;
  }

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int cvBandSetup(CVodeMem cv_mem, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3)
{
  sunindextype ier;
  CVDlsMem cvdls_mem;
  
  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  /* setup system matrix */
  ier = CVDlsSetupMatrix(cv_mem, vtemp1, vtemp2, vtemp3);
  if (ier < 0)  return(-1);
  if (ier > 0)  return(1);

  /* Do LU factorization of M */
  ier = BandGBTRF(cvdls_mem->d_M, cvdls_mem->d_lpivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    cvdls_mem->d_last_flag = (long int) ier;
    return(1);
  }
  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int cvBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  CVDlsMem cvdls_mem;
  realtype *bd;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(cvdls_mem->d_M, cvdls_mem->d_lpivots, bd);

  /* If CV_BDF, scale the correction to account for change in gamma */
  if ((cv_mem->cv_lmm == CV_BDF) && (cv_mem->cv_gamrat != ONE)) {
    N_VScale(TWO/(ONE + cv_mem->cv_gamrat), b, b);
  }

  cvdls_mem->d_last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static int cvBandFree(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) cv_mem->cv_lmem;

  DestroyMat(cvdls_mem->d_M);
  DestroyMat(cvdls_mem->d_savedJ);
  DestroyArray(cvdls_mem->d_lpivots);
  free(cvdls_mem);
  cv_mem->cv_lmem = NULL;

  return(0);
}

