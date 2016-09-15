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
 * This is the implementation file for the ARKBAND linear solver.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_band.h>
#include "arkode_direct_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_math.h>

/* Constants */
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* ARKBAND linit, lsetup, lsolve, and lfree routines */
static int arkBandInit(ARKodeMem ark_mem);
static int arkBandSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);
static int arkBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);
static int arkBandFree(ARKodeMem ark_mem);

/* ARKBAND minit, msetup, msolve, mfree and mtimes routines */
static int arkMassBandInit(ARKodeMem ark_mem);
static int arkMassBandSetup(ARKodeMem ark_mem, N_Vector vtemp1, 
			    N_Vector vtemp2, N_Vector vtemp3);
static int arkMassBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight);
static int arkMassBandFree(ARKodeMem ark_mem);
static int arkMassBandMultiply(N_Vector v, N_Vector Mv, 
			       realtype t, void *user_data);


/*---------------------------------------------------------------
 ARKBand:

 This routine initializes the memory record and sets various function
 fields specific to the band linear solver module.  ARKBand first calls
 the existing lfree routine if this is not NULL.  It then sets the
 ark_linit, ark_lsetup, ark_lsolve, and ark_lfree fields in (*arkode_mem)
 to be arkBandInit, arkBandSetup, arkBandSolve, and arkBandFree,
 respectively.  It allocates memory for a structure of type
 ARKDlsMemRec and sets the ark_lmem field in (*arkode_mem) to the
 address of this structure.  It sets setupNonNull in (*arkode_mem) to be
 TRUE, d_mu to be mupper, d_ml to be mlower, and the d_jac field to be 
 arkDlsBandDQJac.  Finally, it allocates memory for M, savedJ, and pivot.  
 The ARKBand return value is SUCCESS = 0, LMEM_FAIL = -1, or 
 LIN_ILL_INPUT = -2.

 NOTE: The band linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKBand will first 
       test for compatible a compatible N_Vector internal
       representation by checking that the function 
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKBand(void *arkode_mem, long int N, long int mupper, long int mlower)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKBAND", "ARKBand", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */  
  ark_mem->ark_linit  = arkBandInit;
  ark_mem->ark_lsetup = arkBandSetup;
  ark_mem->ark_lsolve = arkBandSolve;
  ark_mem->ark_lfree  = arkBandFree;
  ark_mem->ark_lsolve_type = 2;
  
  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_BAND;
  
  /* Initialize Jacobian-related data */
  arkdls_mem->d_jacDQ = TRUE;
  arkdls_mem->d_bjac = NULL;
  arkdls_mem->d_J_data = NULL;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  ark_mem->ark_setupNonNull = TRUE;

  /* Initialize counters */
  arkDlsInitializeCounters(arkdls_mem);
  
  /* Load problem dimension */
  arkdls_mem->d_n = N;

  /* Load half-bandwiths in arkdls_mem */
  arkdls_mem->d_ml = mlower;
  arkdls_mem->d_mu = mupper;

  /* Test ml and mu for legality */
  if ((arkdls_mem->d_ml < 0) || (arkdls_mem->d_mu < 0) || (arkdls_mem->d_ml >= N) || (arkdls_mem->d_mu >= N)) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  arkdls_mem->d_smu = SUNMIN(N-1, arkdls_mem->d_mu + arkdls_mem->d_ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_M = NewBandMat(N, arkdls_mem->d_mu, arkdls_mem->d_ml, arkdls_mem->d_smu);
  if (arkdls_mem->d_M == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_savedJ = NULL;
  arkdls_mem->d_savedJ = NewBandMat(N, arkdls_mem->d_mu, arkdls_mem->d_ml, arkdls_mem->d_mu);
  if (arkdls_mem->d_savedJ == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_lpivots = NULL;
  arkdls_mem->d_lpivots = NewLintArray(N);
  if (arkdls_mem->d_lpivots == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    DestroyMat(arkdls_mem->d_savedJ);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkBandInit:

 This routine does remaining initializations specific to the band
 linear solver.
---------------------------------------------------------------*/
static int arkBandInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  arkDlsInitializeCounters(arkdls_mem);

  /* Set Jacobian function and data, depending on jacDQ */
  if (arkdls_mem->d_jacDQ) {
    arkdls_mem->d_bjac = arkDlsBandDQJac;
    arkdls_mem->d_J_data = ark_mem;
  } else {
    arkdls_mem->d_J_data = ark_mem->ark_user_data;
  }

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkBandSetup:

 This routine does the setup operations for the band linear 
 solver. It makes a decision whether or not to call the Jacobian 
 evaluation routine based on various state variables, and if not 
 it uses the saved copy.  In any case, it constructs the Newton 
 matrix  A = M - gamma*J, updates counters, and calls the band 
 LU factorization routine.
---------------------------------------------------------------*/
static int arkBandSetup(ARKodeMem ark_mem, int convfail, 
			N_Vector ypred, N_Vector fpred, 
			booleantype *jcurPtr, N_Vector vtemp1,
			N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int i, j, ier, ml, mu, N, M, is, ie;
  DlsMat A, Mass;
  ARKDlsMem arkdls_mem;
  ARKDlsMassMem arkdls_mass_mem;
  int retval;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkdls_mem->d_nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    BandCopy(arkdls_mem->d_savedJ, arkdls_mem->d_M, arkdls_mem->d_mu, arkdls_mem->d_ml);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    arkdls_mem->d_nje++;
    arkdls_mem->d_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SetToZero(arkdls_mem->d_M); 

    retval = arkdls_mem->d_bjac(arkdls_mem->d_n, arkdls_mem->d_mu, arkdls_mem->d_ml, 
				ark_mem->ark_tn, ypred, fpred, arkdls_mem->d_M, 
				arkdls_mem->d_J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKBAND", "arkBandSetup", MSGD_JACFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }

    BandCopy(arkdls_mem->d_M, arkdls_mem->d_savedJ, arkdls_mem->d_mu, arkdls_mem->d_ml);

  }

  /* Scale J by -gamma */
  BandScale(-ark_mem->ark_gamma, arkdls_mem->d_M);

  /* Add mass matrix to get A = M-gamma*J*/
  if (ark_mem->ark_mass_matrix) {

    /* Compute mass matrix */
    arkdls_mass_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
    SetToZero(arkdls_mass_mem->d_M);
    retval = arkdls_mass_mem->d_bmass(arkdls_mass_mem->d_n, arkdls_mass_mem->d_mu, 
				      arkdls_mass_mem->d_ml, ark_mem->ark_tn, 
				      arkdls_mass_mem->d_M, arkdls_mass_mem->d_M_data, 
				      vtemp1, vtemp2, vtemp3);
    arkdls_mass_mem->d_nme++;
    if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKBAND", 
		      "arkBandSetup",  MSGD_MASSFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_RECVR;
      return(1);
    }

    /* perform matrix sum */
    ml = arkdls_mem->d_M->ml;
    mu = arkdls_mem->d_M->mu;
    N = arkdls_mem->d_M->N;
    M = arkdls_mem->d_M->M;
    A = arkdls_mem->d_M;
    Mass = arkdls_mass_mem->d_M;
    for (j=0; j<N; j++) {                /* loop over columns */
      is = (0 > j-mu) ? 0 : j-mu;        /* colum nonzero bounds */
      ie = (M-1 < j+ml) ? M-1 : j+ml;
      for (i=is; i<=ie; i++) {           /* loop over rows */
	BAND_ELEM(A,i,j) += BAND_ELEM(Mass,i,j);
      }
    }

  } else {
    AddIdentity(arkdls_mem->d_M);
  }

  /* Do LU factorization of M */
  ier = BandGBTRF(arkdls_mem->d_M, arkdls_mem->d_lpivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    arkdls_mem->d_last_flag = ier;
    return(1);
  }
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkBandSolve:

 This routine handles the solve operation for the band linear solver
 by calling the band backsolve routine.  The return value is 0.
---------------------------------------------------------------*/
static int arkBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  ARKDlsMem arkdls_mem;
  realtype *bd;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(arkdls_mem->d_M, arkdls_mem->d_lpivots, bd);

  /* scale the correction to account for change in gamma */
  if (ark_mem->ark_gamrat != ONE) 
    N_VScale(TWO/(ONE + ark_mem->ark_gamrat), b, b);

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkBandFree:

 This routine frees memory specific to the band linear solver.
---------------------------------------------------------------*/
static int arkBandFree(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  DestroyMat(arkdls_mem->d_M);
  DestroyMat(arkdls_mem->d_savedJ);
  DestroyArray(arkdls_mem->d_lpivots);
  free(arkdls_mem);
  ark_mem->ark_lmem = NULL;

  return(0);
}




/*---------------------------------------------------------------
 ARKMassBand:

 This routine initializes the memory record and sets various 
 function fields specific to the band mass matrix linear solver 
 module.  ARKMassBand first calls the existing mfree routine if 
 this is not NULL.  It then sets the ark_minit, ark_msetup, 
 ark_msolve, and ark_mfree fields in (*arkode_mem) to be 
 arkMassBandInit, arkMassBandSetup, arkMassBandSolve, and 
 arkMassBandFree, respectively.  It allocates memory for a 
 structure of type ARKDlsMassMemRec and sets the ark_mass_mem 
 field in (*arkode_mem) to the address of this structure.  It 
 sets MassSetupNonNull in (*arkode_mem) to be TRUE, d_mu to be 
 mupper and d_ml to be mlower. Finally, it allocates memory for 
 M and pivot.  The ARKMassBand return value is SUCCESS = 0, 
 LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.

 NOTE: The band linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKMassBand will first 
       test for compatible a compatible N_Vector internal
       representation by checking that the function 
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKMassBand(void *arkode_mem, long int N, long int mupper, 
		long int mlower, ARKDlsBandMassFn bmass)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKBAND", "ARKMassBand", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKMassBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_mfree != NULL) ark_mem->ark_mfree(ark_mem);

  /* Set four main function fields in ark_mem, enable mass matrix */
  ark_mem->ark_mass_matrix = TRUE;
  ark_mem->ark_minit  = arkMassBandInit;
  ark_mem->ark_msetup = arkMassBandSetup;
  ark_mem->ark_msolve = arkMassBandSolve;
  ark_mem->ark_mfree  = arkMassBandFree;
  ark_mem->ark_mtimes = arkMassBandMultiply;
  ark_mem->ark_mtimes_data = (void *) ark_mem;
  ark_mem->ark_msolve_type = 2;
  
  /* Get memory for ARKDlsMassMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMassMem) malloc(sizeof(struct ARKDlsMassMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKMassBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_BAND;
  
  /* Initialize mass-matrix-related data */
  arkdls_mem->d_nme = 0;
  arkdls_mem->d_bmass = bmass;
  arkdls_mem->d_M_data = NULL;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  ark_mem->ark_MassSetupNonNull = TRUE;

  /* Load problem dimension */
  arkdls_mem->d_n = N;

  /* Load half-bandwiths in arkdls_mem */
  arkdls_mem->d_ml = mlower;
  arkdls_mem->d_mu = mupper;

  /* Test ml and mu for legality */
  if ((arkdls_mem->d_ml < 0) || (arkdls_mem->d_mu < 0) || (arkdls_mem->d_ml >= N) || (arkdls_mem->d_mu >= N)) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKMassBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  arkdls_mem->d_smu = SUNMIN(N-1, arkdls_mem->d_mu + arkdls_mem->d_ml);

  /* Allocate memory for M, M_lu, and pivot arrays */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_M = NewBandMat(N, arkdls_mem->d_mu, arkdls_mem->d_ml, arkdls_mem->d_mu);
  if (arkdls_mem->d_M == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKMassBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_M_lu = NULL;
  arkdls_mem->d_M_lu = NewBandMat(N, arkdls_mem->d_mu, arkdls_mem->d_ml, arkdls_mem->d_smu);
  if (arkdls_mem->d_M_lu == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKMassBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_lpivots = NULL;
  arkdls_mem->d_lpivots = NewLintArray(N);
  if (arkdls_mem->d_lpivots == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKMassBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    DestroyMat(arkdls_mem->d_M_lu);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_mass_mem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkMassBandInit:

 This routine does remaining initializations specific to the band
 mass matrix linear solver.
---------------------------------------------------------------*/
static int arkMassBandInit(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  arkdls_mem->d_nme = 0;

  /* Set mass matrix data */
  arkdls_mem->d_M_data = ark_mem->ark_user_data;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassBandSetup:

 This routine does the setup operations for the band mass matrix
 solver. It constructs the mass matrix, M, updates counters, 
 and calls the band LU factorization routine.
---------------------------------------------------------------*/
static int arkMassBandSetup(ARKodeMem ark_mem, N_Vector vtemp1, 
			    N_Vector vtemp2, N_Vector vtemp3)
{
  long int ier;
  ARKDlsMassMem arkdls_mem;
  int retval;

  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  SetToZero(arkdls_mem->d_M); 
  retval = arkdls_mem->d_bmass(arkdls_mem->d_n, arkdls_mem->d_mu, 
			       arkdls_mem->d_ml, ark_mem->ark_tn, 
			       arkdls_mem->d_M, arkdls_mem->d_M_data, 
			       vtemp1, vtemp2, vtemp3);
  arkdls_mem->d_nme++;
  if (retval < 0) {
    arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKBAND", "arkMassBandSetup", 
		    MSGD_MASSFUNC_FAILED);
    arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_RECVR;
    return(1);
  }

  /* Copy M into M_lu for LU decomposition */
  BandCopy(arkdls_mem->d_M, arkdls_mem->d_M_lu, arkdls_mem->d_mu, arkdls_mem->d_ml);

  /* Do LU factorization of M */
  ier = BandGBTRF(arkdls_mem->d_M_lu, arkdls_mem->d_lpivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    arkdls_mem->d_last_flag = ier;
    return(1);
  }
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassBandSolve:

 This routine handles the solve operation for the band mass 
 matrix solver by calling the band backsolve routine.  The 
 return value is 0.
---------------------------------------------------------------*/
static int arkMassBandSolve(ARKodeMem ark_mem, N_Vector b, 
			    N_Vector weight)
{
  ARKDlsMassMem arkdls_mem;
  realtype *bd;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  bd = N_VGetArrayPointer(b);
  BandGBTRS(arkdls_mem->d_M_lu, arkdls_mem->d_lpivots, bd);
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassBandFree:

 This routine frees memory specific to the band linear solver.
---------------------------------------------------------------*/
static int arkMassBandFree(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;

  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  DestroyMat(arkdls_mem->d_M);
  DestroyMat(arkdls_mem->d_M_lu);
  DestroyArray(arkdls_mem->d_lpivots);
  free(arkdls_mem);
  ark_mem->ark_mass_mem = NULL;

  return(0);
}


/*---------------------------------------------------------------
 arkMassBandMultiply performs a matrix-vector product, 
 multiplying the current mass matrix by a given vector.
---------------------------------------------------------------*/                  
static int arkMassBandMultiply(N_Vector v, N_Vector Mv, 
			       realtype t, void *arkode_mem)
{
  /* extract the DlsMassMem structure from the user_data pointer */
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;
  realtype *vdata=NULL, *Mvdata=NULL;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKBAND", 
		    "arkMassBandMultiply", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* access the vector arrays (since they must be serial vectors) */
  vdata = N_VGetArrayPointer(v);
  Mvdata = N_VGetArrayPointer(Mv);
  if (vdata == NULL || Mvdata == NULL)
    return(1);

  /* perform matrix-vector product and return */
  BandMatvec(arkdls_mem->d_M, vdata, Mvdata);
  return(0);
}


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
