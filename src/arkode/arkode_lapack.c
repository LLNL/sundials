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
 * This is the implementation file for a ARKODE dense linear solver
 * using BLAS and LAPACK functions.
 *--------------------------------------------------------------*/

/* NOTE: the only operation that does not use Blas/Lapack functions
   is matrix plus mass matrix (in calculating M-gamma*J in lsetup) */

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_lapack.h>
#include "arkode_direct_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_math.h>

/* Constants */
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)


/*===============================================================
 PROTOTYPES FOR PRIVATE FUNCTIONS
===============================================================*/

/* ARKLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int arkLapackDenseInit(ARKodeMem ark_mem);
static int arkLapackDenseSetup(ARKodeMem ark_mem, int convfail, 
			       N_Vector yP, N_Vector fctP, 
			       booleantype *jcurPtr, N_Vector tmp1, 
			       N_Vector tmp2, N_Vector tmp3);
static int arkLapackDenseSolve(ARKodeMem ark_mem, N_Vector b, 
			       N_Vector weight, N_Vector yC, 
			       N_Vector fctC);
static void arkLapackDenseFree(ARKodeMem ark_mem);

/* ARKLAPACK DENSE minit, msetup, msolve, and mfree routines */ 
static int arkMassLapackDenseInit(ARKodeMem ark_mem);
static int arkMassLapackDenseSetup(ARKodeMem ark_mem, N_Vector tmp1, 
				   N_Vector tmp2, N_Vector tmp3);
static int arkMassLapackDenseSolve(ARKodeMem ark_mem, N_Vector b, 
				   N_Vector weight);
static void arkMassLapackDenseFree(ARKodeMem ark_mem);
static int arkMassLapackDenseMultiply(N_Vector v, N_Vector Mv, 
				      realtype t, void *user_data);

/* ARKLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int arkLapackBandInit(ARKodeMem ark_mem);
static int arkLapackBandSetup(ARKodeMem ark_mem, int convfail, 
			      N_Vector yP, N_Vector fctP, 
			      booleantype *jcurPtr, N_Vector tmp1, 
			      N_Vector tmp2, N_Vector tmp3);
static int arkLapackBandSolve(ARKodeMem ark_mem, N_Vector b, 
			      N_Vector weight, N_Vector yC, 
			      N_Vector fctC);
static void arkLapackBandFree(ARKodeMem ark_mem);

/* ARKLAPACK BAND minit, msetup, msolve, and mfree routines */ 
static int arkMassLapackBandInit(ARKodeMem ark_mem);
static int arkMassLapackBandSetup(ARKodeMem ark_mem, N_Vector tmp1, 
				  N_Vector tmp2, N_Vector tmp3);
static int arkMassLapackBandSolve(ARKodeMem ark_mem, N_Vector b, 
				  N_Vector weight);
static void arkMassLapackBandFree(ARKodeMem ark_mem);
static int arkMassLapackBandMultiply(N_Vector v, N_Vector Mv, 
				     realtype t, void *user_data);

/*===============================================================
 EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
===============================================================*/
              
/*---------------------------------------------------------------
 ARKLapackDense:

 This routine initializes the memory record and sets various 
 function fields specific to the linear solver module.  
 ARKLapackDense first calls the existing lfree routine if this is 
 not NULL.  Then it sets the ark_linit, ark_lsetup, ark_lsolve, 
 ark_lfree fields in (*arkode_mem) to be arkLapackDenseInit, 
 arkLapackDenseSetup, arkLapackDenseSolve, and arkLapackDenseFree, 
 respectively.  It allocates memory for a structure of type 
 ARKDlsMemRec and sets the ark_lmem field in (*arkode_mem) to the
 address of this structure.  It sets setupNonNull in (*arkode_mem) 
 to TRUE, and the d_jac field to the default arkDlsDenseDQJac. 
 Finally, it allocates memory for M, pivots, and savedJ.
 The return value is SUCCESS = 0, or LMEM_FAIL = -1.

 NOTE: The dense linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKLapackDense will 
       first test for a compatible N_Vector internal 
       representation by checking that N_VGetArrayPointer and 
       N_VSetArrayPointer exist.  Of course, other vector 
       implementations may also have these functions set, so 
       this test is not sufficient to guarantee use of the 
       serial NVECTOR package.
---------------------------------------------------------------*/
int ARKLapackDense(void *arkode_mem, int N)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", 
		    "ARKLapackDense", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL ||
      ark_mem->ark_tempv->ops->nvsetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", 
		    "ARKLapackDense", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree !=NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  ark_mem->ark_linit  = arkLapackDenseInit;
  ark_mem->ark_lsetup = arkLapackDenseSetup;
  ark_mem->ark_lsolve = arkLapackDenseSolve;
  ark_mem->ark_lfree  = arkLapackDenseFree;
  ark_mem->ark_lsolve_type = 1;

  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackDense", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  arkdls_mem->d_jacDQ  = TRUE;
  arkdls_mem->d_djac   = NULL;
  arkdls_mem->d_J_data = NULL;

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  ark_mem->ark_setupNonNull = TRUE;

  /* Set problem dimension */
  arkdls_mem->d_n = (long int) N;

  /* Allocate memory for M, pivot array, and savedJ */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_pivots = NULL;
  arkdls_mem->d_savedJ = NULL;

  arkdls_mem->d_M = NewDenseMat(arkdls_mem->d_n, arkdls_mem->d_n);
  if (arkdls_mem->d_M == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackDense", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_pivots = NewIntArray(N);
  if (arkdls_mem->d_pivots == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackDense", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_savedJ = NewDenseMat(arkdls_mem->d_n, arkdls_mem->d_n);
  if (arkdls_mem->d_savedJ == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackDense", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    DestroyArray(arkdls_mem->d_pivots);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKLapackBand:

 This routine initializes the memory record and sets various 
 function fields specific to the band linear solver module. It 
 first calls the existing lfree routine if this is not NULL.  It 
 then sets the ark_linit, ark_lsetup, ark_lsolve, and ark_lfree 
 fields in (*arkode_mem) to be arkLapackBandInit, 
 arkLapackBandSetup, arkLapackBandSolve, and arkLapackBandFree, 
 respectively.  It allocates memory for a structure of type 
 ARKLapackBandMemRec and sets the ark_lmem field in (*arkode_mem) 
 to the address of this structure.  It sets setupNonNull in 
 (*arkode_mem) to be TRUE, mu to be mupper, and ml to be mlower.  
 Finally, it allocates memory for M, pivots, and savedJ.  The 
 ARKLapackBand return value is ARKDLS_SUCCESS=0, 
 ARKDLS_MEM_FAIL=-1, or ARKDLS_ILL_INPUT=-2.

 NOTE: The ARKLAPACK linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKLapackBand will first 
       test for compatible a compatible N_Vector internal
       representation by checking that the function 
       N_VGetArrayPointer exists.  Again, this test is insufficient
       to guarantee the serial NVECTOR package, but it's a start.
---------------------------------------------------------------*/                  
int ARKLapackBand(void *arkode_mem, int N, int mupper, int mlower)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */  
  ark_mem->ark_linit  = arkLapackBandInit;
  ark_mem->ark_lsetup = arkLapackBandSetup;
  ark_mem->ark_lsolve = arkLapackBandSolve;
  ark_mem->ark_lfree  = arkLapackBandFree;
  ark_mem->ark_lsolve_type = 2;
  
  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_BAND;

  /* Initialize Jacobian-related data */
  arkdls_mem->d_jacDQ  = TRUE;
  arkdls_mem->d_bjac   = NULL;
  arkdls_mem->d_J_data = NULL;

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  ark_mem->ark_setupNonNull = TRUE;
  
  /* Load problem dimension */
  arkdls_mem->d_n = (long int) N;

  /* Load half-bandwiths in arkdls_mem */
  arkdls_mem->d_ml = (long int) mlower;
  arkdls_mem->d_mu = (long int) mupper;

  /* Test ml and mu for legality */
  if ((arkdls_mem->d_ml < 0) || (arkdls_mem->d_mu < 0) || 
      (arkdls_mem->d_ml >= arkdls_mem->d_n) || 
      (arkdls_mem->d_mu >= arkdls_mem->d_n)) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  arkdls_mem->d_smu = SUNMIN(arkdls_mem->d_n-1,
			  arkdls_mem->d_mu + arkdls_mem->d_ml);

  /* Allocate memory for M, pivot array, and savedJ */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_pivots = NULL;
  arkdls_mem->d_savedJ = NULL;

  arkdls_mem->d_M = NewBandMat(arkdls_mem->d_n, arkdls_mem->d_mu, 
			       arkdls_mem->d_ml, arkdls_mem->d_smu);
  if (arkdls_mem->d_M == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }  
  arkdls_mem->d_pivots = NewIntArray(N);
  if (arkdls_mem->d_pivots == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_savedJ = NewBandMat(arkdls_mem->d_n, arkdls_mem->d_mu, 
				    arkdls_mem->d_ml, arkdls_mem->d_smu);
  if (arkdls_mem->d_savedJ == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKLapackBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    DestroyArray(arkdls_mem->d_pivots);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKMassLapackDense:

 This routine initializes the memory record and sets various 
 function fields specific to the mass matrix solver module.  
 ARKMassLapackDense first calls the existing mfree routine if 
 this is not NULL.  Then it sets the ark_minit, ark_msetup, 
 ark_msolve, ark_mfree fields in (*arkode_mem) to be 
 arkMassLapackDenseInit, arkMassLapackDenseSetup, 
 arkMassLapackDenseSolve, and arkMassLapackDenseFree, 
 respectively.  It allocates memory for a structure of type 
 ARKDlsMassMemRec and sets the ark_mass_mem field in 
 (*arkode_mem) to the address of this structure.  It sets 
 MassSetupNonNull in (*arkode_mem) to TRUE.  Finally, it 
 allocates memory for M and pivots. The return value is 
 SUCCESS = 0, or LMEM_FAIL = -1.

 NOTE: The dense linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKMassLapackDense will 
       first test for a compatible N_Vector internal 
       representation by checking that N_VGetArrayPointer and 
       N_VSetArrayPointer exist.  Of course, other vector 
       implementations may also have these functions set, so 
       this test is not sufficient to guarantee use of the 
       serial NVECTOR package.
---------------------------------------------------------------*/
int ARKMassLapackDense(void *arkode_mem, int N, 
		       ARKDlsDenseMassFn dmass)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", 
		    "ARKMassLapackDense", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL ||
      ark_mem->ark_tempv->ops->nvsetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", 
		    "ARKMassLapackDense", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_mfree !=NULL) ark_mem->ark_mfree(ark_mem);

  /* Set related function fields in ark_mem, enable mass matrix */
  ark_mem->ark_mass_matrix = TRUE;
  ark_mem->ark_minit  = arkMassLapackDenseInit;
  ark_mem->ark_msetup = arkMassLapackDenseSetup;
  ark_mem->ark_msolve = arkMassLapackDenseSolve;
  ark_mem->ark_mfree  = arkMassLapackDenseFree;
  ark_mem->ark_mtimes = arkMassLapackDenseMultiply;
  ark_mem->ark_mtimes_data = (void *) ark_mem;
  ark_mem->ark_msolve_type = 1;

  /* Get memory for ARKDlsMassMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMassMem) malloc(sizeof(struct ARKDlsMassMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKMassLapackDense", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_DENSE;

  /* Initialize mass-matrix-related data */
  arkdls_mem->d_dmass = dmass;
  arkdls_mem->d_M_data = NULL;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  ark_mem->ark_MassSetupNonNull = TRUE;

  /* Set problem dimension */
  arkdls_mem->d_n = (long int) N;

  /* Allocate memory for M and pivot array */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_pivots = NULL;

  arkdls_mem->d_M = NewDenseMat(arkdls_mem->d_n, arkdls_mem->d_n);
  if (arkdls_mem->d_M == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKMassLapackDense", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_pivots = NewIntArray(N);
  if (arkdls_mem->d_pivots == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKMassLapackDense", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_mass_mem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKMassLapackBand:

 This routine initializes the memory record and sets various 
 function fields specific to the band mass matrix solver module. 
 It first calls the existing mfree routine if this is not NULL.  
 It then sets the ark_minit, ark_msetup, ark_msolve, and ark_mfree 
 fields in (*arkode_mem) to be arkMassLapackBandInit, 
 arkMassLapackBandSetup, arkMassLapackBandSolve, and 
 arkMassLapackBandFree, respectively.  It allocates memory for a
 structure of type ARKMassLapackBandMemRec and sets the 
 ark_mass_mem field in (*arkode_mem) to the address of this 
 structure.  It sets MassSetupNonNull in (*arkode_mem) to be TRUE, 
 mu to be mupper, and ml to be mlower.  Finally, it allocates 
 memory for M and pivots.  The ARKMassLapackBand return value is 
 ARKDLS_SUCCESS=0, ARKDLS_MEM_FAIL=-1, or ARKDLS_ILL_INPUT=-2.

 NOTE: The ARKLAPACK linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKMassLapackBand will first 
       test for compatible a compatible N_Vector internal
       representation by checking that the function 
       N_VGetArrayPointer exists.  Again, this test is insufficient
       to guarantee the serial NVECTOR package, but it's a start.
---------------------------------------------------------------*/                  
int ARKMassLapackBand(void *arkode_mem, int N, int mupper, 
		      int mlower, ARKDlsBandMassFn bmass)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", 
		    "ARKMassLapackBand", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", 
		    "ARKMassLapackBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_mfree != NULL) ark_mem->ark_mfree(ark_mem);

  /* Set four main function fields in ark_mem, enable mass matrix */
  ark_mem->ark_mass_matrix = TRUE;
  ark_mem->ark_minit  = arkMassLapackBandInit;
  ark_mem->ark_msetup = arkMassLapackBandSetup;
  ark_mem->ark_msolve = arkMassLapackBandSolve;
  ark_mem->ark_mfree  = arkMassLapackBandFree;
  ark_mem->ark_mtimes = arkMassLapackBandMultiply;
  ark_mem->ark_mtimes_data = (void *) ark_mem;
  ark_mem->ark_msolve_type = 2;
  
  /* Get memory for ARKDlsMassMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMassMem) malloc(sizeof(struct ARKDlsMassMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKMassLapackBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_BAND;

  /* Initialize mass-matrix-related data */
  arkdls_mem->d_bmass = bmass;
  arkdls_mem->d_M_data = NULL;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  ark_mem->ark_MassSetupNonNull = TRUE;
  
  /* Load problem dimension */
  arkdls_mem->d_n = (long int) N;

  /* Load half-bandwiths in arkdls_mem */
  arkdls_mem->d_ml = (long int) mlower;
  arkdls_mem->d_mu = (long int) mupper;

  /* Test ml and mu for legality */
  if ((arkdls_mem->d_ml < 0) || (arkdls_mem->d_mu < 0) || 
      (arkdls_mem->d_ml >= arkdls_mem->d_n) || 
      (arkdls_mem->d_mu >= arkdls_mem->d_n)) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", 
		    "ARKMassLapackBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  arkdls_mem->d_smu = SUNMIN(arkdls_mem->d_n-1,
			  arkdls_mem->d_mu + arkdls_mem->d_ml);

  /* Allocate memory for M and pivot array */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_pivots = NULL;

  arkdls_mem->d_M = NewBandMat(arkdls_mem->d_n, arkdls_mem->d_mu, 
			       arkdls_mem->d_ml, arkdls_mem->d_smu);
  if (arkdls_mem->d_M == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKMassLapackBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }  
  arkdls_mem->d_pivots = NewIntArray(N);
  if (arkdls_mem->d_pivots == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", 
		    "ARKMassLapackBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_mass_mem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*===============================================================
  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
===============================================================*/

/*---------------------------------------------------------------
 arkLapackDenseInit does remaining initializations specific to 
 the dense linear solver.
---------------------------------------------------------------*/                  
static int arkLapackDenseInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  
  arkdls_mem->d_nje   = 0;
  arkdls_mem->d_nfeDQ = 0;
  arkdls_mem->d_nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (arkdls_mem->d_jacDQ) {
    arkdls_mem->d_djac = arkDlsDenseDQJac;
    arkdls_mem->d_J_data = ark_mem;
  } else {
    arkdls_mem->d_J_data = ark_mem->ark_user_data;
  }

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkLapackDenseSetup does the setup operations for the dense 
 linear solver. It makes a decision whether or not to call the 
 Jacobian evaluation routine based on various state variables, 
 and if not it uses the saved copy. In any case, it constructs 
 the Newton matrix A = M - gamma*J updates counters, and calls 
 the dense LU factorization routine.
---------------------------------------------------------------*/                  
static int arkLapackDenseSetup(ARKodeMem ark_mem, int convfail,
			       N_Vector yP, N_Vector fctP,
			       booleantype *jcurPtr, N_Vector tmp1, 
			       N_Vector tmp2, N_Vector tmp3)
{
  ARKDlsMem arkdls_mem;
  ARKDlsMassMem arkdls_mass_mem;
  realtype dgamma, fact, *Acol_j, *Mcol_j;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, lenmat, i, j;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  intn = (int) arkdls_mem->d_n;
  lenmat = arkdls_mem->d_M->ldata ;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkdls_mem->d_nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  /* If jok = TRUE, use saved copy of J */
  if (jok) {
    *jcurPtr = FALSE;
    dcopy_f77(&lenmat, arkdls_mem->d_savedJ->data, &one, 
	      arkdls_mem->d_M->data, &one);
    
  /* If jok = FALSE, call jac routine for new J value */
  } else {
    arkdls_mem->d_nje++;
    arkdls_mem->d_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SetToZero(arkdls_mem->d_M);

    retval = arkdls_mem->d_djac(arkdls_mem->d_n, ark_mem->ark_tn, 
				yP, fctP, arkdls_mem->d_M, 
				arkdls_mem->d_J_data, tmp1, tmp2, tmp3);

    if (retval == 0) {
      dcopy_f77(&lenmat, arkdls_mem->d_M->data, &one, 
		arkdls_mem->d_savedJ->data, &one);
    } else if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKLAPACK", 
		      "arkLapackDenseSetup", MSGD_JACFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }
  }

  /* Scale J by -gamma */
  fact = -ark_mem->ark_gamma;
  dscal_f77(&lenmat, &fact, arkdls_mem->d_M->data, &one);
  
  /* Add mass matrix to get A = M-gamma*J*/
  if (ark_mem->ark_mass_matrix) {

    /* Compute mass matrix */
    arkdls_mass_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
    SetToZero(arkdls_mass_mem->d_M);
    retval = arkdls_mass_mem->d_dmass(arkdls_mass_mem->d_n, 
				      ark_mem->ark_tn, 
				      arkdls_mass_mem->d_M, 
				      arkdls_mass_mem->d_M_data, 
				      tmp1, tmp2, tmp3);
    arkdls_mass_mem->d_nme++;
    if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKDENSE", 
		      "arkLapackDenseSetup",  MSGD_MASSFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_RECVR;
      return(1);
    }

    /* Add to A */
    for (j=0; j<arkdls_mem->d_M->N; j++) {
      Acol_j = arkdls_mem->d_M->cols[j];
      Mcol_j = arkdls_mass_mem->d_M->cols[j];
      for (i=0; i<arkdls_mem->d_M->M; i++) 
	Acol_j[i] += Mcol_j[i];
    }
  } else {
    AddIdentity(arkdls_mem->d_M);
  }

  /* Do LU factorization of M */
  dgetrf_f77(&intn, &intn, arkdls_mem->d_M->data, 
	     &intn, arkdls_mem->d_pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  arkdls_mem->d_last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}


/*---------------------------------------------------------------
 arkLapackDenseSolve handles the solve operation for the dense 
 linear solver by calling the dense backsolve routine.
---------------------------------------------------------------*/                  
static int arkLapackDenseSolve(ARKodeMem ark_mem, N_Vector b, 
			       N_Vector weight, N_Vector yC, 
			       N_Vector fctC)
{
  ARKDlsMem arkdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  intn = (int) arkdls_mem->d_n;

  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &intn, &one, arkdls_mem->d_M->data, &intn, 
	     arkdls_mem->d_pivots, bd, &intn, &ier, 1); 

  if (ier > 0) return(1);

  /* scale the correction to account for change in gamma */
  if (ark_mem->ark_gamrat != ONE) {
    fact = TWO/(ONE + ark_mem->ark_gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }
  
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkLapackDenseFree frees memory specific to the dense linear 
 solver.
---------------------------------------------------------------*/                  
static void arkLapackDenseFree(ARKodeMem ark_mem)
{
  ARKDlsMem  arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  
  DestroyMat(arkdls_mem->d_M);
  DestroyArray(arkdls_mem->d_pivots);
  DestroyMat(arkdls_mem->d_savedJ);
  free(arkdls_mem); 
  arkdls_mem = NULL;
}


/*===============================================================
  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE MASS MATRIX
===============================================================*/

/*---------------------------------------------------------------
 arkMassLapackDenseInit does remaining initializations specific to 
 the dense linear solver.
---------------------------------------------------------------*/                  
static int arkMassLapackDenseInit(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  arkdls_mem->d_nme = 0;
  arkdls_mem->d_M_data = ark_mem->ark_user_data;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassLapackDenseSetup does the setup operations for the dense 
 mass matrix solver. It calls the mass matrix evaluation routine,
 updates counters, and calls the dense LU factorization routine.
---------------------------------------------------------------*/                  
static int arkMassLapackDenseSetup(ARKodeMem ark_mem, N_Vector tmp1, 
				   N_Vector tmp2, N_Vector tmp3)
{
  ARKDlsMassMem arkdls_mem;
  int ier, retval;
  int intn;

  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  intn = (int) arkdls_mem->d_n;

  SetToZero(arkdls_mem->d_M);
  retval = arkdls_mem->d_dmass(arkdls_mem->d_n, ark_mem->ark_tn, 
			       arkdls_mem->d_M, arkdls_mem->d_M_data, 
			       tmp1, tmp2, tmp3);
  arkdls_mem->d_nme++;
  if (retval < 0) {
    arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKLAPACK", 
		    "arkMassLapackDenseSetup", MSGD_MASSFUNC_FAILED);
    arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_UNRECVR;
    return(-1);
  } else if (retval > 0) {
    arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_RECVR;
    return(1);
  }

  /* Do LU factorization of M */
  dgetrf_f77(&intn, &intn, arkdls_mem->d_M->data, 
	     &intn, arkdls_mem->d_pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  arkdls_mem->d_last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}


/*---------------------------------------------------------------
 arkMassLapackDenseSolve handles the solve operation for the
 dense mass matrix solver by calling the dense backsolve routine.
---------------------------------------------------------------*/                  
static int arkMassLapackDenseSolve(ARKodeMem ark_mem, N_Vector b, 
				   N_Vector weight)
{
  ARKDlsMassMem arkdls_mem;
  realtype *bd;
  int ier, one = 1;
  int intn;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  intn = (int) arkdls_mem->d_n;
  bd = N_VGetArrayPointer(b);
  dgetrs_f77("N", &intn, &one, arkdls_mem->d_M->data, &intn, 
	     arkdls_mem->d_pivots, bd, &intn, &ier, 1); 
  if (ier > 0) return(1);
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassLapackDenseFree frees memory specific to the dense mass
 matrix solver.
---------------------------------------------------------------*/                  
static void arkMassLapackDenseFree(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  DestroyMat(arkdls_mem->d_M);
  DestroyArray(arkdls_mem->d_pivots);
  free(arkdls_mem); 
  arkdls_mem = NULL;
}


/*---------------------------------------------------------------
 arkMassLapackDenseMultiply performs a matrix-vector product, 
 multiplying the current mass matrix by a given vector.
---------------------------------------------------------------*/                  
static int arkMassLapackDenseMultiply(N_Vector v, N_Vector Mv, 
				      realtype t, void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;
  realtype *vdata=NULL, *Mvdata=NULL;
  realtype *Mcol_j;
  int i, j;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", 
		    "arkMassLapackDenseMultiply", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }

  /* extract the DlsMassMem structure from the user_data pointer */
  ark_mem = (ARKodeMem) arkode_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* zero out the result */
  N_VConst(0.0, Mv);

  /* access the vector arrays (since they must be serial vectors) */
  vdata = N_VGetArrayPointer(v);
  Mvdata = N_VGetArrayPointer(Mv);
  if (vdata == NULL || Mvdata == NULL)
    return(1);

  /* perform matrix-vector product and return */
  for (j=0; j<arkdls_mem->d_M->N; j++) {
    Mcol_j = arkdls_mem->d_M->cols[j];
    for (i=0; i<arkdls_mem->d_M->M; i++) 
      Mvdata[i] += Mcol_j[i]*vdata[j];
  }
  return(0);
}


/*===============================================================
  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
===============================================================*/

/*---------------------------------------------------------------
 arkLapackBandInit does remaining initializations specific to 
 the band linear solver.
---------------------------------------------------------------*/                  
static int arkLapackBandInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  arkdls_mem->d_nje   = 0;
  arkdls_mem->d_nfeDQ = 0;
  arkdls_mem->d_nstlj = 0;

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
 arkLapackBandSetup does the setup operations for the band linear 
 solver. It makes a decision whether or not to call the Jacobian 
 evaluation routine based on various state variables, and if not 
 it uses the saved copy. In any case, it constructs the Newton 
 matrix A = M - gamma*J, updates counters, and calls the band LU
 factorization routine.
---------------------------------------------------------------*/                  
static int arkLapackBandSetup(ARKodeMem ark_mem, int convfail, 
			      N_Vector yP, N_Vector fctP, 
			      booleantype *jcurPtr, N_Vector tmp1, 
			      N_Vector tmp2, N_Vector tmp3)
{
  ARKDlsMem arkdls_mem;
  ARKDlsMassMem arkdls_mass_mem;
  realtype dgamma, fact, *Acol_j, *Mcol_j;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, iml, imu, lenmat, ldmat, i, j, colSize;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  intn = (int) arkdls_mem->d_n;
  iml = (int) arkdls_mem->d_ml;
  imu = (int) arkdls_mem->d_mu;
  lenmat = arkdls_mem->d_M->ldata;
  ldmat = arkdls_mem->d_M->ldim;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkdls_mem->d_nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  /* If jok = TRUE, use saved copy of J */
  if (jok) {
    *jcurPtr = FALSE;
    dcopy_f77(&lenmat, arkdls_mem->d_savedJ->data, 
	      &one, arkdls_mem->d_M->data, &one);
    
  /* If jok = FALSE, call jac routine for new J value */
  } else {
    arkdls_mem->d_nje++;
    arkdls_mem->d_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SetToZero(arkdls_mem->d_M);

    retval = arkdls_mem->d_bjac(arkdls_mem->d_n, arkdls_mem->d_mu, 
				arkdls_mem->d_ml, ark_mem->ark_tn, 
				yP, fctP, arkdls_mem->d_M, 
				arkdls_mem->d_J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&lenmat, arkdls_mem->d_M->data, &one, 
		arkdls_mem->d_savedJ->data, &one);
    } else if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKLAPACK", 
		      "arkLapackBandSetup", MSGD_JACFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }
  }
  
  /* Scale J by -gamma */
  fact = -ark_mem->ark_gamma;
  dscal_f77(&lenmat, &fact, arkdls_mem->d_M->data, &one);
  
  /* Add mass matrix to get A = M-gamma*J*/
  if (ark_mem->ark_mass_matrix) {

    /* Compute mass matrix */
    arkdls_mass_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
    SetToZero(arkdls_mass_mem->d_M);
    retval = arkdls_mass_mem->d_bmass(arkdls_mass_mem->d_n, arkdls_mass_mem->d_mu, 
				      arkdls_mass_mem->d_ml, ark_mem->ark_tn, 
				      arkdls_mass_mem->d_M, arkdls_mass_mem->d_M_data, 
				      tmp1, tmp2, tmp3);
    arkdls_mass_mem->d_nme++;
    if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKLAPACK", 
		      "arkLapackBandSetup",  MSGD_MASSFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_RECVR;
      return(1);
    }

    /* Add to A -- CURRENTLY ASSUMES THAT BOTH MATRICES HAVE 
                   THE SAME BAND STRUCTURE AND COLUMN SIZE */
    colSize = arkdls_mem->d_M->mu + arkdls_mem->d_M->ml + 1;
    for (j=0; j<arkdls_mem->d_M->M; j++) {
      Acol_j = arkdls_mem->d_M->cols[j] + arkdls_mem->d_M->s_mu - arkdls_mem->d_M->mu;
      Mcol_j = arkdls_mass_mem->d_M->cols[j] + arkdls_mass_mem->d_M->s_mu 
	     - arkdls_mass_mem->d_M->mu;
      for (i=0; i<colSize; i++) 
	Acol_j[i] += Mcol_j[i];
    }
  } else {
    AddIdentity(arkdls_mem->d_M);
  }
  
  /* Do LU factorization of M */
  dgbtrf_f77(&intn, &intn, &iml, &imu, arkdls_mem->d_M->data, 
	     &ldmat, arkdls_mem->d_pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  arkdls_mem->d_last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);

}


/*---------------------------------------------------------------
 arkLapackBandSolve handles the solve operation for the band 
 linear solver by calling the band backsolve routine.
---------------------------------------------------------------*/                  
static int arkLapackBandSolve(ARKodeMem ark_mem, N_Vector b, 
			      N_Vector weight, N_Vector yC, 
			      N_Vector fctC)
{
  ARKDlsMem arkdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  intn = (int) arkdls_mem->d_n;
  iml = (int) arkdls_mem->d_ml;
  imu = (int) arkdls_mem->d_mu;
  ldmat = arkdls_mem->d_M->ldim;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &intn, &iml, &imu, &one, arkdls_mem->d_M->data, 
	     &ldmat, arkdls_mem->d_pivots, bd, &intn, &ier, 1);
  if (ier > 0) return(1);

  /* scale the correction to account for change in gamma */
  if (ark_mem->ark_gamrat != ONE) {
    fact = TWO/(ONE + ark_mem->ark_gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkLapackBandFree frees memory specific to the band linear 
 solver.
---------------------------------------------------------------*/                  
static void arkLapackBandFree(ARKodeMem ark_mem)
{
  ARKDlsMem  arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  
  DestroyMat(arkdls_mem->d_M);
  DestroyArray(arkdls_mem->d_pivots);
  DestroyMat(arkdls_mem->d_savedJ);
  free(arkdls_mem); 
  arkdls_mem = NULL;
}


/*===============================================================
  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND MASS MATRIX
===============================================================*/

/*---------------------------------------------------------------
 arkMassLapackBandInit does remaining initializations specific to 
 the band linear solver.
---------------------------------------------------------------*/                  
static int arkMassLapackBandInit(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  arkdls_mem->d_nme = 0;
  arkdls_mem->d_M_data = ark_mem->ark_user_data;
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassLapackBandSetup does the setup operations for the band 
 mass matrix solver. It constructs the mass matrix M, updates 
 counters, and calls the band LU factorization routine.
---------------------------------------------------------------*/                  
static int arkMassLapackBandSetup(ARKodeMem ark_mem, N_Vector tmp1, 
				  N_Vector tmp2, N_Vector tmp3)
{
  ARKDlsMassMem arkdls_mem;
  int ier, retval;
  int intn, iml, imu, ldmat;

  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  intn = (int) arkdls_mem->d_n;
  iml = (int) arkdls_mem->d_ml;
  imu = (int) arkdls_mem->d_mu;
  ldmat = arkdls_mem->d_M->ldim;

  SetToZero(arkdls_mem->d_M);
  retval = arkdls_mem->d_bmass(arkdls_mem->d_n, arkdls_mem->d_mu, 
			       arkdls_mem->d_ml, ark_mem->ark_tn, 
			       arkdls_mem->d_M, arkdls_mem->d_M_data, 
			       tmp1, tmp2, tmp3);
  arkdls_mem->d_nme++;
  if (retval < 0) {
    arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKLAPACK", 
		    "arkMassLapackBandSetup", MSGD_MASSFUNC_FAILED);
    arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_UNRECVR;
    return(-1);
  } else if (retval > 0) {
    arkdls_mem->d_last_flag = ARKDLS_MASSFUNC_RECVR;
    return(1);
  }
  
  /* Do LU factorization of M */
  dgbtrf_f77(&intn, &intn, &iml, &imu, arkdls_mem->d_M->data, 
	     &ldmat, arkdls_mem->d_pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  arkdls_mem->d_last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}


/*---------------------------------------------------------------
 arkMassLapackBandSolve handles the solve operation for the band 
 mass matrix solver by calling the band backsolve routine.
---------------------------------------------------------------*/                  
static int arkMassLapackBandSolve(ARKodeMem ark_mem, N_Vector b, 
				  N_Vector weight)
{
  ARKDlsMassMem arkdls_mem;
  realtype *bd;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  intn = (int) arkdls_mem->d_n;
  iml = (int) arkdls_mem->d_ml;
  imu = (int) arkdls_mem->d_mu;
  ldmat = arkdls_mem->d_M->ldim;
  bd = N_VGetArrayPointer(b);
  dgbtrs_f77("N", &intn, &iml, &imu, &one, arkdls_mem->d_M->data, 
	     &ldmat, arkdls_mem->d_pivots, bd, &intn, &ier, 1);
  if (ier > 0) return(1);
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkMassLapackBandFree frees memory specific to the band mass
 matrix solver.
---------------------------------------------------------------*/                  
static void arkMassLapackBandFree(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  DestroyMat(arkdls_mem->d_M);
  DestroyArray(arkdls_mem->d_pivots);
  free(arkdls_mem); 
  arkdls_mem = NULL;
}


/*---------------------------------------------------------------
 arkMassLapackBandMultiply performs a matrix-vector product, 
 multiplying the current mass matrix by a given vector.
---------------------------------------------------------------*/                  
static int arkMassLapackBandMultiply(N_Vector v, N_Vector Mv, 
				     realtype t, void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;
  realtype *vdata=NULL, *Mvdata=NULL;
  realtype *Mcol_j;
  int colSize, s_mu, i, j;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", 
		    "arkMassLapackBandMultiply", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }

  /* extract the DlsMassMem structure from the user_data pointer */
  ark_mem = (ARKodeMem) arkode_mem;
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* zero out the result */
  N_VConst(0.0, Mv);

  /* access the vector arrays (since they must be serial vectors) */
  vdata = N_VGetArrayPointer(v);
  Mvdata = N_VGetArrayPointer(Mv);
  if (vdata == NULL || Mvdata == NULL)
    return(1);

  /* perform matrix-vector product and return */
  colSize = arkdls_mem->d_M->mu + arkdls_mem->d_M->ml + 1;
  s_mu = arkdls_mem->d_M->s_mu;
  for (j=0; j<arkdls_mem->d_M->M; j++) {
    Mcol_j = arkdls_mem->d_M->cols[j] + arkdls_mem->d_M->s_mu - arkdls_mem->d_M->mu;
    for (i=0; i<colSize; i++) 
      Mvdata[i+j-s_mu] += Mcol_j[i]*vdata[j];
  }
  return(0);
}


/*---------------------------------------------------------------
     EOF
---------------------------------------------------------------*/                  
