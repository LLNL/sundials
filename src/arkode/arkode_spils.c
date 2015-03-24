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
 * This is the implementation file for the ARKSPILS linear solvers.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_spils_impl.h"

/* constants */
#define MAX_DQITERS  3  /* max. # of attempts to recover in DQ J*v */
#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define ONE    RCONST(1.0)


/*===============================================================
   OPTIONAL INPUT and OUTPUT
===============================================================*/

/*---------------------------------------------------------------
 ARKSpilsSetPrecType
---------------------------------------------------------------*/
int ARKSpilsSetPrecType(void *arkode_mem, int pretype)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetPrecType", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetPrecType", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetPrecType", MSGS_BAD_PRETYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  arkspils_mem->s_pretype = pretype;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassPrecType
---------------------------------------------------------------*/
int ARKSpilsSetMassPrecType(void *arkode_mem, int pretype)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassPrecType", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassPrecType", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE)  && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassPrecType", MSGS_BAD_PRETYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  arkspils_mem->s_pretype = pretype;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetGSType
---------------------------------------------------------------*/
int ARKSpilsSetGSType(void *arkode_mem, int gstype)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetGSType", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetGSType", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if ((arkspils_mem->s_type != SPILS_SPGMR) &&
      (arkspils_mem->s_type != SPILS_SPFGMR)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetGSType", MSGS_BAD_LSTYPE);
    fprintf(stderr,"solver type = %i\n",arkspils_mem->s_type);
    return(ARKSPILS_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetGSType", MSGS_BAD_GSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }
  arkspils_mem->s_gstype = gstype;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassGSType
---------------------------------------------------------------*/
int ARKSpilsSetMassGSType(void *arkode_mem, int gstype)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassGSType", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassGSType", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  if ((arkspils_mem->s_type != SPILS_SPGMR) &&
      (arkspils_mem->s_type != SPILS_SPFGMR)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassGSType", MSGS_BAD_LSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassGSType", MSGS_BAD_GSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }
  arkspils_mem->s_gstype = gstype;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 Function : ARKSpilsSetMaxl
---------------------------------------------------------------*/
int ARKSpilsSetMaxl(void *arkode_mem, int maxl)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMaxl", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(NULL, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMaxl", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if ((arkspils_mem->s_type == SPILS_SPGMR) &&
      (arkspils_mem->s_type == SPILS_SPFGMR)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMaxl", MSGS_BAD_LSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  mxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;
  arkspils_mem->s_maxl = mxl;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 Function : ARKSpilsSetMassMaxl
---------------------------------------------------------------*/
int ARKSpilsSetMassMaxl(void *arkode_mem, int maxl)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassMaxl", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassMaxl", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  if ((arkspils_mem->s_type == SPILS_SPGMR) &&
      (arkspils_mem->s_type == SPILS_SPFGMR)) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassMaxl", MSGS_BAD_LSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  mxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;
  arkspils_mem->s_maxl = mxl;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetEpsLin
---------------------------------------------------------------*/
int ARKSpilsSetEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetEpsLin", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetEpsLin", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetEpsLin", MSGS_BAD_EPLIN);
    return(ARKSPILS_ILL_INPUT);
  }
  arkspils_mem->s_eplifac = (eplifac == ZERO) ? ARKSPILS_EPLIN : eplifac;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassEpsLin
---------------------------------------------------------------*/
int ARKSpilsSetMassEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassEpsLin", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassEpsLin", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassEpsLin", MSGS_BAD_EPLIN);
    return(ARKSPILS_ILL_INPUT);
  }
  arkspils_mem->s_eplifac = (eplifac == ZERO) ? ARKSPILS_EPLIN : eplifac;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetPreconditioner
---------------------------------------------------------------*/
int ARKSpilsSetPreconditioner(void *arkode_mem, 
			      ARKSpilsPrecSetupFn pset, 
			      ARKSpilsPrecSolveFn psolve)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetPreconditioner", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  arkspils_mem->s_pset   = pset;
  arkspils_mem->s_psolve = psolve;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassPreconditioner
---------------------------------------------------------------*/
int ARKSpilsSetMassPreconditioner(void *arkode_mem, 
				  ARKSpilsMassPrecSetupFn pset, 
				  ARKSpilsMassPrecSolveFn psolve)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassPreconditioner", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassPreconditioner", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  arkspils_mem->s_pset   = pset;
  arkspils_mem->s_psolve = psolve;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetJacTimesVecFn
---------------------------------------------------------------*/
int ARKSpilsSetJacTimesVecFn(void *arkode_mem, 
			     ARKSpilsJacTimesVecFn jtv)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetJacTimesVecFn", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (jtv != NULL) {
    arkspils_mem->s_jtimesDQ = FALSE;
    arkspils_mem->s_jtimes   = jtv;
  } else {
    arkspils_mem->s_jtimesDQ = TRUE;
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassTimesVecFn
---------------------------------------------------------------*/
int ARKSpilsSetMassTimesVecFn(void *arkode_mem, 
			      ARKSpilsMassTimesVecFn mtv,
			      void *mtimes_data)
{
  ARKodeMem ark_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetJacTimesVecFn", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (mtv == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassTimesVecFn", "non-NULL function must be supplied");
    return(ARKSPILS_LMEM_NULL);
  }

  /* set arguments into ark_mem data structure */
  ark_mem->ark_mtimes = mtv;
  ark_mem->ark_mtimes_data = mtimes_data;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetWorkSpace
---------------------------------------------------------------*/
int ARKSpilsGetWorkSpace(void *arkode_mem, long int *lenrwLS, 
			 long int *leniwLS)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  int maxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetWorkSpace", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  switch(arkspils_mem->s_type) {
  case SPILS_SPGMR:
    maxl = arkspils_mem->s_maxl;
    *lenrwLS = ark_mem->ark_lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;
    *leniwLS = ark_mem->ark_liw1*(maxl + 5);
    break;
  case SPILS_SPBCG:
    *lenrwLS = ark_mem->ark_lrw1 * 9;
    *leniwLS = ark_mem->ark_liw1 * 9;
    break;
  case SPILS_SPTFQMR:
    *lenrwLS = ark_mem->ark_lrw1*11;
    *leniwLS = ark_mem->ark_liw1*11;
    break;
  case SPILS_PCG:
    *lenrwLS = ark_mem->ark_lrw1 * 4;
    *leniwLS = ark_mem->ark_liw1 * 4 + 1;
    break;
  case SPILS_SPFGMR:
    maxl = arkspils_mem->s_maxl;
    *lenrwLS = ark_mem->ark_lrw1*(2*maxl + 4) + maxl*(maxl + 4) + 1;
    *leniwLS = ark_mem->ark_liw1*(2*maxl + 4);
    break;
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetMassWorkSpace
---------------------------------------------------------------*/
int ARKSpilsGetMassWorkSpace(void *arkode_mem, long int *lenrwMLS, 
			     long int *leniwMLS)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;
  int maxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetMassWorkSpace", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetMassWorkSpace", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  switch(arkspils_mem->s_type) {
  case SPILS_SPGMR:
    maxl = arkspils_mem->s_maxl;
    *lenrwMLS = ark_mem->ark_lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;
    *leniwMLS = ark_mem->ark_liw1*(maxl + 5);
    break;
  case SPILS_SPBCG:
    *lenrwMLS = ark_mem->ark_lrw1 * 9;
    *leniwMLS = ark_mem->ark_liw1 * 9;
    break;
  case SPILS_SPTFQMR:
    *lenrwMLS = ark_mem->ark_lrw1*11;
    *leniwMLS = ark_mem->ark_liw1*11;
    break;
  case SPILS_PCG:
    *lenrwMLS = ark_mem->ark_lrw1 * 4;
    *leniwMLS = ark_mem->ark_liw1 * 4 + 1;
    break;
  case SPILS_SPFGMR:
    maxl = arkspils_mem->s_maxl;
    *lenrwMLS = ark_mem->ark_lrw1*(2*maxl + 4) + maxl*(maxl + 4) + 1;
    *leniwMLS = ark_mem->ark_liw1*(2*maxl + 4);
    break;
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumPrecEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumPrecEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *npevals = arkspils_mem->s_npe;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassPrecEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumMassPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecEvals", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *npevals = arkspils_mem->s_npe;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumPrecSolves
---------------------------------------------------------------*/
int ARKSpilsGetNumPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumPrecSolves", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *npsolves = arkspils_mem->s_nps;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassPrecSolves
---------------------------------------------------------------*/
int ARKSpilsGetNumMassPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecSolves", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecSolves", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *npsolves = arkspils_mem->s_nps;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumLinIters
---------------------------------------------------------------*/
int ARKSpilsGetNumLinIters(void *arkode_mem, long int *nliters)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumLinIters", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *nliters = arkspils_mem->s_nli;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassIters
---------------------------------------------------------------*/
int ARKSpilsGetNumMassIters(void *arkode_mem, long int *nmiters)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassLinIters", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassLinIters", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *nmiters = arkspils_mem->s_nli;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumConvFails
---------------------------------------------------------------*/
int ARKSpilsGetNumConvFails(void *arkode_mem, long int *nlcfails)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumConvFails", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *nlcfails = arkspils_mem->s_ncfl;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassConvFails
---------------------------------------------------------------*/
int ARKSpilsGetNumMassConvFails(void *arkode_mem, long int *nmcfails)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassConvFails", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassConvFails", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *nmcfails = arkspils_mem->s_ncfl;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumJtimesEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumJtimesEvals(void *arkode_mem, long int *njvevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumJtimesEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *njvevals = arkspils_mem->s_njtimes;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumRhsEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumRhsEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumRhsEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *nfevalsLS = arkspils_mem->s_nfes;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetLastFlag
---------------------------------------------------------------*/
int ARKSpilsGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetLastFlag", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *flag = arkspils_mem->s_last_flag;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetLastMassFlag
---------------------------------------------------------------*/
int ARKSpilsGetLastMassFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetLastMassFlag", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetLastMassFlag", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *flag = arkspils_mem->s_last_flag;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetReturnFlagName
---------------------------------------------------------------*/
char *ARKSpilsGetReturnFlagName(long int flag)
{
  char *name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKSPILS_SUCCESS:
    sprintf(name,"ARKSPILS_SUCCESS");
    break; 
  case ARKSPILS_MEM_NULL:
    sprintf(name,"ARKSPILS_MEM_NULL");
    break;
  case ARKSPILS_LMEM_NULL:
    sprintf(name,"ARKSPILS_LMEM_NULL");
    break;
  case ARKSPILS_MASSMEM_NULL:
    sprintf(name,"ARKSPILS_MASSMEM_NULL");
    break;
  case ARKSPILS_ILL_INPUT:
    sprintf(name,"ARKSPILS_ILL_INPUT");
    break;
  case ARKSPILS_MEM_FAIL:
    sprintf(name,"ARKSPILS_MEM_FAIL");
    break;
  case ARKSPILS_PMEM_NULL:
    sprintf(name,"ARKSPILS_PMEM_NULL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*===============================================================
   ARKSPILS private functions
===============================================================*/

/*---------------------------------------------------------------
 ARKSpilsAtimes:

 This routine generates the matrix-vector product z = Av, where
 A = M - gamma*J. The product M*v is obtained either by calling 
 the mtimes routine or by just using v (if M=I).  The product 
 J*v is obtained by calling the jtimes routine. It is then scaled 
 by -gamma and added to M*v to obtain A*v. The return value is 
 the same as the values returned by jtimes and mtimes -- 
 0 if successful, nonzero otherwise.
---------------------------------------------------------------*/
int ARKSpilsAtimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;
  int jtflag, mtflag;

  ark_mem = (ARKodeMem) arkode_mem;
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  jtflag = arkspils_mem->s_jtimes(v, z, ark_mem->ark_tn, 
				  arkspils_mem->s_ycur, 
				  arkspils_mem->s_fcur, 
				  arkspils_mem->s_j_data, 
				  arkspils_mem->s_ytemp);
  arkspils_mem->s_njtimes++;
  if (jtflag != 0) return(jtflag);

  /* Compute mass matrix vector product and add to result */
  if (ark_mem->ark_mass_matrix) {
    mtflag = ARKSpilsMtimes(arkode_mem, v, arkspils_mem->s_ytemp);
    if (mtflag != 0) return(mtflag);
    N_VLinearSum(ONE, arkspils_mem->s_ytemp, -ark_mem->ark_gamma, z, z);
  } else {
    N_VLinearSum(ONE, v, -ark_mem->ark_gamma, z, z);
  }

  return(0);
}


/*---------------------------------------------------------------
 ARKSpilsPSolve:

 This routine interfaces between the generic Sp***Solve routine
 (within the SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG solver) and the 
 user's psolve routine.  It passes to psolve all required state 
 information from arkode_mem.  Its return value is the same as 
 that returned by psolve. Note that the generic SP*** solver 
 guarantees that ARKSpilsPSolve will not be called in the case 
 in which preconditioning is not done. This is the only case in 
 which the user's psolve routine is allowed to be NULL.
---------------------------------------------------------------*/
int ARKSpilsPSolve(void *arkode_mem, N_Vector r, N_Vector z, int lr)
{
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;
  int retval;

  ark_mem = (ARKodeMem) arkode_mem;
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* This call is counted in nps within the ARKSp***Solve routine */
  retval = arkspils_mem->s_psolve(ark_mem->ark_tn, 
				  arkspils_mem->s_ycur, 
				  arkspils_mem->s_fcur, r, z, 
				  ark_mem->ark_gamma, 
				  arkspils_mem->s_delta, lr, 
				  arkspils_mem->s_P_data, 
				  arkspils_mem->s_ytemp);

  return(retval);     
}


/*---------------------------------------------------------------
 ARKSpilsMtimes:

 This routine generates the matrix-vector product z = Mv, where
 M is the system mass matrix, by calling the user-supplied mtimes
 routine.. The return value is the same as the value returned 
 by mtimes -- 0 if successful, nonzero otherwise.
---------------------------------------------------------------*/
int ARKSpilsMtimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem ark_mem;
  int retval;
  ark_mem = (ARKodeMem) arkode_mem;
  retval = ark_mem->ark_mtimes(v, z, ark_mem->ark_tn, 
			       ark_mem->ark_mtimes_data);
  ark_mem->ark_mass_mult++;

  return(retval);
}


/*---------------------------------------------------------------
 ARKSpilsMPSolve:

 This routine interfaces between the generic Sp***Solve routine
 (within the SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG solver) and the 
 user's mass matrix psolve routine.  It passes to psolve all 
 required state information from arkode_mem.  Its return value is
 the same as that returned by psolve. Note that the generic SP*** 
 solver guarantees that ARKSpilsMPSolve will not be called in the 
 case in which preconditioning is not done. This is the only case 
 in which the user's psolve routine is allowed to be NULL.
---------------------------------------------------------------*/
int ARKSpilsMPSolve(void *arkode_mem, N_Vector r, N_Vector z, int lr)
{
  ARKodeMem       ark_mem;
  ARKSpilsMassMem arkspils_mem;
  int retval;

  ark_mem = (ARKodeMem) arkode_mem;
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* This call is counted in nps within the ARKSp***Solve routine */
  retval = arkspils_mem->s_psolve(ark_mem->ark_tn, r, z, 
				  arkspils_mem->s_delta, lr, 
				  arkspils_mem->s_P_data, 
				  arkspils_mem->s_ytemp);
  return(retval);     
}


/*---------------------------------------------------------------
 ARKSpilsDQJtimes:

 This routine generates a difference quotient approximation to
 the Jacobian times vector f_y(t,y) * v. The approximation is 
 Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) 
 is input, i.e. the WRMS norm of v/vnrm is 1.
---------------------------------------------------------------*/
int ARKSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
		     N_Vector y, N_Vector fy,
		     void *data, N_Vector work)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  realtype sig, siginv;
  int iter, retval;

  /* data is arkode_mem */
  ark_mem = (ARKodeMem) data;
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, ark_mem->ark_ewt);

  for (iter=0; iter<MAX_DQITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = ark_mem->ark_fi(t, work, Jv, ark_mem->ark_user_data); 
    arkspils_mem->s_nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    /* If fi failed recoverably, shrink sig and retry */
    sig *= PT25;

  }

  /* If retval still isn't 0, return with a recoverable failure */
  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
