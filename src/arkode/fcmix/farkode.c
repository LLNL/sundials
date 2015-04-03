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
 * This is the implementation file for the Fortran interface to
 * the ARKODE package.  See farkode.h for usage.
 * NOTE: some routines are necessarily stored elsewhere to avoid
 * linking problems.  Therefore, see also farkpreco.c, farkpsol.c,
 * farkjtimes.c, farkadapt.c and farkexpstab.c for all the 
 * available options.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_band.h>
#include <arkode/arkode_dense.h>
#include <arkode/arkode_klu.h>
#include <arkode/arkode_superlumt.h>
#include <arkode/arkode_spgmr.h>
#include <arkode/arkode_spbcgs.h>
#include <arkode/arkode_sptfqmr.h>
#include <arkode/arkode_spfgmr.h>
#include <arkode/arkode_pcg.h>

/*=============================================================*/

/* Definitions for global variables shared between Fortran/C 
   interface routines */
void     *ARK_arkodemem;
long int *ARK_iout;
realtype *ARK_rout;
int       ARK_nrtfn;
int       ARK_ls;
int       ARK_mass_ls;

/*=============================================================*/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_IMP_FUN(realtype *T, realtype *Y, realtype *YDOT,
			   long int *IPAR, realtype *RPAR, int *IER);
  extern void FARK_EXP_FUN(realtype *T, realtype *Y, realtype *YDOT,
			   long int *IPAR, realtype *RPAR, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to initialize ARKode memory 
   structure; functions as an all-in-one interface to the C 
   routines ARKodeCreate, ARKodeSetUserData, ARKodeInit, and 
   ARKodeSStolerances (or ARKodeSVtolerances); see farkode.h 
   for further details */
void FARK_MALLOC(realtype *t0, realtype *y0, int *imex, 
		 int *iatol, realtype *rtol, realtype *atol, 
		 long int *iout, realtype *rout, 
		 long int *ipar, realtype *rpar, int *ier) {

  N_Vector Vatol;
  FARKUserData ARK_userdata;
  realtype reltol, abstol;

  *ier = 0;

  /* Check for required vector operations */
  if(F2C_ARKODE_vec->ops->nvgetarraypointer == NULL) {
    *ier = -1;
    printf("Error: getarraypointer vector operation is not implemented.\n\n");
    return;
  }
  if(F2C_ARKODE_vec->ops->nvsetarraypointer == NULL) {
    *ier = -1;
    printf("Error: setarraypointer vector operation is not implemented.\n\n");
    return;
  }
  if(F2C_ARKODE_vec->ops->nvcloneempty == NULL) {
    *ier = -1;
    printf("Error: cloneempty vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  ARK_arkodemem = NULL;
  Vatol = NULL;

  /* initialize global constants to zero */
  ARK_nrtfn = 0;
  ARK_ls = 0;
  ARK_mass_ls = 0;

  /* Create ARKODE object */
  ARK_arkodemem = ARKodeCreate();
  if (ARK_arkodemem == NULL) {
    *ier = -1;
    return;
  }

  /* Set and attach user data */
  ARK_userdata = NULL;
  ARK_userdata = (FARKUserData) malloc(sizeof *ARK_userdata);
  if (ARK_userdata == NULL) {
    *ier = -1;
    return;
  }
  ARK_userdata->rpar = rpar;
  ARK_userdata->ipar = ipar;
  *ier = ARKodeSetUserData(ARK_arkodemem, ARK_userdata);
  if(*ier != ARK_SUCCESS) {
    free(ARK_userdata); ARK_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Set data in F2C_ARKODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_ARKODE_vec);

  /* Call ARKodeInit based on imex argument */
  switch (*imex) {
  case 0:  /* purely implicit */
    *ier = ARKodeInit(ARK_arkodemem, NULL, FARKfi, 
		      *t0, F2C_ARKODE_vec);
    break;
  case 1:  /* purely explicit */
    *ier = ARKodeInit(ARK_arkodemem, FARKfe, NULL, 
		      *t0, F2C_ARKODE_vec);
    break;
  case 2:  /* imex */
    *ier = ARKodeInit(ARK_arkodemem, FARKfe, FARKfi, 
		      *t0, F2C_ARKODE_vec);
    break;
  }
    
  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* On failure, exit */
  if(*ier != ARK_SUCCESS) {
    free(ARK_userdata);
    ARK_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Set tolerances -- if <= 0, keep as defaults */
  reltol=1.e-4;
  abstol=1.e-9;
  if (*rtol > 0.0)  reltol = *rtol;
  switch (*iatol) {
  case 1:
    if (*atol > 0.0)  abstol = *atol;
    *ier = ARKodeSStolerances(ARK_arkodemem, reltol, abstol); 
    break;
  case 2:
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      free(ARK_userdata);
      ARK_userdata = NULL;
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    if (N_VMin(Vatol) <= 0.0)  N_VConst(abstol, Vatol);
    *ier = ARKodeSVtolerances(ARK_arkodemem, reltol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if(*ier != ARK_SUCCESS) {
    free(ARK_userdata); 
    ARK_userdata = NULL;
    *ier = -1;
    return;
  }

  /* store pointers to optional output arrays in global vars */
  ARK_iout = iout;
  ARK_rout = rout;

  /* Store the unit roundoff in rout for user access */
  ARK_rout[5] = UNIT_ROUNDOFF;

  return;
}

/*=============================================================*/

/* Fortran interface routine to re-initialize ARKode memory 
   structure; functions as an all-in-one interface to the C 
   routines ARKodeReInit and ARKodeSStolerances (or 
   ARKodeSVtolerances); see farkode.h for further details */
void FARK_REINIT(realtype *t0, realtype *y0, int *imex, int *iatol, 
		 realtype *rtol, realtype *atol, int *ier) {

  N_Vector Vatol;
  realtype reltol, abstol;
  *ier = 0;

  /* Initialize all pointers to NULL */
  Vatol = NULL;

  /* Set data in F2C_ARKODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_ARKODE_vec);

  /* Call ARKodeReInit based on imex argument */
  switch (*imex) {
  case 0:  /* purely implicit */
    *ier = ARKodeReInit(ARK_arkodemem, NULL, FARKfi, 
			*t0, F2C_ARKODE_vec);
    break;
  case 1:  /* purely explicit */
    *ier = ARKodeReInit(ARK_arkodemem, FARKfe, NULL, 
			*t0, F2C_ARKODE_vec);
    break;
  case 2:  /* imex */
    *ier = ARKodeReInit(ARK_arkodemem, FARKfe, FARKfi, 
			*t0, F2C_ARKODE_vec);
    break;
  }

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* On failure, exit */
  if (*ier != ARK_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Set tolerances */
  reltol=1.e-4;
  abstol=1.e-9;
  if (*rtol > 0.0)  reltol = *rtol;
  switch (*iatol) {
  case 1:
    if (*atol > 0.0)  abstol = *atol;
    *ier = ARKodeSStolerances(ARK_arkodemem, reltol, abstol); 
    break;
  case 2:
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    if (N_VMin(Vatol) <= 0.0)  N_VConst(abstol, Vatol);
    *ier = ARKodeSVtolerances(ARK_arkodemem, reltol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if (*ier != ARK_SUCCESS) {
    *ier = -1;
    return;
  }

  return;
}

/*=============================================================*/

/* Fortran interface routine to re-initialize ARKode memory 
   structure for a problem with a new size but similar time 
   scale; functions as an all-in-one interface to the C 
   routines ARKodeResize (and potentially ARKodeSVtolerances); 
   see farkode.h for further details */
void FARK_RESIZE(realtype *t0, realtype *y0, realtype *hscale,
		 int *itol, realtype *rtol, realtype *atol, int *ier) {

  *ier = 0;

  /* Set data in F2C_ARKODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_ARKODE_vec);
  
  /* Call ARKodeResize (currently does not allow Fortran 
     user-supplied vector resize function) */
  *ier = ARKodeResize(ARK_arkodemem, F2C_ARKODE_vec, *hscale,
		      *t0, NULL, NULL);

  /* Reset data pointer */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* On failure, exit */
  if (*ier != ARK_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Set tolerances, based on itol argument */
  if (*itol) {
    N_Vector Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = ARKodeSVtolerances(ARK_arkodemem, *rtol, Vatol);
    N_VDestroy(Vatol);
  }

  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeSetDefaults; see 
   farkode.h for further details */
void FARK_SETDEFAULTS(int *ier) {
  *ier = ARKodeSetDefaults(ARK_arkodemem);
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines having integer 
   arguments; see farkode.h for further details */
void FARK_SETIIN(char key_name[], long int *ival, int *ier) {
  if (!strncmp(key_name, "ORDER", 5)) 
    *ier = ARKodeSetOrder(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "DENSE_ORDER", 11)) 
    *ier = ARKodeSetDenseOrder(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "LINEAR", 6)) 
    *ier = ARKodeSetLinear(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "NONLINEAR", 9)) 
    *ier = ARKodeSetNonlinear(ARK_arkodemem);
  else if (!strncmp(key_name, "FIXEDPOINT", 10)) 
    *ier = ARKodeSetFixedPoint(ARK_arkodemem, (long int) *ival);
  else if (!strncmp(key_name, "NEWTON", 6)) 
    *ier = ARKodeSetNewton(ARK_arkodemem);
  else if (!strncmp(key_name, "EXPLICIT", 8)) 
    *ier = ARKodeSetExplicit(ARK_arkodemem);
  else if (!strncmp(key_name, "IMPLICIT", 8)) 
    *ier = ARKodeSetImplicit(ARK_arkodemem);
  else if (!strncmp(key_name, "IMEX", 4)) 
    *ier = ARKodeSetImEx(ARK_arkodemem);
  else if (!strncmp(key_name, "IRK_TABLE_NUM", 13)) 
    *ier = ARKodeSetIRKTableNum(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "ERK_TABLE_NUM", 13)) 
    *ier = ARKodeSetERKTableNum(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "ARK_TABLE_NUM", 13)) 
    *ier = ARKodeSetARKTableNum(ARK_arkodemem, (int) ival[0], (int) ival[1]);
  else if (!strncmp(key_name, "MAX_NSTEPS", 10)) 
    *ier = ARKodeSetMaxNumSteps(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "HNIL_WARNS", 10)) 
    *ier = ARKodeSetMaxHnilWarns(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "PREDICT_METHOD", 14)) 
    *ier = ARKodeSetPredictorMethod(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "MAX_ERRFAIL", 11)) 
    *ier = ARKodeSetMaxErrTestFails(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "MAX_CONVFAIL", 12)) 
    *ier = ARKodeSetMaxConvFails(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "MAX_NITERS", 10)) 
    *ier = ARKodeSetMaxNonlinIters(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "ADAPT_SMALL_NEF", 15))
    *ier = ARKodeSetSmallNumEFails(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "LSETUP_MSBP", 11))
    *ier = ARKodeSetMaxStepsBetweenLSet(ARK_arkodemem, (int) *ival);
  else {
    *ier = -99;
    printf("FARKSETIIN: Unrecognized key.\n\n");
  }
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines having real
   arguments; see farkode.h for further details */
void FARK_SETRIN(char key_name[], realtype *rval, int *ier) {
  if (!strncmp(key_name, "INIT_STEP", 9)) 
    *ier = ARKodeSetInitStep(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "MAX_STEP", 8)) 
    *ier = ARKodeSetMaxStep(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "MIN_STEP", 8)) 
    *ier = ARKodeSetMinStep(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "STOP_TIME", 9)) 
    *ier = ARKodeSetStopTime(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "NLCONV_COEF", 11)) 
    *ier = ARKodeSetNonlinConvCoef(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_CFL", 9))
    *ier = ARKodeSetCFLFraction(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_SAFETY", 12)) 
    *ier = ARKodeSetSafetyFactor(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_BIAS", 10))
    *ier = ARKodeSetErrorBias(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_GROWTH", 12)) 
    *ier = ARKodeSetMaxGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_BOUNDS", 12))
    *ier = ARKodeSetFixedStepBounds(ARK_arkodemem, rval[0], rval[1]);
  else if (!strncmp(key_name, "ADAPT_ETAMX1", 12))
    *ier = ARKodeSetMaxFirstGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_ETAMXF", 12))
    *ier = ARKodeSetMaxEFailGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_ETACF", 11))
    *ier = ARKodeSetMaxCFailGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "NONLIN_CRDOWN", 11))
    *ier = ARKodeSetNonlinCRDown(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "NONLIN_RDIV", 9))
    *ier = ARKodeSetNonlinRDiv(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "LSETUP_DGMAX", 12))
    *ier = ARKodeSetDeltaGammaMax(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "FIXED_STEP", 10))
    *ier = ARKodeSetFixedStep(ARK_arkodemem, *rval);
  else {
    *ier = -99;
    printf("FARKSETRIN: Unrecognized key: %s\n\n",key_name);
  }
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeSetAdaptivityMethod; 
   see farkode.h for further details */
void FARK_SETADAPTMETHOD(int *imethod, int *idefault, int *ipq, 
			 realtype *params, int *ier) {

  *ier = ARKodeSetAdaptivityMethod(ARK_arkodemem, *imethod, 
				   *idefault, *ipq, params);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeSetERKTable; see 
   farkode.h for further details */
void FARK_SETERKTABLE(int *s, int *q, int *p, realtype *c, realtype *A, 
		      realtype *b, realtype *b2, int *ier) {
  *ier = ARKodeSetERKTable(ARK_arkodemem, *s, *q, *p, c, A, b, b2);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeSetIRKTable; see 
   farkode.h for further details */
void FARK_SETIRKTABLE(int *s, int *q, int *p, realtype *c, realtype *A, 
		      realtype *b, realtype *b2, int *ier) {
  *ier = ARKodeSetIRKTable(ARK_arkodemem, *s, *q, *p, c, A, b, b2);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeSetARKTables; see 
   farkode.h for further details */
void FARK_SETARKTABLES(int *s, int *q, int *p, realtype *c, 
		       realtype *Ai, realtype *Ae, realtype *b, 
		       realtype *b2, int *ier) {
  *ier = ARKodeSetARKTables(ARK_arkodemem, *s, *q, 
			    *p, c, Ai, Ae, b, b2);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeSetDiagnostics; see 
   farkode.h for further details */
void FARK_SETDIAGNOSTICS(char fname[], int *flen, int *ier) {
  char *filename=NULL;
  FILE *DFID=NULL;
  int i;

  /* copy fname into array of specified length */
  filename = (char *) malloc((*flen)*sizeof(char));
  for (i=0; i<*flen; i++)  filename[i] = fname[i];

  /* open diagnostics output file */
  DFID = fopen(filename,"w");
  if (DFID == NULL) {
    *ier = 1;
    return;
  }
  *ier = ARKodeSetDiagnostics(ARK_arkodemem, DFID);
  free(filename);
  return;
}

/*=============================================================*/

/* Fortran routine to close diagnostics output file; see farkode.h 
   for further details */
void FARK_STOPDIAGNOSTICS(int *ier) {
  ARKodeMem ark_mem;
  if (ARK_arkodemem == NULL) {
    *ier = 1;
    return;
  }
  ark_mem = (ARKodeMem) ARK_arkodemem;

  if (ark_mem->ark_diagfp == NULL) {
    *ier = 1;
    return;
  }
  *ier = fclose(ark_mem->ark_diagfp);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKDense; see farkode.h for 
   further details */
void FARK_DENSE(long int *neq, int *ier) {
  *ier = ARKDense(ARK_arkodemem, *neq);
  ARK_ls = ARK_LS_DENSE;
  return;
}

/* Fortran interface to C routine ARKMassDense; see farkode.h 
   for further details */
void FARK_MASSDENSE(long int *neq, int *ier) {
  *ier = ARKMassDense(ARK_arkodemem, *neq, NULL);
  ARK_mass_ls = ARK_LS_DENSE;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKBand; see farkode.h for 
   further details */
void FARK_BAND(long int *neq, long int *mupper, 
	       long int *mlower, int *ier) {
  *ier = ARKBand(ARK_arkodemem, *neq, *mupper, *mlower);
  ARK_ls = ARK_LS_BAND;
  return;
}

/* Fortran interface to C routine ARKMassBand; see farkode.h 
   for further details */
void FARK_MASSBAND(long int *neq, long int *mupper, 
		   long int *mlower, int *ier) {
  *ier = ARKMassBand(ARK_arkodemem, *neq, *mupper, 
		     *mlower, NULL);
  ARK_mass_ls = ARK_LS_BAND;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKSpgmr and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_SPGMR(int *pretype, int *gstype, int *maxl, 
		realtype *delt, int *ier) {

  *ier = ARKSpgmr(ARK_arkodemem, *pretype, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPGMR;
  return;
}

/* Fortran interface to C routine ARKMassSpgmr and it's 
   associated "set" routines; see farkode.h for further details */
void FARK_MASSSPGMR(int *pretype, int *gstype, int *maxl, 
		    realtype *delt, int *ier) {

  *ier = ARKMassSpgmr(ARK_arkodemem, *pretype, *maxl, NULL, NULL);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPGMR;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKSpbcg and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_SPBCG(int *pretype, int *maxl, realtype *delt, int *ier) {

  *ier = ARKSpbcg(ARK_arkodemem, *pretype, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPBCG;
  return;
}

/* Fortran interface to C routine ARKMassSpbcg and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_MASSSPBCG(int *pretype, int *maxl, 
		    realtype *delt, int *ier) {

  *ier = ARKMassSpbcg(ARK_arkodemem, *pretype, *maxl, NULL, NULL);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPBCG;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKSptfqmr and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_SPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier) {

  *ier = ARKSptfqmr(ARK_arkodemem, *pretype, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPTFQMR;
  return;
}

/* Fortran interface to C routine ARKMassSptfqmr and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_MASSSPTFQMR(int *pretype, int *maxl, 
		      realtype *delt, int *ier) {

  *ier = ARKMassSptfqmr(ARK_arkodemem, *pretype, *maxl, NULL, NULL);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPTFQMR;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKSpfgmr and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_SPFGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier) {

  *ier = ARKSpfgmr(ARK_arkodemem, *pretype, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPFGMR;
  return;
}

/* Fortran interface to C routine ARKMassSpfgmr and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_MASSSPFGMR(int *pretype, int *gstype, int *maxl, 
		     realtype *delt, int *ier) {

  *ier = ARKMassSpfgmr(ARK_arkodemem, *pretype, *maxl, NULL, NULL);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPFGMR;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKPcg and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_PCG(int *pretype, int *maxl, realtype *delt, int *ier) {

  *ier = ARKPcg(ARK_arkodemem, *pretype, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_PCG;
  return;
}

/* Fortran interface to C routine ARKMassPcg and it's associated 
   "set" routines; see farkode.h for further details */
void FARK_MASSPCG(int *pretype, int *maxl, 
		  realtype *delt, int *ier) {

  *ier = ARKMassPcg(ARK_arkodemem, *pretype, *maxl, NULL, NULL);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_PCG;
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines for the ARKSpgmr solver; 
   see farkode.h for further details */
void FARK_SPGMRREINIT(int *pretype, int *gstype, 
		      realtype *delt, int *ier) {

  *ier = ARKSpilsSetPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPGMR;
  return;
}

/* Fortran interface to C "set" routines for the ARKMassSpgmr solver; 
   see farkode.h for further details */
void FARK_MASSSPGMRREINIT(int *pretype, int *gstype, 
			  realtype *delt, int *ier) {

  *ier = ARKSpilsSetMassPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPGMR;
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines for the ARKSpbcg solver; 
   see farkode.h for further details */
void FARK_SPBCGREINIT(int *pretype, int *maxl, 
		      realtype *delt, int *ier) {

  *ier = ARKSpilsSetPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMaxl(ARK_arkodemem, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPBCG;
  return;
}

/* Fortran interface to C "set" routines for the ARKMassSpbcg solver; 
   see farkode.h for further details */
void FARK_MASSSPBCGREINIT(int *pretype, int *maxl, 
			  realtype *delt, int *ier) {

  *ier = ARKSpilsSetMassPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassMaxl(ARK_arkodemem, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPBCG;
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines for the ARKSptfqmr
   solver; see farkode.h for further details */
void FARK_SPTFQMRREINIT(int *pretype, int *maxl, 
			realtype *delt, int *ier) {

  *ier = ARKSpilsSetPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMaxl(ARK_arkodemem, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPTFQMR;
  return;
}

/* Fortran interface to C "set" routines for the ARKMassSptfqmr
   solver; see farkode.h for further details */
void FARK_MASSSPTFQMRREINIT(int *pretype, int *maxl, 
			    realtype *delt, int *ier) {

  *ier = ARKSpilsSetMassPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassMaxl(ARK_arkodemem, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPTFQMR;
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines for the ARKSpfgmr solver; 
   see farkode.h for further details */
void FARK_SPFGMRREINIT(int *pretype, int *gstype, 
		       realtype *delt, int *ier) {

  *ier = ARKSpilsSetPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_SPFGMR;
  return;
}

/* Fortran interface to C "set" routines for the ARKMassSpfgmr solver; 
   see farkode.h for further details */
void FARK_MASSSPFGMRREINIT(int *pretype, int *gstype, 
			   realtype *delt, int *ier) {

  *ier = ARKSpilsSetMassPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassGSType(ARK_arkodemem, *gstype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_SPFGMR;
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines for the ARKPcg solver; 
   see farkode.h for further details */
void FARK_PCGREINIT(int *pretype, int *maxl, 
		    realtype *delt, int *ier) {

  *ier = ARKSpilsSetPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMaxl(ARK_arkodemem, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_ls = ARK_LS_PCG;
  return;
}

/* Fortran interface to C "set" routines for the ARKMassPcg solver; 
   see farkode.h for further details */
void FARK_MASSPCGREINIT(int *pretype, int *maxl, 
			realtype *delt, int *ier) {

  *ier = ARKSpilsSetMassPrecType(ARK_arkodemem, *pretype);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassMaxl(ARK_arkodemem, *maxl);
  if (*ier != ARKSPILS_SUCCESS) return;

  *ier = ARKSpilsSetMassEpsLin(ARK_arkodemem, *delt);
  if (*ier != ARKSPILS_SUCCESS) return;

  ARK_mass_ls = ARK_LS_PCG;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKode (the main integrator); 
   see farkode.h for further details */
void FARK_ARKODE(realtype *tout, realtype *t, realtype *y, 
		 int *itask, int *ier) {

  /* attach user solution array to solver memory */
  N_VSetArrayPointer(y, F2C_ARKODE_vec);

  /* call ARKode solver */
  *ier = ARKode(ARK_arkodemem, *tout, F2C_ARKODE_vec, t, *itask);

  /* detach user solution array from solver memory */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* Load optional outputs in iout & rout */
  ARKodeGetWorkSpace(ARK_arkodemem,
		     &ARK_iout[0],          /* LENRW   */
		     &ARK_iout[1]);         /* LENIW   */
  ARKodeGetIntegratorStats(ARK_arkodemem, 
			   &ARK_iout[2],    /* NST     */
			   &ARK_iout[3],    /* NST_STB */
			   &ARK_iout[4],    /* NST_ACC */ 
			   &ARK_iout[5],    /* NST_ATT */ 
			   &ARK_iout[6],    /* NFE     */ 
			   &ARK_iout[7],    /* NFI     */ 
			   &ARK_iout[8],    /* NSETUPS */ 
			   &ARK_iout[9],    /* NETF    */ 
			   &ARK_rout[0],    /* H0U     */
			   &ARK_rout[1],    /* HU      */ 
			   &ARK_rout[2],    /* HCUR    */ 
			   &ARK_rout[3]);   /* TCUR    */ 
  ARKodeGetTolScaleFactor(ARK_arkodemem, 
			  &ARK_rout[4]);    /* TOLSFAC */
  ARKodeGetNonlinSolvStats(ARK_arkodemem,
                          &ARK_iout[10],    /* NNI     */
                          &ARK_iout[11]);   /* NCFN    */
  
  /* If root finding is on, load those outputs as well */
  if (ARK_nrtfn != 0)
    ARKodeGetNumGEvals(ARK_arkodemem, &ARK_iout[12]);  /* NGE */

  /* Attach linear solver outputs */
  switch(ARK_ls) {
  case ARK_LS_DENSE:
  case ARK_LS_BAND:
  case ARK_LS_LAPACKDENSE:
  case ARK_LS_LAPACKBAND:
    ARKDlsGetWorkSpace(ARK_arkodemem, &ARK_iout[13], &ARK_iout[14]);  /* LENRWLS, LENIWLS */
    ARKDlsGetLastFlag(ARK_arkodemem, &ARK_iout[15]);                  /* LSTF  */
    ARKDlsGetNumRhsEvals(ARK_arkodemem, &ARK_iout[16]);               /* NFELS */
    ARKDlsGetNumJacEvals(ARK_arkodemem, &ARK_iout[17]);               /* NJE   */
    break;
  case ARK_LS_KLU:
  case ARK_LS_SUPERLUMT:
    ARKSlsGetLastFlag(ARK_arkodemem, &ARK_iout[15]);                  /* LSTF  */
    ARKSlsGetNumJacEvals(ARK_arkodemem, &ARK_iout[17]);               /* NJE   */
    break;
  case ARK_LS_SPGMR:
  case ARK_LS_SPBCG:
  case ARK_LS_SPTFQMR:
  case ARK_LS_SPFGMR:
  case ARK_LS_PCG:
    ARKSpilsGetWorkSpace(ARK_arkodemem, &ARK_iout[13], &ARK_iout[14]); /* LENRWLS, LENIWLS */
    ARKSpilsGetLastFlag(ARK_arkodemem, &ARK_iout[15]);                 /* LSTF  */
    ARKSpilsGetNumRhsEvals(ARK_arkodemem, &ARK_iout[16]);              /* NFELS */
    ARKSpilsGetNumJtimesEvals(ARK_arkodemem, &ARK_iout[17]);           /* NJTV  */
    ARKSpilsGetNumPrecEvals(ARK_arkodemem, &ARK_iout[18]);             /* NPE   */
    ARKSpilsGetNumPrecSolves(ARK_arkodemem, &ARK_iout[19]);            /* NPS   */
    ARKSpilsGetNumLinIters(ARK_arkodemem, &ARK_iout[20]);              /* NLI   */
    ARKSpilsGetNumConvFails(ARK_arkodemem, &ARK_iout[21]);             /* NCFL  */
  }

  /* Attach mass matrix linear solver outputs */
  switch(ARK_mass_ls) {
  case ARK_LS_DENSE:
  case ARK_LS_BAND:
  case ARK_LS_LAPACKDENSE:
  case ARK_LS_LAPACKBAND:
    ARKDlsGetMassWorkSpace(ARK_arkodemem, &ARK_iout[22], &ARK_iout[23]);  /* LENRWMS, LENIWMS */
    ARKDlsGetLastMassFlag(ARK_arkodemem, &ARK_iout[24]);                  /* LSTMF */
    ARKDlsGetNumMassEvals(ARK_arkodemem, &ARK_iout[25]);                  /* NME   */
    break;
  case ARK_LS_KLU:
  case ARK_LS_SUPERLUMT:
    ARKSlsGetLastMassFlag(ARK_arkodemem, &ARK_iout[24]);                  /* LSTMF */
    ARKSlsGetNumMassEvals(ARK_arkodemem, &ARK_iout[25]);                  /* NME   */
    break;
  case ARK_LS_SPGMR:
  case ARK_LS_SPBCG:
  case ARK_LS_SPTFQMR:
  case ARK_LS_SPFGMR:
  case ARK_LS_PCG:
    ARKSpilsGetMassWorkSpace(ARK_arkodemem, &ARK_iout[22], &ARK_iout[23]); /* LENRWMS, LENIWMS */
    ARKSpilsGetLastMassFlag(ARK_arkodemem, &ARK_iout[24]);                 /* LSTMF */
    ARKSpilsGetNumMassPrecEvals(ARK_arkodemem, &ARK_iout[25]);             /* NMPE  */
    ARKSpilsGetNumMassPrecSolves(ARK_arkodemem, &ARK_iout[26]);            /* NMPS  */
    ARKSpilsGetNumMassIters(ARK_arkodemem, &ARK_iout[27]);                 /* NMLI  */
    ARKSpilsGetNumMassConvFails(ARK_arkodemem, &ARK_iout[28]);             /* NMCFL */
  }
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeGetDky; see farkode.h 
   for further details */
void FARK_DKY(realtype *t, int *k, realtype *dky, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(dky, F2C_ARKODE_vec);

  /* call ARKodeGetDky */
  *ier = 0;
  *ier = ARKodeGetDky(ARK_arkodemem, *t, *k, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeGetErrWeights; see 
   farkode.h for further details */
void FARK_GETERRWEIGHTS(realtype *eweight, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(eweight, F2C_ARKODE_vec);

  /* call ARKodeGetErrWeights */
  *ier = 0;
  *ier = ARKodeGetErrWeights(ARK_arkodemem, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeGetEstLocalErrors; see 
   farkode.h for further details */
void FARK_GETESTLOCALERR(realtype *ele, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(ele, F2C_ARKODE_vec);

  /* call ARKodeGetEstLocalErrors */
  *ier = 0;
  *ier = ARKodeGetEstLocalErrors(ARK_arkodemem, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeFree; see farkode.h for 
   further details */
void FARK_FREE() {

  ARKodeMem ark_mem;
  ark_mem = (ARKodeMem) ARK_arkodemem;
  free(ark_mem->ark_user_data); ark_mem->ark_user_data = NULL;

  ARKodeFree(&ARK_arkodemem);

  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);
  N_VDestroy(F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKodeWriteParameters; see 
   farkode.h for further details */
void FARK_WRITEPARAMETERS(int *ier) {
  *ier = ARKodeWriteParameters(ARK_arkodemem, stdout);
  return;
}

/*=============================================================*/

/* C interface to user-supplied FORTRAN function FARKEFUN; see 
   farkode.h for further details */
int FARKfe(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  int ier;
  realtype *ydata, *dydata;
  FARKUserData ARK_userdata;
  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);
  ARK_userdata = (FARKUserData) user_data;

  FARK_EXP_FUN(&t, ydata, dydata, ARK_userdata->ipar, 
	       ARK_userdata->rpar, &ier);
  return(ier);
}

/*=============================================================*/

/* C interface to user-supplied FORTRAN function FARKIFUN; see 
   farkode.h for further details */
int FARKfi(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  int ier;
  realtype *ydata, *dydata;
  FARKUserData ARK_userdata;
  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);
  ARK_userdata = (FARKUserData) user_data;

  FARK_IMP_FUN(&t, ydata, dydata, ARK_userdata->ipar, 
	       ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
