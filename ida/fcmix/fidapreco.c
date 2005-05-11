/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-05-11 23:10:54 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * The C function FIDAPSet is to interface between the IDASPGMR/IDASPBCG
 * module and the user-supplied preconditioner setup routine FIDAPSET.
 * Note the use of the generic name FIDA_PSET below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "idaspgmr.h"       /* IDASPGMR prototypes                            */
#include "idaspbcg.h"       /* IDASPBCG prototypes                            */
#include "fida.h"           /* actual function names, prototypes and global
			       variables                                      */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FIDA_PSET(realtype*, realtype*, realtype*, realtype*,
		      realtype*, realtype*, realtype*, realtype*,
		      realtype*, realtype*, int*);

extern void FIDA_PSOL(realtype*, realtype*, realtype*, realtype*,
		      realtype*, realtype*, realtype*, realtype*,
		      realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_SPGMRSETPREC(int *flag, int *ier)
{
  *ier = 0;

  if (*flag == 0) IDASpgmrSetPreconditioner(IDA_idamem, NULL, NULL, NULL);
  else {
    IDASpgmrSetPreconditioner(IDA_idamem, (IDASpilsPrecSetupFn) FIDAPSet,
			      (IDASpilsPrecSolveFn) FIDAPSol, NULL);
    if (F2C_IDA_ewtvec == NULL) F2C_IDA_ewtvec = N_VClone(F2C_IDA_ewtvec);
  }

  return;
}

/*************************************************/

void FIDA_SPBCGSETPREC(int *flag, int *ier)
{
  *ier = 0;

  if (*flag == 0) IDASpbcgSetPreconditioner(IDA_idamem, NULL, NULL, NULL);
  else {
    IDASpbcgSetPreconditioner(IDA_idamem, (IDASpilsPrecSetupFn) FIDAPSet,
			      (IDASpilsPrecSolveFn) FIDAPSol, NULL);
    if (F2C_IDA_ewtvec == NULL) F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
  }

  return;
}

/*************************************************/

int FIDAPSet(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	     realtype c_j, void *prec_data,
	     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *yy_data, *yp_data, *rr_data, *ewtdata, *v1data, *v2data, *v3data;
  realtype h;
  int ier;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = ewtdata = NULL;
  v1data = v2data = v3data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, yy, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  v1data = N_VGetArrayPointer(vtemp1);
  v2data = N_VGetArrayPointer(vtemp2);
  v3data = N_VGetArrayPointer(vtemp3);

  /* Call user-supplied routine */
  FIDA_PSET(&t, yy_data, yp_data, rr_data, &c_j, ewtdata, &h,
	    v1data, v2data, v3data, &ier);

  return(ier);
}

/*************************************************/

int FIDAPSol(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	     N_Vector rvec, N_Vector zvec,
	     realtype c_j, realtype delta, void *prec_data,
	     N_Vector vtemp1)
{
  realtype *yy_data, *yp_data, *rr_data, *ewtdata, *rdata, *zdata, *v1data;
  int ier;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = ewtdata = rdata = zdata = v1data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, yy, F2C_IDA_ewtvec);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  rdata   = N_VGetArrayPointer(rvec);
  zdata   = N_VGetArrayPointer(zvec);
  v1data  = N_VGetArrayPointer(vtemp1);

  /* Call user-supplied routine */
  FIDA_PSOL(&t, yy_data, yp_data, rr_data, rdata, zdata,
	    &c_j, &delta, ewtdata, v1data, &ier);

  return(ier);
}
