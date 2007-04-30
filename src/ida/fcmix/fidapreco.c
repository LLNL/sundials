/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2007-04-30 19:29:00 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * The C function FIDAPSet is to interface between the IDASPILS
 * modules and the user-supplied preconditioner setup routine FIDAPSET.
 * Note the use of the generic name FIDA_PSET below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* actual fn. names, prototypes and global vars.*/
#include "ida_impl.h" /* definition of IDAMem type                    */

#include <ida/ida_spils.h>

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_PSET(realtype*, realtype*, realtype*, realtype*,
                        realtype*, realtype*, realtype*, 
                        long int*, realtype*,
                        realtype*, realtype*, realtype*, 
                        int*);
  
  extern void FIDA_PSOL(realtype*, realtype*, realtype*, realtype*,
                        realtype*, realtype*, realtype*, realtype*,
                        realtype*, 
                        long int*, realtype*,
                        realtype*, int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_SPILSSETPREC(int *flag, int *ier)
{
  *ier = 0;

  if (*flag == 0) {

    *ier = IDASpilsSetPreconditioner(IDA_idamem, NULL, NULL);

  } else {

    if (F2C_IDA_ewtvec == NULL) {
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
      if (F2C_IDA_ewtvec == NULL) {
        *ier = -1;
        return;
      }
    }

    *ier = IDASpilsSetPreconditioner(IDA_idamem, (IDASpilsPrecSetupFn) FIDAPSet,
				     (IDASpilsPrecSolveFn) FIDAPSol);
  }

  return;
}

/*************************************************/

int FIDAPSet(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	     realtype c_j, void *user_data,
	     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *yy_data, *yp_data, *rr_data, *ewtdata, *v1data, *v2data, *v3data;
  realtype h;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = ewtdata = NULL;
  v1data = v2data = v3data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  v1data = N_VGetArrayPointer(vtemp1);
  v2data = N_VGetArrayPointer(vtemp2);
  v3data = N_VGetArrayPointer(vtemp3);

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_PSET(&t, yy_data, yp_data, rr_data, &c_j, ewtdata, &h,
            IDA_userdata->ipar, IDA_userdata->rpar,
	    v1data, v2data, v3data, &ier);

  return(ier);
}

/*************************************************/

int FIDAPSol(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	     N_Vector rvec, N_Vector zvec,
	     realtype c_j, realtype delta, void *user_data,
	     N_Vector vtemp1)
{
  realtype *yy_data, *yp_data, *rr_data, *ewtdata, *rdata, *zdata, *v1data;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = ewtdata = user_data = zdata = v1data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  rdata   = N_VGetArrayPointer(rvec);
  zdata   = N_VGetArrayPointer(zvec);
  v1data  = N_VGetArrayPointer(vtemp1);

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_PSOL(&t, yy_data, yp_data, rr_data, rdata, zdata,
	    &c_j, &delta, ewtdata, 
            IDA_userdata->ipar, IDA_userdata->rpar,
            v1data, &ier);

  return(ier);
}
