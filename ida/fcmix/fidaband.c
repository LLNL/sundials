/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2006-01-24 22:17:29 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for IDA/IDABAND, for the case of
 * a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_band.h"         /* IDABAND prototypes                             */
#include "fida.h"             /* function names, prototypes, global variables   */
#include "ida_impl.h"         /* definition of IDAMem type                      */
#include "sundials_nvector.h" /* definitions of type N_Vector and vector macros */
#include "sundials_types.h"   /* definition of type realtype                    */

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_BJAC(long int*, long int*, long int*, long int*,
                        realtype*, realtype*, realtype*, realtype*,
                        realtype*, realtype*, realtype*, realtype*,
                        long int*, realtype*,
                        realtype*, realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_BANDSETJAC(int *flag, int *ier)
{
  IDAMem ida_mem;

  *ier = 0;

  if (*flag == 0) *ier = IDABandSetJacFn(IDA_idamem, NULL, NULL);
  else {
    ida_mem = (IDAMem) IDA_idamem;    
    *ier = IDABandSetJacFn(IDA_idamem, (IDABandJacFn) FIDABandJac, NULL);
    if (F2C_IDA_ewtvec == NULL) 
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
  }

  return;
}

/*************************************************/

int FIDABandJac(long int N, long int mupper, long int mlower,
		BandMat J, realtype t,
		N_Vector yy, N_Vector yp, N_Vector rr,
		realtype c_j, void *jac_data, BandMat Jac,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *yy_data, *yp_data, *rr_data, *jacdata, *ewtdata, *v1data, *v2data, *v3data;
  realtype h;
  long int eband;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = jacdata = ewtdata = NULL;
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
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  eband = (J->smu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  IDA_userdata = (FIDAUserData) jac_data;

  /* Call user-supplied routine */
  FIDA_BJAC(&N, &mupper, &mlower, &eband, &t, yy_data, yp_data, rr_data,
	    &c_j, jacdata, ewtdata, &h, 
            IDA_userdata->ipar, IDA_userdata->rpar,
            v1data, v2data, v3data, &ier);

  return(ier);
}
