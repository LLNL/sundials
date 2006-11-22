/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-22 00:12:50 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for IDA/IDALAPACK, for the case of
 * a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* function names, prototypes, global vars.*/
#include "ida_impl.h" /* definition of IDAMem type               */

#include <ida/ida_band.h>

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_BJAC(int*, int*, int*, int*,
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
  *ier = 0;

  if (*flag == 0) {
    *ier = IDADlsSetJacFn(IDA_idamem, NULL, NULL);
  } else {
    if (F2C_IDA_ewtvec == NULL) {
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
      if (F2C_IDA_ewtvec == NULL) {
        *ier = -1;
        return;
      }
    }
    *ier = IDADlsSetJacFn(IDA_idamem, (IDADlsBandJacFn) FIDABandJac, NULL);
  }

  return;
}

/*************************************************/

int FIDALapackBandJac(int N, int mupper, int mlower,
                      realtype t, realtype c_j, 
                      N_Vector yy, N_Vector yp, N_Vector rr,
                      DlsMat J, void *jac_data,
                      N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *yy_data, *yp_data, *rr_data, *jacdata, *ewtdata, *v1data, *v2data, *v3data;
  realtype h;
  int eband;
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

  eband = (J->s_mu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  IDA_userdata = (FIDAUserData) jac_data;

  /* Call user-supplied routine */
  FIDA_BJAC(&N, &mupper, &mlower, &eband, &t, yy_data, yp_data, rr_data,
            &c_j, jacdata, ewtdata, &h, 
            IDA_userdata->ipar, IDA_userdata->rpar,
            v1data, v2data, v3data, &ier);

  return(ier);
}
