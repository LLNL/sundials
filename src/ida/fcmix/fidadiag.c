/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Aaron Collier @ LLNL
 *-----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 *-----------------------------------------------------------------
 * Fortran/C interface routines for IDA/IDADLS, for the case
 * of a user-supplied Jacobian approximation routine.
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* actual function names, prototypes and global vars.*/
#include "ida_impl.h" /* definition of IDAMem type                         */

#include <ida/ida_direct.h>
#include <sunmatrix/sunmatrix_diagonal.h>

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_DIAGJAC(realtype* T, realtype* Y, realtype* YP, 
                           realtype* R, realtype* J, realtype* CJ, 
                           realtype* EWT, realtype* H, long int* IPAR, 
                           realtype* RPAR, realtype* V1, 
                           realtype* V2, realtype* V3, int* IER);

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_DIAGSETJAC(int *ier)
{
  if (F2C_IDA_ewtvec == NULL) {
    F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
    if (F2C_IDA_ewtvec == NULL) {
      *ier = -1;
      return;
    }
  }
  *ier = IDADlsSetJacFn(IDA_idamem, FIDADiagJac);
  return;
}

/*************************************************/

int FIDADiagJac(realtype t, realtype c_j, N_Vector yy, N_Vector yp,
                N_Vector rr, SUNMatrix J, void *user_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *ypdata, *rdata, *jacdata, *ewtdata, *v1data, *v2data, *v3data;
  realtype h;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  ydata = ypdata = rdata = jacdata = ewtdata = NULL;
  v1data = v2data = v3data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);

  /* Get pointers to vector data */
  ydata   = N_VGetArrayPointer(yy);
  ypdata  = N_VGetArrayPointer(yp);
  rdata   = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  jacdata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(J));

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine*/
  FIDA_DIAGJAC(&t, ydata, ypdata, rdata,
            jacdata, &c_j, ewtdata, &h, 
            IDA_userdata->ipar, IDA_userdata->rpar,
            v1data, v2data, v3data, &ier);

  return(ier);
}
