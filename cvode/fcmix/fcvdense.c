/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2004-06-21 23:07:12 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVDENSE, for the case   
 * of a user-supplied Jacobian approximation routine.            
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype                    */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvdense.h"       /* CVDense prototype, type DenseMat                */

/***************************************************************************/

/* Prototype of the Fortran routine */
void FCV_DJAC(long int*, realtype*, realtype*, realtype*, realtype*, 
              realtype*, realtype*, realtype*, realtype*);

/***************************************************************************/

void FCV_DENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) CVDenseSetJacFn(CV_cvodemem, NULL);
  else            CVDenseSetJacFn(CV_cvodemem, FCVDenseJac);
}

/***************************************************************************/

/* C function CVDenseJac interfaces between CVODE and a Fortran subroutine
   FCVDJAC for solution of a linear system with dense Jacobian approximation.
   Addresses of arguments are passed to FCVDJAC, using the macro 
   DENSE_COL from DENSE and the routine N_VGetData from NVECTOR.
   Auxiliary data is assumed to be communicated by Common. */

void FCVDenseJac(long int N, DenseMat J, realtype t, 
                 N_Vector y, N_Vector fy, void *jac_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  N_Vector ewt;
  realtype *ydata, *fydata, *jacdata, *ewtdata, *v1data, *v2data, *v3data;

  CVodeGetErrWeights(CV_cvodemem, &ewt);

  ydata   = (realtype *) N_VGetData(y);
  fydata  = (realtype *) N_VGetData(fy);
  v1data  = (realtype *) N_VGetData(vtemp1);
  v2data  = (realtype *) N_VGetData(vtemp2);
  v3data  = (realtype *) N_VGetData(vtemp3);
  ewtdata = (realtype *) N_VGetData(ewt);

  jacdata = DENSE_COL(J,0);


  FCV_DJAC(&N, &t, ydata, fydata, jacdata, ewtdata, v1data, v2data, v3data); 

}
