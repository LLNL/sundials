/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2004-06-18 21:33:49 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * The C function FCVPSet is to interface between the CVSPGMR module  
 * and the user-supplied preconditioner setup routine FCVPSET.        
 * Note the use of the generic name FCV_PSET below.                   
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype                    */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvspgmr.h"       /* CVSpgmr prototype                               */

/*********************************************************************/

/* Prototype of the Fortran routine */
void FCV_PSET(realtype*, realtype*, realtype*, booleantype*, 
              booleantype*, realtype*, realtype*, realtype*, realtype*, 
              realtype*, realtype*, realtype*, int*);

/***************************************************************************/

void FCV_SPGMRSETPSET(int *flag, int *ier)
{
  if (*flag == 0) CVSpgmrSetPrecSetupFn(CV_cvodemem, NULL);
  else            CVSpgmrSetPrecSetupFn(CV_cvodemem, FCVPSet);
}


/***************************************************************************/

/* C function FCVPSet to interface between CVODE and a Fortran subroutine
   FCVPSET for setup of a Krylov preconditioner.
   Addresses of Nlocal, t, jok, gamma, h, uround, y, fy, ewt, vtemp1, vtemp2, 
   vtemp3, and the address jcurPtr are passed to FCVPSET, using
   the routine N_VGetData from NVECTOR.  A return flag ier from FCVPSET
   is returned by FCVPSet.
   Auxiliary data is assumed to be communicated by common blocks. */

int FCVPSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma,
            void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  N_Vector ewt;
  realtype h, uround;
  realtype *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  int ier = 0;

  CVodeGetErrWeights(CV_cvodemem, &ewt);
  CVodeGetLastStep(CV_cvodemem, &h);
  uround = UNIT_ROUNDOFF;

  ydata   = (realtype *) N_VGetData(y);
  fydata  = (realtype *) N_VGetData(fy);
  ewtdata = (realtype *) N_VGetData(ewt);
  v1data  = (realtype *) N_VGetData(vtemp1);
  v2data  = (realtype *) N_VGetData(vtemp2);
  v3data  = (realtype *) N_VGetData(vtemp3);

  FCV_PSET(&t, ydata, fydata, &jok, jcurPtr, &gamma, ewtdata,
           &h, &uround, v1data, v2data, v3data, &ier);

  return(ier);
}
