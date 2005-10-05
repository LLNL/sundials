/*
 * -----------------------------------------------------------------
 * $Revision: 1.21 $
 * $Date: 2005-10-05 20:31:20 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * The C function FCVPSet is to interface between the CVSP*
 * module and the user-supplied preconditioner setup routine FCVPSET.
 * Note the use of the generic name FCV_PSET below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvsptfqmr.h"      /* CVSptfqmr prototype                            */
#include "cvspbcg.h"        /* CVSpbcg prototype                              */
#include "cvspgmr.h"        /* CVSpgmr prototype                              */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                               */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/*********************************************************************/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FCV_PSET(realtype*, realtype*, realtype*, booleantype*, 
		     booleantype*, realtype*, realtype*, realtype*,
		     realtype*, realtype*, int*);
extern void FCV_PSOL(realtype*, realtype*, realtype*, realtype*,  realtype*,
		     realtype*, realtype*, int*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_SPGMRSETPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSpgmrSetPreconditioner(CV_cvodemem, NULL, NULL, NULL);
  } else {
    *ier = CVSpgmrSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol, NULL);
  }
}

/***************************************************************************/

void FCV_SPBCGSETPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSpbcgSetPreconditioner(CV_cvodemem, NULL, NULL, NULL);
  } else {
    *ier = CVSpbcgSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol, NULL);
  }
}

/***************************************************************************/

void FCV_SPTFQMRSETPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSptfqmrSetPreconditioner(CV_cvodemem, NULL, NULL, NULL);
  } else {
    *ier = CVSptfqmrSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol, NULL);
  }
}

/***************************************************************************/

/* C function FCVPSet to interface between CVODE and a Fortran subroutine
   FCVPSET for setup of a Krylov preconditioner.
   Addresses of t, y, fy, jok, gamma, h, vtemp1, vtemp2, vtemp3, and the 
   address jcurPtr are passed to FCVPSET, using the routine
   N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVPSET is returned by FCVPSet.
   Auxiliary data is assumed to be communicated by common blocks. */

int FCVPSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma,
            void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype h;
  realtype *ydata, *fydata, *v1data, *v2data, *v3data;

  int ier = 0;
  
  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  FCV_PSET(&t, ydata, fydata, &jok, jcurPtr, &gamma, &h,
           v1data, v2data, v3data, &ier);

  return(ier);
}

/***************************************************************************/

/* C function FCVPSol to interface between CVODE and a Fortran subroutine
   FCVPSOL for solution of a Krylov preconditioner.
   Addresses of t, y, fy, gamma, delta, lr, vtemp, r, and z are
   passed to FCVPSOL, using the routine N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVPSOL is returned by FCVPSol.
   Auxiliary data is assumed to be communicated by Common blocks. */

int FCVPSol(realtype t, N_Vector y, N_Vector fy, 
            N_Vector r, N_Vector z,
            realtype gamma, realtype delta,
            int lr, void *P_data, N_Vector vtemp)
{
  realtype *ydata, *fydata, *vtdata, *rdata, *zdata;

  int ier = 0;

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  vtdata  = N_VGetArrayPointer(vtemp);
  rdata   = N_VGetArrayPointer(r);
  zdata   = N_VGetArrayPointer(z);

  FCV_PSOL(&t, ydata, fydata, vtdata, &gamma, &delta, rdata, &lr, zdata,
           &ier);

  return(ier);
}
