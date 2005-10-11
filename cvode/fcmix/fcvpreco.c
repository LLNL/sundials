/*
 * -----------------------------------------------------------------
 * $Revision: 1.22 $
 * $Date: 2005-10-11 16:02:39 $
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
#include "fcvode.h"         /* actual fn. names, prototypes and global vars.  */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */
#include "cvode_impl.h"     /* definition of CVodeMem type                    */

/*********************************************************************/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_PSET(realtype*, realtype*, realtype*,  /* T, Y, FY */
                       booleantype*, booleantype*,       /* JOK, JCUR */
                       realtype*, realtype*,             /* GAMMA, H */
                       long int*, realtype*,             /* IPAR, RPAR */
                       realtype*, realtype*, realtype*,  /* W1, W2, W3 */
                       int*);                            /* IER */

  extern void FCV_PSOL(realtype*, realtype*, realtype*,  /* T, Y, FY */
                       realtype*, realtype*,             /* R, Z */
                       realtype*, realtype*,             /* GAMMA, DELTA */
                       int*,                             /* LR */
                       long int*, realtype*,             /* IPAR, RPAR */
                       realtype*,                        /* WRK */
                       int*);                            /* IER */

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_SPGMRSETPREC(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag == 0) {
    *ier = CVSpgmrSetPreconditioner(CV_cvodemem, NULL, NULL, NULL);
  } else {
    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVSpgmrSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol, cv_mem->cv_f_data);
  }
}

/***************************************************************************/

void FCV_SPBCGSETPREC(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag == 0) {
    *ier = CVSpbcgSetPreconditioner(CV_cvodemem, NULL, NULL, NULL);
  } else {
    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVSpbcgSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol, cv_mem->cv_f_data);
  }
}

/***************************************************************************/

void FCV_SPTFQMRSETPREC(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag == 0) {
    *ier = CVSptfqmrSetPreconditioner(CV_cvodemem, NULL, NULL, NULL);
  } else {
    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVSptfqmrSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol, cv_mem->cv_f_data);
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
  int ier = 0;
  realtype *ydata, *fydata, *v1data, *v2data, *v3data;
  realtype h;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  CV_userdata = (FCVUserData) P_data;

  FCV_PSET(&t, ydata, fydata, &jok, jcurPtr, &gamma, &h,
           CV_userdata->ipar, CV_userdata->rpar,
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
  int ier = 0;
  realtype *ydata, *fydata, *vtdata, *rdata, *zdata;
  FCVUserData CV_userdata;

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  vtdata  = N_VGetArrayPointer(vtemp);
  rdata   = N_VGetArrayPointer(r);
  zdata   = N_VGetArrayPointer(z);

  CV_userdata = (FCVUserData) P_data;

  FCV_PSOL(&t, ydata, fydata, rdata, zdata, &gamma, &delta, &lr, 
           CV_userdata->ipar, CV_userdata->rpar, vtdata, &ier);

  return(ier);
}
