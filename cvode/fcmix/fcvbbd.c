/*
 * -----------------------------------------------------------------
 * $Revision: 1.30 $
 * $Date: 2006-02-10 00:03:09 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with the
 * CVBBDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide
 * a standard interface to the C code of the CVBBDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"           /* actual function names, prototypes, global vars.*/
#include "fcvbbd.h"           /* prototypes of interfaces to CVBBDPRE           */

#include "cvode.h"            /* CVODE constants and prototypes                 */
#include "cvode_bbdpre.h"     /* prototypes of CVBBDPRE functions and macros    */

#include "cvode_sptfqmr.h"    /* prototypes of CVSPTFQMR interface routines     */
#include "cvode_spbcgs.h"     /* prototypes of CVSPBCG interface routines       */
#include "cvode_spgmr.h"      /* prototypes of CVSPGMR interface routines       */

#include "sundials_nvector.h" /* definition of type N_Vector                    */
#include "sundials_types.h"   /* definition of type realtype                    */

/***************************************************************************/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_GLOCFN(long int*,                        /* NLOC          */
                         realtype*, realtype*, realtype*,  /* T, YLOC, GLOC */
                         long int*, realtype*,             /* IPAR, RPAR    */
                         int *ier);                        /* IER           */

  extern void FCV_COMMFN(long int*,                        /* NLOC          */
                         realtype*, realtype*,             /* T, Y          */
                         long int*, realtype*,             /* IPAR, RPAR    */
                         int *ier);                        /* IER           */

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_BBDINIT(long int *Nloc, long int *mudq, long int *mldq, 
                 long int *mu, long int *ml, realtype* dqrely, int *ier)
{

  /* 
     First call CVBBDPrecAlloc to initialize CVBBDPRE module:
     Nloc       is the local vector size
     mudq,mldq  are the half-bandwidths for computing preconditioner blocks
     mu, ml     are the half-bandwidths of the retained preconditioner blocks
     dqrely     is the difference quotient relative increment factor
     FCVgloc    is a pointer to the CVLocalFn function
     FCVcfn     is a pointer to the CVCommFn function 
  */

  CVBBD_Data = CVBBDPrecAlloc(CV_cvodemem, *Nloc, *mudq, *mldq, *mu, *ml, 
                              *dqrely, FCVgloc, FCVcfn);
  if (CVBBD_Data == NULL) *ier = -1; 
  else                    *ier = 0;

  return; 
}

/***************************************************************************/

void FCV_BBDSPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     Call CVBBDSptfqmr to specify the SPTFQMR linear solver:
     pretype    is the preconditioner type
     maxl       is the maximum Krylov dimension
     delt       is the linear convergence tolerance factor 
  */

  *ier = CVBBDSptfqmr(CV_cvodemem, *pretype, *maxl, CVBBD_Data);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPTFQMR;
}

/***************************************************************************/

void FCV_BBDSPBCG(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     Call CVBBDSpbcg to specify the SPBCG linear solver:
     pretype    is the preconditioner type
     maxl       is the maximum Krylov dimension
     delt       is the linear convergence tolerance factor 
  */

  *ier = CVBBDSpbcg(CV_cvodemem, *pretype, *maxl, CVBBD_Data);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPBCG;
}

/***************************************************************************/

void FCV_BBDSPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{
  /* 
     Call CVBBDSpgmr to specify the SPGMR linear solver:
     pretype    is the preconditioner type
     gstype     is the Gram-Schmidt process type
     maxl       is the maximum Krylov dimension
     delt       is the linear convergence tolerance factor 
  */

  *ier = CVBBDSpgmr(CV_cvodemem, *pretype, *maxl, CVBBD_Data);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPGMR;
}

/***************************************************************************/

void FCV_BBDREINIT(long int *Nloc, long int *mudq, long int *mldq, 
                   realtype* dqrely, int *ier)
{
  /* 
     First call CVReInitBBD to re-initialize CVBBDPRE module:
     CVBBD_Data  is the pointer to P_data
     mudq,mldq   are the half-bandwidths for computing preconditioner blocks
     dqrely      is the difference quotient relative increment factor
     FCVgloc     is a pointer to the CVLocalFn function
     FCVcfn      is a pointer to the CVCommFn function 
  */

  *ier = CVBBDPrecReInit(CVBBD_Data, *mudq, *mldq, *dqrely, FCVgloc, FCVcfn);
}

/***************************************************************************/

/* C function FCVgloc to interface between CVBBDPRE module and a Fortran 
   subroutine FCVLOCFN. */

int FCVgloc(long int Nloc, realtype t, N_Vector yloc, N_Vector gloc,
            void *f_data)
{
  int ier;
  realtype *yloc_data, *gloc_data;
  FCVUserData CV_userdata;

  yloc_data = N_VGetArrayPointer(yloc);
  gloc_data = N_VGetArrayPointer(gloc);

  CV_userdata = (FCVUserData) f_data;

  FCV_GLOCFN(&Nloc, &t, yloc_data, gloc_data, 
             CV_userdata->ipar, CV_userdata->rpar,
             &ier);
  return(ier);
}

/***************************************************************************/

/* C function FCVcfn to interface between CVBBDPRE module and a Fortran 
   subroutine FCVCOMMF. */

int FCVcfn(long int Nloc, realtype t, N_Vector y, void *f_data)
{
  int ier;
  realtype *yloc;
  FCVUserData CV_userdata;

  yloc = N_VGetArrayPointer(y);

  CV_userdata = (FCVUserData) f_data;

  FCV_COMMFN(&Nloc, &t, yloc, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}

/***************************************************************************/

/* C function FCVBBDOPT to access optional outputs from CVBBD_Data */

void FCV_BBDOPT(long int *lenrwbbd, long int *leniwbbd, long int *ngebbd)
{
  CVBBDPrecGetWorkSpace(CVBBD_Data, lenrwbbd, leniwbbd);
  CVBBDPrecGetNumGfnEvals(CVBBD_Data, ngebbd);
}


/***************************************************************************/

/* C function FCVBBDFREE to interface to CVBBDPrecFree, to free memory 
   created by CVBBDPrecAlloc */

void FCV_BBDFREE(void)
{
  CVBBDPrecFree(&CVBBD_Data);
}
