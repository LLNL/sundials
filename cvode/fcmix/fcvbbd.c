/*
 * -----------------------------------------------------------------
 * $Revision: 1.17 $
 * $Date: 2004-10-21 20:55:05 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
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

#include "cvbbdpre.h"      /* prototypes of CVBBDPRE functions and macros */
#include "cvode.h"         /* CVODE constants and prototypes              */
#include "cvspgmr.h"       /* prototypes of CVSPGMR interface routines    */
#include "fcvbbd.h"        /* prototypes of interfaces to CVBBDPRE        */
#include "fcvode.h"        /* actual function names, prototypes and
			      global variables                            */
#include "nvector.h"       /* definition of type N_Vector                 */
#include "sundialstypes.h" /* definition of type realtype                 */

/***************************************************************************/

/* Prototypes of the Fortran routines */

extern void FCV_GLOCFN(long int*, realtype*, realtype*, realtype*);
extern void FCV_COMMFN(long int*, realtype*, realtype*);

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
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPGMR_SUCCESS) return;

  CV_ls = 4;
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

  *ier = CVBBDPrecReInit(CVBBD_Data, *mudq, *mldq,
                         *dqrely, FCVgloc, FCVcfn);
}

/***************************************************************************/

/* C function FCVgloc to interface between CVBBDPRE module and a Fortran 
   subroutine FCVLOCFN. */

void FCVgloc(long int Nloc, realtype t, N_Vector yloc, N_Vector gloc,
             void *f_data)
{
  realtype *yloc_data, *gloc_data;
  
  yloc_data = N_VGetArrayPointer(yloc);
  gloc_data = N_VGetArrayPointer(gloc);

  FCV_GLOCFN(&Nloc, &t, yloc_data, gloc_data);
}

/***************************************************************************/

/* C function FCVcfn to interface between CVBBDPRE module and a Fortran 
   subroutine FCVCOMMF. */


void FCVcfn(long int Nloc, realtype t, N_Vector y, void *f_data)
{
  realtype *yloc;

  yloc = N_VGetArrayPointer(y);

  FCV_COMMFN(&Nloc, &t, yloc);

}

/***************************************************************************/

/* C function FCVBBDOPT to access optional outputs from CVBBD_Data */

void FCV_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  CVBBDPrecGetWorkSpace(CVBBD_Data, lenrpw, lenipw);
  CVBBDPrecGetNumGfnEvals(CVBBD_Data, nge);
}


/***************************************************************************/

/* C function FCVBBDFREE to interface to CVBBDPrecFree, to free memory 
   created by CVBBDPrecAlloc */

void FCV_BBDFREE()
{
  CVBBDPrecFree(CVBBD_Data);
}
