/******************************************************************
 *                                                                *
 * File          : sensfpvbbd.c                                   *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 25 August 2000                                 *
 *----------------------------------------------------------------*
 * Fortran/C interface routines for use of sensitivity variant of *
 * PVODE with PVBBDPRE preconditioner module                      *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "cvode.h"    /* CVODE constants and prototypes                    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "fpvbbd.h"   /* prototypes of interfaces to PVBBDPRE              */
#include "pvbbdpre.h" /* prototypes of PVBBDPRE functions, macros          */
#include "mpi.h"      /* MPI types and constants                           */
#include "sensfpvbbd.h" /* sensitivity version of PVBBDPRE interface       */
#include "sensitivity.h" /* sensitivity data types and prototypes          */

/***************************************************************************/

void SFPV_BBDIN (integer *mudq, integer *mldq, integer *mu, integer *ml,
                real* dqrely, int *pretype, int *gstype, int *maxl,
                real *delt, int *ier)
{
  integer Nloc;
  CVodeMem cv_mem;
  void *s_data;

  cv_mem = (CVodeMem) CV_cvodemem;
  s_data = cv_mem->cv_f_data;

  /* First call PVBBDAlloc to initialize PVBBDPRE module:
     *Nloc       is the local vector size
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     SPVgloc     is a pointer to the SPVLocalFn function
     PVcfn       is a pointer to the PVCommFn function
     s_data      is the pointer to s_data                             */

  Nloc = ((machEnvType)PV_machEnv)->local_vec_length;
  PVBBD_Data = PVBBDAlloc (Nloc, *mudq, *mldq, *mu, *ml, *dqrely, 
                           SPVgloc, PVcfn, s_data);
  if (PVBBD_Data == NULL) {
    *ier = -1;
    return;
  }

  /* Call SensCVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretyp     is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     PVBBDPrecon is a pointer to the preconditioner setup routine
     PVBBDPSol   is a pointer to the preconditioner solve routine
     PVBBD_Data  is the pointer to P_data                             */

  SensCVSpgmr(CV_cvodemem, *pretype, *gstype, *maxl, *delt, 
           PVBBDPrecon, PVBBDPSol, PVBBD_Data);

  *ier = 0;
}

/***************************************************************************/

/* C function SPVgloc to interface between PVBBDPRE module and a Fortran 
   subroutine PVLOCFN. */

void SPVgloc(integer Nloc, real t, real *yloc, real *gloc, void *s_data)
{
  real *p;
  SensData sdata;

  sdata = (SensData) s_data;
  p = sdata->p;

  SFPV_GLOCFN (&Nloc, &t, yloc, gloc, p);
}
