/* File fpvbbdin1.c:  Fortran/C and C/Fortran interface routine for use
                      of PVODE with PVBBDPRE preconditioner module, for
                      the case of a user-supplied Jtimes routine.
   Version of 11 January 2002 */

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "cvode.h"    /* CVODE constants and prototypes */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "fpvbbd.h"   /* prototypes of interfaces to PVBBDPRE  */
#include "pvbbdpre.h" /* prototypes of PVBBDPRE functions, macros  */
#include "mpi.h"      /* MPI types and constants  */

/***************************************************************************/

void FPV_BBDIN1 (integer *mudq, integer *mldq, integer *mu, integer *ml,
                real* dqrely, int *pretype, int *gstype, int *maxl,
                real *delt, int *ier)
{
  integer Nloc;

  /* First call PVBBDAlloc to initialize PVBBDPRE module:
     *Nloc       is the local vector size
     *mudq,*mldq are the half-bandwidths for computing preconditioner blocks
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     *dqrely     is the difference quotient relative increment factor
     PVgloc      is a pointer to the PVLocalFn function
     PVcfn       is a pointer to the PVCommFn function
     NULL        is the pointer to f_data                             */

  Nloc = ((machEnvType)PV_machEnv)->local_vec_length;
  PVBBD_Data = PVBBDAlloc (Nloc, *mudq, *mldq, *mu, *ml, *dqrely, 
                           PVgloc, PVcfn, NULL);
  if (PVBBD_Data == NULL) {
    *ier = -1;
    return;
  }

  /* Call CVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretyp     is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     PVBBDPrecon is a pointer to the preconditioner setup routine
     PVBBDPSol   is a pointer to the preconditioner solve routine
     PVBBD_Data  is the pointer to P_data
     CVjtimes    is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                               */

  *ier = CVSpgmr (CV_cvodemem, *pretype, *gstype, *maxl, *delt, 
           PVBBDPrecon, PVBBDPSol, PVBBD_Data, CVJtimes, NULL);

}
