/*
 * -----------------------------------------------------------------
 * $Revision: 1.3.2.1 $
 * $Date: 2005-01-26 22:05:10 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Michael Wittman, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the CVBBDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CVBBDPRE_IMPL_H
#define _CVBBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvbbdpre.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type: CVBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* passed by user to CVBBDPrecAlloc, used by PrecSetup/PrecSolve */
  long int mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  CVLocalFn gloc;
  CVCommFn cfn;

  /* set by CVBBDPrecSetup and used by CVBBDPrecSolve */
  BandMat savedJ;
  BandMat savedP;
  long int *pivots;

  /* set by CVBBDPrecAlloc and used by CVBBDPrecSetup */
  long int n_local;

  /* available for optional output: */
  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* Pointer to cvode_mem */
  void *cvode_mem;

} *CVBBDPrecData;


/* Error Messages */

#define _CVBBDALLOC_        "CVBBDAlloc-- "
#define MSGBBDP_CVMEM_NULL  _CVBBDALLOC_ "Integrator memory is NULL.\n\n"
#define MSGBBDP_BAD_NVECTOR _CVBBDALLOC_ "A required vector operation is not implemented.\n\n"

#define MSGBBDP_PDATA_NULL "CVBBDPrecGet*-- BBDPrecData is NULL.\n\n"

#define MSGBBDP_NO_PDATA "CVBBDSpgmr-- BBDPrecData is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
