/*
 * -----------------------------------------------------------------
 * $Revision: 1.2.2.1 $
 * $Date: 2005-01-26 22:05:06 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Michael Wittman, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the CVBANDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CVBANDPRE_IMPL_H
#define _CVBANDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvbandpre.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type: CVBandPrecData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* Data set by user in CVBandPrecAlloc: */
  long int N;
  long int ml, mu;

  /* Data set by CVBandPrecSetup: */
  BandMat savedJ;
  BandMat savedP;
  long int *pivots;

  /* Rhs calls */
  long int nfeBP;

  /* Pointer to cvode_mem */
  void *cvode_mem;

} *CVBandPrecData;

/* Error Messages */

#define _CVBALLOC_        "CVBandPreAlloc-- "
#define MSGBP_CVMEM_NULL  _CVBALLOC_ "Integrator memory is NULL.\n\n"
#define MSGBP_BAD_NVECTOR _CVBALLOC_ "A required vector operation is not implemented.\n\n"

#define MSGBP_PDATA_NULL "CVBandPrecGet*-- BandPrecData is NULL.\n\n"

#define MSGBP_NO_PDATA "CVBPSpgmr-- BandPrecData is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
