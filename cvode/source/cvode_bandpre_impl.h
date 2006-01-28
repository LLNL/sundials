/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-28 00:47:27 $
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

#include "cvode_bandpre.h"
#include "sundials_band.h"

/*
 * -----------------------------------------------------------------
 * Type: CVBandPrecData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* Data set by user in CVBandPrecAlloc */

  long int N;
  long int ml, mu;

  /* Data set by CVBandPrecSetup */

  BandMat savedJ;
  BandMat savedP;
  long int *pivots;

  /* Rhs calls */

  long int nfeBP;

  /* Pointer to cvode_mem */

  void *cvode_mem;

} *CVBandPrecData;

/*
 * -----------------------------------------------------------------
 * CVBANDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBP_CVMEM_NULL "Integrator memory is NULL."
#define MSGBP_MEM_FAIL "A memory request failed."
#define MSGBP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBP_PDATA_NULL "BandPrecData is NULL."


#ifdef __cplusplus
}
#endif

#endif
