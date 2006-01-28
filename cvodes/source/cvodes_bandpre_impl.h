/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-28 00:47:17 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the CVBANDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CVSBANDPRE_IMPL_H
#define _CVSBANDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_bandpre.h"
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
#define MSGBP_PDATA_NULL "CVBANDPRE memory is NULL."

#define MSGBP_CAMEM_NULL "cvadj_mem = NULL illegal."

#ifdef __cplusplus
}
#endif

#endif
