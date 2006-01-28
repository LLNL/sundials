/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-01-28 00:47:17 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the CVBBDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CVSBBDPRE_IMPL_H
#define _CVSBBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_bbdpre.h"
#include "sundials_band.h"

  /*
   * -----------------------------------------------------------------
   * Type: CVBBDPrecData
   * -----------------------------------------------------------------
   */

  typedef struct {

    /* passed by user to CVBBDPrecAlloc and used by PrecSetup/PrecSolve */

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

    /* available for optional output */

    long int rpwsize;
    long int ipwsize;
    long int nge;

    /* pointer to cvode_mem */

    void *cvode_mem;

  } *CVBBDPrecData;


  /*
   * -----------------------------------------------------------------
   * Type: CVBBDPrecDataB
   * -----------------------------------------------------------------
   */

  typedef struct {

    /* BBD user functions (glocB and cfnB) for backward run */
    CVLocalFnB glocB;
    CVCommFnB  cfnB;
    
    /* BBD prec data */
    void *bbd_dataB;

  } *CVBBDPrecDataB;

  /*
   * -----------------------------------------------------------------
   * CVBBDPRE error messages
   * -----------------------------------------------------------------
   */

#define MSGBBDP_CVMEM_NULL "Integrator memory is NULL."
#define MSGBBDP_MEM_FAIL "A memory request failed."
#define MSGBBDP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBDP_PDATA_NULL "CVBBDPRE memory is NULL."

#define MSGBBDP_CAMEM_NULL "cvadj_mem = NULL illegal."
#define MSGBBDP_PMEMB_NULL "CVBBDPRE memory is NULL for the backward integration."
#define MSGBBDP_BAD_T "Bad t for interpolation. Abort!"


#ifdef __cplusplus
}
#endif

#endif
