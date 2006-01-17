/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-17 23:30:21 $
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

  /* CVBBDAlloc error messages */

#define _CVBBDALLOC_        "CVBBDAlloc-- "
#define MSGBBDP_CVMEM_NULL  _CVBBDALLOC_ "Integrator memory is NULL.\n\n"
#define MSGBBDP_BAD_NVECTOR _CVBBDALLOC_ "A required vector operation is not implemented.\n\n"

  /* CVBBDPrecGet* error message */

#define MSGBBDP_PDATA_NULL "CVBBDPrecReInit/CVBBDPrecGet*-- BBDPrecData is NULL.\n\n"

  /* CVBBDSp* error message */

#define MSGBBDP_NO_PDATA "CVBBDSp*-- BBDPrecData is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
