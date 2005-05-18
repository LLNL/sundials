/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2005-05-18 18:17:13 $
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
 * CVBBDPRE error messages
 * -----------------------------------------------------------------
 */

/* CVBBDAlloc error messages */

#define _CVBBDALLOC_        "CVBBDAlloc-- "
#define MSGBBDP_CVMEM_NULL  _CVBBDALLOC_ "Integrator memory is NULL.\n\n"
#define MSGBBDP_BAD_NVECTOR _CVBBDALLOC_ "A required vector operation is not implemented.\n\n"

/* CVBBDPrecGet* error message */

#define MSGBBDP_PDATA_NULL "CVBBDPrecGet*-- BBDPrecData is NULL.\n\n"

/* CVBBDSp* error message */

#define MSGBBDP_NO_PDATA "CVBBDSp*-- BBDPrecData is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
