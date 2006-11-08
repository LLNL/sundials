/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:05 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the CPBBDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CPBBDPRE_IMPL_H
#define _CPBBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_bbdpre.h>
#include <sundials/sundials_band.h>

/*
 * -----------------------------------------------------------------
 * Type: CPBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* passed by user to CPBBDPrecAlloc and used by PrecSetup/PrecSolve */

  long int mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  CPLocalFn gloc;
  CPCommFn cfn;

  /* set by CPBBDPrecSetup and used by CPBBDPrecSolve */

  BandMat savedJ;
  BandMat savedP;
  long int *pivots;

  /* set by CPBBDPrecAlloc and used by CPBBDPrecSetup */

  long int n_local;

  /* available for optional output */

  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to cpode_mem */

  void *cpode_mem;

} *CPBBDPrecData;

/*
 * -----------------------------------------------------------------
 * CPBBDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBBDP_CPMEM_NULL "Integrator memory is NULL."
#define MSGBBDP_MEM_FAIL "A memory request failed."
#define MSGBBDP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBDP_PDATA_NULL "BBDPrecData is NULL."
#define MSGBBDP_FUNC_FAILED "The gloc or cfn routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
