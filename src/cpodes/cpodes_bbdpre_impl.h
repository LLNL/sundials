/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Implementation header file for the CPBBDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CPBBDPRE_IMPL_H
#define _CPBBDPRE_IMPL_H

#include <cpodes/cpodes_bbdpre.h>
#include <sundials/sundials_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Type: CPBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct CPBBDPrecDataRec {

  /* passed by user to CPBBDPrecAlloc and used by PrecSetup/PrecSolve */
  int mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  CPBBDLocalRhsFn glocE;
  CPBBDLocalResFn glocI;
  CPBBDCommFn cfn;

  /* additional N_Vector needed by cpBBDPrecSetupImpl */
  N_Vector tmp4;

  /* set by CPBBDPrecSetup and used by CPBBDPrecSolve */
  DlsMat savedJ;
  DlsMat savedP;
  long int *pivots;

  /* set by CPBBDPrecAlloc and used by CPBBDPrecSetup */
  int n_local;

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

#define MSGBBD_CPMEM_NULL "Integrator memory is NULL."
#define MSGBBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBBD_MEM_FAIL "A memory request failed."
#define MSGBBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBD_PMEM_NULL   "BBD peconditioner memory is NULL. CPBBDPrecInit must be called."
#define MSGBBD_FUNC_FAILED "The gloc or cfn routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
