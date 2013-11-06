/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2007-12-19 20:26:42 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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
