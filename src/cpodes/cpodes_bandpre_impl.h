/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Implementation header file for the CPBANDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CPBANDPRE_IMPL_H
#define _CPBANDPRE_IMPL_H

#include <cpodes/cpodes_bandpre.h>
#include <sundials/sundials_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Type: CPBandPrecData
 * -----------------------------------------------------------------
 */

typedef struct CPBandPrecDataRec {

  /* Data set by user in CPBandPrecInit */
  int N;
  int ml, mu;

  /* Data set by CPBandPrecSetup */
  DlsMat savedJ;
  DlsMat savedP;
  long int *pivots;

  /* Function evaluations for DQ Jacobian approximation */
  long int nfeBP;

  /* Pointer to cpode_mem */
  void *cpode_mem;

} *CPBandPrecData;

/*
 * -----------------------------------------------------------------
 * CPBANDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBP_CPMEM_NULL "Integrator memory is NULL."
#define MSGBP_LMEM_NULL "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBP_MEM_FAIL "A memory request failed."
#define MSGBP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBP_PMEM_NULL "Band preconditioner memory is NULL. CPBandPrecInit must be called."
#define MSGBP_FUNC_FAILED "The ODE function failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
