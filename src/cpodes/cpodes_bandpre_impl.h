/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2007-12-19 20:26:42 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * begincopyright(llns)
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * endcopyright(llns)
 * -----------------------------------------------------------------
 * Implementation header file for the CPBANDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CPBANDPRE_IMPL_H
#define _CPBANDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_bandpre.h>
#include <sundials/sundials_band.h>

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
