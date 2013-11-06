/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * Implementation header file for the IDAS SuperLUMT linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _IDASSLUMT_IMPL_H
#define _IDASSLUMT_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _SLUMT_H
#define _SLUMT_H
#include "pdsp_defs.h"
#endif
/*
 * -----------------------------------------------------------------
 * Definition of SLUMTData
 * -----------------------------------------------------------------
 */
 
typedef struct SLUMTDataRec {
 
  /* Structure for SuperLUMT-specific data */
 
  SuperMatrix *s_A, *s_AC, *s_L, *s_U, *s_B;
  Gstat_t *Gstat;
  int *perm_r, *perm_c;
  int num_threads; 
  double diag_pivot_thresh; 
  superlumt_options_t *superlumt_options;
 
} *SLUMTData;
 
#ifdef __cplusplus
} 
#endif 
 
#endif 
