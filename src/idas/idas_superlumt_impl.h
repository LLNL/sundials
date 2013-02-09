/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2013, Lawrence Livermore National Security
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
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
#include "/g/g16/carol/OEISU/SuperLU_MT_2.0/SRC/pdsp_defs.h"
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
