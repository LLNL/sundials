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
 * Implementation header file for the IDAS KLU linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _IDASKLU_IMPL_H
#define _IDASKLU_IMPL_H

#ifndef _S_KLU_H
#define _S_KLU_H
#include "/home/carol/KLU/KLU/Include/klu.h"
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Definition of KLUData
 * -----------------------------------------------------------------
 */
 
typedef struct KLUDataRec {
 
  /* Structure for KLU-specific data */
 
  klu_symbolic *s_Symbolic;
  klu_numeric  *s_Numeric;
  klu_common    s_Common;
 
} *KLUData;
 
#ifdef __cplusplus
} 
#endif 
 
#endif 
