/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-05-26 18:37:06 $
 * ----------------------------------------------------------------- 
 * Programmers: Michael Wittman, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the CVBANDPRE module.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvbandpre_impl_h
#define _cvbandpre_impl_h

#include "cvbandpre.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type: CVBandPrecData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* Data set by user in CVBandPrecAlloc: */
  long int N;
  long int ml, mu;

  /* Data set by CVBandPrecSetup: */
  BandMat savedJ;
  BandMat savedP;
  long int *pivots;

  /* Rhs calls */
  long int nfeBP;

  /* Pointer to cvode_mem */
  void *cvode_mem;

} *CVBandPrecData;

#endif

#ifdef __cplusplus
}
#endif
