/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2004-07-22 21:21:16 $
 * ----------------------------------------------------------------- 
 * Programmers: Michael Wittman, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the CVBBDPRE module.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvbbdpre_impl_h
#define _cvbbdpre_impl_h

#include "cvbbdpre.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type: CVBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* passed by user to CVBBDPrecAlloc, used by PrecSetup/PrecSolve */
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

  /* available for optional output: */
  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* Pointer to cvode_mem */
  void *cvode_mem;

} *CVBBDPrecData;

#endif

#ifdef __cplusplus
}
#endif
