/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2004-10-08 23:24:44 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 *  -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINBBDPRE module header file (private version)
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _KBBDPRE_IMPL_H
#define _KBBDPRE_IMPL_H

#include "band.h"
#include "kinbbdpre.h"
#include "kinsol_impl.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Definition of KBBDData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* passed by user to KINBBDPrecAlloc, used by pset/psolve functions */

  long int ml, mu;
  KINLocalFn gloc;
  KINCommFn gcomm;

  /* relative error for the Jacobian DQ routine */

  realtype rel_uu;

  /* allocated for use by KINBBDPrecSetup */

  N_Vector vtemp3;

  /* set by KINBBDPrecSetup and used by KINBBDPrecSolve */

  BandMat PP;
  long int *pivots;

  /* set by KINBBDPrecAlloc and used by KINBBDPrecSetup */

  long int n_local;

  /* available for optional output */

  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to KINSol memory */

  KINMem kin_mem;

} *KBBDPrecData;

#endif

#ifdef __cplusplus
}
#endif
