/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-06-02 23:22:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 *  -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/kinsol/LICENSE
 * -----------------------------------------------------------------
 * KINBBDPRE module header file (private version)
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _kbbdpre_impl_h
#define _kbbdpre_impl_h

#include "kinsol_impl.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "band.h"
#include "kinbbdpre.h"

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
