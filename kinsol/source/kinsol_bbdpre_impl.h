/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-25 22:18:39 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 *  -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINBBDPRE module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _KBBDPRE_IMPL_H
#define _KBBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinsol_bbdpre.h"
#include "sundials_band.h"
#include "kinsol_impl.h"

/*
 * -----------------------------------------------------------------
 * Definition of KBBDData
 * -----------------------------------------------------------------
 */

typedef struct {

  /* passed by user to KINBBDPrecAlloc, used by pset/psolve functions */

  long int mudq, mldq, mukeep, mlkeep;
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

/*
 *-----------------------------------------------------------------
 * KINBBDPRE error messages
 *-----------------------------------------------------------------
 */

/* KINBBDPrecAlloc error messages */

#define KINBBDALLOC        "KINBBDPrecAlloc-- "
#define MSGBBD_KINMEM_NULL KINBBDALLOC "KINSOL Memory is NULL.\n\n"
#define MSGBBD_BAD_NVECTOR KINBBDALLOC "A required vector operation is not implemented.\n\n"

/* KINBBDPrecGet* error message */

#define MSGBBD_PDATA_NULL "KINBBDPrecGet*-- KBBDPrecData is NULL. \n\n"

/* KINBBDSpgmr/KINBBDSpbcg error message */

#define MSGBBD_NO_PDATA "KINBBDSp*-- KBBDPrecData is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
