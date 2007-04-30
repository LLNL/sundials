/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007-04-30 17:43:09 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDABBDPRE
 * module, for a band-block-diagonal preconditioner, i.e. a
 * block-diagonal matrix with banded blocks, for use with IDA
 * and an IDASPILS linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _IDABBDPRE_IMPL_H
#define _IDABBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <ida/ida_bbdpre.h>
#include <sundials/sundials_band.h>

/*
 * -----------------------------------------------------------------
 * Definition of IBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct IBBDPrecDataRec {

  /* passed by user to IDABBDPrecAlloc and used by
     IDABBDPrecSetup/IDABBDPrecSolve functions */

  int mudq, mldq, mukeep, mlkeep;
  realtype rel_yy;
  IDABBDLocalFn glocal;
  IDABBDCommFn gcomm;

  /* allocated for use by IDABBDPrecSetup */

  N_Vector tempv4;

  /* set by IDABBDPrecon and used by IDABBDPrecSolve */

  DlsMat PP;
  int *pivots;

  /* set by IDABBDPrecAlloc and used by IDABBDPrecSetup */

  int n_local;

  /* available for optional output */

  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to ida_mem */

  void *ida_mem;

} *IBBDPrecData;

/*
 * -----------------------------------------------------------------
 * IDABBDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBBD_MEM_NULL    "Integrator memory is NULL."
#define MSGBBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBBD_MEM_FAIL    "A memory request failed."
#define MSGBBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBD_PMEM_NULL   "BBD peconditioner memory is NULL. IDABBDPrecInit must be called."
#define MSGBBD_FUNC_FAILED "The Glocal or Gcomm routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
