/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2004-06-18 21:36:13 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/ida/LICENSE
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDABBDPRE
 * module, for a band-block-diagonal preconditioner, i.e. a
 * block-diagonal matrix with banded blocks, for use with IDA/IDAS
 * and IDASpgmr.  
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _ibbdpre_impl_h
#define _ibbdpre_impl_h

#include "idabbdpre.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

/*********************** Definition of IBBDPrecData *****************/

typedef struct {

  /* passed by user to IDABBDPrecAlloc, used by 
     IDABBDPrecSetup/IDABBDPrecSolve functions: */
  long int mudq, mldq, mukeep, mlkeep;
  realtype rel_yy;
  IDALocalFn glocal;
  IDACommFn gcomm;

 /* allocated for use by IDABBDPrecSetup */
  N_Vector tempv4;

  /* set by IDABBDPrecon and used by IDABBDPrecSolve: */
  BandMat PP;
  long int *pivots;

  /* set by IDABBDPrecAlloc and used by IDABBDPrecSetup */
  long int n_local;

  /* available for optional output: */
  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* Pointer to ida_mem */
  IDAMem IDA_mem;

} *IBBDPrecData;

#endif

#ifdef __cplusplus
}
#endif
