/*******************************************************************
 *                                                                 *
 * File          : fnvector_parallel.h                             *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 29 March 2002                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file (companion of nvector_parallel.h) contains the        *
 * definitions needed for the Fortran callable wrappers to         *
 * M_EnvInit_Parallel and M_EnvFree_Parallel (these definitions    *
 * are based on the machine specific information for Fortran       *
 * externals given in the header file fcmixpar.h).                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_parallel_h
#define included_fnvector_parallel_h

#include "fcmixpar.h" /* Machine specific definitions for Fortran externals */

/* Fortran callable wrappers to M_EnvInit_Parallel and M_EnvFree_Parallel */ 

#if (CRAY)
  
#define F_MENVINITP  FMENVINITP
#define F_MENVFREEP  FMENVFREEP

#elif  (UNDERSCORE)

#define F_MENVINITP  fmenvinitp_
#define F_MENVFREEP  fmenvfreep_

#else

#define F_MENVINITP  fmenvinitp
#define F_MENVFREEP  fmenvfreep

#endif


#endif
#ifdef __cplusplus
}
#endif
