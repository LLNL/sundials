/*******************************************************************
 *                                                                 *
 * File          : fnvector_parallel.h                             *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 06 June 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file (companion of nvector_parallel.h) contains the        *
 * definitions needed for the Fortran callable wrappers to         *
 * NV_SpecInit_Parallel and NV_SpecFree_Parallel (these definitions*
 * are based on the machine specific information for Fortran       *
 * externals given in the header file fcmixpar.h).                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_parallel_h
#define included_fnvector_parallel_h

#include "fcmixpar.h" /* Machine specific definitions for Fortran externals */

/* Fortran callable wrappers to NV_SpecInit_Parallel and NV_SpecFree_Parallel */ 

#if (CRAY)
  
#define F_NVSPECINITP  FNVSPECINITP
#define F_NVSPECFREEP  FNVSPECFREEP

#elif  (UNDERSCORE)

#define F_NVSPECINITP  fnvspecinitp_
#define F_NVSPECFREEP  fnvspecfreep_

#else

#define F_NVSPECINITP  fnvspecinitp
#define F_NVSPECFREEP  fnvspecfreep

#endif


#endif
#ifdef __cplusplus
}
#endif
