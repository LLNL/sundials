/*******************************************************************
 *                                                                 *
 * File          : fnvector_serial.h                               *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 06 June 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file (companion of nvector_serial.h) contains the          *
 * definitions needed for the Fortran callable wrappers to         *
 * NV_SpecInit_Serial and NV_SpecFree_Serial (these definitions    *
 * are based on the machine specific information for Fortran       *
 * externals given in the header file fcmixpar.h).                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_serial_h
#define included_fnvector_serial_h

#include "fcmixpar.h" /* Machine specific definitions for Fortran externals */

/* Fortran callable wrappers to NV_SpecInit_Serial and NV_SpecFree_Serial */ 

#if (CRAY)
  
#define F_NVSPECINITS  FNVSPECINITS
#define F_NVSPECFREES  FNVSPECFREES

#elif  (UNDERSCORE)

#define F_NVSPECINITS  fnvspecinits_
#define F_NVSPECFREES  fnvspecfrees_

#else

#define F_NVSPECINITS  fnvspecinits
#define F_NVSPECFREES  fnvspecfrees

#endif

#endif
#ifdef __cplusplus
}
#endif
