/*******************************************************************
 *                                                                 *
 * File          : fnvector_serial.h                               *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 29 March 2002                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file (companion of nvector_serial.h) contains the          *
 * definitions needed for the Fortran callable wrappers to         *
 * M_EnvInit_Serial and M_EnvFree_Serial (these definitions are    *
 * based on the machine specific information for Fortran           *
 * externals given in the header file fcmixpar.h).                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_serial_h
#define included_fnvector_serial_h

#include "fcmixpar.h" /* Machine specific definitions for Fortran externals */

/* Fortran callable wrappers to M_EnvInit_Serial and M_EnvFree_Serial */ 

#if (CRAY)
  
#define F_MENVINITS  FMENVINITS
#define F_MENVFREES  FMENVFREES

#elif  (UNDERSCORE)

#define F_MENVINITS  fmenvinits_
#define F_MENVFREES  fmenvfrees_

#else

#define F_MENVINITS  fmenvinits
#define F_MENVFREES  fmenvfrees

#endif

#endif
#ifdef __cplusplus
}
#endif
