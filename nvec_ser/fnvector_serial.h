/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2004-07-22 21:12:43 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban, LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/shared/LICENSE
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the 
 * definitions needed for the Fortran initialization of serial 
 * vector operations
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_serial_h
#define included_fnvector_serial_h

/* Fortran callable wrappers to NV_SpecInit_Serial and NV_SpecFree_Serial */ 

#if defined(SUNDIALS_UNDERSCORE_NONE)

#define FNV_INITS  fnvinits
#define FNV_FREES  fnvfrees

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#define FNV_INITS  fnvinits__
#define FNV_FREES  fnvfrees__

#else

#define FNV_INITS  fnvinits_
#define FNV_FREES  fnvfrees_

#endif

#endif
#ifdef __cplusplus
}
#endif
