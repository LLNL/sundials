/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2004-07-26 17:26:53 $
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

#if defined(SUNDIALS_CASE_LOWER)

#define FNV_INITS  fnvinits
#define FNV_FREES  fnvfrees

#elif defined(SUNDIALS_CASE_UPPER)

#define FNV_INITS  FNVINITS
#define FNV_FREES  FNVFREES

#endif

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#if defined(SUNDIALS_CASE_LOWER)

#define FNV_INITS  fnvinits__
#define FNV_FREES  fnvfrees__

#elif defined(SUNDIALS_CASE_UPPER)

#define FNV_INITS  FNVINITS__
#define FNV_FREES  FNVFREES__

#endif

#else

#if defined(SUNDIALS_CASE_LOWER)

#define FNV_INITS  fnvinits_
#define FNV_FREES  fnvfrees_

#elif defined(SUNDIALS_CASE_UPPER)

#define FNV_INITS  FNVINITS_
#define FNV_FREES  FNVFREES_

#endif

#endif

#endif
#ifdef __cplusplus
}
#endif
