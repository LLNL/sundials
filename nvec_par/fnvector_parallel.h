/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2004-07-26 17:27:04 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban, LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/shared/LICENSE
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the 
 * definitions needed for the Fortran initialization of parallel
 * vector operations
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_parallel_h
#define included_fnvector_parallel_h

/* Fortran callable wrappers to NV_SpecInit_Parallel and NV_SpecFree_Parallel */ 

#if defined(SUNDIALS_UNDERSCORE_NONE)

#if defined(SUNDIALS_CASE_LOWER)

#define FNV_INITP  fnvinitp
#define FNV_FREEP  fnvfreep

#elif defined(SUNDIALS_CASE_UPPER)

#define FNV_INITP  FNVINITP
#define FNV_FREEP  FNVFREEP

#endif

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#if defined(SUNDIALS_CASE_LOWER)

#define FNV_INITP  fnvinitp__
#define FNV_FREEP  fnvfreep__

#elif defined(SUNDIALS_CASE_UPPER)

#define FNV_INITP  FNVINITP__
#define FNV_FREEP  FNVFREEP__

#endif

#else

#if defined(SUNDIALS_CASE_LOWER)

#define FNV_INITP  fnvinitp_
#define FNV_FREEP  fnvfreep_

#elif defined(SUNDIALS_CASE_UPPER)

#define FNV_INITP  FNVINITP_
#define FNV_FREEP  FNVFREEP_

#endif

#endif


#endif
#ifdef __cplusplus
}
#endif
