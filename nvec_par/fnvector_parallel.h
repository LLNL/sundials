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
 * NV_SpecInit_Parallel and NV_SpecFree_Parallel                   *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef included_fnvector_parallel_h
#define included_fnvector_parallel_h

/* Fortran callable wrappers to NV_SpecInit_Parallel and NV_SpecFree_Parallel */ 

#if defined(SUNDIALS_UNDERSCORE_NONE)

#define FNV_INITP  fnvinitp
#define FNV_FREEP  fnvfreep

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#define FNV_INITP  fnvinitp__
#define FNV_FREEP  fnvfreep__

#else

#define FNV_INITP  fnvinitp_
#define FNV_FREEP  fnvfreep_

#endif


#endif
#ifdef __cplusplus
}
#endif
