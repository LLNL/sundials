/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2005-03-19 00:10:51 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the
 * definitions needed for the initialization of parallel
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_PARALLEL_H
#define _FNVECTOR_PARALLEL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include "sundials_config.h"
#endif

#if defined(F77_FUNC)

#define FNV_INITP F77_FUNC(fnvinitp, FNVINITP)
#define FNV_FREEP F77_FUNC(fnvfreep, FNVFREEP)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITP fnvinitp
#define FNV_FREEP fnvfreep

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITP FNVINITP
#define FNV_FREEP FNVFREEP

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITP fnvinitp_
#define FNV_FREEP fnvfreep_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITP FNVINITP_
#define FNV_FREEP FNVFREEP_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITP fnvinitp__
#define FNV_FREEP fnvfreep__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITP FNVINITP__
#define FNV_FREEP FNVFREEP__

#endif

/* Prototypes of exported funcions */

void FNV_INITP(long int *nlocal, long int *nglobal, int *ier);
void FNV_FREEP(void);

#ifdef __cplusplus
}
#endif

#endif
